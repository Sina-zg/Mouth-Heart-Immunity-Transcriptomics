#----Loading Package----
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

#----File of TCR Paths for NICM and PD----
paths <- list(
  PD_High = c(
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_PD1.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_PD5.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_PD6.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_PD7.csv"
  ),
  PD_Low = c(
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_PD2.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_PD3.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_PD4.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_PD8.csv"
  ),
  Prog = c(
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_CV5.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_CV12.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_CV109.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_CV110.csv"
  ),
  Surv = c(
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_CV39.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_CV87.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_CV103.csv",
    "/Users/sina_zg/Desktop/Sina_air/Sequencing Analysis/filtered_contig_annotations_CV106.csv"
  )
)

#----Load and label----
file_list  <- unlist(paths, use.names = FALSE)
sample_ids <- sub("^.*filtered_contig_annotations_(.*)\\.csv$", "\\1", basename(file_list))
conditions <- rep(names(paths), lengths(paths))

all_contigs <- bind_rows(lapply(seq_along(file_list), function(i){
  df <- read.csv(file_list[i], stringsAsFactors = FALSE, check.names = FALSE)
  df$sample_id <- sample_ids[i]
  df$condition <- conditions[i]
  df
}))


contigs_prod <- all_contigs %>% filter(tolower(productive) == "true")

#----Audit----
cat("Columns in your contig files:\n"); print(sort(names(all_contigs)))


std_10x_cols <- function(df) {
  
  pick <- function(cands) {
    nm <- cands[cands %in% names(df)][1]
    if (is.na(nm)) NA_character_ else nm
  }
  
  
  nm_barcode  <- pick(c("barcode","cell_id","Cell","cell"))
  nm_chain    <- pick(c("chain","locus"))
  nm_cdr3aa   <- pick(c("cdr3_aa","cdr3","amino_acid","aa_sequence"))
  nm_cdr3nt   <- pick(c("cdr3_nt","cdr3.nucleotide","nucleotide"))
  nm_v        <- pick(c("v_gene","v_call","TRAV","TRBV","V_gene","bestVGene"))
  nm_j        <- pick(c("j_gene","j_call","TRAJ","TRBJ","J_gene","bestJGene"))
  nm_product  <- pick(c("productive","productive.reads","is_productive"))
  nm_reads    <- pick(c("reads","read_count"))
  nm_umis     <- pick(c("umis","umi_count","umi"))
  
  
  rn <- c(
    setNames(nm_barcode, "barcode"),
    setNames(nm_chain,   "chain"),
    setNames(nm_cdr3aa,  "cdr3_aa"),
    setNames(nm_cdr3nt,  "cdr3_nt"),
    setNames(nm_v,       "v_gene"),
    setNames(nm_j,       "j_gene"),
    setNames(nm_product, "productive"),
    setNames(nm_reads,   "reads"),
    setNames(nm_umis,    "umis")
  )
  rn <- rn[!is.na(rn)]            
  df  <- dplyr::rename(df, !!!rn)
  
  
  if (!"cdr3_aa" %in% names(df) && "cdr3" %in% names(df)) {
    df <- dplyr::rename(df, cdr3_aa = cdr3)
  }
  
  df
}

all_contigs <- std_10x_cols(all_contigs)

#----Audit----
needed <- c("barcode","chain","cdr3_aa","v_gene","j_gene")
missing <- setdiff(needed, names(all_contigs))
if (length(missing)) {
  stop("Missing required columns after standardization: ",
       paste(missing, collapse = ", "),
       "\nPlease print names(all_contigs) and share; we’ll map them.")
}


contigs_prod <- all_contigs %>%
  dplyr::mutate(productive = tolower(as.character(productive))) %>%
  dplyr::filter(productive %in% c("true","t","1"))



#----Sanity check----
stopifnot(all(c("barcode","chain","cdr3_aa","v_gene","j_gene") %in% names(contigs_prod)))

#----α and one β per cell----
pick_first_present <- function(cands, nms) {
  nm <- cands[cands %in% nms][1]; if (is.na(nm)) return(NA_character_); nm
}
umi_col  <- pick_first_present(c("umis","umi_count","umi"), names(contigs_prod))
read_col <- pick_first_present(c("reads","read_count"), names(contigs_prod))
if (is.na(umi_col) || is.na(read_col)) {
  message("UMI/reads columns not found — will resolve using row order within (sample, barcode, chain).")
  contigs_prod$umi_val  <- 0
  contigs_prod$read_val <- 0
} else {
  contigs_prod$umi_val  <- suppressWarnings(as.numeric(contigs_prod[[umi_col]]))
  contigs_prod$read_val <- suppressWarnings(as.numeric(contigs_prod[[read_col]]))
}

ranked <- contigs_prod %>%
  group_by(sample_id, condition, barcode, chain) %>%
  arrange(desc(umi_val), desc(read_val), .by_group = TRUE) %>%
  slice(1) %>%    
  ungroup()


resolved <- ranked %>%
  filter(chain %in% c("TRA","TRB")) %>%
  select(sample_id, condition, barcode, chain, cdr3_aa, v_gene, j_gene) %>%
  tidyr::pivot_wider(
    id_cols = c(sample_id, condition, barcode),
    names_from = chain,
    values_from = c(cdr3_aa, v_gene, j_gene),
    names_sep = "."
  )

#----Clonotypes per subject----
strip_allele <- function(x) sub("\\*.*$", "", ifelse(is.na(x), "", x))

clono_tbl <- resolved %>%
  transmute(
    sample_id,
    condition,
    CDR3a = cdr3_aa.TRA %||% "",
    CDR3b = cdr3_aa.TRB %||% "",
    Va    = strip_allele(v_gene.TRA),
    Ja    = strip_allele(j_gene.TRA),
    Vb    = strip_allele(v_gene.TRB),
    Jb    = strip_allele(j_gene.TRB)
  ) %>%
  filter(CDR3b != "") 


clono_counts <- clono_tbl %>%
  group_by(sample_id, condition, CDR3a, CDR3b, Vb, Jb) %>%  
  summarise(Frequency = n(), .groups = "drop")

#----Draft GLIPH2 columns----
gliph2_ready <- clono_counts %>%
  transmute(
    CDR3b = CDR3b,
    Vb    = Vb,
    Jb    = Jb,
    CDR3a = CDR3a,                                
    `Subject:condition` = paste0(sample_id, ":", condition),
    Frequency = Frequency
  )


outdir <- "gliph2_inputs"
dir.create(outdir, showWarnings = FALSE)

run1 <- gliph2_ready %>% filter(grepl(":(Prog|PD_High)$", `Subject:condition`))
run2 <- gliph2_ready %>% filter(grepl(":(Surv|PD_Low)$",  `Subject:condition`))
all_in_one <- gliph2_ready

write_tsv(run1,      file.path(outdir, "GLIPH2_Prog_vs_PDHigh.txt"))
write_tsv(run2,      file.path(outdir, "GLIPH2_Surv_vs_PDLow.txt"))
write_tsv(all_in_one,file.path(outdir, "GLIPH2_all.txt"))

#----Simple summary for Audit----
summary_tbl <- gliph2_ready %>%
  tidyr::separate(`Subject:condition`, into = c("Subject","Cond"), sep = ":", remove = FALSE) %>%
  count(Cond, name = "n_clonotypes") %>%
  arrange(Cond)

print(summary_tbl)
cat("\nWrote:\n",
    file.path(outdir, "GLIPH2_Prog_vs_PDHigh.txt"), "\n",
    file.path(outdir, "GLIPH2_Surv_vs_PDLow.txt"), "\n",
    file.path(outdir, "GLIPH2_all.txt"), "\n")


#----Saving & Running GLIPH2----
desktop_outdir <- file.path(path.expand("~/Desktop"), "GLIPH2_inputs")
dir.create(desktop_outdir, showWarnings = FALSE, recursive = TRUE)

run1      <- gliph2_ready %>% dplyr::filter(grepl(":(Prog|PD_High)$", `Subject:condition`))
run2      <- gliph2_ready %>% dplyr::filter(grepl(":(Surv|PD_Low)$",  `Subject:condition`))
all_in_one <- gliph2_ready

readr::write_tsv(run1,      file.path(desktop_outdir, "GLIPH2_Prog_vs_PDHigh.txt"))
readr::write_tsv(run2,      file.path(desktop_outdir, "GLIPH2_Surv_vs_PDLow.txt"))
readr::write_tsv(all_in_one,file.path(desktop_outdir, "GLIPH2_all.txt"))

cat("Wrote to Desktop/GLIPH2_inputs:\n",
    file.path(desktop_outdir, "GLIPH2_Prog_vs_PDHigh.txt"), "\n",
    file.path(desktop_outdir, "GLIPH2_Surv_vs_PDLow.txt"), "\n",
    file.path(desktop_outdir, "GLIPH2_all.txt"), "\n")




#----Export GLIPH2 and Motif Analysis----
library(readr); library(dplyr); library(stringr)

run1_file <- "~/Desktop/GLIPH2_inputs/gliph2_run1.csv"  
x1 <- read_csv(run1_file, show_col_types = FALSE)

#----Safety---
stopifnot(all(c("pattern","Sample") %in% names(x1)))   
x1 <- x1 %>% filter(!is.na(pattern), !is.na(Sample))

#----classify groups----
x1 <- x1 %>%
  mutate(class = case_when(
    str_starts(Sample, "PD") ~ "PD",
    str_starts(Sample, "CV") ~ "CV",
    TRUE ~ "Other"
  ))

#----level presence per motif----
m1 <- x1 %>%
  filter(class %in% c("PD","CV")) %>%
  group_by(pattern, class, Sample) %>%
  summarise(present = TRUE, .groups = "drop") %>%
  group_by(pattern, class) %>%
  summarise(n_donors = n_distinct(Sample), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = class, values_from = n_donors, values_fill = 0)

#----keep motifs present in ≥2 PD  AND ≥2 Progressor----
m1_kept <- m1 %>% filter(PD >= 2, CV >= 2)

cat("Run1 (Prog + PD High): motifs with ≥2 PD donors AND ≥2 Progressor donors = ",
    nrow(m1_kept), "\n", sep = "")


out1 <- "~/Desktop/GLIPH2_results/run1_motifs_PD>=2_AND_Prog>=2.csv"
dir.create(dirname(out1), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(m1_kept %>% arrange(desc(PD + CV)), out1)
cat("Saved details to: ", out1, "\n", sep = "")


#----Survivor vs Progressors----
run2_file <- "~/Desktop/GLIPH2_inputs/gliph2_run2.csv"  
x2 <- read_csv(run2_file, show_col_types = FALSE)

stopifnot(all(c("pattern","Sample") %in% names(x2)))
x2 <- x2 %>% filter(!is.na(pattern), !is.na(Sample))


x2 <- x2 %>%
  mutate(class = case_when(
    str_starts(Sample, "PD") ~ "PD",
    str_starts(Sample, "CV") ~ "CV",
    TRUE ~ "Other"
  ))

m2 <- x2 %>%
  filter(class %in% c("PD","CV")) %>%
  group_by(pattern, class, Sample) %>%
  summarise(present = TRUE, .groups = "drop") %>%
  group_by(pattern, class) %>%
  summarise(n_donors = n_distinct(Sample), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = class, values_from = n_donors, values_fill = 0)

#----keep motifs present in ≥2 PDLow AND ≥2 Survivor----
m2_kept <- m2 %>% filter(PD >= 2, CV >= 2)

cat("Run2 (Survivor + PD Low): motifs with ≥2 PD donors AND ≥2 Survivor donors = ",
    nrow(m2_kept), "\n", sep = "")

out2 <- "~/Desktop/GLIPH2_results/run2_motifs_PD>=2_AND_Surv>=2.csv"
dir.create(dirname(out2), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(m2_kept %>% arrange(desc(PD + CV)), out2)
cat("Saved details to: ", out2, "\n", sep = "")

suppressPackageStartupMessages({library(dplyr); library(readr)})

res_dir <- "~/Desktop/GLIPH2_results"
run1_csv <- file.path(res_dir, "run1_motifs_PD>=2_AND_Prog>=2.csv")
run2_csv <- file.path(res_dir, "run2_motifs_PD>=2_AND_Surv>=2.csv")

#----Loading----
if (exists("m1_kept")) {
  s1 <- unique(m1_kept$pattern)
} else {
  stopifnot(file.exists(run1_csv))
  s1 <- unique(readr::read_csv(run1_csv, show_col_types = FALSE)$pattern)
}

if (exists("m2_kept")) {
  s2 <- unique(m2_kept$pattern)
} else {
  stopifnot(file.exists(run2_csv))
  s2 <- unique(readr::read_csv(run2_csv, show_col_types = FALSE)$pattern)
}

#----Sets----
shared      <- intersect(s1, s2)
unique_g1   <- setdiff(s1, s2)   
unique_g2   <- setdiff(s2, s1)   

#----Preview----
cat("\n=== Motif overlap (robust: ≥2 PD donors & ≥2 CV donors per run) ===\n")
cat("Group1 only (Prog + PD High): ", length(unique_g1), "\n")
cat("Group2 only (Surv + PD Low) : ", length(unique_g2), "\n")
cat("Shared between groups       : ", length(shared),   "\n")

#----Saving----
readr::write_tsv(tibble(pattern = sort(unique_g1)), file.path(res_dir, "motifs_unique_to_Group1.txt"))
readr::write_tsv(tibble(pattern = sort(unique_g2)), file.path(res_dir, "motifs_unique_to_Group2.txt"))
readr::write_tsv(tibble(pattern = sort(shared)),    file.path(res_dir, "motifs_shared_Group1_and_Group2.txt"))

cat("\nSaved lists to:\n", res_dir, "\n")


#----Inputs----
run3_file <- "~/Desktop/GLIPH2_inputs/gliph2_run3.csv"    
res_dir   <- "~/Desktop/GLIPH2_results"
fig_dir   <- "~/Desktop/GLIPH2_figures"; dir.create(fig_dir, FALSE, TRUE)

#----Load Motifs----
if (!exists("m1_kept")) m1_kept <- read_csv(file.path(res_dir, "run1_motifs_PD>=2_AND_Prog>=2.csv"), show_col_types = FALSE)
if (!exists("m2_kept")) m2_kept <- read_csv(file.path(res_dir, "run2_motifs_PD>=2_AND_Surv>=2.csv"), show_col_types = FALSE)


#----Applying Fisher and Length Score----
pick_score_col <- function(df, which = c("fisher","length")) {
  which <- match.arg(which)
  pat <- if (which == "fisher") "(?i)^fisher[_\\.]?score$" else "(?i)^length[_\\.]?score$"
  nm  <- names(df)[grepl(pat, names(df), perl = TRUE)]
  if (length(nm)) nm[1] else NA_character_
}


extract_scores <- function(df) {
  cf <- pick_score_col(df, "fisher")
  cl <- pick_score_col(df, "length")
  if (is.na(cf) || is.na(cl)) {
    stop("Could not find Fisher_Score and/or Length_Score columns in this run. ",
         "Found columns: ", paste(names(df), collapse=", "))
  }
  df %>%
    dplyr::filter(!is.na(pattern)) %>%
    dplyr::transmute(pattern,
                     Fisher_Score = suppressWarnings(as.numeric(.data[[cf]])),
                     Length_Score = suppressWarnings(as.numeric(.data[[cl]]))) %>%
    dplyr::group_by(pattern) %>%
    dplyr::summarise(Fisher_Score = min(Fisher_Score, na.rm = TRUE),
                     Length_Score = min(Length_Score, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::filter(is.finite(Fisher_Score), is.finite(Length_Score))
}

#----Score Table----
scores_run1 <- extract_scores(x1)   
scores_run2 <- extract_scores(x2)  

#----Apply the Score and Individual Counts----
m1_kept <- m1 %>%
  dplyr::filter(PD >= 2, CV >= 2) %>%
  dplyr::inner_join(scores_run1 %>%
                      dplyr::filter(Fisher_Score < 0.05, Length_Score < 0.05),
                    by = "pattern")

m2_kept <- m2 %>%
  dplyr::filter(PD >= 2, CV >= 2) %>%
  dplyr::inner_join(scores_run2 %>%
                      dplyr::filter(Fisher_Score < 0.05, Length_Score < 0.05),
                    by = "pattern")

#----Sanity Check----
cat("Run1 kept after scores: ", nrow(m1_kept), " motifs\n", sep = "")
cat("Run2 kept after scores: ", nrow(m2_kept), " motifs\n", sep = "")

#----Saving----
readr::write_csv(m1_kept %>% dplyr::arrange(Fisher_Score, Length_Score),
                 "~/Desktop/GLIPH2_results/run1_motifs_PD>=2_AND_Prog>=2_scores<0.05.csv")
readr::write_csv(m2_kept %>% dplyr::arrange(Fisher_Score, Length_Score),
                 "~/Desktop/GLIPH2_results/run2_motifs_PD>=2_AND_Surv>=2_scores<0.05.csv")


s1 <- unique(m1_kept$pattern)  
s2 <- unique(m2_kept$pattern)  

shared_patterns    <- intersect(s1, s2)
unique_g1_patterns <- setdiff(s1, s2)
unique_g2_patterns <- setdiff(s2, s1)


s1 <- unique(m1_kept$pattern)  
s2 <- unique(m2_kept$pattern)  
shared_patterns      <- intersect(s1, s2)     
unique_g1_patterns   <- setdiff(s1, s2)
unique_g2_patterns   <- setdiff(s2, s1)

pick_freq_col <- function(df) {
  cand <- c("Frequency","Freq","TcRb","n","N","count","Counts","CloneCount")
  nm <- intersect(cand, names(df))
  if (length(nm)) nm[1] else NA_character_
}

#----Shared Motifs Progressor and PD High----
freq1 <- pick_freq_col(x1)
x1_std <- x1 %>%
  filter(!is.na(pattern), !is.na(Sample)) %>%
  mutate(
    condition = ifelse(str_starts(Sample, "CV"), "Prog", "PD_High"),
    Subject   = Sample,
    Frequency = if (!is.na(freq1)) suppressWarnings(as.numeric(.data[[freq1]])) else 1
  ) %>%
  select(Subject, condition, pattern, Frequency)

# ----Shared Motifs Survivor and PD Low----
freq2 <- pick_freq_col(x2)
x2_std <- x2 %>%
  filter(!is.na(pattern), !is.na(Sample)) %>%
  mutate(
    condition = ifelse(str_starts(Sample, "CV"), "Surv", "PD_Low"),
    Subject   = Sample,
    Frequency = if (!is.na(freq2)) suppressWarnings(as.numeric(.data[[freq2]])) else 1
  ) %>%
  select(Subject, condition, pattern, Frequency)

#----Combine----
joined_with_pattern <- bind_rows(x1_std, x2_std) %>%
  mutate(
    Panel = dplyr::case_when(
      condition == "Prog"    ~ "Progressor",
      condition == "PD_High" ~ "PD High CTL",
      condition == "PD_Low"  ~ "PD Low CTL",
      condition == "Surv"    ~ "Survivor",
      TRUE ~ condition
    ),
    SuperGroup = dplyr::case_when(
      Panel %in% c("Progressor","PD High CTL") ~ "Group1 (Prog + PD High CTL)",
      Panel %in% c("PD Low CTL","Survivor")    ~ "Group2 (Surv + PD Low CTL)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Panel), !is.na(SuperGroup), !is.na(pattern))

#----Frequencies per supergroup and motif----
counts_by_pattern <- joined_with_pattern %>%
  group_by(SuperGroup, pattern) %>%
  summarise(Frequency = sum(Frequency, na.rm = TRUE), .groups = "drop")


g2_freq <- joined_with_pattern %>%
  filter(Panel %in% c("PD Low CTL","Survivor")) %>%   
  group_by(pattern) %>% summarise(freq_g2 = sum(Frequency), .groups="drop")

right_patterns <- g2_freq %>%
  filter(pattern %in% unique_g2_patterns) %>%
  arrange(desc(freq_g2)) %>% pull(pattern) %>%
  unique() %>% { .[seq_len(min(20, length(.)))] }

left_patterns   <- unique_g1_patterns
middle_patterns <- shared_patterns

#----Ordering----
left_order <- joined_with_pattern %>%
  filter(pattern %in% left_patterns, Panel %in% c("Progressor","PD High CTL")) %>%
  group_by(pattern) %>% summarise(freq = sum(Frequency), .groups="drop") %>%
  arrange(desc(freq)) %>% pull(pattern)

mid_order <- joined_with_pattern %>%
  filter(pattern %in% middle_patterns) %>%
  group_by(pattern) %>% summarise(freq = sum(Frequency), .groups="drop") %>%
  arrange(desc(freq)) %>% pull(pattern)

right_order <- joined_with_pattern %>%
  filter(pattern %in% right_patterns, Panel %in% c("PD Low CTL","Survivor")) %>%
  group_by(pattern) %>% summarise(freq = sum(Frequency), .groups="drop") %>%
  arrange(desc(freq)) %>% pull(pattern)

ordered_patterns <- c(left_order, mid_order, right_order)

#----Comprehensive Table----
plot_freq <- joined_with_pattern %>%
  filter(pattern %in% ordered_patterns, Panel %in% c("Progressor","PD High CTL","PD Low CTL","Survivor")) %>%
  group_by(Panel, pattern) %>%
  summarise(Frequency = sum(Frequency), .groups="drop") %>%
  complete(Panel = factor(c("Progressor","PD High CTL","PD Low CTL","Survivor"),
                          levels = c("Progressor","PD High CTL","PD Low CTL","Survivor")),
           pattern = ordered_patterns,
           fill = list(Frequency = 0)) %>%
  mutate(pattern = factor(pattern, levels = ordered_patterns))


b_left   <- length(left_order)
b_middle <- length(left_order) + length(mid_order)


y_lim <- max(plot_freq$Frequency, na.rm = TRUE)


cols <- c(
  "Progressor"   = "darkred",
  "PD High CTL"  = "#f03e3e",
  "PD Low CTL"   = "#1e3a8a",
  "Survivor"     = "#3b5bdb"
)

# Make nice pseudo-log breaks up to y_lim
make_breaks <- function(ymax){
  bases <- c(0, 1, 2, 5)
  mags  <- 10^(0:6)
  br <- sort(unique(as.vector(outer(bases, mags, `*`))))
  br[br <= ymax]
}
y_breaks <- make_breaks(y_lim)


mk_panel <- function(df, panel_name, show_x = FALSE, y_lab = NULL) {
  ggplot(dplyr::filter(df, Panel == panel_name),
         aes(x = pattern, y = Frequency)) +
    geom_col(fill = cols[[panel_name]], width = 0.9) +
    scale_y_continuous(
      trans  = scales::pseudo_log_trans(base = 2, sigma = 1),
      breaks = y_breaks,
      limits = c(0, y_lim),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(x = if (show_x) "Motif (pattern)" else NULL,
         y = if (is.null(y_lab)) "Clone frequency (pseudo-log2)" else y_lab) +
    theme_bw(base_size = 10) +
    theme(
      # no titles
      plot.title        = element_blank(),
      axis.title.x      = element_text(face = "bold", size = 12),
      axis.title.y      = element_text(face = "bold", size = 12),
      axis.text.x       = if (show_x)
        element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15, face = "bold")
      else
        element_blank(),
      axis.text.y       = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    ) +
    geom_vline(xintercept = c(b_left + 0.5, b_middle + 0.5),
               linetype = 2, linewidth = 0.3)
}

# Rebuild and save
p1 <- mk_panel(plot_freq, "Progressor",   show_x = FALSE, y_lab = "Clone frequency (pseudo-log2)")
p2 <- mk_panel(plot_freq, "PD High CTL",  show_x = FALSE)
p3 <- mk_panel(plot_freq, "PD Low CTL",   show_x = FALSE)
p4 <- mk_panel(plot_freq, "Survivor",     show_x = TRUE)

final_plot <- p1 / p2 / p3 / p4 + patchwork::plot_layout(heights = c(1,1,1,1))
out_file <- file.path(fig_dir, "GLIPH2_4layers_FreqLOG2_noTitles.png")
ggsave(out_file, final_plot, width = 22, height = 12, dpi = 600, bg = "white")
cat("Saved figure:\n", out_file, "\n")


#----Motif Heatmap----

donor_map <- joined_with_pattern %>%
  mutate(Group = case_when(
    Panel %in% c("Progressor","PD High CTL") ~ "Group1",
    Panel %in% c("PD Low CTL","Survivor")    ~ "Group2",
    TRUE ~ NA_character_
  )) %>% filter(!is.na(Group), pattern %in% ordered_patterns)

#----per-individual counts----
heat_df <- donor_map %>%
  group_by(Group, Subject, pattern) %>%
  summarise(n_tcr = sum(Frequency), .groups = "drop") %>%
  mutate(pattern = factor(pattern, levels = ordered_patterns))


all_subjects_g1 <- donor_map %>% filter(Group=="Group1") %>% distinct(Subject) %>% arrange(Subject) %>% pull(Subject)
all_subjects_g2 <- donor_map %>% filter(Group=="Group2") %>% distinct(Subject) %>% arrange(Subject) %>% pull(Subject)

heat_df_full <- bind_rows(
  heat_df %>% filter(Group=="Group1") %>%
    complete(Subject = all_subjects_g1,
             pattern  = factor(ordered_patterns, levels = ordered_patterns),
             fill = list(n_tcr = 0)) %>% mutate(Group="Group1"),
  heat_df %>% filter(Group=="Group2") %>%
    complete(Subject = all_subjects_g2,
             pattern  = factor(ordered_patterns, levels = ordered_patterns),
             fill = list(n_tcr = 0)) %>% mutate(Group="Group2")
) %>%
  mutate(pattern = factor(pattern, levels = ordered_patterns),
         Subject = forcats::fct_rev(factor(Subject)))

FILL_CAP <- 4

p_heat_g1 <- ggplot(dplyr::filter(heat_df_full, Group == "Group1"),
                    aes(x = pattern, y = Subject, fill = pmin(n_tcr, FILL_CAP))) +
  geom_tile() +
  scale_x_discrete(drop = FALSE) +
  scale_fill_gradient(
    name   = "Number of\nTCRs per donor",
    limits = c(0, FILL_CAP),
    breaks = 0:FILL_CAP,
    low    = "white",
    high   = "#3b5bdb"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid  = element_blank(),
        legend.position = "right") +
  geom_vline(xintercept = c(b_left + 0.5, b_middle + 0.5),
             linetype = 2, linewidth = 0.3)

p_heat_g2 <- ggplot(dplyr::filter(heat_df_full, Group == "Group2"),
                    aes(x = pattern, y = Subject, fill = pmin(n_tcr, FILL_CAP))) +
  geom_tile() +
  scale_x_discrete(drop = FALSE) +
  scale_fill_gradient(
    name   = "Number of\nTCRs per donor",
    limits = c(0, FILL_CAP),
    breaks = 0:FILL_CAP,
    low    = "white",
    high   = "#3b5bdb"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid  = element_blank(),
        legend.position = "right") +
  geom_vline(xintercept = c(b_left + 0.5, b_middle + 0.5),
             linetype = 2, linewidth = 0.3)


#----Assemple the Table with Heatmap----
p_heat_g2_noleg <- p_heat_g2 + guides(fill = "none")


final_with_heat <- (p_heat_g1 / p_heat_g2_noleg) / (p1 / p2 / p3 / p4) +
  plot_layout(heights = c(0.9, 0.9, 4), guides = "collect")

fig_dir <- "~/Desktop/GLIPH2_figures"; dir.create(fig_dir, FALSE, TRUE)
outfile <- file.path(fig_dir, "GLIPH2_heatmaps_on_top_of_4layers_NO_MASK_right20.png")
ggsave(outfile, final_with_heat, width = 22, height = 14, dpi = 600, bg = "white")
cat("Saved:", outfile, "\n")











