
#-----Libraries-----
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scRepertoire)
  library(data.table)   
  library(scales)
})

#----Path to Matrix----
rna_paths <- list(
  CV5_Unstim   = "Path to sample1 unstim",
  CV5_MYH6     = "Path to sample1 MYH6",
  CV5_DUF1002  = "Path to sample1 DUF1002",
  CV12_Unstim  = "Path to sample2 unstim",
  CV12_MYH6    = "Path to sample2 MYH6",
  CV12_DUF1002 = "Path to sample2 DUF1002",
  CV109_Unstim = "Path to sample3 unstim",
  CV109_MYH6   = "Path to sample3 MYH6",
  CV109_DUF1002= "Path to sample3 DUF1002",
  CV110_Unstim = "Path to sample4 unstim",
  CV110_MYH6   = "Path to sample4 MYH6",
  CV110_DUF1002= "Path to sample4 DUF1002"
)

#----Seurat Object Per Sample----
make_seu <- function(sample_id, path){
  mtx <- Read10X(data.dir = path)
  seu <- CreateSeuratObject(mtx, project = sample_id, min.cells = 3, min.features = 200)
  seu$sample_id <- sample_id
  seu$Patient   <- sub("^(CV\\d+).*", "\\1", sample_id)
  raw_cond      <- sub(".*_([Uu]nstim(?:ulated)?|MYH6|DUF(?:1002)?)$", "\\1", sample_id)
  seu$Condition <- ifelse(grepl("^unstim", raw_cond, ignore.case = TRUE), "Unstim",
                          ifelse(grepl("^myh6$", raw_cond, ignore.case = TRUE), "MYH6", "DUF1002"))
  
  seu <- RenameCells(seu, add.cell.id = sample_id)
  return(seu)
}

objs <- mapply(make_seu, names(rna_paths), rna_paths, SIMPLIFY = FALSE)

objs <- lapply(objs, function(s) {
  DefaultAssay(s) <- "RNA"
  
  s <- SCTransform(
    s,
    vst.flavor = "v2",              
    variable.features.n = 3000,
    verbose = FALSE
  )
  s
})

#----Patients for Refrence----
ref_ids <- c("CV5_Unstim", "CV12_Unstim", "CV109_Unstim", "CV110_Unstim")
ref_objs <- objs[names(objs) %in% ref_ids]

#----SCT per sample----
objs <- lapply(objs, function(s) {
  DefaultAssay(s) <- "RNA"
  SCTransform(s, vst.flavor = "v2", variable.features.n = 3000, verbose = FALSE)
})

features <- SelectIntegrationFeatures(objs, nfeatures = 3000)
objs     <- PrepSCTIntegration(objs, anchor.features = features)

#----anchors----
anchors <- FindIntegrationAnchors(object.list = objs, normalization.method = "SCT",
                                  anchor.features = features, reference = match(ref_ids, names(objs)))

 seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(seu) <- "integrated"
seu <- RunPCA(seu, npcs = 40, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.6)


#----Path to TCR Congis----
contig_paths <- list(
  CV5_Unstim   = "Path to contig sample 1 unstim.csv",
  CV5_MYH6     = "Path to contig sample 1 MYH6.csv",
  CV5_DUF1002  = "Path to contig sample 1 DUF1002.csv",
  CV12_Unstim  = "Path to contig sample 2 unstim.csv",
  CV12_MYH6    = "Path to contig sample 2 MYH6.csv",
  CV12_DUF1002 = "Path to contig sample 2 DUF1002.csv",
  CV109_Unstim = "Path to contig sample 3 unstim.csv",
  CV109_MYH6   = "Path to contig sample 3 MYH6.csv",
  CV109_DUF1002= "Path to contig sample 3 DUF1002.csv",
  CV110_Unstim = "Path to contig sample 4 unstim.csv",
  CV110_MYH6   = "Path to contig sample 4 MYH6.csv",
  CV110_DUF1002= "Path to contig sample 4 DUF1002.csv"
)

contigs_raw <- lapply(contig_paths, read.csv, stringsAsFactors = FALSE)

#----Barode Matching----
contigs_list <- Map(function(df, sid){
  if(!"barcode" %in% colnames(df)) stop(paste("Missing 'barcode' in", sid))
  df$barcode   <- paste0(sid, "_", df$barcode)
  df$sample_id <- sid
  df
}, contigs_raw, names(contigs_raw))

#----Combine TCRs----
combined.TCR <- combineTCR(contigs_list, removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
for(i in seq_along(combined.TCR)) combined.TCR[[i]]$sample_id <- names(combined.TCR)[i]

#----Map CTaa----
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(data.table); library(scRepertoire)
})

#----Sanity----
stopifnot(exists("seu"), exists("contigs_list"))
stopifnot("sample_id" %in% colnames(seu@meta.data))

sids <- sort(unique(seu$sample_id))

#----Match TCR Barcode----
seurat_prefix_map <- setNames(paste0(sids, "_"), sids)

contigs_rebuilt <- lapply(sids, function(sid){
  df <- fix_contig_cols(contigs_list[[sid]])
  
  
  core <- sub("^.*_", "", df$barcode)  
 
  df$barcode   <- paste0(seurat_prefix_map[sid], core) 
  df$sample_id <- sid
  df
})
names(contigs_rebuilt) <- sids

combined.TCR <- scRepertoire::combineTCR(
  contigs_rebuilt, removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE
)
names(combined.TCR) <- names(contigs_rebuilt)
for (i in seq_along(combined.TCR)) combined.TCR[[i]]$sample_id <- names(combined.TCR)[i]

suppressPackageStartupMessages({ library(dplyr); library(data.table) })

stopifnot(exists("seu"), exists("combined.TCR"))

#----Sanity----
expected_sids <- sort(unique(seu$sample_id))
present_sids  <- sort(names(combined.TCR))

cat("\n# Expected sample_ids from Seurat (n=", length(expected_sids), "):\n", sep="")
print(expected_sids)
cat("\n# Present sample_ids in combined.TCR (n=", length(present_sids), "):\n", sep="")
print(present_sids)

missing_in_TCR <- setdiff(expected_sids, present_sids)
extra_in_TCR   <- setdiff(present_sids, expected_sids)

cat("\n# Missing in TCR:", if(length(missing_in_TCR)) paste(missing_in_TCR, collapse=", ") else "NONE", "\n")
cat("# Extra in TCR:  ", if(length(extra_in_TCR))   paste(extra_in_TCR,   collapse=", ") else "NONE", "\n")

#----Summary----
tcr_summary <- lapply(present_sids, function(sid){
  df <- combined.TCR[[sid]]
  data.frame(
    sample_id     = sid,
    tcr_rows      = nrow(df),
    tcr_unique_bc = length(unique(df$barcode)),
    tcr_CTaa_nonNA= sum(!is.na(df$CTaa)),
    stringsAsFactors = FALSE
  )
}) |> bind_rows()

#----Exact Matching----
seu_counts <- as.data.frame(table(seu$sample_id)) |>
  dplyr::rename(sample_id = Var1, seu_cells = Freq)

match_cells <- sapply(present_sids, function(sid){
  sum(combined.TCR[[sid]]$barcode %in% colnames(seu)[seu$sample_id == sid])
})

coverage <- data.frame(sample_id = present_sids, match_cells = as.integer(match_cells))
summary_tbl <- tcr_summary |>
  left_join(seu_counts, by="sample_id") |>
  left_join(coverage,  by="sample_id") |>
  mutate(
    pct_cells_with_TCR = round(100 * match_cells / seu_cells, 1),
    pct_CTaa_in_TCR    = ifelse(tcr_rows>0, round(100 * tcr_CTaa_nonNA / tcr_rows, 1), NA)
  ) |>
  arrange(sample_id)

cat("\n# Per-sample mapping summary:\n")
print(summary_tbl)

cat("\n# Sanity verdict:\n")
if (length(missing_in_TCR)==0 && length(extra_in_TCR)==0 &&
    all(summary_tbl$match_cells > 0) && all(summary_tbl$seu_cells > 0)) {
  cat("All expected samples are present, barcodes align, and every sample has mapped TCRs.\n")
} else {
  cat("wrong.\n")
}


#----Expansion evaluation----
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(data.table)
  library(ggplot2); library(scales); library(patchwork)
  library(Seurat)
})

#----Define expansion----
present_all3   <- per_cond_wide %>% filter(Unstim > 0, MYH6 > 0, DUF1002 > 0)
expanded_both  <- present_all3 %>% filter(MYH6 > Unstim, DUF1002 > Unstim) %>%
  mutate(TotalStim = MYH6 + DUF1002) %>% arrange(desc(TotalStim))
cat("# expanded-in-both clones:", nrow(expanded_both), "\n")

TOP_N <- 30
expanded_clones <- head(expanded_both$Clone, TOP_N)
assign("expanded_clones", expanded_clones, envir = .GlobalEnv)


plot_df <- per_cond %>% filter(Clone %in% expanded_clones)
plot_df$Condition <- factor(plot_df$Condition, levels = c("Unstim","MYH6","DUF1002"))

#----Reordering the legend----
clone_levels <- per_cond %>% filter(Condition=="DUF1002", Clone %in% expanded_clones) %>%
  arrange(desc(Count)) %>% pull(Clone)
plot_df$Clone <- factor(plot_df$Clone, levels = clone_levels)
clone_colors  <- setNames(hue_pal()(length(clone_levels)), clone_levels)

#----Extract TRA/TRB----
extract_TRA_TRB <- function(ctaa, ctgene = NULL) {
  n <- length(ctaa)
  if (is.null(ctgene)) ctgene <- rep(NA_character_, n)
  TRA <- rep(NA_character_, n)
  TRB <- rep(NA_character_, n)
  
  parts_list <- strsplit(ifelse(is.na(ctaa), "", ctaa), "_", fixed = TRUE)
  genes_list <- strsplit(ifelse(is.na(ctgene), "", ctgene), "_", fixed = TRUE)
  
  for (i in seq_len(n)) {
    p <- parts_list[[i]]
    g <- genes_list[[i]]
    if (length(p) == 0) next
    
    if (length(g) == length(p) && all(g %in% c("TRA","TRB"))) {
      
      for (j in seq_along(g)) {
        if (g[j] == "TRA") TRA[i] <- p[j]
        if (g[j] == "TRB") TRB[i] <- p[j]
      }
    } else {
      
      if (length(p) == 2) {
        if (grepl("^CASS", p[1])) { TRB[i] <- p[1]; TRA[i] <- p[2] }
        else if (grepl("^CASS", p[2])) { TRB[i] <- p[2]; TRA[i] <- p[1] }
        else { TRA[i] <- p[1]; TRB[i] <- p[2] }
      } else if (length(p) == 1) {
        if (!is.na(ctgene[i]) && grepl("TRB", ctgene[i])) TRB[i] <- p[1] else TRA[i] <- p[1]
      }
    }
  }
  data.frame(TRA = TRA, TRB = TRB, stringsAsFactors = FALSE)
}

#----Build per-row TRA/TRB from combined.TCR----
DT_raw <- rbindlist(lapply(names(combined.TCR), function(sid){
  df <- as.data.frame(combined.TCR[[sid]])
  n  <- nrow(df)
  ctgene <- if ("CTgene" %in% names(df)) df$CTgene else rep(NA_character_, n)
  ctaa   <- if ("CTaa"   %in% names(df)) df$CTaa   else rep(NA_character_, n)
  ch <- extract_TRA_TRB(ctaa, ctgene)
  data.table(
    sample_id = rep(sid, n),
    patient   = rep(sub("_.*$", "", sid), n),
    TRA = ch$TRA,
    TRB = ch$TRB
  )
}), fill = TRUE)

cat("# rows:", nrow(DT_raw),
    " TRA NAs:", sum(is.na(DT_raw$TRA)),
    " TRB NAs:", sum(is.na(DT_raw$TRB)), "\n")


cons_TRB_by_TRA <- DT_raw[!is.na(TRA) & !is.na(TRB),
                          .N, by = .(patient, TRA, TRB)][order(patient, TRA, -N)
                          ][, .SD[1], by = .(patient, TRA)
                          ][, .(patient, TRA, TRB_cons = TRB)]

cons_TRA_by_TRB <- DT_raw[!is.na(TRA) & !is.na(TRB),
                          .N, by = .(patient, TRB, TRA)][order(patient, TRB, -N)
                          ][, .SD[1], by = .(patient, TRB)
                          ][, .(patient, TRB, TRA_cons = TRA)]

DT_imp <- DT_raw |>
  left_join(cons_TRB_by_TRA, by = c("patient","TRA")) |>
  left_join(cons_TRA_by_TRB, by = c("patient","TRB")) |>
  mutate(
    TRB = ifelse(is.na(TRB) & !is.na(TRB_cons), TRB_cons, TRB),
    TRA = ifelse(is.na(TRA) & !is.na(TRA_cons), TRA_cons, TRA)
  ) |>
  select(-TRB_cons, -TRA_cons)

cat("# still NA after fill — TRA:", sum(is.na(DT_imp$TRA)),
    " TRB:", sum(is.na(DT_imp$TRB)), "\n")

#----Normalized clonotype label----
DT_imp[, Clone := paste0(ifelse(is.na(TRA) | TRA=="","NA",TRA), "_",
                         ifelse(is.na(TRB) | TRB=="","NA",TRB))]

DT_counts <- DT_imp[, .N, by = .(Clone, sample_id)]
setnames(DT_counts, "N", "Count")

sample_to_condition <- unique(DT_counts[, .(sample_id)])[
  , Condition := sub(".*_(Unstim|MYH6|DUF1002)$", "\\1", sample_id)]

per_cond <- merge(DT_counts, sample_to_condition, by = "sample_id")[
  , .(Count = sum(Count)), by = .(Clone, Condition)]
per_cond_wide <- dcast(as.data.table(per_cond), Clone ~ Condition, value.var = "Count", fill = 0)

#----Expanded Clones Shared in BOTH MYH6 & DUF1002 vs Unstim----
expanded_clones <- per_cond_wide %>%
  filter(Unstim > 0, MYH6 > 0, DUF1002 > 0) %>%
  filter(MYH6 > Unstim, DUF1002 > Unstim) %>%
  arrange(desc(DUF1002 + MYH6)) %>%
  pull(Clone)

cat("# expanded clones:", length(expanded_clones), "\n")
if (length(expanded_clones) == 0) {
  print(head(per_cond_wide))   
  stop("No expanded clones detected after proper parsing.")
}

#----Plotting----
df_long <- per_cond %>%
  filter(Clone %in% expanded_clones) %>%
  mutate(Condition = factor(Condition, levels = c("Unstim","MYH6","DUF1002")))

clone_levels <- per_cond_wide %>%
  filter(Clone %in% expanded_clones) %>%
  arrange(desc(DUF1002)) %>%
  pull(Clone)

df_long$Clone <- factor(df_long$Clone, levels = clone_levels)

p_alluvial <- ggplot(
  df_long,
  aes(x = Condition, stratum = Clone, alluvium = Clone, y = Count, fill = Clone)
) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", color = "black", alpha = 0.25) +
  geom_stratum(color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = 16, color = "black", face = "bold"),
    axis.title.y = element_text(size = 16, color = "black", face = "bold"),
    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black")
  ) +
  xlab("Condition") +
  ylab("Clone Count") +
  guides(fill = guide_legend(title = "Clones"))

p_alluvial
ggsave("Path to Save/TCR_alluvial_plot.png", p_alluvial, width = 18, height = 8, dpi = 600)
ggsave("/Path to Save/TCR_alluvial_plot.pdf", p_alluvial, width = 18, height = 8, dpi = 600)


#--------UMAP Projection for Expanded Cells--------
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(data.table); library(Seurat); library(ggplot2)
})

stopifnot(exists("seu"), exists("combined.TCR"), exists("expanded_clones"))
if (!"Condition" %in% colnames(seu@meta.data)) {
  seu$Condition <- sub(".*_(Unstim|MYH6|DUF1002)$", "\\1", seu$sample_id)
}
seu$Condition <- factor(seu$Condition, levels = c("Unstim","MYH6","DUF1002"))


extract_core <- function(x) {
  
  last <- sub("^.*_", "", x)
  
  last <- sub("-1$", "", last)
  
  last <- toupper(last)
  last <- sub("^([^ACGT]*)([ACGT]+).*?$", "\\2", last, perl = TRUE)
  last
}

#----Seurat sanity table----
seu_cells <- colnames(seu)
seu_df <- data.frame(
  seucell   = seu_cells,
  sample_id = seu$sample_id,
  Condition = seu$Condition,
  stringsAsFactors = FALSE
)
seu_df$core <- extract_core(seu_df$seucell)

cat("\n# CHECK A1 — Seurat colname examples:\n")
print(head(seu_df, 5))


seu_keys <- seu_df %>%
  select(sample_id, core, seucell)
dup_keys <- seu_keys %>% count(sample_id, core) %>% filter(n > 1)
if (nrow(dup_keys) > 0) {
  warning("Some (sample_id, core) map to multiple Seurat colnames (unexpected). Showing first few:")
  print(head(dup_keys, 5))
}


#----TCR side-----
tcr_list <- lapply(names(combined.TCR), function(sid){
  df <- as.data.frame(combined.TCR[[sid]])
  if (!("barcode" %in% names(df))) stop(paste("combined.TCR element missing 'barcode' for", sid))
  data.frame(
    sample_id = sid,
    barcode   = df$barcode,
    stringsAsFactors = FALSE
  )
})
tcr_df <- bind_rows(tcr_list)
tcr_df$core <- extract_core(tcr_df$barcode)

cat("\n# CHECK B1 — TCR barcode examples:\n")
print(head(tcr_df, 5))

#----Overlap by sample_id----
overlap_tbl <- tcr_df %>%
  mutate(key = paste(sample_id, core)) %>%
  distinct(sample_id, core) %>%
  mutate(in_seurat = paste(sample_id, core) %in% paste(seu_keys$sample_id, seu_keys$core)) %>%
  group_by(sample_id) %>%
  summarise(
    n_tcr_keys = n(),
    n_match    = sum(in_seurat),
    overlap_pct = round(100 * n_match / pmax(n_tcr_keys, 1), 1),
    .groups = "drop"
  ) %>% arrange(sample_id)

cat("\n# CHECK B2 — Overlap by (sample_id, core) with Seurat:\n")
print(overlap_tbl)



#----Expansion bins----
counts <- tcr_cells %>% count(CloneNorm, Condition, name = "N")
cut_size <- function(n) {
  if (is.na(n)) return(NA_character_)
  if (n == 1)        "singleton"
  else if (n <= 5)   "small"
  else if (n <= 20)  "medium"
  else               "large"
}
counts$ExpCat <- vapply(counts$N, cut_size, character(1))

cells27 <- tcr_cells %>%
  inner_join(counts, by = c("CloneNorm","Condition")) %>%
  rename(cell = seucell) %>%
  select(cell, Condition, ExpCat)

cat("# Per-condition 27-clone cell counts by expansion bin:\n")
print(cells27 %>% count(Condition, ExpCat) %>% arrange(Condition, ExpCat))

#----UMAP frames----
emb <- Embeddings(seu, "umap"); colnames(emb) <- c("UMAP_1","UMAP_2")
md  <- seu@meta.data[, c("sample_id","Condition")]
df_bg <- data.frame(cell = colnames(seu), emb, md, row.names = NULL, check.names = FALSE)

#----CD4 background----
get_cd4_cells <- function(seu){
  m <- seu@meta.data
  ann_cols <- intersect(c("CellType","celltype","cell_type","major_celltype","fine_celltype","annot","annotation"), colnames(m))
  if (length(ann_cols)) {
    for (c in ann_cols) {
      if (any(grepl("CD4", m[[c]], ignore.case = TRUE)))
        return(rownames(m)[grepl("CD4", m[[c]], ignore.case = TRUE)])
    }
  }
  assay_to_use <- if ("SCT" %in% names(seu@assays)) "SCT" else "RNA"
  DefaultAssay(seu) <- assay_to_use
  vars <- intersect(c("CD3D","CD4","CD8A"), rownames(seu[[assay_to_use]]))
  if (length(vars) >= 2) {
    df <- FetchData(seu, vars = vars, assay = assay_to_use, layer = "data")
    rownames(df) <- colnames(seu)
    keep <- rep(TRUE, nrow(df))
    if ("CD3D" %in% vars) keep <- keep & (df$CD3D > 0)
    if (all(c("CD4","CD8A") %in% vars)) keep <- keep & (df$CD4 >= df$CD8A)
    return(rownames(df)[keep])
  }
  rownames(seu@meta.data)  
}
cd4_cells <- get_cd4_cells(seu)
df_bg$isCD4 <- df_bg$cell %in% cd4_cells

#----All 27 expanded clones----
df_fg <- cells27 %>%
  inner_join(df_bg[, c("cell","UMAP_1","UMAP_2","Condition")], by = c("cell","Condition"))

bg_col   <- "#d9d9d9"
exp_cols <- c(
  singleton = "lightblue",  
  small     = "orange",  
  medium    = "red",  
  large     = "brown"   
)

plot_one <- function(cond) {
  sub_bg <- df_bg %>% filter(Condition == cond, isCD4)
  sub_fg <- df_fg %>% filter(Condition == cond)
  
  cat(sprintf("[%s] bg=%d CD4 cells, fg=%d 27-clone cells\n", cond, nrow(sub_bg), nrow(sub_fg)))
  
  ggplot() +
    geom_point(data = sub_bg, aes(UMAP_1, UMAP_2),
               color = bg_col, fill = bg_col, size = 0.6, alpha = 0.6, shape = 16) +
    geom_point(data = sub_fg, aes(UMAP_1, UMAP_2, fill = ExpCat),
               color = "black", size = 1.9, stroke = 0.28, shape = 21, alpha = 0.95) +
    scale_fill_manual(values = exp_cols, drop = FALSE) +
    coord_equal() +
    labs(title = paste0(cond, " — CD4 T cells; 27 expanded clones highlighted"),
         x = "UMAP 1", y = "UMAP 2", fill = "Expansion") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

p_unstim <- plot_one("Unstim")
p_myh6   <- plot_one("MYH6")
p_duf    <- plot_one("DUF1002")

p_unstim; p_myh6; p_duf
ggsave("UMAP_CD4_27clones_Unstim.png",  p_unstim, width=5.5, height=5, dpi=600)
ggsave("UMAP_CD4_27clones_MYH6.png",    p_myh6,   width=5.5, height=5, dpi=600)
ggsave("UMAP_CD4_27clones_DUF1002.png", p_duf,    width=5.5, height=5, dpi=600)
ggsave("UMAP_CD4_27clones_Unstim.pdf",  p_unstim, width=5.5, height=5, dpi=600)
ggsave("UMAP_CD4_27clones_MYH6.pdf",    p_myh6,   width=5.5, height=5, dpi=600)
ggsave("UMAP_CD4_27clones_DUF1002.pdf", p_duf,    width=5.5, height=5, dpi=600)



#--------Define which Clusters Contain Expanded Cells and What TCR Gene enriched there--------

#----Build Seurat keys----
seu_df <- data.frame(
  seucell   = colnames(seu),
  sample_id = seu$sample_id,
  stringsAsFactors = FALSE, check.names = FALSE
) %>%
  mutate(core = sub("^.*_", "", seucell),    
         core = sub("-1$", "", core))

seu_keys <- seu_df %>% select(sample_id, core, seucell)

#---- 27-clone cells in Seurat-----
stopifnot(exists("cells27"))

DT27_cells <- unique(cells27$cell)
cat("# matched 27-clone cells in Seurat:", length(DT27_cells), "\n")


#----Which clusters contain ≥1 of those cells-----
cl_with27 <- seu$seurat_clusters[colnames(seu) %in% DT27_cells] %>%
  as.character() %>% unique()
cat("# clusters with ≥1 expanded-clone cell:", ifelse(length(cl_with27)>20, paste(sort(cl_with27), collapse=", "), "NONE"), "\n")

cat("# matched 27-clone cells in Seurat:", length(DT27_cells), "\n")
cat("# clusters with ≥1 expanded-clone cell: ",
    if (length(cl_with27) > 0) paste(sort(cl_with27), collapse = ", ") else "NONE",
    "\n")

#----Prep points & centers----
emb <- Embeddings(seu, "umap")
plot_pts <- data.frame(
  cell    = colnames(seu),
  UMAP_1  = emb[,1],
  UMAP_2  = emb[,2],
  cluster = as.character(seu$seurat_clusters),
  stringsAsFactors = FALSE, check.names = FALSE
)

plot_pts$highlight <- plot_pts$cluster %in% cl_with27

centers <- plot_pts %>%
  filter(highlight) %>%
  group_by(cluster) %>%
  summarize(cx = median(UMAP_1),
            cy = median(UMAP_2),
            .groups = "drop")

#----Colors & theme----
bg_col <- "#d9d9d9"

if (length(cl_with27) > 0) {
  cl_hl <- sort(unique(centers$cluster))
  cluster_colors <- setNames(hue_pal()(length(cl_hl)), cl_hl)
} else {
  cluster_colors <- character(0)
}

sexy_white <- theme_classic(base_size = 14) +
  theme(
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    plot.title  = element_text(face = "bold")
  )

#----Plotting----
p <- ggplot() +
  
  geom_point(data = plot_pts,
             aes(UMAP_1, UMAP_2),
             color = bg_col, alpha = 0.35, size = 0.45, shape = 16)

if (length(cl_with27) > 0) {
  p <- p +
    geom_point(
      data = subset(plot_pts, highlight),
      aes(UMAP_1, UMAP_2, fill = cluster),
      shape = 21, color = "black", stroke = 0.20, alpha = 0.95, size = 1.2
    ) +
    scale_fill_manual(
      values = cluster_colors,
      guide = guide_legend(title = "Clusters w/ 27-clone cells")
    )
}

p <- p +
  coord_equal() +
  labs(
    title = "Merged UMAP (RNA)\nClusters containing ≥1 of the 27 expanded clones are highlighted",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  sexy_white

p


#----cluster subsets with more than 5 cells------
target_clusters <- c("19","4","9","7","16","5","6","26","27","23","12")

plot_pts$cluster   <- as.character(plot_pts$cluster)
plot_pts$highlight <- plot_pts$cluster %in% target_clusters

#----Colors for the selected clusters----
if (!exists("cluster_colors") || !length(intersect(names(cluster_colors), target_clusters))) {
  cluster_colors <- setNames(hue_pal()(length(target_clusters)), target_clusters)
} else {
  cluster_colors <- cluster_colors[names(cluster_colors) %in% target_clusters]
  missing_cols <- setdiff(target_clusters, names(cluster_colors))
  if (length(missing_cols)) {
    cluster_colors <- c(cluster_colors,
                        setNames(hue_pal()(length(missing_cols)), missing_cols))
  }
}

bg_col <- "#d9d9d9"

p <- ggplot() +
  
  geom_point(
    data = plot_pts, aes(UMAP_1, UMAP_2),
    color = bg_col, alpha = 0.35, size = 0.45, shape = 16
  ) +
  
  geom_point(
    data = subset(plot_pts, highlight),
    aes(UMAP_1, UMAP_2, fill = cluster),
    shape = 21, color = "black", stroke = 0.22, alpha = 0.95, size = 1.25
  ) +
  scale_fill_manual(values = cluster_colors,
                    guide = guide_legend(title = "Selected clusters")) +
  coord_equal(clip = "off") +
  labs(
    title = "Merged UMAP (RNA)\nSelected clusters highlighted",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 40, 10, 20)
  )

p
ggsave("UMAP_selected_clusters_no_annotations.png", p, width = 10, height = 7, dpi = 600)
ggsave("UMAP_selected_clusters_no_annotations.pdf", p, width = 10, height = 7, dpi = 600)



#--------Defining CD4 CTL Clusters----

#----Config of cytotoxic genes----
cyto_genes  <- c("GZMB","PRF1","NKG7")
assay_to_use <- if ("SCT" %in% names(seu@assays)) "SCT" else "RNA"
expr_layer   <- "data"   

#----CD4 selection----
get_cd4_cells <- function(seu_obj, assay = assay_to_use, layer = expr_layer){
  md <- seu_obj@meta.data
  ann_cols <- intersect(c("CellType","celltype","cell_type","major_celltype",
                          "fine_celltype","annot","annotation"), colnames(md))
  if (length(ann_cols)) {
    for (c in ann_cols) {
      hit <- rownames(md)[grepl("CD4", md[[c]], ignore.case = TRUE)]
      if (length(hit)) return(hit)
    }
  }
  
  needed <- intersect(c("CD3D","CD4","CD8A"), rownames(seu_obj[[assay]]))
  if (length(needed) >= 2) {
    df <- FetchData(seu_obj, vars = needed, assay = assay, layer = layer)
    rownames(df) <- colnames(seu_obj)
    keep <- rep(TRUE, nrow(df))
    if ("CD3D" %in% colnames(df)) keep <- keep & (df$CD3D > 0)
    if (all(c("CD4","CD8A") %in% colnames(df))) keep <- keep & (df$CD4 >= df$CD8A)
    return(rownames(df)[keep])
  }
  warning("CD4 fallback not possible — returning all cells.")
  colnames(seu_obj)
}

cd4_cells <- get_cd4_cells(seu, assay = assay_to_use, layer = expr_layer)

#----Computing cytotoxic enrichment per CD4 cluster----
genes_present <- cyto_genes[cyto_genes %in% rownames(seu[[assay_to_use]])]
stopifnot(length(genes_present) > 0)

dat <- FetchData(
  seu, vars = c("seurat_clusters", genes_present),
  cells = cd4_cells, assay = assay_to_use, layer = expr_layer
)
colnames(dat)[1] <- "cluster"
dat$cluster <- as.character(dat$cluster)


frac_df <- dat %>%
  mutate(across(all_of(genes_present), ~ .x > 0, .names = "det_{col}")) %>%
  group_by(cluster) %>%
  summarise(across(starts_with("det_"), ~ mean(.x, na.rm = TRUE)*100, .names = "{.col}"),
            .groups = "drop")

#----average expression per cluster----
avg_df <- dat %>%
  group_by(cluster) %>%
  summarise(across(all_of(genes_present), ~ mean(.x, na.rm = TRUE)), .groups = "drop")


avg_long <- avg_df %>%
  pivot_longer(-cluster, names_to = "gene", values_to = "avg_expr") %>%
  group_by(gene) %>%
  mutate(z = (avg_expr - mean(avg_expr, na.rm=TRUE)) / (sd(avg_expr, na.rm=TRUE)+1e-9)) %>%
  ungroup()

cyto_score <- avg_long %>% group_by(cluster) %>%
  summarise(cytotoxic_score = mean(z, na.rm = TRUE), .groups = "drop")


thr_score <- quantile(cyto_score$cytotoxic_score, 0.75, na.rm = TRUE)
thr_frac  <- 50  

frac_cols <- grep("^det_", names(frac_df), value = TRUE)
summary_tbl <- cyto_score %>%
  left_join(frac_df, by = "cluster") %>%
  mutate(any_marker_high = apply(select(., all_of(frac_cols)), 1, function(x) any(x >= thr_frac, na.rm=TRUE)),
         enriched_cyto_cd4 = cytotoxic_score >= thr_score & any_marker_high)

enriched_clusters <- summary_tbl %>% filter(enriched_cyto_cd4) %>% pull(cluster)
message("Enriched (navy) clusters: ", paste(enriched_clusters, collapse = ", "))

#----Build UMAP frame----
emb <- Embeddings(seu, "umap")
plot_df <- data.frame(
  cell    = cd4_cells,
  UMAP_1  = emb[cd4_cells, 1],
  UMAP_2  = emb[cd4_cells, 2],
  cluster = as.character(seu$seurat_clusters[cd4_cells]),
  stringsAsFactors = FALSE, check.names = FALSE
)

plot_df$highlight <- plot_df$cluster %in% enriched_clusters

navy_col <- "#08306b"
bg_col   <- "#d9d9d9"


centers <- plot_df %>%
  group_by(cluster) %>%
  summarise(cx = median(UMAP_1), cy = median(UMAP_2), .groups = "drop") %>%
  filter(cluster %in% enriched_clusters)


p_umap_cyto <- ggplot() +
  
  geom_point(data = subset(plot_df, !highlight),
             aes(UMAP_1, UMAP_2),
             color = bg_col, alpha = 0.35, size = 0.45, shape = 16) +
  
  geom_point(data = subset(plot_df, highlight),
             aes(UMAP_1, UMAP_2),
             color = "black", fill = navy_col, size = 1.4, alpha = 0.95, shape = 21, stroke = 0.25) +
  coord_equal() +
  labs(title = "Merged UMAP\nClusters enriched for GZMB/PRF1/NKG7",
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "bold"))

p_umap_cyto
ggsave("Path to Save/UMAP_CD4_cytotoxic_clusters_navy.png", p_umap_cyto, width=8, height=7, dpi=400)
ggsave("Path to Save/UMAP_CD4_cytotoxic_clusters_navy.pdf", p_umap_cyto, width=8, height=7, dpi=400)



#----DE within 27-clone cells----
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(tibble); library(ggplot2) })

#----Subset to the 27-clone cells----
stopifnot(exists("cell_map"))
cells_27 <- unique(cell_map$cell)
seu_27   <- subset(seu, cells = cells_27)

if (!"Condition" %in% colnames(seu_27@meta.data)) {
  seu_27$Condition <- sub(".*_(Unstim|MYH6|DUF1002)$", "\\1", seu_27$sample_id)
}
seu_27$Condition <- factor(seu_27$Condition, levels = c("Unstim","MYH6","DUF1002"))
Idents(seu_27)   <- seu_27$Condition

#----RNA assay----
DefaultAssay(seu_27) <- "RNA"


seu_27 <- tryCatch({
  JoinLayers(seu_27, assay = "RNA")
}, error = function(e) {
  message("JoinLayers(RNA) failed: ", e$message)
  stop("Could not join layers for RNA assay; aborting DE.")
})


rna_layers <- tryCatch(Layers(seu_27[["RNA"]]), error = function(e) character(0))
if (!"data" %in% rna_layers) {
  message("Normalizing RNA (LogNormalize) to create 'data' layer")
  seu_27 <- NormalizeData(seu_27, assay = "RNA",
                          normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = FALSE)
}

#----Wilcoxon test comparing----
run_markers_safe <- function(obj, ident1, ident2) {
  FindMarkers(
    obj, ident.1 = ident1, ident.2 = ident2,
    group.by = "Condition",
    test.use = "wilcox",
    min.pct = 0.05,
    logfc.threshold = 0.10,
    only.pos = TRUE
  )
}

markers_MYH6 <- tryCatch(run_markers_safe(seu_27, "MYH6", "Unstim"),
                         error = function(e){ message("FindMarkers MYH6 error: ", e$message); NULL })
markers_DUF  <- tryCatch(run_markers_safe(seu_27, "DUF1002", "Unstim"),
                         error = function(e){ message("FindMarkers DUF error: ", e$message); NULL })

stopifnot(!is.null(markers_MYH6), !is.null(markers_DUF))

#----Shared upregulated genes----
markers_MYH6 <- markers_MYH6 %>% rownames_to_column("gene") %>% mutate(comp = "MYH6_vs_Unstim")
markers_DUF  <- markers_DUF  %>% rownames_to_column("gene") %>% mutate(comp = "DUF_vs_Unstim")

shared_up <- inner_join(
  markers_MYH6 %>% filter(p_val_adj < 0.05, avg_log2FC > 0),
  markers_DUF  %>% filter(p_val_adj < 0.05, avg_log2FC > 0),
  by = "gene",
  suffix = c("_MYH6","_DUF")
) %>%
  mutate(mean_log2FC = (avg_log2FC_MYH6 + avg_log2FC_DUF)/2) %>%
  arrange(desc(mean_log2FC))

cat("\nTop shared UP genes in 27-clone cells (MYH6 & DUF1002 > Unstim):\n")
print(shared_up %>% select(gene, mean_log2FC, p_val_adj_MYH6, p_val_adj_DUF) %>% head(100))

#----Preview bar of top hits----
if (nrow(shared_up)) {
  top_n_show <- min(30, nrow(shared_up))
  ggplot(shared_up %>% head(top_n_show),
         aes(x = reorder(gene, mean_log2FC), y = mean_log2FC)) +
    geom_col(alpha = 0.9) +
    coord_flip() +
    theme_classic(base_size = 13) +
    labs(x = NULL, y = "Mean log2FC (MYH6 & DUF > Unstim)",
         title = "Shared upregulated genes in 27-clone cells (RNA, joined layers)")
}


#----One Plot Per Gene----

#----10 TCR-direct genes-----
tcr_direct_genes <- c("JUN","DUSP4","DUSP16","TNFAIP3","RASGEF1B",
                      "BTG1","PDCD4","LDHA","HIF1A","STK17A","IL7R","NFKBIA")

stopifnot(exists("cell_map"), exists("expanded_clones"))
cells_use <- unique(cell_map$cell)


pick_expr_source <- function(seu) {
  
  if ("SCT" %in% names(seu@assays)) {
    lyr <- tryCatch(Layers(seu[["SCT"]]), error = function(e) character(0))
    if ("data" %in% lyr) return(list(assay="SCT", layer="data", log1p_needed=FALSE))
  }
  
  if ("RNA" %in% names(seu@assays)) {
    lyr <- tryCatch(Layers(seu[["RNA"]]), error = function(e) character(0))
    if ("data" %in% lyr)   return(list(assay="RNA", layer="data",   log1p_needed=FALSE))
    if ("counts" %in% lyr) return(list(assay="RNA", layer="counts", log1p_needed=TRUE))
  }
  
  fa <- names(seu@assays)[1]
  lyr <- tryCatch(Layers(seu[[fa]]), error = function(e) character(0))
  if (length(lyr) == 0) stop("No layers found in assay: ", fa)
  list(assay=fa, layer=lyr[1], log1p_needed=TRUE)
}

src <- pick_expr_source(seu)
assay_use <- src$assay
layer_use <- src$layer
message(sprintf("Using %s/%s (log1p_needed=%s) for expression.",
                assay_use, layer_use, src$log1p_needed))

#----cells x genes matrix----
get_cells_by_genes <- function(seu, assay, layer, genes, cells, log1p_needed=FALSE) {
  DefaultAssay(seu) <- assay
  if (!layer %in% Layers(seu[[assay]])) {
    stop("Requested layer '", layer, "' not present in assay '", assay, "'. Available: ",
         paste(Layers(seu[[assay]]), collapse=", "))
  }
  M <- GetAssayData(seu, assay = assay, layer = layer)  
  genes_present <- intersect(genes, rownames(M))
  cells_present <- intersect(cells, colnames(M))
  if (!length(genes_present)) stop("None of the requested genes are present in ", assay, "/", layer)
  if (!length(cells_present)) stop("None of the requested cells are present in ", assay, "/", layer)
  
  sub <- M[genes_present, cells_present, drop = FALSE]  
  mat <- t(as.matrix(sub))                               
  
  
  out <- matrix(NA_real_, nrow = length(cells), ncol = length(genes),
                dimnames = list(cells, genes))
  out[match(cells_present, cells), match(genes_present, genes)] <- mat
  
  if (isTRUE(log1p_needed)) out <- log1p(out)
  out
}

genes_use <- intersect(tcr_direct_genes, rownames(seu[[assay_use]]))
if (!length(genes_use)) stop("None of the requested genes found in assay: ", assay_use)

mat_cxg <- get_cells_by_genes(
  seu, assay = assay_use, layer = layer_use,
  genes = genes_use, cells = cells_use, log1p_needed = src$log1p_needed
)


#----Build clone-level----
expr_df <- cell_map %>%
  select(cell, CloneNorm, Condition) %>%
  mutate(Condition = factor(Condition, levels = c("Unstim","MYH6","DUF1002"))) %>%
  
  left_join(
    as_tibble(mat_cxg, rownames = "cell") %>%
      pivot_longer(-cell, names_to = "gene", values_to = "log1p_exp"),
    by = "cell"
  ) %>%
  filter(!is.na(log1p_exp), gene %in% genes_use) %>%
  group_by(gene, CloneNorm, Condition) %>%
  summarise(mean_log1p = mean(log1p_exp, na.rm = TRUE), .groups = "drop")

#----BarPlot summary----
summarize_bars <- function(df_gene) {
  df_gene %>%
    group_by(Condition) %>%
    summarise(
      mu  = mean(mean_log1p, na.rm = TRUE),
      sem = sd(mean_log1p,  na.rm = TRUE) / sqrt(sum(is.finite(mean_log1p))),
      .groups = "drop"
    )
}


dot_colors <- c("Unstim"="#6e6e6e", "MYH6"="#8b0000", "DUF1002"="#001f5b")
bar_fills  <- c("Unstim"="#e4e4e4", "MYH6"="#f2c8c8", "DUF1002"="#c9d3f6")

#----Single-gene plot function----
plot_gene_single <- function(gene, out_dir = NULL,
                             width = 3.8, height = 3.4, dpi = 600) {
  df_g  <- expr_df %>% filter(gene == !!gene)
  if (!nrow(df_g)) stop("No data for gene: ", gene)
  sum_g <- summarize_bars(df_g)
  
  p <- ggplot(sum_g, aes(x = Condition, y = mu, fill = Condition)) +
    geom_col(width = 0.58, alpha = 0.95, color = NA) +
    geom_errorbar(aes(ymin = mu - sem, ymax = mu + sem),
                  width = 0.16, linewidth = 0.55) +
    geom_point(
      data = df_g,
      aes(x = Condition, y = mean_log1p, color = Condition),
      position = position_jitter(width = 0.10, height = 0),
      size = 2.2, alpha = 0.9, stroke = 0
    ) +
    scale_fill_manual(values = bar_fills) +
    scale_color_manual(values = dot_colors) +
    labs(title = gene, y = "Mean log1p(count) per clone", x = NULL) +
    theme_classic(base_size = 13) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.x  = element_text(angle = 18, hjust = 1),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(out_dir, paste0("BAR_", gene, ".png")), p,
           width = width, height = height, dpi = dpi)
    ggsave(file.path(out_dir, paste0("BAR_", gene, ".pdf")), p,
           width = width, height = height, dpi = dpi)
  }
  p
}

#----Plots with paired Wilcoxon----
safe_paired_wilcox <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  n <- length(x)
  if (n < 2) return(NA_real_)
  d <- y - x
  if (all(d == 0)) return(1) 
  suppressWarnings(
    tryCatch(
      wilcox.test(y, x, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
  )
}

#----paired rank-biserial effect size----
rank_biserial <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  n <- length(x); if (n < 2) return(NA_real_)
  d <- y - x
  if (all(d == 0)) return(0)
  r <- rank(abs(d))
  Wpos <- sum(r[d > 0]); Wneg <- sum(r[d < 0])
  (Wpos - Wneg) / (n * (n + 1) / 2)
}

#----P value instead of Stars helper----
p_to_stars <- function(p) {
  if (is.na(p)) return("n/a")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  "ns"
}

#----Compute per-gene stats----
compute_gene_stats <- function(df_g, sum_g) {
  wide <- tidyr::pivot_wider(df_g, names_from = Condition, values_from = mean_log1p)
  p_MYH6 <- safe_paired_wilcox(wide$Unstim, wide$MYH6)
  p_DUF  <- safe_paired_wilcox(wide$Unstim, wide$DUF1002)
  r_MYH6 <- rank_biserial(wide$Unstim, wide$MYH6)
  r_DUF  <- rank_biserial(wide$Unstim, wide$DUF1002)
  
  
  y_top <- max(sum_g$mu + sum_g$sem, na.rm = TRUE)
  pad   <- 0.10 * y_top
  list(
    p_MYH6 = p_MYH6, p_DUF = p_DUF,
    r_MYH6 = r_MYH6, r_DUF = r_DUF,
    star_MYH6 = p_to_stars(p_MYH6),
    star_DUF  = p_to_stars(p_DUF),
    y_MYH6 = y_top + 0.8*pad,
    y_DUF  = y_top + 1.7*pad
  )
}

#----Plot per gene WITH stats----
plot_gene_single <- function(gene, out_dir = NULL,
                             width = 3.8, height = 3.4, dpi = 600) {
  df_g  <- expr_df %>% dplyr::filter(gene == !!gene)
  if (!nrow(df_g)) stop("No data for gene: ", gene)
  sum_g <- summarize_bars(df_g)
  st    <- compute_gene_stats(df_g, sum_g)
  
  p <- ggplot(sum_g, aes(x = Condition, y = mu, fill = Condition)) +
    geom_col(width = 0.58, alpha = 0.95, color = NA) +
    geom_errorbar(aes(ymin = mu - sem, ymax = mu + sem),
                  width = 0.16, linewidth = 0.55) +
    geom_point(
      data = df_g,
      aes(x = Condition, y = mean_log1p, color = Condition),
      position = position_jitter(width = 0.10, height = 0),
      size = 2.2, alpha = 0.9, stroke = 0
    ) +
    
    geom_segment(data = data.frame(y = st$y_MYH6),
                 aes(x = 1, xend = 2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.5, color = "black") +
    geom_text(data = data.frame(y = st$y_MYH6*1.03, lab = st$star_MYH6),
              aes(x = 1.5, y = y, label = lab),
              inherit.aes = FALSE, size = 3.6) +
    geom_segment(data = data.frame(y = st$y_DUF),
                 aes(x = 1, xend = 3, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.5, color = "black") +
    geom_text(data = data.frame(y = st$y_DUF*1.03, lab = st$star_DUF),
              aes(x = 2, y = y, label = lab),
              inherit.aes = FALSE, size = 3.6) +
    scale_fill_manual(values = c("Unstim"="#e4e4e4","MYH6"="#f2c8c8","DUF1002"="#c9d3f6")) +
    scale_color_manual(values = c("Unstim"="#6e6e6e","MYH6"="#8b0000","DUF1002"="#001f5b")) +
    labs(
      title = gene,
      subtitle = sprintf("MYH6 vs Unstim: %s  |  DUF vs Unstim: %s", st$star_MYH6, st$star_DUF),
      y = "Mean log1p(count) per clone", x = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.x  = element_text(angle = 18, hjust = 1),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(out_dir, paste0("BAR_", gene, ".png")), p, width = width, height = height, dpi = dpi)
    ggsave(file.path(out_dir, paste0("BAR_", gene, ".pdf")), p, width = width, height = height, dpi = dpi)
  }
  p
}

#----Supplementary Stat Tables----
export_gene_stats <- function(out_csv = NULL) {
  stats_all <- lapply(unique(expr_df$gene), function(g) {
    df_g  <- expr_df %>% dplyr::filter(gene == !!g)
    sum_g <- summarize_bars(df_g)
    st    <- compute_gene_stats(df_g, sum_g)
    data.frame(
      gene = g,
      n_pairs_MYH6 = sum(stats::complete.cases(tidyr::pivot_wider(df_g, names_from = Condition, values_from = mean_log1p))[c("Unstim","MYH6")]),
      n_pairs_DUF  = sum(stats::complete.cases(tidyr::pivot_wider(df_g, names_from = Condition, values_from = mean_log1p))[c("Unstim","DUF1002")]),
      p_MYH6 = st$p_MYH6, p_DUF = st$p_DUF,
      rrb_MYH6 = st$r_MYH6, rrb_DUF = st$r_DUF,
      stringsAsFactors = FALSE
    )
  }) %>% dplyr::bind_rows() %>% dplyr::arrange(p_MYH6)
  
  if (!is.null(out_csv)) utils::write.csv(stats_all, out_csv, row.names = FALSE)
  stats_all
}


plot_gene_single <- function(gene, out_dir = NULL,
                             width = 3.8, height = 3.4, dpi = 600) {
  df_g  <- expr_df %>% dplyr::filter(gene == !!gene)
  if (!nrow(df_g)) stop("No data for gene: ", gene)
  sum_g <- summarize_bars(df_g)
  st    <- compute_gene_stats(df_g, sum_g)
  
  
  y_max_dot <- max(df_g$mean_log1p, na.rm = TRUE)
  y_top <- max(sum_g$mu + sum_g$sem, y_max_dot, na.rm = TRUE)
  pad   <- 0.08 * y_top
  st$y_MYH6 <- y_top + 0.8 * pad
  st$y_DUF  <- y_top + 1.7 * pad
  
  p <- ggplot(sum_g, aes(x = Condition, y = mu, fill = Condition)) +
    geom_col(width = 0.58, alpha = 0.95, color = NA) +
    geom_errorbar(aes(ymin = mu - sem, ymax = mu + sem),
                  width = 0.16, linewidth = 0.55) +
    geom_point(
      data = df_g,
      aes(x = Condition, y = mean_log1p, color = Condition),
      position = position_jitter(width = 0.10, height = 0),
      size = 3.2, alpha = 0.9, stroke = 0
    ) +
    
    geom_segment(data = data.frame(y = st$y_MYH6),
                 aes(x = 1, xend = 2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.6, color = "black") +
    geom_text(data = data.frame(y = st$y_MYH6 * 1.03, lab = st$star_MYH6),
              aes(x = 1.5, y = y, label = lab),
              inherit.aes = FALSE, size = 8) +
    geom_segment(data = data.frame(y = st$y_DUF),
                 aes(x = 1, xend = 3, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 1, color = "black") +
    geom_text(data = data.frame(y = st$y_DUF * 1.03, lab = st$star_DUF),
              aes(x = 2, y = y, label = lab),
              inherit.aes = FALSE, size = 5.5) +
    scale_fill_manual(values = c("Unstim"="#e4e4e4","MYH6"="#f2c8c8","DUF1002"="#c9d3f6")) +
    scale_color_manual(values = c("Unstim"="#6e6e6e","MYH6"="#8b0000","DUF1002"="#001f5b")) +
    labs(
      title = gene,
      subtitle = sprintf("MYH6 vs Unstim: %s  |  DUF vs Unstim: %s",
                         st$star_MYH6, st$star_DUF),
      y = "Mean log1p(count) per clone", x = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.text.x  = element_text(angle = 18, hjust = 1),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(out_dir, paste0("BAR_", gene, ".png")), p, width = width, height = height, dpi = dpi)
    ggsave(file.path(out_dir, paste0("BAR_", gene, ".pdf")), p, width = width, height = height, dpi = dpi)
  }
  p
}
plot_gene_single("JUN")



format_p <- function(p) {
  if (is.na(p)) return("p = n/a")
  if (p < 1e-4) return(sprintf("p = %.1e", p))
  sprintf("p = %.4f", p)
}

plot_gene_single <- function(gene, out_dir = NULL,
                             width = 3.8, height = 3.6, dpi = 600) {
  df_g  <- expr_df %>% dplyr::filter(gene == !!gene)
  if (!nrow(df_g)) stop("No data for gene: ", gene)
  sum_g <- summarize_bars(df_g)
  st    <- compute_gene_stats(df_g, sum_g)
  
  
  y_max_dot <- max(df_g$mean_log1p, na.rm = TRUE)
  y_top     <- max(sum_g$mu + sum_g$sem, y_max_dot, na.rm = TRUE)
  pad       <- 0.10 * y_top
  y_myh6    <- y_top + 0.9 * pad
  y_duf     <- y_top + 1.9 * pad
  
  
  lab_M <- format_p(st$p_MYH6)
  lab_D <- format_p(st$p_DUF)
  
  p <- ggplot(sum_g, aes(x = Condition, y = mu, fill = Condition)) +
    geom_col(width = 0.58, alpha = 0.95, color = NA) +
    geom_errorbar(aes(ymin = mu - sem, ymax = mu + sem),
                  width = 0.16, linewidth = 0.6) +
    geom_point(
      data = df_g,
      aes(x = Condition, y = mean_log1p, color = Condition),
      position = position_jitter(width = 0.10, height = 0),
      size = 4.0, alpha = 0.95, stroke = 0
    ) +
    
    geom_segment(
      data = data.frame(y = y_myh6),
      aes(x = 1, xend = 2, y = y, yend = y),
      inherit.aes = FALSE, linewidth = 0.8, color = "black"
    ) +
    geom_segment(
      data = data.frame(y = y_myh6),
      aes(x = 1, xend = 1, y = y - 0.03*y, yend = y),
      inherit.aes = FALSE, linewidth = 0.8, color = "black"
    ) +
    geom_segment(
      data = data.frame(y = y_myh6),
      aes(x = 2, xend = 2, y = y - 0.03*y, yend = y),
      inherit.aes = FALSE, linewidth = 0.8, color = "black"
    ) +
    geom_text(
      data = data.frame(y = y_myh6 * 1.04, lab = lab_M),
      aes(x = 1.5, y = y, label = lab),
      inherit.aes = FALSE, size = 5.6
    ) +
    
    geom_segment(
      data = data.frame(y = y_duf),
      aes(x = 1, xend = 3, y = y, yend = y),
      inherit.aes = FALSE, linewidth = 0.8, color = "black"
    ) +
    geom_segment(
      data = data.frame(y = y_duf),
      aes(x = 1, xend = 1, y = y - 0.03*y, yend = y),
      inherit.aes = FALSE, linewidth = 0.8, color = "black"
    ) +
    geom_segment(
      data = data.frame(y = y_duf),
      aes(x = 3, xend = 3, y = y - 0.03*y, yend = y),
      inherit.aes = FALSE, linewidth = 0.8, color = "black"
    ) +
    geom_text(
      data = data.frame(y = y_duf * 1.04, lab = lab_D),
      aes(x = 2, y = y, label = lab),
      inherit.aes = FALSE, size = 5.6
    ) +
    scale_fill_manual(values = c("Unstim"="#e4e4e4","MYH6"="#f2c8c8","DUF1002"="#c9d3f6")) +
    scale_color_manual(values = c("Unstim"="#6e6e6e","MYH6"="#8b0000","DUF1002"="#001f5b")) +
    labs(title = gene, y = "Mean log1p(count) per clone", x = NULL) +
    theme_classic(base_size = 13) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.y  = element_text(size = 14),
      axis.text.x  = element_text(size = 14, angle = 18, hjust = 1),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      panel.grid = element_blank()
    )
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(out_dir, paste0("BAR_", gene, ".png")), p,
           width = width, height = height, dpi = dpi)
    ggsave(file.path(out_dir, paste0("BAR_", gene, ".pdf")), p,
           width = width, height = height, dpi = dpi)
  }
  p
}

plot_gene_single("NFKBIA")
cell_map %>%
  count(Condition) %>%
  arrange(Condition)


format_pnum <- function(p) {
  if (is.na(p)) return("")
  if (p < 1e-4) return(sprintf("%.1e", p))
  sprintf("%.4f", p)
}


bar_fills  <- c("Unstim"="#e4e4e4", "MYH6"="#f2c8c8", "DUF1002"="#c9d3f6")  
dot_fills  <- c("Unstim"="#6e6e6e", "MYH6"="#8b0000", "DUF1002"="navy")  

plot_gene_single <- function(gene, out_dir = NULL,
                             width = 2.8, height = 5.0, dpi = 600,
                             dot_size = 3.2, dot_stroke = 0.5) {
  
  df_g  <- expr_df %>% dplyr::filter(gene == !!gene)
  if (!nrow(df_g)) stop("No data for gene: ", gene)
  sum_g <- summarize_bars(df_g)
  st    <- compute_gene_stats(df_g, sum_g)
  
  
  y_max_dot <- max(df_g$mean_log1p, na.rm = TRUE)
  y_top     <- max(sum_g$mu + sum_g$sem, y_max_dot, na.rm = TRUE)
  pad       <- 0.12 * y_top
  y_myh6    <- y_top + 1.0 * pad
  y_duf     <- y_top + 2.0 * pad
  
  lab_M <- format_pnum(st$p_MYH6)
  lab_D <- format_pnum(st$p_DUF)
  
  p <- ggplot(sum_g, aes(x = Condition, y = mu, fill = Condition)) +
    
    geom_col(width = 0.58, alpha = 0.95, color = NA) +
    geom_errorbar(aes(ymin = mu - sem, ymax = mu + sem),
                  width = 0.16, linewidth = 0.6) +
    scale_fill_manual(values = bar_fills, drop = FALSE) +
    
    
    ggnewscale::new_scale_fill() +
    
    
    geom_point(
      data = df_g,
      aes(x = Condition, y = mean_log1p, fill = Condition),
      position = position_jitter(width = 0.08, height = 0),
      shape = 21, size = dot_size, alpha = 0.95, color = "black", stroke = dot_stroke
    ) +
    scale_fill_manual(values = dot_fills, drop = FALSE) +
    
    
    geom_segment(data = data.frame(y = y_myh6),
                 aes(x = 1, xend = 2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    geom_segment(data = data.frame(y = y_myh6),
                 aes(x = 1, xend = 1, y = y - 0.03*y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    geom_segment(data = data.frame(y = y_myh6),
                 aes(x = 2, xend = 2, y = y - 0.03*y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    geom_text(data = data.frame(y = y_myh6 * 1.03, lab = lab_M),
              aes(x = 1.5, y = y, label = lab),
              inherit.aes = FALSE, size = 4.0) +
    
    geom_segment(data = data.frame(y = y_duf),
                 aes(x = 1, xend = 3, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    geom_segment(data = data.frame(y = y_duf),
                 aes(x = 1, xend = 1, y = y - 0.03*y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    geom_segment(data = data.frame(y = y_duf),
                 aes(x = 3, xend = 3, y = y - 0.03*y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.7, color = "black") +
    geom_text(data = data.frame(y = y_duf * 1.03, lab = lab_D),
              aes(x = 2, y = y, label = lab),
              inherit.aes = FALSE, size = 4.0) +
    
    labs(title = gene, y = "Mean log1p(count) per clone", x = NULL) +
    theme_classic(base_size = 14) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 17),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.y  = element_text(size = 14),
      axis.text.x  = element_text(size = 14, angle = 20, hjust = 1),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = margin(20, 20, 10, 20)
    )
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    file_base <- file.path(out_dir, paste0("BAR_", gene))
    ggsave(paste0(file_base, ".png"), p, width = width, height = height, dpi = dpi)
    ggsave(paste0(file_base, ".pdf"), p, width = width, height = height, dpi = dpi)
  }
  p
}
plot_gene_single("JUN", width = 2.5, height = 5.5)   
plot_gene_single("TNFAIP3", width = 3, height = 5)   
plot_gene_single("IFNG", width = 2.8, height = 5.2)  

main_fig_genes <- c("TNFAIP3","NFKBIA","DUSP4","DUSP16","RASGEF1B","HIF1A","JUN")
for (g in intersect(main_fig_genes, unique(expr_df$gene))) {
  print(plot_gene_single(g,
                         out_dir = "Path to Save/Figure_bars",
                         width = 2.8, height = 5.0))
}



#----Return to RNA After Confirming the Expansion for Heatmap on Total CD4----
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject)
  library(dplyr);  library(tidyr)
  library(ComplexHeatmap); library(circlize)
})

#----TCR activation gene set----
tcr_genes <- c(
  
  "JUN","FOS","FOSB","JUNB","DUSP4","DUSP1","DUSP2","DUSP16","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3","IER2","IER3",
  
  "NFKBIA","NFKBIZ","TNFAIP3","NFKB1","NFKB2","RELA","REL",
  
  "RASGEF1B","CD69","IL2RA","CXCR3","CXCR4","CCL5",
  
  "HIF1A","LDHA",
  
  "TNF","IFNG"
)

#----CD4 selection----
get_cd4_cells <- function(seu_obj){
  md <- seu_obj@meta.data
  ann_cols <- intersect(c("CellType","celltype","cell_type","major_celltype",
                          "fine_celltype","annot","annotation"), colnames(md))
  if (length(ann_cols)) {
    for (c in ann_cols) {
      hit <- rownames(md)[grepl("CD4", md[[c]], ignore.case = TRUE)]
      if (length(hit)) return(hit)
    }
  }
  
  assay_pref <- if ("SCT" %in% names(seu_obj@assays)) "SCT" else "RNA"
  DefaultAssay(seu_obj) <- assay_pref
  feats <- intersect(c("CD3D","CD4","CD8A"), rownames(seu_obj[[assay_pref]]))
  if (length(feats) >= 2) {
    df <- FetchData(seu_obj, vars = feats, assay = assay_pref, layer = "data")
    rownames(df) <- colnames(seu_obj)
    keep <- rep(TRUE, nrow(df))
    if ("CD3D" %in% feats) keep <- keep & (df$CD3D > 0)
    if (all(c("CD4","CD8A") %in% feats)) keep <- keep & (df$CD4 >= df$CD8A)
    return(rownames(df)[keep])
  }
  warning("CD4 fallback not possible — returning all cells.")
  colnames(seu_obj)
}

cd4_cells <- get_cd4_cells(seu)
seu_cd4 <- subset(seu, cells = cd4_cells)


if (!"Condition" %in% colnames(seu_cd4@meta.data)) {
  seu_cd4$Condition <- sub(".*_(Unstim|MYH6|DUF1002)$", "\\1", seu_cd4$sample_id)
}
seu_cd4$Condition <- factor(seu_cd4$Condition, levels = c("Unstim","MYH6","DUF1002"))

#-----Expression matrix (genes × cells)----
pick_expr_source <- function(seu){
  if ("SCT" %in% names(seu@assays)) {
    lyr <- tryCatch(Layers(seu[["SCT"]]), error = function(e) character(0))
    if ("data" %in% lyr) return(list(assay="SCT", layer="data", log1p=FALSE))
  }
  if ("RNA" %in% names(seu@assays)) {
    lyr <- tryCatch(Layers(seu[["RNA"]]), error = function(e) character(0))
    if ("data" %in% lyr)   return(list(assay="RNA", layer="data",   log1p=FALSE))
    if ("counts" %in% lyr) return(list(assay="RNA", layer="counts", log1p=TRUE))
  }
  fa <- names(seu@assays)[1]
  lyr <- tryCatch(Layers(seu[[fa]]), error = function(e) character(0))
  if (length(lyr)==0) stop("No layers available in assay: ", fa)
  list(assay=fa, layer=lyr[1], log1p=TRUE)
}

src <- pick_expr_source(seu_cd4)
DefaultAssay(seu_cd4) <- src$assay
if (!src$layer %in% Layers(seu_cd4[[src$assay]]))
  stop("Layer '", src$layer, "' not found in assay '", src$assay, "'.")

#----Intersect genes with object----
genes_present <- intersect(tcr_genes, rownames(seu_cd4[[src$assay]]))
if (!length(genes_present)) stop("None of the TCR genes present in assay: ", src$assay)

#----features x cells----
M <- GetAssayData(seu_cd4, assay = src$assay, layer = src$layer)
M <- M[genes_present, colnames(seu_cd4), drop = FALSE]
if (isTRUE(src$log1p)) M <- log1p(M)

#----Expression summary by condition (genes × conditions)----
cond <- seu_cd4$Condition
expr_means <- sapply(levels(cond), function(k) {
  cols <- names(cond)[which(cond == k)]
  if (length(cols)==0) return(rep(NA_real_, nrow(M)))
  rowMeans(M[, cols, drop = FALSE], na.rm = TRUE)
})
rownames(expr_means) <- rownames(M)  

#----Z-score per gene across the 3 conditions----
z_expr <- t(scale(t(expr_means)))   
z_expr[is.na(z_expr)] <- 0

#----Differential expression----
suppressWarnings(try({ seu_cd4 <- JoinLayers(seu_cd4) }, silent = TRUE))
DefaultAssay(seu_cd4) <- if ("SCT" %in% names(seu_cd4@assays)) "SCT" else "RNA"
if (DefaultAssay(seu_cd4) == "SCT") {
  suppressMessages({ seu_cd4 <- PrepSCTFindMarkers(seu_cd4, verbose = FALSE) })
}
Idents(seu_cd4) <- seu_cd4$Condition

.safe_join <- function(obj) {
  suppressWarnings(try({ obj <- JoinLayers(obj) }, silent = TRUE))
  obj
}
.ensure_condition_idents <- function(obj) {
  if (!"Condition" %in% colnames(obj@meta.data)) {
    obj$Condition <- sub(".*_(Unstim|MYH6|DUF1002)$", "\\1", obj$sample_id)
  }
  obj$Condition <- factor(obj$Condition, levels = c("Unstim","MYH6","DUF1002"))
  Idents(obj) <- obj$Condition
  obj
}
.ensure_rna_norm_and_joined <- function(obj) {
  if (!"RNA" %in% names(obj@assays)) stop("No RNA assay to fall back to.")
  DefaultAssay(obj) <- "RNA"
  
  if (!"data" %in% Layers(obj[["RNA"]])) {
    message("Normalizing RNA for DE fallback...")
    obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", verbose = FALSE)
  }
  
  obj <- .safe_join(obj)
  DefaultAssay(obj) <- "RNA"
  obj
}

run_markers_pair <- function(obj, ident1, ident2,
                             min.pct = 0.05,
                             logfc.threshold = 0) {
  
  if (!"RNA" %in% names(obj@assays)) {
    stop("RNA assay not found in object.")
  }
  
  DefaultAssay(obj) <- "RNA"
  
  
  if (!"data" %in% Layers(obj[["RNA"]])) {
    message("Normalizing RNA for DE...")
    obj <- NormalizeData(
      obj,
      assay = "RNA",
      normalization.method = "LogNormalize",
      verbose = FALSE
    )
  }
  
  
  suppressWarnings({
    obj <- JoinLayers(obj, assay = "RNA")
  })
  
  res <- FindMarkers(
    obj,
    ident.1         = ident1,
    ident.2         = ident2,
    group.by        = "Condition",
    min.pct         = min.pct,
    logfc.threshold = logfc.threshold,
    test.use        = "wilcox",
    assay           = "RNA"
  )
  
  message("DE ran on RNA (Wilcoxon, joined layers).")
  res
}

duf_out <- run_markers_pair(seu_cd4, "DUF1002", "Unstim")
m6_out  <- run_markers_pair(seu_cd4, "MYH6",    "Unstim")

de_duf <- duf_out$res
de_m6  <- m6_out$res



#----DotPlot for Manuscript----
#----Genes----
genes_TCR_IE <- c("NR4A1","NR4A3","EGR2","EGR3","RGS1","SOCS1")
genes_MAPK   <- c("DUSP2","DUSP4","DUSP16")
genes_NFKB   <- c("TRAF2","NFKBIA","TNFAIP3","RELB","NFKB2","RASGEF1B")
genes_ISG    <- c("IFI6","IF I44L","ISG15","MX1","IFIT1","IFIT3","RSAD2","OAS1","OASL","CMPK2","HERC5","PLSCR1","IRF7")

genes_ISG[genes_ISG == "IF I44L"] <- "IFI44L"
genes_AP     <- c("CTSL","LAMP3")
genes_MET    <- c("PFKFB3","LDHA","HIF1A")


cohorts <- list(
  "TCR / immediate-early"    = genes_TCR_IE,
  "MAPK feedback"            = genes_MAPK,
  "NF-κB pathway"            = genes_NFKB,
  "Interferon / ISG program" = genes_ISG,
  "Antigen processing"       = genes_AP,
  "Metabolic activation"     = genes_MET
)


get_cd4_cells <- function(seu_obj){
  md <- seu_obj@meta.data
  ann_cols <- intersect(c("CellType","celltype","cell_type","major_celltype",
                          "fine_celltype","annot","annotation"), colnames(md))
  if (length(ann_cols)) {
    for (c in ann_cols) {
      hit <- rownames(md)[grepl("CD4", md[[c]], ignore.case = TRUE)]
      if (length(hit)) return(hit)
    }
  }
  
  assay_pref <- if ("SCT" %in% names(seu_obj@assays)) "SCT" else "RNA"
  DefaultAssay(seu_obj) <- assay_pref
  feats <- intersect(c("CD3D","CD4","CD8A"), rownames(seu_obj[[assay_pref]]))
  if (length(feats) >= 2) {
    df <- FetchData(seu_obj, vars = feats, assay = assay_pref, layer = "data")
    rownames(df) <- colnames(seu_obj)
    keep <- rep(TRUE, nrow(df))
    if ("CD3D" %in% feats) keep <- keep & (df$CD3D > 0)
    if (all(c("CD4","CD8A") %in% feats)) keep <- keep & (df$CD4 >= df$CD8A)
    return(rownames(df)[keep])
  }
  colnames(seu_obj)
}

stopifnot(exists("seu"))
if (!"Condition" %in% colnames(seu@meta.data)) {
  seu$Condition <- sub(".*_(Unstim|MYH6|DUF1002)$", "\\1", seu$sample_id)
}
seu$Condition <- factor(seu$Condition, levels = c("Unstim","MYH6","DUF1002"))

cd4_cells <- get_cd4_cells(seu)
seu_use   <- subset(seu, cells = cd4_cells)

#----Rebuilding Robust assay/layer----
pick_expr_source <- function(seu_obj){
  if ("SCT" %in% names(seu_obj@assays) && "data" %in% Layers(seu_obj[["SCT"]]))
    return(list(assay="SCT", layer="data"))
  if ("RNA" %in% names(seu_obj@assays)) {
    lyr <- Layers(seu_obj[["RNA"]])
    if ("data" %in% lyr)   return(list(assay="RNA", layer="data"))
    if ("counts" %in% lyr) return(list(assay="RNA", layer="counts"))
  }
  fa <- names(seu_obj@assays)[1]
  lyr <- Layers(seu_obj[[fa]])
  list(assay=fa, layer=ifelse(length(lyr), lyr[1], "data"))
}
src <- pick_expr_source(seu_use)
DefaultAssay(seu_use) <- src$assay

#----Supp Table (Condition × gene) with mean & % detected-----
genes_all <- unique(unlist(cohorts))
genes_present <- intersect(genes_all, rownames(seu_use[[src$assay]]))
dropped <- setdiff(genes_all, genes_present)
if (length(dropped)) message("Dropped (not in assay): ", paste(dropped, collapse = ", "))

mat <- FetchData(
  seu_use, vars = genes_present, cells = colnames(seu_use),
  assay = src$assay, layer = src$layer
) |> as.matrix()

long_df <- as.data.frame(mat) %>%
  tibble::rownames_to_column("cell") %>%
  mutate(Condition = seu_use$Condition[cell]) %>%
  pivot_longer(-c(cell, Condition), names_to = "gene", values_to = "expr")

dotdf <- long_df %>%
  group_by(Condition, gene) %>%
  summarise(
    mean = mean(expr, na.rm = TRUE),
    pct  = 100 * mean(expr > 0, na.rm = TRUE),
    .groups = "drop"
  )

#----Reorder genes ASCENDING within each cohort----

map_tbl <- lapply(names(cohorts), function(nm){
  data.frame(gene = intersect(cohorts[[nm]], genes_present),
             Pathway = nm, stringsAsFactors = FALSE)
}) |> bind_rows()

dotdf <- dotdf %>% inner_join(map_tbl, by = "gene")


rank_tbl <- dotdf %>%
  filter(Condition %in% c("MYH6","DUF1002")) %>%
  group_by(Pathway, gene) %>%
  summarise(score = max(mean), .groups = "drop") %>%
  arrange(Pathway, score)  


dotdf$gene    <- factor(dotdf$gene, levels = rank_tbl$gene)
dotdf$Pathway <- factor(dotdf$Pathway, levels = names(cohorts))  

#----Plotting----
p_dot <- ggplot(dotdf, aes(x = Condition, y = gene)) +
  geom_point(aes(size = pct, color = mean), stroke = 0.25, shape = 16) +
  scale_size_area(max_size = 7, breaks = c(10,25,50,75,100), name = "% cells") +
  scale_color_viridis_c(name = "Mean expr", option = "C") +
  facet_grid(Pathway ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Pathway-resolved activation — CD4 cells", x = NULL, y = NULL) +
  theme_classic(base_size = 13) +  
  theme(
    
    strip.text       = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey95", color = NA),
    
    
    axis.text.y      = element_text(face = "italic", size = 14),   
    axis.text.x      = element_text(angle = 90, hjust = 1, size = 14),  
    
    axis.ticks.y     = element_blank(),
    panel.grid       = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 15)
  )


#----Preview----
p_dot

#----Saving----
out_dir <- "Path to save/Dotplot_Cohorts_ASC_CD4"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_dir, "Dotplot_Cohorts_ASC_CD4.png"), p_dot, width = 4, height = 9, dpi = 600)
ggsave(file.path(out_dir, "Dotplot_Cohorts_ASC_CD4.pdf"),  p_dot, width = 4, height = 9)
