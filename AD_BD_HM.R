#====== Set up working directory and packages ======
setwd("/path/to/working/directory/")

library(readxl)      # To read Excel metadata
library(vegan)       # For alpha diversity
library(tidyverse)   # Data manipulation
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
#install MaAsLin2 if needed 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("Maaslin2")

library(Maaslin2)
library(readr)
library(stringr)
library(forcats)
library(scales)
library(patchwork)




#====== Load in data and preprocessing Species ======

#load in taxon data from humann processing 
buglist <- read_tsv("/path/to/buglist.tsv")

# Define ending patterns
samples_to_remove <- c("69$", "73$", "66$", "78$", "43$", "37$","90$")
pattern <- paste(samples_to_remove, collapse = "|")


# View data structure
head(buglist)

# Retain only sample ID 
colnames(buglist) <- colnames(buglist) %>%
  str_replace("^([^_]+)_.*", "\\1")

# Keep clades and columns starting with "c"
buglist_filtered <- buglist %>%
  select(clade, starts_with("c"))

# Extract species-level (s__) names explicitly, removing t__ if present
buglist_filtered <- buglist_filtered %>%
  mutate(clade = str_extract(clade, "s__[^|]+")) %>%
  mutate(clade = str_remove(clade, "s__")) %>%
  drop_na(clade)

#rename incorrecrly classsified Porphyromonas bobii to Porphyromonas gingivalis 
buglist_filtered <- buglist_filtered %>%
  mutate(clade = str_replace(clade, "Porphyromonas_bobii", "Porphyromonas_gingivalis"))


# remove duplicate species, keep first occurrence
buglist_filtered <- buglist_filtered %>%
  distinct(clade, .keep_all = TRUE)

# Keep clade column and columns not matching the pattern
buglist_filtered <- buglist_filtered %>%
  select(clade, !matches(pattern))


# Verify
head(buglist_filtered)

# Inspect resulting column
head(buglist_filtered$clade)

# load in metadata file
metadata <- read_excel("/path/to/metadata.xlsx")

# Ensure that the metadata has a column named 'Sample' matching HUMAnN3 data column names
head(metadata)

# Simplify Sample column to keep only text before the first underscore
metadata <- metadata %>%
  mutate(Sample = str_replace(Sample, "^([^_]+)_.*", "\\1"))

# Remove 'ICM' from the SampleType column
metadata$SampleType <- gsub("ICM", "", metadata$SampleType)

# Remove rows where SampleType is now an empty string
metadata <- metadata[metadata$SampleType != "", ]

metadata <- metadata %>%
  filter(!str_detect(Sample, pattern))

# Verify the change
head(metadata)


#====== PCoA ======

# recreate numeric data
buglist_numeric <- buglist_filtered %>%
  column_to_rownames("clade") %>%
  t() %>%
  as.data.frame()

# Create the 'Sample' column from the old column names which are row names
buglist_numeric$Sample <- rownames(buglist_numeric)

# Remove the row names now as they are column names
rownames(buglist_numeric) <- NULL

# Reorder the columns so 'Sample' is the first column
buglist_numeric <- buglist_numeric[, c("Sample", setdiff(names(buglist_numeric), "Sample"))]

# View the result
head(buglist_numeric)

# Check for NAs and fix
na_count <- sum(is.na(buglist_numeric))
if (na_count > 0) {
  buglist_numeric <- buglist_numeric %>% select(where(~ !any(is.na(.))))
}

# Verify final clean data
dim(buglist_numeric)
str(buglist_numeric)


# Merging the SampleType from metadata into buglist_filtered by the 'Sample' column
buglist <- merge(buglist_numeric, metadata[, c("Sample", "SampleType")], by = "Sample", all.x = TRUE)

# Reordering columns so that SampleType is the second column
buglist <- buglist %>%
  select(Sample, SampleType, everything())

# Reorder columns so SampleType is second
buglist <- buglist %>%
  select(Sample, SampleType, everything()) %>%
  drop_na()  # remove rows with any NA

# Print the result
print(buglist)

# Set rownames only if Sample column exists
rownames(buglist) <- buglist$Sample
buglist <- buglist[, -1]  # Drop 'Sample' column if needed
buglist <- buglist[, -1] #drop SampleType column

# Match metadata order to abundance data
metadata <- metadata[match(rownames(buglist), metadata$Sample), ]

# ---------- Normalize to relative abundance ----------
abund_rel <- buglist / rowSums(buglist)

# Write to CSV
write.csv(abund_rel, file = "/path/to/data/files/abund_rel_species.csv", row.names = T)

# ---------- Bray-Curtis distance ----------
bray_dist <- vegdist(abund_rel, method = "bray")

# ---------- PCoA ----------
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 2)  # k=2 for 2D
pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample <- rownames(pcoa_df)

# Merge with metadata
pcoa_df <- left_join(pcoa_df, metadata, by = "Sample")


# PERMANOVA to test for significant differences in composition
permanova_result <- adonis2(bray_dist ~ SampleType, data = metadata)
permanova_result <- adonis2(bray_dist ~ SampleType, data = metadata)
pval <- permanova_result$`Pr(>F)`[1]
r2   <- permanova_result$R2[1]


label_text <- paste0("PERMANOVA: R² = ", round(r2, 3), 
                     ", p = ", format.pval(pval, digits = 3, eps = .001))



eig_vals <- pcoa_res$eig
var_explained <- eig_vals / sum(eig_vals)
x_lab <- paste0("PCoA 1 (", round(var_explained[1] * 100, 1), "%)")
y_lab <- paste0("PCoA 2 (", round(var_explained[2] * 100, 1), "%)")

PCoA_S <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = SampleType)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linetype = 2) +
  theme_minimal() +
  labs(
    title = "Microbial Community at Species Level PCoA based on Bray-Curtis Distance",
    x = x_lab,
    y = y_lab,
    color = "Sample Type",
    caption = label_text
  )

ggsave("/path/to/plots/PCoA_S.svg", plot = PCoA_S, device = "svg", width = 10, height = 8)



#======heatmap======

# ---------- HUMAnN buglist ----------
buglist <- read_tsv("/path/to/buglist.tsv")

# samples to remove
samples_to_remove <- c("69$", "73$", "66$", "78$", "43$", "37$", "90$")
pattern <- paste(samples_to_remove, collapse = "|")

# keep only simplified sample IDs
colnames(buglist) <- colnames(buglist) %>%
  str_replace("^([^_]+)_.*", "\\1")

# keep clade + sample columns
buglist_filtered <- buglist %>%
  select(clade, starts_with("c"))

# keep species-level clades (s__) and drop t__
buglist_filtered <- buglist_filtered %>%
  mutate(clade = str_extract(clade, "s__[^|]+")) %>%
  mutate(clade = str_remove(clade, "s__")) %>%
  drop_na(clade)

# fix Porphyromonas bobii → Porphyromonas gingivalis
buglist_filtered <- buglist_filtered %>%
  mutate(clade = str_replace(clade, "Porphyromonas_bobii", "Porphyromonas_gingivalis"))

# remove duplicate species
buglist_filtered <- buglist_filtered %>%
  distinct(clade, .keep_all = TRUE)

# remove unwanted samples by suffix
buglist_filtered <- buglist_filtered %>%
  select(clade, !matches(pattern))

# ---------- convert to numeric matrix (samples × species) ----------
buglist_numeric <- buglist_filtered %>%
  column_to_rownames("clade") %>%
  t() %>%
  as.data.frame()

buglist_numeric$Sample <- rownames(buglist_numeric)
rownames(buglist_numeric) <- NULL

# reorder so Sample first
buglist_numeric <- buglist_numeric[, c("Sample", setdiff(names(buglist_numeric), "Sample"))]

# remove any columns with NAs (if present)
na_count <- sum(is.na(buglist_numeric))
if (na_count > 0) {
  buglist_numeric <- buglist_numeric %>% select(where(~ !any(is.na(.))))
}

#load and clean metadata

metadata <- read_excel("/path/to/metadata.xlsx")

metadata <- metadata %>%
  mutate(Sample = str_replace(Sample, "^([^_]+)_.*", "\\1"))

# remove "ICM" from SampleType and empty rows
metadata$SampleType <- gsub("ICM", "", metadata$SampleType)
metadata <- metadata[metadata$SampleType != "", ]

# drop removed samples
metadata <- metadata %>%
  filter(!str_detect(Sample, pattern))


#Merge abundance with metadata

buglist <- merge(
  buglist_numeric,
  metadata[, c("Sample", "SampleType")],
  by = "Sample",
  all.x = TRUE
)

buglist <- buglist %>%
  select(Sample, SampleType, everything()) %>%
  drop_na()

# for phyloseq: make sample IDs rownames and drop Sample/SampleType
rownames(buglist) <- buglist$Sample
abundance_mat <- buglist %>%
  select(-Sample, -SampleType) %>%
  as.matrix()


#Build taxonomy table from HUMAnN clade strings

full_tax_table <- read_tsv("/path/to/buglist.tsv")

colnames(full_tax_table) <- colnames(full_tax_table) %>%
  str_replace("^([^_]+)_.*", "\\1")

full_tax_table <- full_tax_table %>%
  select(clade, starts_with("c"))

tax_strings <- full_tax_table %>% select(clade)

# keep full taxonomy strings (k__ ... s__ ... t__)
tax_strings <- tax_strings %>%
  filter(grepl("^k__.*\\|p__.*\\|c__.*\\|o__.*\\|f__.*\\|g__.*\\|s__.*\\|t__", clade))

# strip t__ tail
tax_strings$clade <- sub("\\|t__.*", "", tax_strings$clade)

# split into levels
tax_df <- tidyr::separate(
  tax_strings,
  col  = clade,
  into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  sep  = "\\|",
  remove = TRUE
)

# drop duplicate taxonomy rows
tax_df_unique <- tax_df %>% distinct()

# clean prefixes (k__ etc.)
tax_df_clean <- tax_df_unique %>%
  mutate(across(everything(), ~ sub("^[a-z]__*", "", .))) %>%
  as.data.frame()

# fix species naming consistency
tax_df_clean <- tax_df_clean %>%
  mutate(Species = str_replace(Species, "Porphyromonas_bobii", "Porphyromonas_gingivalis"))

# rownames = species (must match colnames of abundance_mat)
rownames(tax_df_clean) <- make.unique(tax_df_clean$Species)

# subset taxonomy to only species present in abundance_mat
tax_df_clean <- tax_df_clean[rownames(tax_df_clean) %in% colnames(abundance_mat), ]

# ensure same order
tax_df_clean <- tax_df_clean[match(colnames(abundance_mat), rownames(tax_df_clean)), ]


#Build phyloseq object

otu_phy <- otu_table(t(abundance_mat), taxa_are_rows = TRUE)

metadata_phylo <- metadata %>%
  filter(Sample %in% rownames(abundance_mat)) %>%
  as.data.frame()
rownames(metadata_phylo) <- metadata_phylo$Sample
metadata_phylo$Sample <- NULL

meta_phy <- sample_data(metadata_phylo)

tax_table_matrix <- as.matrix(tax_df_clean)
tax_phy <- tax_table(tax_table_matrix)

ps <- phyloseq(otu_phy, meta_phy, tax_phy)

#====== phyloseq alpha diversity ======
AD_S_SHAN_SIMP <-plot_richness(ps, x = "SampleType", measures = c("observed", "Shannon","Simpson")) +
  geom_boxplot(aes(fill = SampleType), alpha = 0.4)+
  theme_minimal()+
  labs(title = "Species Level Alpha Diversity by Sample Type ")

ggsave("/path/to/plots/AD_S_SHAN_SIMP.svg", plot = AD_S_SHAN_SIMP, device = "svg", width = 10, height = 8)



#6) Final abundance heatmap (top 30 taxa, log1p, split by SampleType)


# extract OTU matrix (samples × taxa)
otu_matrix <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) {
  otu_matrix <- t(otu_matrix)
}

# log1p transform raw counts
otu_matrix <- log1p(otu_matrix)

# transpose for ComplexHeatmap: taxa = rows, samples = cols
otu_matrix <- t(otu_matrix)

# top 30 taxa by total abundance
taxa_sums <- rowSums(otu_matrix)
top_30_taxa <- names(sort(taxa_sums, decreasing = TRUE)[1:30])
otu_matrix_top30 <- otu_matrix[top_30_taxa, ]

# sample annotations (SampleType)
sample_ann <- data.frame(SampleType = sample_data(ps)$SampleType)
rownames(sample_ann) <- colnames(otu_matrix_top30)

sample_type_levels <- unique(sample_ann$SampleType)
sample_type_colors <- setNames(
  viridis(length(sample_type_levels), option = "D"),
  sample_type_levels
)
sample_colors <- list(SampleType = sample_type_colors)

# color scale for log1p counts
color_scale <- colorRamp2(
  c(0, 1.5, 3),
  inferno(256)[c(1, 128, 256)]
)

# draw heatmap (on screen)
ht <- Heatmap(
  otu_matrix_top30,
  name = "Log1p Count",
  col = color_scale,
  top_annotation = HeatmapAnnotation(df = sample_ann, col = sample_colors),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = sample_ann$SampleType,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 8),
  column_title = "Samples",
  row_title = "Top 30 Taxa"
)
draw(ht)




#====== Sanguinis plot ======


# Boxplots with MaAsLin2 p & q only (Ss)


features <- buglist
metadata_Masslin2 <- metadata

# Convert to a base R data frame to allow row names
metadata_Masslin2 <- as.data.frame(metadata_Masslin2)
rownames(metadata_Masslin2) <- metadata_Masslin2$Sample
metadata_Masslin2$Sample <- NULL


# Check how many sample names match between metadata and features
length(intersect(rownames(metadata_Masslin2), rownames(features))) == nrow(metadata)

fit_data <- Maaslin2(
  input_data = features,
  input_metadata = metadata_Masslin2,
  output = "maaslin2_output",        # output folder will be created
  fixed_effects = c("SampleType")    # test this variable
)


# --- Params ---
out_dir <- "/path/to/CVD_S_Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sp_ss <- "Streptococcus_sanguinis"

stopifnot(exists("buglist"), exists("metadata"))

# ========= 1) Align metadata & abundance =========
normalize_id <- function(x) {
  x %>%
    as.character() %>%
    trimws() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[/\\\\]", "_") %>%
    str_replace_all("\\.fastq\\.gz$|\\.fastq$|\\.fq\\.gz$|\\.fq$", "")
}

# Align IDs
buglist2 <- buglist
rownames(buglist2) <- normalize_id(rownames(buglist2))

metadata2 <- metadata %>%
  mutate(Sample = normalize_id(Sample))

common_ids <- intersect(rownames(buglist2), metadata2$Sample)
if (length(common_ids) == 0) stop("No overlapping sample IDs after normalization.")

buglist2  <- buglist2[common_ids, , drop = FALSE]
metadata2 <- filter(metadata2, Sample %in% common_ids)

# Ensure species exists
if (!sp_ss %in% colnames(buglist2)) {
  stop(paste0("Missing species in buglist: ", sp_ss))
}

# Build df for Streptococcus sanguinis only

buglist2 <- buglist2 %>%
  as.data.frame()

# If there is already a Sample column, drop it – we want normalized rownames instead
if ("Sample" %in% colnames(buglist2)) {
  buglist2$Sample <- NULL
}

df_ss <- buglist2 %>%
  tibble::rownames_to_column("Sample") %>%
  select(Sample, !!sp_ss) %>%
  left_join(metadata2, by = "Sample") %>%
  rename(RelativeAbundance = !!sp_ss) %>%
  mutate(
    RelativeAbundance = suppressWarnings(as.numeric(RelativeAbundance)),
    # map sample type labels for plotting
    SampleType = case_when(
      str_detect(SampleType, regex("^surviv",  ignore_case = TRUE)) ~ "Survivors",
      str_detect(SampleType, regex("^progres", ignore_case = TRUE)) ~ "Progressors",
      TRUE ~ as.character(SampleType)
    ),
    SampleType = fct_relevel(factor(SampleType), "Survivors", "Progressors")
  ) %>%
  filter(!is.na(SampleType), !is.na(RelativeAbundance))

stopifnot(nrow(df_ss) > 0)

# ========= 2) MaAsLin2 subtitle (p, q, stars) =========
get_subtitle_pq <- function(species) {
  res_path <- file.path("maaslin2_output", "all_results.tsv")
  if (!file.exists(res_path)) return(NULL)
  
  maaslin_res <- read_tsv(res_path, show_col_types = FALSE)
  row <- maaslin_res %>%
    filter(metadata == "SampleType") %>%
    mutate(Species = str_replace_all(feature, " ", "_")) %>%
    filter(Species == species) %>%
    slice_min(qval, n = 1, with_ties = FALSE) %>%
    mutate(
      stars = case_when(
        qval < 0.001 ~ "***",
        qval < 0.01  ~ "**",
        qval < 0.05  ~ "*",
        TRUE         ~ "ns"
      ),
      p_str = ifelse(is.finite(pval), formatC(pval, format = "g", digits = 3), "NA"),
      q_str = ifelse(is.finite(qval), formatC(qval, format = "g", digits = 3), "NA")
    )
  
  if (nrow(row) == 0) return(NULL)
  sprintf("p=%s • q=%s • %s", row$p_str, row$q_str, row$stars)
}

subtitle_text <- get_subtitle_pq(sp_ss)

# ========= 3) Summary stats for mean ± SE =========
summ <- df_ss %>%
  group_by(SampleType) %>%
  summarise(
    mean = mean(RelativeAbundance, na.rm = TRUE),
    se   = sd(RelativeAbundance, na.rm = TRUE) / sqrt(sum(!is.na(RelativeAbundance))),
    .groups = "drop"
  )

# ========= 4) Y-axis formatting =========
y_is_fraction <- max(df_ss$RelativeAbundance, na.rm = TRUE) <= 1
y_scale <- if (y_is_fraction) scale_y_continuous(labels = percent_format(accuracy = 1)) else scale_y_continuous()
y_lab   <- if (y_is_fraction) "Relative Abundance (%)" else "Relative Abundance"

# ========= 5) Final plot =========
p_ss <- ggplot(df_ss, aes(x = SampleType, y = RelativeAbundance)) +
  # jittered points by group
  geom_point(
    aes(fill = SampleType),
    shape = 21, size = 2.8, stroke = 0.8, color = "black",
    position = position_jitter(width = 0.08, height = 0)
  ) +
  # mean ± SE
  geom_errorbar(
    data = summ,
    aes(x = SampleType, ymin = mean - se, ymax = mean + se),
    inherit.aes = FALSE,
    width = 0.2, size = 0.7, color = "black"
  ) +
  geom_segment(
    data = summ,
    aes(
      x = as.numeric(SampleType) - 0.1,
      xend = as.numeric(SampleType) + 0.1,
      y = mean,
      yend = mean
    ),
    inherit.aes = FALSE,
    linewidth = 0.7,
    color = "black",
    lineend = "round"
  ) +
  scale_fill_manual(values = c("Survivors" = "#22c55e", "Progressors" = "#ef4444")) +
  y_scale +
  labs(
    title    = "Streptococcus sanguinis",
    subtitle = subtitle_text,
    x = NULL,
    y = y_lab
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

p_ss


ggsave(
  "/path/to/plots/Boxplot_Ss_custom.svg",
  p_ss, width = 6.5, height = 5.2, device = "svg"
)







#====== Data Genus ======
#load in taxon data 
buglist_genus <- read_tsv("/path/to/buglist.tsv")

samples_to_remove_genus <- c("69$", "73$", "66$", "78$", "43$", "37$")
pattern_genus <- paste(samples_to_remove_genus, collapse = "|")

colnames(buglist_genus) <- colnames(buglist_genus) %>%
  str_replace("^([^_]+)_.*", "\\1")

buglist_filtered_genus <- buglist_genus %>%
  select(clade, starts_with("c")) %>%
  mutate(
    clade = str_extract(clade, "g__[^|]+"),
    clade = str_remove(clade, "g__")
  ) %>%
  drop_na(clade) %>%
  distinct(clade, .keep_all = TRUE) %>%
  select(clade, !matches(pattern_genus))

metadata_genus <- read_excel("/path/to/metadata.xlsx") %>%
  mutate(Sample = str_replace(Sample, "^([^_]+)_.*", "\\1")) %>%
  mutate(SampleType = gsub("ICM", "", SampleType)) %>%
  filter(SampleType != "") %>%
  filter(!str_detect(Sample, pattern_genus))

buglist_numeric_genus <- buglist_filtered_genus %>%
  column_to_rownames("clade") %>%
  t() %>%
  as.data.frame()

buglist_numeric_genus$Sample <- rownames(buglist_numeric_genus)
rownames(buglist_numeric_genus) <- NULL

buglist_numeric_genus <- buglist_numeric_genus %>%
  select(Sample, everything()) %>%
  select(where(~ !any(is.na(.))))

buglist_genus <- merge(
  buglist_numeric_genus,
  metadata_genus[, c("Sample", "SampleType")],
  by = "Sample",
  all.x = TRUE
) %>%
  select(Sample, SampleType, everything()) %>%
  drop_na()

rownames(buglist_genus) <- buglist_genus$Sample
buglist_genus_numeric <- buglist_genus %>% select(-Sample, -SampleType)
abund_rel_genus <- buglist_genus_numeric / rowSums(buglist_genus_numeric)

# Boxplot data
selected_genus <- c("Porphyromonas","GGB4936","GGB12441",
                    "Streptococcus","Rothia","Lachnoanaerobaculum","GGB1202")

plot_box_genus <- abund_rel_genus[, selected_genus]
plot_box_genus$Sample <- rownames(plot_box_genus)

plot_box_genus <- left_join(
  plot_box_genus,
  metadata_genus[, c("Sample","SampleType")],
  by = "Sample"
)

plot_box_long_genus <- pivot_longer(
  plot_box_genus,
  cols = all_of(selected_genus),
  names_to = "Genus",
  values_to = "RelativeAbundance"
)

plot_box_long_genus$Genus <- factor(plot_box_long_genus$Genus, levels = selected_genus)


Boxplots_G <-ggplot(plot_box_long_genus, aes(x = SampleType, y = RelativeAbundance, fill = SampleType)) +
  geom_boxplot(outlier.shape = NA, outlier.size = 2, width = 0.6) +
  geom_jitter(width = 0.2, size = 3.5, alpha = 0.6, shape = 21, color = "black")+
  facet_wrap(~ Genus, scales = "free_y", nrow = 1) +
  labs(
    title = "Relative Abundance of Selected Genus by Sample Type",
    x = "Sample Type",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

ggsave("/path/to/plots/Boxplots_G.svg", plot = Boxplots_G , device = "svg", width = 10, height = 6)


#====== Multi linear regression with covariant adjustment, MaAsLin2 Genera ======


features_g <- buglist_genus
metadata_g_Masslin2 <- metadata_genus

# Convert to a base R data frame to allow row names
metadata_g_Masslin2 <- as.data.frame(metadata_g_Masslin2)
rownames(metadata_g_Masslin2) <- metadata_g_Masslin2$Sample
metadata_g_Masslin2$Sample <- NULL


# Check how many sample names match between metadata and features
length(intersect(rownames(metadata_g_Masslin2), rownames(features_g))) == nrow(metadata_genus)


fit_data <- Maaslin2(
  input_data = features_g,
  input_metadata = metadata_g_Masslin2,
  output = "maaslin2_g_output",        # output folder will be created
  fixed_effects = c("SampleType")    # test this variable
)























