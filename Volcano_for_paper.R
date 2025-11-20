#set working directory
setwd("/path/to/working/directory")

# Required packages
library(tidyverse)
library(ggrepel)

#-------------------------
# 1. Load and merge data
#-------------------------

taxa_codes   <- read_csv("/path/to/taxa_codes.csv")   # Taxonomy + Code
stat_results <- read_csv("/path/to/rnaseq_de.csv")    # Code/log2FC/Pvalues, etc.

# Make sure column names line up
colnames(taxa_codes)       <- c("Taxonomy", "Code")   # adjust if needed
colnames(stat_results)[1]  <- "Code"                  # first column = Code

merged_results <- stat_results %>%
  left_join(taxa_codes, by = "Code")

#-------------------------
# 2. Clean taxonomy and create Genus_Species
#-------------------------

merged_cleaned <- merged_results %>%
  # keep rows that have both genus and species
  filter(str_detect(Taxonomy, "g__") & str_detect(Taxonomy, "s__")) %>%
  mutate(
    Genus   = str_extract(Taxonomy, "g__[^;]+") %>% str_remove("g__"),
    Species = str_extract(Taxonomy, "s__[^;]+") %>% str_remove("s__")
  ) %>%
  # remove Genus prefix if duplicated in Species (e.g., "salivarius_salivarius")
  mutate(
    Species = str_replace(Species, paste0("^", Genus, "_"), "")
  ) %>%
  mutate(
    Genus_Species = paste(Genus, Species)
  ) %>%
  distinct(Genus_Species, .keep_all = TRUE)

#-------------------------
# 3. Create stats for volcano
#-------------------------

merged_cleaned <- merged_cleaned %>%
  mutate(
    neg_log10_Pval = -log10(Pvalues),
    significance = case_when(
      Pvalues < 0.05 & log2FC >  1 ~ "High CD4",
      Pvalues < 0.05 & log2FC < -1 ~ "Low CD4",
      TRUE                         ~ "Not Significant"
    )
  )

# Top 5 upregulated and top 5 downregulated (for labeling)
upregulated <- merged_cleaned %>%
  filter(significance == "High CD4") %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = 5)

downregulated <- merged_cleaned %>%
  filter(significance == "Low CD4") %>%
  arrange(log2FC) %>%   # most negative first
  slice_head(n = 5)

#-------------------------
# 4. volcano_top10 (with labels)
#-------------------------

volcano_top10 <- ggplot(merged_cleaned, aes(x = log2FC, y = neg_log10_Pval)) +
  # base layer: all points
  geom_point(
    aes(fill = significance),
    shape = 21,
    stroke = 0,
    size = 3,
    alpha = 0.8
  ) +
  # outline for significant points
  geom_point(
    data = merged_cleaned %>% filter(significance != "Not Significant"),
    aes(fill = significance),
    shape = 21,
    color = "black",
    stroke = 0.8,
    size = 3,
    alpha = 0.8
  ) +
  # labels for top 5 upregulated
  geom_text_repel(
    data = upregulated,
    aes(label = Genus_Species),
    size = 5.5,
    color = "black",
    nudge_x = 0.5,
    nudge_y = 0.5,
    max.overlaps = 20
  ) +
  # labels for top 5 downregulated
  geom_text_repel(
    data = downregulated,
    aes(label = Genus_Species),
    size = 5.5,
    color = "black",
    nudge_x = -0.5,
    nudge_y = 0.5,
    max.overlaps = 20
  ) +
  scale_fill_manual(values = c(
    "High CD4"       = "red",
    "Low CD4"        = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(
    x    = "Log2 Fold Change",
    y    = "-Log10 p-value",
    fill = "Enriched in"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1),      linetype = "dashed", color = "black") +
  theme(
    axis.line        = element_line(color = "black"),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#-------------------------
# 5. volcano_top10_unlab (no labels)
#-------------------------

volcano_top10_unlab <- ggplot(merged_cleaned, aes(x = log2FC, y = neg_log10_Pval)) +
  geom_point(
    aes(fill = significance),
    shape = 21,
    stroke = 0,
    size = 3,
    alpha = 0.8
  ) +
  geom_point(
    data = merged_cleaned %>% filter(significance != "Not Significant"),
    aes(fill = significance),
    shape = 21,
    color = "black",
    stroke = 0.8,
    size = 3,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c(
    "High CD4"       = "red",
    "Low CD4"        = "blue",
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(
    x    = "Log2 Fold Change",
    y    = "-Log10 p-value",
    fill = "Enriched in"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1),      linetype = "dashed", color = "black") +
  theme(
    axis.line        = element_line(color = "black"),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#save
ggsave("/path/to/plots/volcano_top10.svg",      plot = volcano_top10,      width = 10, height = 6)
ggsave("/path/to/plots/volcano_top10_unlab.svg", plot = volcano_top10_unlab, width = 10, height = 6)

