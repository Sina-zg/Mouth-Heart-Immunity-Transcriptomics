#----------Complex Heatmap--------
#----Libraries----
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

#----Sample Direction----
sample_info <- data.frame(
  sample_id = c("CV5","CV109","CV12","CV110","CV39","CV87","CV103","CV106"),
  group     = c("Progressor","Progressor","Progressor","Progressor",
                "Survivor","Survivor","Survivor","Survivor"),
  stringsAsFactors = FALSE
)

#----Load & QC of Samples----
seurat_list <- lapply(sample_info$sample_id, function(id) {
  mtx <- Read10X(paste0("/Users/sina_zg/Desktop/Sina_Air/Matrix/sample_filtered_feature_bc_matrix_", id))
  obj <- CreateSeuratObject(counts = mtx, project = id)
  obj$sample_id <- id
  obj$group <- sample_info$group[sample_info$sample_id == id]
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")
  obj <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 6500 & percent.mt < 10)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj
})
names(seurat_list) <- sample_info$sample_id

#----Merge counts----
counts_list <- list(); meta_list <- list()
for (id in names(seurat_list)) {
  s <- seurat_list[[id]]
  counts <- GetAssayData(s, layer = "counts")
  colnames(counts) <- paste0(id, "_", colnames(counts))
  meta <- s@meta.data
  rownames(meta) <- paste0(id, "_", rownames(meta))
  counts_list[[id]] <- counts
  meta_list[[id]]   <- meta
}
merged_counts <- do.call(cbind, counts_list)
merged_meta   <- do.call(rbind,  meta_list)

seurat_combined_clean <- CreateSeuratObject(counts = merged_counts, meta.data = merged_meta)
DefaultAssay(seurat_combined_clean) <- "RNA"
seurat_combined_clean <- NormalizeData(seurat_combined_clean, verbose = FALSE)

#----TCR: clone size----

#----Load TCR Samples----
tcr_list <- lapply(sample_info$sample_id, function(id) {
  tcr <- read.csv(paste0("/Users/sina_zg/Desktop/Sina_Air/Sequencing Analysis/filtered_contig_annotations_", id, ".csv"))
  tcr <- tcr[tcr$productive == "true", ]
  tcr$barcode <- paste0(id, "_", tcr$barcode)
  tcr$CTaa <- tcr$cdr3
  tcr[!duplicated(tcr$barcode), ]
})
tcr_all <- do.call(rbind, tcr_list)
tcr_all <- tcr_all[!duplicated(tcr_all$barcode), ]
tcr_all$CTaa <- tcr_all$cdr3

clone_sizes <- tcr_all %>%
  group_by(CTaa) %>%
  summarise(clone_size = n(), .groups = "drop")

tcr_meta <- tcr_all %>%
  left_join(clone_sizes, by = "CTaa") %>%
  mutate(
    tcr_freq = clone_size,
    clone_category = case_when(
      is.na(clone_size) ~ "Singleton",
      clone_size == 1   ~ "Singleton",
      clone_size <= 10   ~ "Small",
      clone_size <= 20  ~ "Medium",
      clone_size > 20   ~ "Large"
    )
  ) %>%
  distinct(barcode, .keep_all = TRUE) %>%
  filter(barcode %in% colnames(seurat_combined_clean)) %>%
  select(barcode, clone_category, tcr_freq)

rownames(tcr_meta) <- tcr_meta$barcode
tcr_meta <- as.data.frame(tcr_meta[, c("clone_category","tcr_freq")])
seurat_combined_clean <- AddMetaData(seurat_combined_clean, metadata = tcr_meta)

#----fill NAs----
seurat_combined_clean$clone_category[is.na(seurat_combined_clean$clone_category)] <- "Singleton"
seurat_combined_clean$tcr_freq[is.na(seurat_combined_clean$tcr_freq)] <- 0

#----Remove no CD4----
seurat_combined_clean <- subset(seurat_combined_clean, subset = CD3D > 0 & CD4 > 0)

#----Genes of interest----
genes_of_interest <- c(
  "SELL","CCR7","IL7R","CD28","XBP1","CD69","FAS","TIGIT","IKZF2","KLRB1","TOX","GZMM","ID2","KLRG1","TBX21",
  "PRF1","GZMA","CCL5","GNLY","GZMB","GZMH","NKG7","CX3CR1","PLEK","CCL4","SLAMF7",
  "EOMES","ZNF683","GZMK","TNF","IFNG","MKI67","FASLG","CD244",
  "PDCD1","TOX2","LAG3","TRDC","TRDV1","TRGV9","CD8A"
)


genes_present <- genes_of_interest[genes_of_interest %in% rownames(seurat_combined_clean)]
missing <- setdiff(genes_of_interest, genes_present)
if (length(missing) > 0) message("Missing genes (not in RNA assay): ", paste(missing, collapse = ", "))

#----Scale & matrix-----
seurat_combined_clean <- ScaleData(seurat_combined_clean, features = genes_present, verbose = FALSE)
expr_matrix_all <- GetAssayData(seurat_combined_clean, layer = "scale.data")[genes_present, , drop = FALSE]

#----Metadata for heatmap----
seurat_combined_clean$sample_id <- sub("_.*", "", colnames(seurat_combined_clean))
group_map <- setNames(sample_info$group, sample_info$sample_id)
seurat_combined_clean$group <- unname(group_map[seurat_combined_clean$sample_id])

meta_all <- seurat_combined_clean@meta.data[colnames(expr_matrix_all), ] %>%
  transmute(
    Group           = group,
    Patient         = sample_id,
    CloneSize       = clone_category,
    `Mito %`        = percent.mt,
    `Read count`    = nCount_RNA,
    `# of features` = nFeature_RNA,
    `TCR frequency` = tcr_freq
  )

#----Column order----

clone_order <- c("Large","Medium","Small","Singleton")
meta_all$CloneSize <- factor(meta_all$CloneSize, levels = clone_order)

sorted_cols <- meta_all %>%
  mutate(cell = rownames(.)) %>%
  arrange(Group, CloneSize, desc(`TCR frequency`), Patient) %>%  # tie-breakers optional
  pull(cell)

expr_matrix_sorted <- expr_matrix_all[, sorted_cols, drop = FALSE]
meta_sorted        <- meta_all[sorted_cols, , drop = FALSE]

stopifnot(identical(colnames(expr_matrix_sorted), rownames(meta_sorted)))

#----Plot Color Drafting----

patient_levels <- c("CV5","CV12","CV109","CV110","CV39","CV87","CV103","CV106")
meta_sorted$Patient <- factor(meta_sorted$Patient, levels = patient_levels)
patient_cols <- setNames(
  c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"),
  patient_levels
)

max_freq <- max(meta_sorted$`TCR frequency`, na.rm = TRUE)
if (!is.finite(max_freq) || max_freq <= 0) max_freq <- 1

ann_colors <- list(
  Group           = c("Progressor" = "firebrick", "Survivor" = "steelblue"),
  Patient         = patient_cols,
  CloneSize       = c("Singleton" = "gray", "Small" = "orange", "Medium" = "tomato", "Large" = "darkred"),
  `Mito %`        = colorRamp2(c(0, 10),    c("black", "yellow")),
  `Read count`    = colorRamp2(c(0, 15000), c("black", "green")),
  `# of features` = colorRamp2(c(0, 3000),  c("black", "lightblue")),
  `TCR frequency` = colorRamp2(c(0, max_freq), c("gray", "navy"))  
)

#----Annotation---
top_anno_full <- HeatmapAnnotation(
  df = meta_sorted,
  col = ann_colors,
  annotation_height = unit(rep(4, ncol(meta_sorted)), "mm"),
  show_legend = TRUE,
  annotation_name_side = "left"
)

#----Heatmap----
ht_full <- Heatmap(
  expr_matrix_sorted,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow")),
  top_annotation = top_anno_full,
  show_column_names = FALSE,
  cluster_columns = FALSE,           
  cluster_rows = FALSE,              
  row_order = rownames(expr_matrix_sorted),  
  column_split = meta_sorted$Group,  
  row_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(
    title = "Scaled Expression",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  ),
  use_raster = TRUE
)

#----Export PNG/PDF----
png("/Users/sina_zg/Desktop/CD4_CTL_Heatmap_CellStyle_FULL_CloneBlocks_NoSubpanels.png",
    width = 4800, height = 2600, res = 300)
draw(ht_full, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

pdf("/Users/sina_zg/Desktop/CD4_CTL_Heatmap_CellStyle_FULL_CloneBlocks_NoSubpanels.pdf",
    width = 16, height = 9)
draw(ht_full, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


