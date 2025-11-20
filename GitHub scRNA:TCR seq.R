
cran <- c("Seurat","dplyr","ggplot2","Matrix","purrr","glmGamPoi")
for (p in cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
bioc <- c("SingleCellExperiment","scDblFinder")
for (p in bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(Matrix); library(purrr)
  library(SingleCellExperiment); library(scDblFinder); library(glmGamPoi)
  library(future)
})
plan(sequential)
options(future.globals.maxSize = 20 * 1024^3)

set.seed(1234)

# ---- Samples & paths ----
samples <- c("CV5","CV12","CV109","CV110","CV39","CV87","CV103","CV106",
             "PD1","PD2","PD3","PD4","PD5","PD6","PD7","PD8")
group_map <- c(
  CV5="NICM_Progressor", CV12="NICM_Progressor", CV109="NICM_Progressor", CV110="NICM_Progressor",
  CV39="NICM_Survivor",  CV87="NICM_Survivor",  CV103="NICM_Survivor",  CV106="NICM_Survivor",
  PD1="PD_HighCTL", PD2="PD_HighCTL", PD3="PD_HighCTL", PD4="PD_HighCTL",
  PD5="PD_LowCTL",  PD6="PD_LowCTL",  PD7="PD_LowCTL",  PD8="PD_LowCTL"
)
rna_base <- "/Path to Samples/"
rna_path <- function(s) file.path(rna_base, paste0("sample_filtered_feature_bc_matrix_", s))
stopifnot(all(dir.exists(vapply(samples, rna_path, character(1)))))

# ---- helpers ----
mad_hi <- function(x, k=4) { m <- median(x, na.rm=TRUE); m + k*mad(x, na.rm=TRUE) }

#----QC plot----
plot_qc <- function(seu, s, outdir="qc_pngs"){
  if (!dir.exists(outdir)) dir.create(outdir, FALSE)
  p1 <- VlnPlot(seu, features=c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3, pt.size=0) +
    ggtitle(paste(s, "— violin"))
  p2 <- FeatureScatter(seu, feature1="nCount_RNA", feature2="nFeature_RNA") + ggtitle("nCount vs nFeature")
  p3 <- FeatureScatter(seu, feature1="percent.mt", feature2="nFeature_RNA") + ggtitle("%MT vs nFeature")
  png(file.path(outdir, paste0(s, "_QC.png")), width=1800, height=600, res=150)
  print(p1); print(p2); print(p3)
  dev.off()
}

# ----loading----
objs <- list(); qc_rows <- list()
for (s in samples) {
  mtx <- Read10X(rna_path(s))
  seu <- CreateSeuratObject(mtx, project=s, min.features=100, min.cells=3)
  seu$sample_id <- s
  seu$group     <- group_map[[s]]
  seu$cohort    <- if (startsWith(s,"CV")) "NICM" else "PD"
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern="^RP[SL]")
  seu[["percent.hb"]] <- PercentageFeatureSet(seu, pattern="^HB[AB]")
  
  raw <- ncol(seu)
  
  hard_keep <- with(seu@meta.data,
                    (nFeature_RNA > 500) & (nFeature_RNA < 6500) & (percent.mt < 10))
  seu <- subset(seu, cells = rownames(seu@meta.data)[hard_keep])
  
  
  nf_hi <- mad_hi(seu$nFeature_RNA, k=4)
  nc_hi <- mad_hi(seu$nCount_RNA,   k=4)
  mt_hi <- min(10, mad_hi(seu$percent.mt, k=3)) # never above 10%, often stricter
  
  keep2 <- with(seu@meta.data,
                (nFeature_RNA <= nf_hi) & (nCount_RNA <= nc_hi) & (percent.mt <= mt_hi))
  seu <- subset(seu, cells = rownames(seu@meta.data)[keep2])
  
  plot_qc(seu, s)
  objs[[s]] <- seu
  qc_rows[[s]] <- data.frame(sample_id=s,
                             raw=raw,
                             after_hard=sum(hard_keep),
                             after_adaptive=ncol(seu),
                             mt_cap_used=mt_hi,
                             stringsAsFactors=FALSE)
}
qc_stage1 <- bind_rows(qc_rows)
write.csv(qc_stage1, "qc_stage1_counts.csv", row.names=FALSE)

# ----doublets----
singlets <- list(); dbl_rows <- list()
for (s in samples) {
  seu <- objs[[s]]
  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce, samples=rep(s, ncol(sce)), dbr=0.04, verbose=FALSE)
  seu$scDblFinder.class <- colData(sce)$scDblFinder.class
  seu <- subset(seu, cells = rownames(seu@meta.data)[seu$scDblFinder.class=="singlet"])
  singlets[[s]] <- seu
  dbl_rows[[s]] <- data.frame(sample_id=s, after_doublet=ncol(seu))
}
qc_stage2 <- left_join(qc_stage1, bind_rows(dbl_rows), by="sample_id")

#----Normalization & SCTransform----
sct_list <- lapply(singlets, function(seu){
  SCTransform(seu, assay="RNA", variable.features.n=3000,
              vars.to.regress="percent.mt", method="glmGamPoi", verbose=FALSE)
})
sct_rows <- data.frame(sample_id = names(sct_list),
                       after_SCT = vapply(sct_list, ncol, integer(1)))
qc_stage3 <- left_join(qc_stage2, sct_rows, by="sample_id")

#----Filtering out non CD4----
gate_cd3d_cd4_drop_cd8a <- function(seu) {
  DefaultAssay(seu) <- "SCT"
  X <- GetAssayData(seu, slot="data")
  g <- rownames(X)
  v <- function(gene) if (gene %in% g) X[gene,] else rep(NA_real_, ncol(seu))
  cd3d <- v("CD3D"); cd4 <- v("CD4"); cd8a <- v("CD8A")
  keep <- (cd3d > 1) & (cd4 > 1) & (is.na(cd8a) | cd8a <= 1)
  subset(seu, cells = colnames(seu)[keep])
}
gated_list <- lapply(sct_list, gate_cd3d_cd4_drop_cd8a)

gate_rows <- data.frame(sample_id = names(gated_list),
                        after_gate = vapply(gated_list, ncol, integer(1)))
qc_final <- qc_stage3 |>
  left_join(gate_rows, by="sample_id") |>
  mutate(removed_by_gate = after_SCT - after_gate) |>
  arrange(sample_id)
write.csv(qc_final, "qc_stage_all_counts.csv", row.names=FALSE)
print(qc_final)

# Save clean, per-sample objects for integration
saveRDS(gated_list, "qc_clean_per_sample_SCT_CD4only.rds")
message("QC completed. Wrote qc_stage1_counts.csv, qc_stage_all_counts.csv, qc_pngs/*.png, and qc_clean_per_sample_SCT_CD4only.rds")


suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })
library(future); plan(sequential); set.seed(1234)

#----Load samples after QC----
gated_list <- readRDS("Path to qc_clean_per_sample_SCT_CD4only.rds")

#----Keep anchors biology-driven----
gated_list <- lapply(gated_list, function(x){
  DefaultAssay(x) <- "SCT"
  v <- VariableFeatures(x)
  if (length(v)) VariableFeatures(x) <- v[!grepl("^(TR[ABDG]|IG[HKL])", v)]
  x
})

features <- SelectIntegrationFeatures(object.list = gated_list, nfeatures = 3000)
gated_list <- PrepSCTIntegration(object.list = gated_list, anchor.features = features, verbose = FALSE)
gated_list <- lapply(gated_list, function(su) RunPCA(su, features = features, npcs = 30, verbose = FALSE))

sizes   <- vapply(gated_list, ncol, integer(1))
ref_idx <- order(sizes, decreasing = TRUE)[1:2]  

anchors <- FindIntegrationAnchors(
  object.list = gated_list, normalization.method = "SCT",
  anchor.features = features, dims = 1:30, reduction = "rpca",
  reference = ref_idx, k.anchor = 3, verbose = TRUE
)

integrated <- IntegrateData(
  anchorset = anchors, normalization.method = "SCT",
  dims = 1:30, k.weight = 15, verbose = TRUE
)
DefaultAssay(integrated) <- "integrated"

#----UMAP----
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30, graph.name = "pca_snn", verbose = FALSE)
integrated <- FindClusters(integrated, graph.name = "pca_snn", resolution = 0.3, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30, n.neighbors = 30, min.dist = 0.3, verbose = FALSE)

#----sanity check----
DefaultAssay(integrated) <- "SCT"
if (!"percent.mt" %in% colnames(integrated@meta.data)) {
  integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^MT-")
}

res_col <- "pca_snn_res.0.3"
if (!res_col %in% colnames(integrated@meta.data)) {
  stop("Expected clustering column ", res_col, " not found. Re-run FindClusters with resolution = 0.3.")
}
Idents(integrated) <- res_col


my_colors <- c("#66ccff","#3C5488","purple","#00a651","red",
               "#ffbf00","#E377C2","#7F7F7F","#BCBD22","#17BECF","yellow","red")
nclu <- nlevels(Idents(integrated))
cols_use <- if (length(my_colors) >= nclu) my_colors[1:nclu] else grDevices::colorRampPalette(my_colors)(nclu)
theme_umap <- theme_minimal(base_size=16) +
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text=element_blank(),
        axis.ticks=element_blank(), panel.grid=element_blank())

print(
  DimPlot(integrated, reduction="umap", label=FALSE, repel=TRUE, pt.size=3, cols=cols_use, raster=TRUE) +
    ggtitle("UMAP — Integrated CD4 T cells") + theme_umap
)


if ("cohort" %in% colnames(integrated@meta.data)) {
  print(DimPlot(integrated, reduction="umap", split.by="cohort", ncol=2,
                label=FALSE, repel=TRUE, pt.size=2.5, cols=cols_use, raster=TRUE) +
          ggtitle("UMAP — split by Cohort") + theme_umap)
}
if ("group" %in% colnames(integrated@meta.data)) {
  print(DimPlot(integrated, reduction="umap", split.by="group", ncol=2,
                label=FALSE, repel=TRUE, pt.size= 3, cols=cols_use, raster=TRUE) +
          ggtitle("UMAP — split by Group") + theme_umap)
}


#----Find all Markers----
integrated <- PrepSCTFindMarkers(integrated)

markers_all <- FindAllMarkers(
  integrated, assay="SCT", only.pos=TRUE,
  min.pct=0.25, logfc.threshold=0.25, test.use="wilcox"
)

markers_all <- markers_all |>
  filter(!grepl("^(TR[ABDG]|IG[HKL])", gene),
         !grepl("^MT-", gene),
         !grepl("^RP[SL]", gene))

#----Top-20 & 50 per cluster----
top20 <- markers_all |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) |>
  ungroup()
top50_c3 <- FindMarkers(integrated, ident.1 = "3", assay="SCT",
                        only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.0, test.use="wilcox") |>
  tibble::rownames_to_column("gene") |>
  arrange(desc(avg_log2FC)) |>
  filter(!grepl("^(TR[ABDG]|IG[HKL])", gene),
         !grepl("^MT-", gene),
         !grepl("^RP[SL]", gene)) |>
  slice_head(n = 50)

#----Summaries----
cat("\nTop 20 per cluster (gene names):\n")
split(top20$gene, top20$cluster) |> lapply(\(x) cat(paste(x, collapse=", "), "\n"))

cat("\nCluster 3 — Top 50 (gene names):\n")
cat(paste(top50_c3$gene, collapse=", "), "\n")

#----Saving----
write.csv(markers_all, "markers_all_clusters_SCT.csv", row.names=FALSE)
write.csv(top20,      "markers_top20_per_cluster.csv", row.names=FALSE)
write.csv(top50_c3,   "cluster3_top50.csv", row.names=FALSE)

#----annotation----
res_col <- "pca_snn_res.0.3"
stopifnot(res_col %in% colnames(integrated@meta.data))
Idents(integrated) <- res_col

annot_map <- c(
  "0" = "CCR4+ memory",
  "1" = "Naive/CM",
  "2" = "Th1 EM",
  "3" = "Treg",
  "4" = "CD4 CTL",
  "5" = "Activated memory"
)

clu_vec <- as.character(integrated@meta.data[[res_col]])
new_lab <- plyr::mapvalues(clu_vec, from = names(annot_map), to = unname(annot_map), warn_missing = FALSE)


unmapped <- !(clu_vec %in% names(annot_map))
new_lab[unmapped] <- paste0("Cluster ", clu_vec[unmapped])

integrated$annotation <- factor(new_lab,
                                levels = c("Naive/CM","CCR4+ memory","Activated memory","Th1 EM","Treg","CD4 CTL",
                                           sort(unique(new_lab[grepl("^Cluster ", new_lab)])))
)
Idents(integrated) <- "annotation"


pal <- c(
  "CCR4+ memory"   = "#66ccff",
  "Naive/CM"       = "#3C5488",
  "Th1 EM"         = "purple",
  "Treg"           = "#00a651",
  "CD4 CTL"        = "red",
  "Activated memory" = "#ffbf00"
)
extra_labs <- setdiff(levels(integrated$annotation), names(pal))
if (length(extra_labs)) {
  extra_cols <- grDevices::colorRampPalette(c("#7F7F7F","#BCBD22","#17BECF","yellow","red"))(length(extra_labs))
  names(extra_cols) <- extra_labs
  pal <- c(pal, extra_cols)
}

theme_umap <- theme_minimal(base_size=16) +
  theme(plot.title=element_text(face="bold",hjust=0.5), axis.text=element_blank(),
        axis.ticks=element_blank(), panel.grid=element_blank())

#----Annotated UMAP----
print(
  DimPlot(integrated, reduction="umap", label=TRUE, repel=TRUE, pt.size=3,
          cols = pal, raster=TRUE) +
    ggtitle("Annotated UMAP — res 0.25") + theme_umap
)


if ("cohort" %in% colnames(integrated@meta.data)) {
  print(DimPlot(integrated, reduction="umap", split.by="cohort", ncol=2,
                label=FALSE, repel=TRUE, pt.size=2.5, cols=pal, raster=TRUE) +
          ggtitle("Annotated UMAP — by Cohort") + theme_umap)
}
if ("group" %in% colnames(integrated@meta.data)) {
  print(DimPlot(integrated, reduction="umap", split.by="group", ncol=2,
                label=FALSE, repel=TRUE, pt.size=2.5, cols=pal, raster=TRUE) +
          ggtitle("Annotated UMAP — by Group") + theme_umap)
}

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

#----PD only----
md <- integrated@meta.data %>%
  dplyr::filter(cohort == "PD", group %in% c("PD_HighCTL","PD_LowCTL"))

#----Reordering----
order_levels <- c("Naive/CM","Th1 EM","CCR4+ memory","Treg","Activated memory","CD4 CTL")
md$annotation <- factor(md$annotation, levels = order_levels)

#----Stacked Bar-plot proportions----
p <- ggplot(md, aes(x = annotation, fill = group)) +
  geom_bar(position = "fill", width = 0.9) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01), limits = c(0,1), expand = c(0,0)) +
  scale_fill_manual(values = c("PD_HighCTL" = "red", "PD_LowCTL" = "#0b1d6b"),
                    labels = c("PD_HighCTL" = "High CTLs", "PD_LowCTL" = "Low CTLs"),
                    name = NULL) +
  labs(title = "Distribution of CD4 T-cells", x = NULL, y = "Proportion") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "right")
print(p)
# If p is your existing plot:
p <- p + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())

print(p)


ggsave(
  filename = "PD_High_vs_LowCTL_CD4pop_proportions.png",
  plot     = p,
  width    = 5, height = 5, units = "in",
  dpi      = 600,
  bg       = "white"
)

ggsave(
  filename = "PD_High_vs_LowCTL_CD4pop_proportions.svg",
  plot     = p,
  width    = 5, height = 5, units = "in",
  device   = "svg",
  bg       = "white"
)


#----NICM only----
md <- integrated@meta.data %>%
  dplyr::filter(cohort == "NICM", group %in% c("NICM_Progressor","NICM_Survivor"))

# Force left→right order for x axis
order_levels <- c("Naive/CM","Th1 EM","CCR4+ memory","Treg","Activated memory","CD4 CTL")
md$annotation <- factor(md$annotation, levels = order_levels)
md$group      <- factor(md$group, levels = c("NICM_Progressor","NICM_Survivor"))

tab_nicm <- as.data.frame(with(md, table(annotation, group))) %>%
  dplyr::group_by(annotation) %>%
  dplyr::mutate(prop = Freq / sum(Freq)) %>%
  dplyr::ungroup()

#----Stack Bar plotting----
p_nicm <- ggplot(tab_nicm, aes(x = annotation, y = prop, fill = group)) +
  geom_col(width = 0.9) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01),
                     limits = c(0,1), expand = c(0,0)) +
  scale_fill_manual(values = c("NICM_Progressor"="red","NICM_Survivor"="#0b1d6b"),
                    labels = c("NICM_Progressor"="Progressor","NICM_Survivor"="Survivor"),
                    name = NULL) +
  labs(title = "Distribution of CD4 T-cells (NICM)", x = NULL, y = "Proportion") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

print(p_nicm)

# Exports
ggsave("NICM_Progressor_vs_Survivor_CD4pop_proportions.png",
       p_nicm, width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("NICM_Progressor_vs_Survivor_CD4pop_proportions.svg",
       p_nicm, width = 5, height = 5, units = "in", device = "svg", bg = "white")


#----Subset CD4 CTL from Progressor vs PD HighCTL----
cd4ctl <- subset(integrated, subset = annotation == "CD4 CTL" & group %in% c("NICM_Progressor","PD_HighCTL"))
cd4ctl$plot_group <- dplyr::recode(cd4ctl$group, "NICM_Progressor"="Progressor", "PD_HighCTL"="PD")
cd4ctl$plot_group <- factor(cd4ctl$plot_group, levels = c("Progressor","PD"))

## 2) Gene list (keep only those present)
gene_order <- c("CD28","ZEB2","TOX","LAG3","PDCD1","SELL","CD38","GZMA","IL2RA",
                "CTLA4","CD52","CD69","GZMK","GNLY","GZMH","GZMM","NKG7","EOMES","IFNG","PRF1","GZMB")
assays <- tryCatch(Seurat::Assays(cd4ctl), error = function(e) names(cd4ctl@assays))
Seurat::DefaultAssay(cd4ctl) <- if ("SCT" %in% assays) "SCT" else if ("RNA" %in% assays) "RNA" else assays[[1]]

present <- intersect(gene_order, rownames(cd4ctl))
stopifnot(length(present) >= 5)

#----Group means and expression----
avg <- AverageExpression(cd4ctl, group.by = "plot_group", assays = DefaultAssay(cd4ctl), slot = "data")[[DefaultAssay(cd4ctl)]]
avg <- avg[present, , drop = FALSE]

expr <- GetAssayData(cd4ctl, slot = "data")[present, , drop = FALSE]
groups <- levels(cd4ctl$plot_group)

pct <- sapply(groups, function(g) {
  cells_g <- colnames(cd4ctl)[cd4ctl$plot_group == g]
  cells_g <- intersect(cells_g, colnames(expr))
  if (length(cells_g) == 0) return(rep(NA_real_, nrow(expr)))
  rowMeans(expr[, cells_g, drop = FALSE] > 0) * 100
})
colnames(pct) <- groups
pct <- pct[present, , drop = FALSE]

#----Row-zscore means----
zs <- t(scale(t(avg)))
zs[!is.finite(zs)] <- 0
zs <- pmax(pmin(zs, 3), -3)

df_z  <- as.data.frame(zs)  |>
  tibble::rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "group", values_to = "z")

df_pc <- as.data.frame(pct) |>
  tibble::rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "group", values_to = "pct")

df <- left_join(df_z, df_pc, by = c("gene","group")) |>
  mutate(group = factor(group, levels = c("Progressor","PD")),
         gene  = factor(gene,  levels = rev(present)))   # reverse so last gene in list appears at top

#----Plotting----
p_dot <- ggplot(df, aes(x = group, y = gene)) +
  geom_point(aes(size = pct, color = z)) +
  scale_size(range = c(2.5, 8), breaks = c(25, 50, 75, 100), name = "Percent Expressed") +
  scale_color_gradientn(colors = c("grey90", "#879BEA", "navy"),
                        limits = c(-3, 3), name = "Average Expression (z)") +
  labs(x = NULL, y = NULL,
       title = "CD4⁺ CTL transcriptional profile",
       subtitle = "Progressor vs PD High CTL (row-zscored per gene)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        legend.position = "right")

print(p_dot)

#----Pearson r similarity----
r <- cor(avg[, "Progressor"], avg[, "PD"], use = "pairwise.complete.obs")
cat(sprintf("\nPearson r across %d genes = %.3f\n", nrow(avg), r))


ggsave("DotPlot_CD4CTL_Progressor_vs_PD_rowZ.png", p_dot, width = 5, height = 7,
       units = "in", dpi = 300, bg = "white")
ggsave("DotPlot_CD4CTL_Progressor_vs_PD_rowZ.svg", p_dot, width = 5, height = 7.0,
       units = "in", device = "svg", bg = "white")




#--------------TCR Analysis-------------------
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(rlang)
})

#----Loading TCR Contigs----
tcr_base  <- "Path to Samples TCRs"
samples   <- c("CV5","CV12","CV109","CV110","CV39","CV87","CV103","CV106",
               "PD1","PD5","PD6","PD7","PD2","PD3","PD4","PD8")
group_map <- c(
  CV5="NICM_Progressor", CV12="NICM_Progressor", CV109="NICM_Progressor", CV110="NICM_Progressor",
  CV39="NICM_Survivor",  CV87="NICM_Survivor",  CV103="NICM_Survivor",  CV106="NICM_Survivor",
  PD1="PD_HighCTL", PD5="PD_HighCTL", PD6="PD_HighCTL", PD7="PD_HighCTL",
  PD2="PD_LowCTL",  PD3="PD_LowCTL",  PD4="PD_LowCTL",  PD8="PD_LowCTL"
)

#-----Barcode Matching-----
extract_10x <- function(x) {
  m <- regexpr("([ACGTN]{14,20}-[0-9]+)$", x, perl = TRUE)
  ifelse(m > 0, regmatches(x, m), NA_character_)
}

old <- colnames(integrated)

#----Audit Barcode Matching----
if (all(grepl("_[0-9]+$", old)) && all(!is.na(extract_10x(sub("_[0-9]+$", "", old))))) {
  message("Cellnames detected as BARCODE-1_<suffix>; stripping suffix.")
  raw_bc <- sub("_[0-9]+$", "", old)
} else if (all(grepl("_", old))) {
  raw_bc <- sub("^.*?_", "", old)
  bad <- is.na(extract_10x(raw_bc))
  if (any(bad)) raw_bc[bad] <- sub("^.*_", "", old[bad])  
} else {
  raw_bc <- old
}

ok_bc <- !is.na(extract_10x(raw_bc))
if (!all(ok_bc)) {
  ex <- head(old[!ok_bc], 10)
  stop("Could not parse 10x barcodes from RNA colnames for ", sum(!ok_bc),
       " cells.\nExamples: ", paste(ex, collapse=" | "),
       "\nAdjust the parsing rules and re-run.")
}

#----Derive sample_id----
meta <- integrated@meta.data
sample_id <- rep(NA_character_, length(old))

has_prefix <- grepl("_", old)
prefix <- ifelse(has_prefix, sub("_.*$", "", old), NA_character_)
keep_prefix <- has_prefix & !is.na(prefix) & prefix %in% samples
sample_id[keep_prefix] <- prefix[keep_prefix]

if ("orig.ident" %in% names(meta)) {
  oi <- as.character(meta$orig.ident)
  use_oi <- !is.na(oi) & nzchar(oi)
  sample_id[use_oi] <- oi[use_oi]
} else if ("sample" %in% names(meta)) {
  sm <- as.character(meta$sample)
  use_sm <- !is.na(sm) & nzchar(sm)
  sample_id[use_sm] <- sm[use_sm]
}

#----Read Contig----
tcr_files <- setNames(file.path(tcr_base, paste0("filtered_contig_annotations_", samples, ".csv")), samples)
tcr_files <- tcr_files[file.exists(tcr_files)]
if (!length(tcr_files)) stop("No contig files found in: ", tcr_base)

contigs_list <- lapply(tcr_files, function(f) read.csv(f, stringsAsFactors = FALSE))
names(contigs_list) <- names(tcr_files)

prep_contigs <- function(df, sid) {
  nm <- tolower(names(df)); names(df) <- nm
  if (!"reads" %in% nm && "read_count" %in% nm) df$reads <- df$read_count
  if (!"umi_count" %in% nm && "umis" %in% nm)   df$umi_count <- df$umis
  if (!"is_cell" %in% nm && "iscell" %in% nm)   df$is_cell <- df$iscell
  if (!"cdr3_aa" %in% nm) {
    if ("cdr3" %in% nm) df$cdr3_aa <- df$cdr3 else if ("amino_acid" %in% nm) df$cdr3_aa <- df$amino_acid
  }
  if (!"reads" %in% names(df))     df$reads <- 0L
  if (!"umi_count" %in% names(df)) df$umi_count <- 0L
  
  req <- c("barcode","chain","cdr3_aa","productive")
  miss <- setdiff(req, names(df))
  if (length(miss)) stop("Sample ", sid, " missing columns: ", paste(miss, collapse=", "),
                         "\nAvailable: ", paste(names(df), collapse=", "))
  
  df2 <- df %>%
    dplyr::filter(is.na(is_cell) | is_cell %in% c(TRUE,"true","True","TRUE")) %>%
    dplyr::filter(productive %in% c(TRUE,"true","True","TRUE")) %>%
    { if ("high_confidence" %in% names(.)) dplyr::filter(., high_confidence %in% c(TRUE,"true","True","TRUE")) else . } %>%
    { if ("full_length" %in% names(.))    dplyr::filter(., full_length    %in% c(TRUE,"true","True","TRUE")) else . } %>%
    dplyr::filter(!is.na(cdr3_aa), cdr3_aa != "") %>%
    dplyr::mutate(chain = toupper(chain)) %>%
    dplyr::filter(chain %in% c("TRA","TRB")) %>%
    dplyr::group_by(barcode, chain) %>%
    dplyr::arrange(dplyr::desc(umi_count), dplyr::desc(reads), .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(barcode, chain, cdr3_aa)
  
  wide <- df2 %>%
    dplyr::mutate(chain = ifelse(chain=="TRA","TRA_aa","TRB_aa")) %>%
    tidyr::pivot_wider(names_from = chain, values_from = cdr3_aa) %>%
    dplyr::filter(!is.na(TRA_aa), !is.na(TRB_aa)) %>%
    dplyr::mutate(sample_id = sid,
                  sample_barcode = paste0(sid, "_", barcode),
                  CTaa = paste0(TRA_aa, "|", TRB_aa)) %>%
    dplyr::select(sample_id, barcode, sample_barcode, TRA_aa, TRB_aa, CTaa)
  wide
}

clono_tbl <- purrr::imap_dfr(contigs_list, prep_contigs)

#----Fill sample_id----
if (any(is.na(sample_id))) {
  uniq_map <- clono_tbl %>%
    distinct(barcode, sample_id) %>%
    count(barcode, name = "n") %>%
    inner_join(clono_tbl %>% distinct(barcode, sample_id), by = "barcode") %>%
    filter(n == 1) %>%
    select(barcode, sample_id)
  
  map_vec <- setNames(uniq_map$sample_id, uniq_map$barcode)
  filled <- unname(map_vec[raw_bc])
  fill_idx <- which(is.na(sample_id) & !is.na(filled))
  if (length(fill_idx)) sample_id[fill_idx] <- filled[fill_idx]
}

if (anyNA(sample_id)) {
  ex <- head(old[is.na(sample_id)], 10)
  stop("Could not determine sample_id for ", sum(is.na(sample_id)),
       " cells.\nExamples: ", paste(ex, collapse=" | "),
       "\ensure 'orig.ident' is present ")
}

canon <- paste0(sample_id, "_", raw_bc)

#----Audit Duplicity----
du <- canon[duplicated(canon)]
if (length(du)) {
  where <- which(canon %in% head(unique(du), 10))
  diag_df <- data.frame(old = old[where],
                        sample_id = sample_id[where],
                        raw_bc = raw_bc[where],
                        canon = canon[where])
  print(diag_df, row.names = FALSE)
  stop("Canonical names not unique")
}

#----Apply rename----
if (!all(colnames(integrated) == canon)) {
  integrated <- Seurat::RenameCells(integrated, new.names = canon)
}
message("RNA colnames normalized to 'sample_BARCODE-1'.")

#----Add sample/group/cohort----
integrated$sample_id <- sample_id
integrated$group  <- unname(group_map[integrated$sample_id])
integrated$cohort <- ifelse(grepl("^CV", integrated$sample_id), "NICM",
                            ifelse(grepl("^PD", integrated$sample_id), "PD", NA))
mm <- setdiff(unique(integrated$sample_id), names(group_map))
if (length(mm)) message("sample_id not in group_map: ", paste(mm, collapse=", "))

#-------JOIN CTaa to Seurat Object------

clono_tbl <- clono_tbl %>% filter(sample_barcode %in% colnames(integrated))

#----Clone sizes----
clono_sizes <- clono_tbl %>% dplyr::count(CTaa)
names(clono_sizes)[names(clono_sizes) == "n"] <- "clone_size"
clono_tbl <- clono_tbl %>% left_join(clono_sizes, by = "CTaa")

#----Expansion bins----
bin_expansion <- function(n){
  if (is.na(n) || n == 0) return("NA")
  if (n == 1) return("Singleton (1)")
  if (n <= 10) return("Small (2-10)")
  if (n <= 20) return("Medium (11-20)")
  "Large (>20)"
}
new_levels <- c("NA","Singleton (1)","Small (2-10)","Medium (11-20)","Large (>20)")
clono_tbl$expansion_category <- factor(vapply(clono_tbl$clone_size, bin_expansion, character(1)),
                                       levels = new_levels)

#----Attach to Seurat metadata----
meta <- integrated@meta.data
for (nm in c("CTaa","TRA_aa","TRB_aa","clone_size","expansion_category")) {
  if (!nm %in% names(meta)) meta[[nm]] <- NA
}
m <- match(clono_tbl$sample_barcode, rownames(meta))
meta$CTaa[m]               <- clono_tbl$CTaa
meta$TRA_aa[m]             <- clono_tbl$TRA_aa
meta$TRB_aa[m]             <- clono_tbl$TRB_aa
meta$clone_size[m]         <- clono_tbl$clone_size
meta$expansion_category[m] <- clono_tbl$expansion_category
integrated@meta.data <- meta

cat("Attached CTaa to", sum(!is.na(integrated$CTaa)),
    "cells out of", ncol(integrated), "total.\n")

#----Audit----
integ_meta <- integrated@meta.data
integ_meta$expansion_category <- factor(
  vapply(integ_meta$clone_size, bin_expansion, character(1)),
  levels = new_levels
)
integrated@meta.data <- integ_meta

#----Plotting----

pal_exp_all <- c("NA"="grey90","Singleton (1)"="#bfe2f5","Small (2-10)"="#ffa31a",
                 "Medium (11-20)"="#d62728","Large (>20)"="#8b1a1a")

#----Labeling----
if (!"group_simple" %in% colnames(integrated@meta.data)) {
  integrated$group_simple <- dplyr::recode(
    integrated$group,
    "NICM_Survivor"   = "Survivor",
    "NICM_Progressor" = "Progressor",
    "PD_HighCTL"      = "PD HighCTL",
    "PD_LowCTL"       = "PD LowCTL",
    .default = as.character(integrated$group)
  )
}

#----Manual Labeling per Cluster----
cluster_col <- if ("annotation" %in% colnames(integrated@meta.data)) "annotation" else "seurat_clusters"
cl_sym <- rlang::sym(cluster_col)
if (cluster_col == "annotation") {
  pref_levels <- c("Naive/CM","Treg","CCR4+ memory","Activated memory","EM","GZMK+","CD4 CTL")
  present <- intersect(pref_levels, unique(integrated$annotation))
  integrated$annotation <- factor(integrated$annotation, levels = c(present, setdiff(unique(integrated$annotation), present)))
}


align_palette <- function(df, aes_col, pal) {
  vals <- unique(as.character(df[[aes_col]]))
  vals <- vals[!is.na(vals)]
  pal[names(pal) %in% vals]
}


mk_df_group <- function(obj, cohort_keep, group_levels_order) {
  out <- obj@meta.data %>%
    filter(cohort == cohort_keep) %>%
    filter(!is.na(group_simple)) %>%
    filter(!is.na(expansion_category), expansion_category != "NA") %>%
    dplyr::count(group_simple, expansion_category, name = "count") %>%
    group_by(group_simple) %>% mutate(prop = count / sum(count)) %>% ungroup() %>%
    mutate(group_simple = factor(group_simple, levels = group_levels_order))
  out
}
mk_df_cluster <- function(obj, cohort_keep, group_levels_order) {
  df <- obj@meta.data %>%
    filter(cohort == cohort_keep) %>%
    filter(!is.na(group_simple)) %>%
    filter(!is.na(expansion_category), expansion_category != "NA") %>%
    dplyr::count(group_simple, !!cl_sym, expansion_category, name = "count") %>%
    group_by(group_simple, !!cl_sym) %>% mutate(prop = count / sum(count)) %>% ungroup()
  
  if (nrow(df) == 0) return(df)
  
  if (cluster_col == "annotation") {
    df[[cluster_col]] <- factor(df[[cluster_col]], levels = levels(integrated$annotation))
  } else {
    df[[cluster_col]] <- factor(df[[cluster_col]])
  }
  df$group_simple <- factor(df$group_simple, levels = group_levels_order)
  df
}
move_cd4ctl_last <- function(x, ref_levels=NULL) {
  levs <- if (!is.null(ref_levels)) ref_levels else levels(x)
  if (is.null(levs)) levs <- unique(as.character(x))
  if ("CD4 CTL" %in% levs) levs <- c(setdiff(levs, "CD4 CTL"), "CD4 CTL")
  factor(x, levels = levs)
}


# -----------Exporting---------
out_dir <- path.expand("Path to Save files")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_both <- function(plot, base, w=6, h=4.5, dpi=600) {
  png_f <- file.path(out_dir, paste0(base, ".png"))
  svg_f <- file.path(out_dir, paste0(base, ".svg"))
  ggsave(png_f, plot, width = w, height = h, dpi = dpi, bg = "white")
  if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
  ggsave(svg_f, plot, width = w, height = h, dpi = dpi, device = "svg", bg = "white")
  message("Saved: ", basename(png_f), " and ", basename(svg_f), " -> ", out_dir)
}
save_vec <- save_both

#----Rebuild group / group_simple / cohort----
meta <- integrated@meta.data

#----Audit Grouping----
if (!"group" %in% names(meta) || all(is.na(meta$group))) {
  if (!exists("group_map")) stop("group_map is missing")
  integrated$group <- unname(group_map[integrated$sample_id])
  meta <- integrated@meta.data
}

#----Labelling----
integrated$group_simple <- dplyr::recode(
  as.character(integrated$group),
  "NICM_Survivor"   = "Survivor",
  "NICM_Progressor" = "Progressor",
  "PD_HighCTL"      = "PD HighCTL",
  "PD_LowCTL"       = "PD LowCTL",
  .default = NA_character_
)

#----Ensure cohort----
if (!"cohort" %in% names(meta) || all(is.na(meta$cohort))) {
  integrated$cohort <- ifelse(grepl("^CV", integrated$sample_id), "NICM",
                              ifelse(grepl("^PD", integrated$sample_id), "PD", NA))
}

#----Sanity Check----
cat("\n--- Diagnostics (after rebuild) ---\n")
print(table(Cohort = integrated$cohort, useNA = "ifany"))
print(table(Group  = integrated$group,  useNA = "ifany"))
print(table(GroupSimple = integrated$group_simple, useNA = "ifany"))
cat("\nExpansion category summary:\n")
print(table(integrated$expansion_category, useNA = "ifany"))

#----Plotting----
pal_exp_all <- c("NA"="grey90","Singleton (1)"="#bfe2f5","Small (2-10)"="#ffa31a",
                 "Medium (11-20)"="#d62728","Large (>20)"="#8b1a1a")

align_palette <- function(df, aes_col, pal) {
  vals <- unique(as.character(df[[aes_col]]))
  vals <- vals[!is.na(vals)]
  pal[names(pal) %in% vals]
}

cluster_col <- if ("annotation" %in% colnames(integrated@meta.data)) "annotation" else "seurat_clusters"
cl_sym <- rlang::sym(cluster_col)
if (cluster_col == "annotation") {
  pref_levels <- c("Naive/CM","Treg","CCR4+ memory","Activated memory","EM","GZMK+","CD4 CTL")
  present <- intersect(pref_levels, unique(integrated$annotation))
  integrated$annotation <- factor(integrated$annotation,
                                  levels = c(present, setdiff(unique(integrated$annotation), present)))
}

#----With NA----
mk_df_group <- function(obj, cohort_keep, group_levels_order) {
  df <- obj@meta.data %>%
    dplyr::filter(cohort == cohort_keep) %>%
    dplyr::filter(!is.na(group_simple)) %>%
    dplyr::filter(!is.na(expansion_category)) %>%
    dplyr::count(group_simple, expansion_category, name = "count") %>%
    dplyr::group_by(group_simple) %>% dplyr::mutate(prop = count / sum(count)) %>% dplyr::ungroup() %>%
    dplyr::mutate(group_simple = factor(group_simple, levels = group_levels_order))
  df
}

mk_df_cluster <- function(obj, cohort_keep, group_levels_order) {
  df <- obj@meta.data %>%
    dplyr::filter(cohort == cohort_keep) %>%
    dplyr::filter(!is.na(group_simple)) %>%
    dplyr::filter(!is.na(expansion_category)) %>%
    dplyr::count(group_simple, !!cl_sym, expansion_category, name = "count") %>%
    dplyr::group_by(group_simple, !!cl_sym) %>% dplyr::mutate(prop = count / sum(count)) %>% dplyr::ungroup()
  if (nrow(df) == 0) return(df)
  if (cluster_col == "annotation") {
    df[[cluster_col]] <- factor(df[[cluster_col]], levels = levels(integrated$annotation))
  } else {
    df[[cluster_col]] <- factor(df[[cluster_col]])
  }
  df$group_simple <- factor(df$group_simple, levels = group_levels_order)
  df
}

move_cd4ctl_last <- function(x, ref_levels=NULL) {
  levs <- if (!is.null(ref_levels)) ref_levels else levels(x)
  if (is.null(levs)) levs <- unique(as.character(x))
  if ("CD4 CTL" %in% levs) levs <- c(setdiff(levs, "CD4 CTL"), "CD4 CTL")
  factor(x, levels = levs)
}

# ----------------- NICM plots -----------------
nicm_levels <- c("Survivor","Progressor")

#----NICM Barplot per Group NA removed----
df_group_nicm <- mk_df_group(integrated, "NICM", nicm_levels) %>%
  dplyr::filter(!is.na(expansion_category), as.character(expansion_category) != "NA")

if (nrow(df_group_nicm)) {
  pal_exp <- align_palette(df_group_nicm, "expansion_category", pal_exp_all)
  p_group_nicm <- ggplot(df_group_nicm, aes(x = group_simple, y = count, fill = expansion_category)) +
    geom_col(width = 0.9, position = "fill") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0)) +
    scale_fill_manual(values = pal_exp, name = "Clonal Expansion", na.translate = FALSE, drop = TRUE) +
    labs(title = "NICM — clonal expansion composition per group", x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face="bold", hjust=0.5), panel.grid = element_blank())
  print(p_group_nicm)
  save_both(p_group_nicm, "NICM_clonality_per_group", 4, 4.5)
}

#----NICM Barplot per Cluster no NA----
df_cluster_nicm <- mk_df_cluster(integrated, "NICM", nicm_levels) %>%
  dplyr::filter(!is.na(expansion_category), as.character(expansion_category) != "NA")

if (nrow(df_cluster_nicm)) {
  if (cluster_col == "annotation") {
    df_cluster_nicm[[cluster_col]] <- move_cd4ctl_last(df_cluster_nicm[[cluster_col]], levels(integrated$annotation))
  }
  pal_exp <- align_palette(df_cluster_nicm, "expansion_category", pal_exp_all)
  p_cluster_nicm <- ggplot(df_cluster_nicm, aes(x = !!cl_sym, y = count, fill = expansion_category)) +
    geom_col(width = 0.9, position = "fill") +
    facet_wrap(~ group_simple, nrow = 1, scales = "free_x") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0)) +
    scale_fill_manual(values = pal_exp, name = "Clonal Expansion", na.translate = FALSE, drop = TRUE) +
    labs(title = paste0("NICM — clonal expansion by ", cluster_col, " (Progressor on right)"),
         x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face="bold", hjust=0.5),
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
          panel.grid = element_blank(),
          strip.text = element_text(face="bold"))
  print(p_cluster_nicm)
  save_both(p_cluster_nicm, paste0("NICM_clonality_by_", cluster_col), 7, 4.8)
}

if (nrow(df_cluster_nicm)) {
  if (cluster_col == "annotation") {
    df_cluster_nicm[[cluster_col]] <- move_cd4ctl_last(df_cluster_nicm[[cluster_col]], levels(integrated$annotation))
  }
  pal_exp <- align_palette(df_cluster_nicm, "expansion_category", pal_exp_all)
  p_cluster_nicm <- ggplot(df_cluster_nicm, aes(x = !!cl_sym, y = prop, fill = expansion_category)) +
    geom_col(width = 0.9, position = "fill") +
    facet_wrap(~ group_simple, nrow = 1, scales = "free_x") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0)) +
    scale_fill_manual(values = pal_exp, name = "Clonal Expansion") +
    labs(title = paste0("NICM — clonal expansion by ", cluster_col, " (Progressor on right)"),
         x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face="bold", hjust=0.5),
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
          panel.grid = element_blank(),
          strip.text = element_text(face="bold"))
  print(p_cluster_nicm)
  save_both(p_cluster_nicm, paste0("NICM_clonality_by_", cluster_col), 7, 4.8)
} else {
  message("No NICM rows per-cluster plot.")
}

#----PD Barplot per Group no NA----
df_group_pd <- mk_df_group(integrated, "PD", pd_levels) %>%
  dplyr::filter(!is.na(expansion_category), as.character(expansion_category) != "NA")

if (nrow(df_group_pd)) {
  pal_exp <- align_palette(df_group_pd, "expansion_category", pal_exp_all)
  p_group_pd <- ggplot(df_group_pd, aes(x = group_simple, y = count, fill = expansion_category)) +
    geom_col(width = 0.9, position = "fill") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0)) +
    scale_fill_manual(values = pal_exp, name = "Clonal Expansion", na.translate = FALSE, drop = TRUE) +
    labs(title = "PD — clonal expansion composition per group", x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face="bold", hjust=0.5), panel.grid = element_blank())
  print(p_group_pd)
  save_both(p_group_pd, "PD_clonality_per_group", 4, 4.5)
}

#----PD Barplot per Cluster no NA----
df_cluster_pd <- mk_df_cluster(integrated, "PD", pd_levels) %>%
  dplyr::filter(!is.na(expansion_category), as.character(expansion_category) != "NA")

if (nrow(df_cluster_pd)) {
  if (cluster_col == "annotation") {
    df_cluster_pd[[cluster_col]] <- move_cd4ctl_last(df_cluster_pd[[cluster_col]], levels(integrated$annotation))
  }
  pal_exp <- align_palette(df_cluster_pd, "expansion_category", pal_exp_all)
  p_cluster_pd <- ggplot(df_cluster_pd, aes(x = !!cl_sym, y = count, fill = expansion_category)) +
    geom_col(width = 0.9, position = "fill") +
    facet_wrap(~ group_simple, nrow = 1, scales = "free_x") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0)) +
    scale_fill_manual(values = pal_exp, name = "Clonal Expansion", na.translate = FALSE, drop = TRUE) +
    labs(title = paste0("PD — clonal expansion by ", cluster_col, " (HighCTL on right)"),
         x = NULL, y = "Proportion") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face="bold", hjust=0.5),
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
          panel.grid = element_blank(),
          strip.text = element_text(face="bold"))
  print(p_cluster_pd)
  save_both(p_cluster_pd, paste0("PD_clonality_by_", cluster_col), 7, 4.8)
}


# ----------------- UMAP per group -----------------
pal_exp_points <- c("Singleton (1)"="#bfe2f5","Small (2-10)"="#ffa31a",
                    "Medium (11-20)"="#d62728","Large (>20)"="#8b1a1a")

um_all <- Seurat::Embeddings(integrated, "umap")
if (is.null(um_all) || ncol(um_all) < 2) stop("No UMAP embeddings found.")
xr <- range(um_all[,1]); yr <- range(um_all[,2])

make_umap_df <- function(obj) {
  emb <- Seurat::Embeddings(obj, reduction = "umap")[, 1:2, drop = FALSE]
  if (is.null(emb) || ncol(emb) < 2) stop("No UMAP embeddings found in this object.")
  colnames(emb) <- c("UMAP_x", "UMAP_y")
  
  meta <- obj@meta.data[, c("expansion_category","group_simple"), drop = FALSE]
  df <- cbind(as.data.frame(emb), meta)
  
  
  keep <- !is.na(df$expansion_category) & as.character(df$expansion_category) != "NA"
  df <- df[keep, , drop = FALSE]
  
  rownames(df) <- NULL
  df$expansion_category <- droplevels(df$expansion_category)
  df
}

plot_umap_pretty <- function(df, title,
                             pt_size = 1.5, alpha = 0.9, stroke = 0.1,
                             border_col = "black", xlim = NULL, ylim = NULL, pal = pal_exp_points) {
  pal2 <- align_palette_pts(df, pal)
  ggplot(df, aes(UMAP_x, UMAP_y)) +
    geom_point(aes(fill = expansion_category),
               shape = 21, size = pt_size, alpha = alpha,
               stroke = stroke, color = border_col) +
    scale_fill_manual(values = pal2, drop = TRUE, name = "Clonal Expansion") +
    coord_fixed(xlim = xlim, ylim = ylim, ratio = 1) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(plot.title = element_text(face="bold", hjust=0.5),
          legend.position = "right",
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

grps <- c("Survivor","Progressor","PD LowCTL","PD HighCTL")
grps <- grps[grps %in% unique(integrated$group_simple)]

for (g in grps) {
  obj_g <- subset(integrated, subset = group_simple == g)
  df_g  <- make_umap_df(obj_g)
  if (!nrow(df_g)) { message("No cells for UMAP: ", g); next }
  p <- plot_umap_pretty(df_g, paste0("UMAP — ", g, " (Clonal Expansion)"),
                        pt_size = 3.5, alpha = 0.9, stroke = 0.1,
                        border_col = "black", xlim = xr, ylim = yr)
  print(p)
  base <- paste0("UMAP_", gsub("[^A-Za-z0-9]+", "_", g), "_expansion")
  save_vec(p, base, 7, 7, 600)
}


