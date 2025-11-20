
#----Packages----

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(Matrix)
  library(matrixStats)
  library(harmony)
  library(dbscan)
  library(tibble)
  library(scales)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

theme_set(theme_bw(base_size = 11))
options(future.globals.maxSize = 20 * 1024^3)


#----Sample Directory----

samples <- tribble(
  ~sample_id, ~group,       ~path_to_outs,
  "P1",       "Progressor", "Path to sample directory",
  "P2",       "Progressor", "Path to sample directory",
  "P3",       "Progressor", "Path to sample directory",
  "P4",       "Progressor", "Path to sample directory",
  "P5",       "Progressor", "Path to sample directory",
  "S1",       "Survivor",   "Path to sample directory",
  "S2",       "Survivor",   "Path to sample directory",
  "S3",       "Survivor",   "Path to sample directory"
)


samples


#----Helpers----

resolve_xenium_outs <- function(p) {
  cands <- unique(c(
    normalizePath(p, mustWork = FALSE),
    normalizePath(file.path(p, "outs"), mustWork = FALSE)
  ))
  for (x in cands) {
    if (is.na(x)) next
    hits <- file.exists(file.path(x, "experiment.xenium")) ||
      file.exists(file.path(x, "cells.parquet")) ||
      file.exists(file.path(x, "cells.csv.gz"))
    if (hits) return(x)
  }
  stop(
    "Couldn't find Xenium 'outs' under: ", p,
    "\nTried: ", paste(cands, collapse = " | "),
    "\Ensuring files exist. 'experiment.xenium' and 'cells.parquet' or 'cells.csv.gz'.",
    call. = FALSE
  )
}

#----Load Xenium samples----
load_one_xenium <- function(path, sample_id, group) {
  message("Loading ", sample_id, " from base path: ", path)
  path <- resolve_xenium_outs(path)
  message(" -> Resolved outs/: ", path)
  
  obj <- LoadXenium(
    data.dir = path,
    fov = sample_id,
    molecule.coordinates = FALSE
  )
  
  obj <- subset(obj, subset = nCount_Xenium > 0)
  obj$sample_id <- sample_id
  obj$group     <- group
  colnames(obj) <- paste(sample_id, colnames(obj), sep = "_")
  obj
}

#----Gene name----
match_genes <- function(genes, obj) {
  rn <- rownames(obj)
  out <- genes[genes %in% rn]
  miss <- setdiff(genes, out)
  if (length(miss)) {
    alt <- gsub("_", "-", miss, fixed = TRUE)
    out <- unique(c(out, alt[alt %in% rn]))
  }
  unique(out)
}

#----Signature scoring and z-transform----
sig_score <- function(obj, genes) {
  g <- intersect(genes, rownames(obj))
  if (!length(g)) return(rep(0, ncol(obj)))
  m <- GetAssayData(obj, assay = "SCT", slot = "data")[g, , drop = FALSE]
  as.numeric(Matrix::colMeans(m))
}

z <- function(v) {
  v2 <- as.numeric(scale(v))
  v2[is.na(v2)] <- 0
  v2
}


#----Load & merge samples----

objs <- purrr::pmap(
  samples,
  function(sample_id, group, path_to_outs) {
    load_one_xenium(path_to_outs, sample_id, group)
  }
)

xen <- Reduce(function(a, b) merge(a, b), objs)
DefaultAssay(xen) <- "Xenium"


xen <- JoinLayers(xen)

#----Sanity----
if (!"sample_id" %in% colnames(xen@meta.data)) {
  xen$sample_id <- sub("^([^_]+)_.*$", "\\1", colnames(xen))
}
if (!"group" %in% colnames(xen@meta.data)) {
  map <- samples %>% distinct(sample_id, group)
  xen$group <- map$group[match(xen$sample_id, map$sample_id)]
}

table(xen$sample_id)
table(xen$group)


#----QC: control probes + nFeature----

ctrl_frac <- rep(0, ncol(xen))
if ("ControlProbe" %in% names(xen@assays)) {
  cp <- Matrix::colSums(GetAssayData(xen[["ControlProbe"]], slot = "counts"))
  xm <- Matrix::colSums(GetAssayData(xen[["Xenium"]],       slot = "counts"))
  ctrl_frac <- cp / pmax(xm, 1)
}
xen$percent_ctrlprobe <- ctrl_frac

q_nf <- xen$nFeature_Xenium
low_g <- as.numeric(quantile(q_nf, 0.02,  na.rm = TRUE))
hi_g  <- as.numeric(quantile(q_nf, 0.999, na.rm = TRUE))
max_ctrl <- 0.15

n_before <- ncol(xen)
xen <- subset(
  xen,
  subset = nFeature_Xenium >= low_g &
    nFeature_Xenium <= hi_g &
    percent_ctrlprobe <= max_ctrl
)
n_after <- ncol(xen)
cat("Cells before QC:", n_before, " after QC:", n_after,
    " (removed:", n_before - n_after, ")\n")


#----SCTransform----

xen <- SCTransform(xen, assay = "Xenium", verbose = FALSE)
DefaultAssay(xen) <- "SCT"


#----Gene panels----

genes_cardio <- c("ACTC1","TNNT2","TNNI3","MYH6","MYH7","MYL2","MYL7","RYR2","PLN","DES","NPPA","NPPB")
genes_endo   <- c("PECAM1","VWF","KDR","FLT1","TEK","CDH5","EMCN","ACKR1","KLF2")
genes_t      <- c("CD3D","CD3E","CD2","TRAC","CD4","CD8A","CD8B","IL7R","CCR7","TCF7")
genes_b      <- c("MS4A1","CD19","CD79A","CD79B","CD74")
genes_plasma <- c("MZB1","JCHAIN","XBP1","IGHG1","IGKC","SDC1")
genes_myelo  <- c("LST1","SPP1","CD68","ITGAM","CSF1R","FCGR3A","CX3CR1","CCR2")
genes_stress <- c("NPPA","NPPB","FOS","JUN","HSPB1","ATF3","DDIT3")

#----Fibroblast / fibrotic genes----
genes_fibro_core <- c(
  "COL1A1","COL1A2","COL3A1","COL6A1","COL6A2",
  "DCN","LUM","DPT","THY1","PDGFRA","FBLN1","FBLN2"
)
genes_myofibro   <- c("ACTA2","TAGLN","MYL9","TPM2","ITGA11","CNN1")
genes_fibrosis   <- c("POSTN","FN1","THBS2","MMP2","TIMP1","COL1A1","COL3A1")

all_feats   <- rownames(xen)
human_syms  <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
feats_human <- intersect(all_feats, human_syms)

rel_panel2 <- unique(c(
  genes_cardio, genes_endo,
  genes_t, genes_b, genes_plasma, genes_myelo, genes_stress,
  genes_fibro_core, genes_myofibro, genes_fibrosis
))

rel_ok2 <- intersect(match_genes(rel_panel2, xen), feats_human)

target_n <- 180
allowlist <- rel_ok2
if (length(allowlist) < target_n) {
  remaining <- setdiff(feats_human, allowlist)
  mat <- as.matrix(GetAssayData(xen, assay = "SCT", slot = "data")[remaining, , drop = FALSE])
  rv  <- matrixStats::rowVars(mat)
  add <- remaining[order(rv, decreasing = TRUE)]
  allowlist <- unique(c(allowlist, head(add, target_n - length(allowlist))))
}
cat("Allowlist :", length(allowlist), "\n")


#----Program category----

set_cardio <- match_genes(genes_cardio, xen)
set_endo   <- match_genes(genes_endo,   xen)
set_t      <- match_genes(genes_t,      xen)
set_b      <- match_genes(genes_b,      xen)
set_plasma <- match_genes(genes_plasma, xen)
set_myelo  <- match_genes(genes_myelo,  xen)
set_stress <- match_genes(genes_stress, xen)

set_fibro  <- match_genes(genes_fibro_core, xen)
set_myofi  <- match_genes(genes_myofibro,   xen)
set_fibrot <- match_genes(genes_fibrosis,   xen)

s_cardio <- sig_score(xen, set_cardio)
s_endo   <- sig_score(xen, set_endo)
s_imm    <- pmax(
  sig_score(xen, set_t),
  sig_score(xen, set_b),
  sig_score(xen, set_plasma),
  sig_score(xen, set_myelo),
  na.rm = TRUE
)
s_stress <- z(sig_score(xen, set_stress))

s_fibro   <- sig_score(xen, set_fibro)
s_myofi   <- sig_score(xen, set_myofi)
s_fibrot  <- sig_score(xen, set_fibrot)

s_cardio_z <- z(s_cardio)
s_myofi_z  <- z(s_myofi)
s_fibrot_z <- z(s_fibrot)

min_score   <- 0.05
thr_stressZ <- 1.0
thr_myofibZ <- 1.0
thr_fibrotZ <- 1.0


base_mat <- cbind(
  cardio = s_cardio,
  endo   = s_endo,
  immune = s_imm,
  fibro  = s_fibro
)
base_max <- apply(base_mat, 1, max, na.rm = TRUE)
base_ix  <- max.col(base_mat, ties.method = "first")
base_lab <- c("Cardiomyocyte","Endothelium","Immune","Fibroblast")[base_ix]
base_lab[base_max < min_score] <- "Other"

focus_category <- base_lab
is_fibro  <- base_lab == "Fibroblast"
is_cardio <- base_lab == "Cardiomyocyte"

# stressed cardiomyocytes
focus_category[is_cardio & s_stress >= thr_stressZ & s_cardio_z >= 0.5] <- "Stressed CM"

#----Fibroblast+Fibrotic----
focus_category[is_fibro & (s_myofi_z  >= thr_myofibZ)] <- "Fibroblast"
focus_category[is_fibro & (s_fibrot_z >= thr_fibrotZ)] <- "Fibroblast"

xen$focus_category <- as.character(focus_category)

#----Categories----
valid_levels <- c("Cardiomyocyte","Stressed CM","Endothelium","Immune","Fibroblast")
xen$focus_category[!(xen$focus_category %in% valid_levels)] <- "Other"

xen$focus_category <- factor(xen$focus_category,
                             levels = c(valid_levels, "Other"))

cat("Category counts :\n")
print(table(xen$focus_category))
cat("\nCategory by group:\n")
print(table(xen$focus_category, xen$group))


keep_levels <- c(
  "Cardiomyocyte","Stressed CM","Endothelium",
  "Immune",
  "Fibroblast","Myofibroblast","Fibrotic Fibroblast"
)

xen_focus <- subset(xen, cells = colnames(xen)[xen$focus_category %in% keep_levels])
DefaultAssay(xen_focus) <- "SCT"

feat_use <- intersect(allowlist, rownames(xen_focus))

set.seed(1337)
xen_focus <- RunPCA(xen_focus, features = feat_use, npcs = 30, verbose = FALSE)

by_var <- if (!"sample_id" %in% colnames(xen_focus@meta.data)) "group" else "sample_id"
message(sprintf("[Harmony] batching by: %s", by_var))

xen_focus <- harmony::RunHarmony(
  object           = xen_focus,
  group.by.vars    = by_var,
  reduction.use    = "pca",
  dims.use         = 1:30,
  assay.use        = DefaultAssay(xen_focus),
  project.dim      = FALSE,
  plot_convergence = FALSE,
  verbose          = FALSE
)

xen_focus <- RunUMAP(xen_focus, reduction = "harmony", dims = 1:30, verbose = FALSE)
xen_focus <- FindNeighbors(xen_focus, reduction = "harmony", dims = 1:30, verbose = FALSE)
xen_focus <- FindClusters(xen_focus, resolution = 0.4, verbose = FALSE)


emb <- Embeddings(xen_focus, "umap")
fc  <- xen_focus$focus_category

#----Low Confident island Removal----
keep <- rep(TRUE, nrow(emb))
min_comp <- 900   
eps_val  <- 0.7   
minPts   <- 20

for (lvl in unique(fc)) {
  idx <- which(fc == lvl)
  if (length(idx) < min_comp) { keep[idx] <- FALSE; next }
  cl  <- dbscan(emb[idx, , drop = FALSE], eps = eps_val, minPts = minPts)$cluster
  sizes <- table(cl)
  small <- names(sizes[sizes < min_comp])
  if (length(small)) keep[idx[cl %in% as.integer(small)]] <- FALSE
}

xen_focus <- subset(xen_focus, cells = Cells(xen_focus)[keep])

cat("Category counts :\n")
print(table(xen_focus$focus_category, xen_focus$group))


#----Harmony UMAP----

my_cols <- c(
  "Cardiomyocyte"       = "#F9B6C1",
  "Stressed CM"         = "#C2185B",
  "Endothelium"         = "#4CAF50",
  "Fibroblast"          = "#F4A582",
  "Myofibroblast"       = "#F4A582",
  "Fibrotic Fibroblast" = "#F4A582",
  "Immune"              = "#4A148C"
)

p_umap_split <- DimPlot(
  xen_focus,
  group.by = "focus_category",
  split.by = "group",
  label = FALSE,
  raster = TRUE,
  ncol = 2,
  cols = my_cols,
  pt.size = 1.5
) + ggtitle("Harmony UMAP — Progressor vs Survivor")

p_umap_split

ggsave("Path to save directory",
       plot = p_umap_split,
       width = 7, height = 3, units = "in",
       dpi = 600, bg = "white")

ggsave("Path to save directory",
       plot = p_umap_split,
       width = 7, height = 3, units = "in",
       dpi = 600, bg = "white", useDingbats = FALSE)

#----Stacked percentage bar plot----

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

keep_levels <- c(
  "Cardiomyocyte","Stressed CM","Endothelium",
  "Fibroblast","Myofibroblast","Fibrotic Fibroblast",
  "Immune"
)

xen_focus_bar <- subset(
  xen_focus,
  cells = colnames(xen_focus)[xen_focus$focus_category %in% keep_levels]
)
xen_focus_bar$focus_category <- factor(xen_focus_bar$focus_category,
                                       levels = keep_levels)

#----Same color as UMAP----
my_cols <- c(
  "Cardiomyocyte"       = "#F9B6C1",
  "Stressed CM"         = "#C2185B",
  "Endothelium"         = "#4CAF50",
  "Fibroblast"          = "#F4A582",
  "Myofibroblast"       = "#F4A582",
  "Fibrotic Fibroblast" = "#F4A582",
  "Immune"              = "#4A148C"
)

#----meta table----
meta_bar <- xen_focus_bar@meta.data %>%
  as.data.frame() %>%
  dplyr::transmute(
    group = factor(group, levels = c("Survivor","Progressor")),
    cat   = focus_category
  )

#----Proportions per group----
prop_tbl <- meta_bar %>%
  dplyr::count(group, cat, name = "n") %>%   
  dplyr::group_by(group) %>%
  dplyr::mutate(p = n / sum(n)) %>%
  dplyr::ungroup()

group_totals <- meta_bar %>%
  dplyr::count(group, name = "N") %>%
  dplyr::mutate(labN = paste0("n=", scales::comma(N)))

#----Stacked % bar plot----
p_stack <- ggplot(prop_tbl, aes(x = group, y = p, fill = cat)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.4) +
  geom_text(data = group_totals,
            aes(x = group, y = 1.02, label = labN),
            inherit.aes = FALSE, size = 4.2, fontface = "bold") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 10),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_fill_manual(values = my_cols, name = NULL) +
  labs(
    x = NULL,
    y = "Cell proportion (%)",
    title    = "Cellular composition by group",
    subtitle = "Survivor vs Progressor"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, margin = margin(b = 8)),
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    panel.grid    = element_blank(),
    legend.position = "right",
    legend.text     = element_text(size = 11)
  )

p_stack

#----Saving----

ggsave("Path to Save directory",
       plot = p_stack,
       width = 4, height = 4.6, units = "in",
       dpi = 600, bg = "white")

ggsave("Path to Save directory",
       plot = p_stack,
       width = 4, height = 4.6, units = "in",
       dpi = 600, bg = "white", useDingbats = FALSE)



#----Heatmap----

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(scales)
})

stopifnot(exists("xen_focus"))
DefaultAssay(xen_focus) <- "SCT"

#----Clusters----

keep_levels_hm <- c(
  "Cardiomyocyte","Stressed CM","Endothelium",
  "Fibroblast","Myofibroblast","Immune"
)

xen_focus_hm <- subset(
  xen_focus,
  cells = colnames(xen_focus)[xen_focus$focus_category %in%
                                c(keep_levels_hm, "Fibrotic Fibroblast")]
)

xen_focus_hm$focus_category <- as.character(xen_focus_hm$focus_category)
xen_focus_hm$focus_category[xen_focus_hm$focus_category == "Fibrotic Fibroblast"] <- "Fibroblast"
xen_focus_hm$focus_category <- factor(xen_focus_hm$focus_category, levels = keep_levels_hm)

stopifnot(all(c("group","focus_category") %in% colnames(xen_focus_hm@meta.data)))


#----Gene panel----

genes_show <- c(
  # cardiomyocyte / stress
  "MYH6","MYH7","NPPA","NPPB","TNNT2","ACTC1",
  # endothelium
  "PECAM1","VWF",
  # fibrosis / ECM
  "COL1A1","COL3A1","POSTN","FN1","MMP2","TIMP1",
  # T-lineage
  "CD4","CD8A","CD2","CD3D",
  # cytotoxic / activation
  "GZMB","GZMH","GZMA","GZMM","PRF1","NKG7","GNLY","MKI67","CD69","LCK",
  # antigen presentation
  "HLA-DQB1","HLA-DQA1","CD74","CD86"
)
genes_show <- unique(genes_show)

#----Sanity----
genes_ok <- genes_show[genes_show %in% rownames(xen_focus_hm)]
if (length(genes_ok) < 5) stop("Fewer than 5 requested genes found in object.")
if (length(genes_ok) < length(genes_show)) {
  message("Missing genes dropped: ", paste(setdiff(genes_show, genes_ok), collapse = ", "))
}


#----Build group and average expression----

xen_focus_hm$grp_cat <- paste(xen_focus_hm$group, xen_focus_hm$focus_category, sep = "|")

cols_order <- c("Survivor","Progressor")
cats_order <- c("Cardiomyocyte","Stressed CM","Endothelium","Fibroblast","Myofibroblast","Immune")
grpcat_levels <- as.vector(outer(cols_order, cats_order, paste, sep = "|"))

avg <- AverageExpression(
  xen_focus_hm,
  assays   = "SCT",
  slot     = "data",
  features = genes_ok,
  group.by = "grp_cat",
  verbose  = FALSE
)$SCT

#----Sanity----
common_cols <- intersect(grpcat_levels, colnames(avg))
if (length(common_cols) == 0) stop("No matching grp_cat columns found. Check group/focus_category values.")
mat <- as.matrix(avg[, common_cols, drop = FALSE])

#----Z-score----
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

gene_order <- genes_ok
mat_z <- mat_z[gene_order, , drop = FALSE]

#----column labels----
colnames(mat_z) <- gsub("\\|", " — ", colnames(mat_z))


#----Horizontal----
mat_h <- t(mat_z)
mat_h <- mat_h[, gene_order, drop = FALSE]

#----Row annotation----
grp_vec_rows <- ifelse(grepl("^Survivor", rownames(mat_h)), "Survivor", "Progressor")
grp_fac_rows <- factor(grp_vec_rows, levels = c("Survivor","Progressor"))
left_anno <- rowAnnotation(
  Group = grp_fac_rows,
  col   = list(Group = c(Survivor = "#2E7D32", Progressor = "#B71C1C")),
  gp    = gpar(col = NA),
  annotation_name_gp  = gpar(fontsize = 10),
  annotation_name_rot = 90
)

#----Gene blocks----
gene_groups <- list(
  Cardiomyocyte = c("MYH6","MYH7","TNNT2","ACTC1"),
  Stressed_CM   = c("NPPA","NPPB"),
  Endothelium   = c("PECAM1","VWF"),
  Fibrosis      = c("COL1A1","COL3A1","POSTN","FN1","MMP2","TIMP1"),
  T_lineage     = c("CD4","CD8A","CD2","CD3D"),
  Cytotoxic     = c("GZMA","GZMB","GZMH","GZMM","PRF1","NKG7","GNLY","MKI67","CD69","LCK"),
  AntigenPres   = c("HLA-DQB1","HLA-DQA1","CD74","CD86","CD80")
)

gene_group_vec <- rep(names(gene_groups), lengths(gene_groups))
names(gene_group_vec) <- unlist(gene_groups)

unmapped <- setdiff(gene_order, names(gene_group_vec))
if (length(unmapped) > 0) {
  gene_groups$Other <- unmapped
  gene_group_vec <- rep(names(gene_groups), lengths(gene_groups))
  names(gene_group_vec) <- unlist(gene_groups)
}
col_split_vec <- factor(gene_group_vec[gene_order], levels = names(gene_groups))

#----Color scale----
col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2C7BB6", "#FFFFFF", "#D7191C"))

#----Plot Size----
n_rows  <- nrow(mat_h)
n_cols  <- ncol(mat_h)
cell_mm <- 5.0
ht_w    <- unit(n_cols * cell_mm, "mm")
ht_h    <- unit(n_rows * cell_mm, "mm")

ht_horiz <- Heatmap(
  mat_h,
  name   = "Z",
  col    = col_fun,
  width  = ht_w,
  height = ht_h,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  column_split    = col_split_vec,
  show_row_names  = TRUE,
  show_column_names = TRUE,
  row_names_gp      = gpar(fontsize = 10),
  row_names_rot     = 0,
  column_names_gp   = gpar(fontsize = 10),
  column_names_rot  = 90,
  row_names_max_width     = unit(140, "mm"),
  column_names_max_height = unit(25, "mm"),
  left_annotation = left_anno,
  heatmap_legend_param = list(title = "Z-score", at = c(-2, 0, 2))
)

#----Preview----
grid.newpage()
draw(
  ht_horiz,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  padding = unit(c(25, 60, 25, 60), "mm")
)


#----Saving----

pad_top_mm    <- 25
pad_right_mm  <- 60
pad_bottom_mm <- 25
pad_left_mm   <- 60
pad_w_in <- (pad_left_mm + pad_right_mm) / 25.4
pad_h_in <- (pad_top_mm  + pad_bottom_mm) / 25.4
w_in <- (n_cols * cell_mm) / 25.4 + pad_w_in
h_in <- (n_rows * cell_mm) / 25.4 + pad_h_in

png_file <- "Path to Save Directory"
pdf_file <- "Path to Save Directory"

png(png_file, width = w_in, height = h_in, units = "in", res = 600,
    bg = "white", type = "cairo")
grid.newpage()
draw(
  ht_horiz,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  padding = unit(c(pad_top_mm, pad_right_mm, pad_bottom_mm, pad_left_mm), "mm")
)
dev.off()

pdf(pdf_file, width = w_in, height = h_in, bg = "white",
    useDingbats = FALSE, family = "Helvetica")
draw(
  ht_horiz,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right",
  padding = unit(c(pad_top_mm, pad_right_mm, pad_bottom_mm, pad_left_mm), "mm")
)
dev.off()

message("Saved heatmap to: ", png_file, " and ", pdf_file)
