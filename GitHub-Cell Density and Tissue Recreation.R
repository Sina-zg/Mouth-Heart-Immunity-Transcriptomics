
#----Density----

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(arrow); library(sf)
})

#----Load sample----
base_dir <- "Path to Spatial Transcriptomics File Per Patient"  
counts   <- readMM(file.path(base_dir, "matrix.mtx.gz"))
barcodes <- fread(file.path(base_dir, "barcodes.tsv.gz"), header = FALSE)$V1
feats    <- fread(file.path(base_dir, "features.tsv.gz"), header = FALSE)  
rownames(counts) <- feats$V2
colnames(counts) <- barcodes

cell_bounds <- read_parquet(file.path(base_dir, "cell_boundaries.parquet"), as_data_frame = TRUE)
stopifnot(all(c("cell_id","vertex_x","vertex_y") %in% names(cell_bounds)))

#----QC, remove zero-transcript cells & align boundaries----
tot_per_cell <- Matrix::colSums(counts)
keep <- tot_per_cell >= 1
counts <- counts[, keep, drop = FALSE]
cell_bounds <- cell_bounds[cell_bounds$cell_id %in% colnames(counts), ]

#---- ROI area µm² to mm²----
pts  <- st_as_sf(cell_bounds[, c("vertex_x","vertex_y")], coords = c("vertex_x","vertex_y"))
hull <- st_convex_hull(st_union(pts))
roi_area_um2 <- as.numeric(st_area(hull))
roi_area_mm2 <- roi_area_um2 / 1e6
if (!is.finite(roi_area_mm2) || roi_area_mm2 <= 0) stop("ROI area failed; check boundary units.")

#----presence = ≥1 molecule----
present_any <- function(genes) {
  sel <- intersect(genes, rownames(counts))
  if (!length(sel)) return(rep(FALSE, ncol(counts)))
  Matrix::colSums(counts[sel, , drop = FALSE] > 0) > 0
}

#----Gene sets----
avail <- rownames(counts)

#--T / macrophage--
cd3_genes <- intersect(c("CD3D","CD3E","CD3G","CD247","CD3Z"), avail)
macrophage_markers <- intersect(c("LYZ","C1QA","C1QB","C1QC","MS4A7","CD68","CSF1R","AIF1","CD163","MRC1","MARCO","CCR2","APOE","TREM2"), avail)

#--Cytotoxic for CD4 CTL--
cyto_genes <- intersect(c("GZMB","PRF1","GNLY","NKG7","GZMA","GZMM","GZMH",
                          "CCL5","CX3CR1","PLEK","ZEB2"), avail)

#--Expanded TCR genes--
tcr_panel <- intersect(c("TRBV7-9","TRBV28","TRBV6-1","TRBV5-1","TRBV2","TRBV5-4","TRBV19","TRBV24-1",
                         "TRAV38-2DV8","TRAV16","TRAV2","TRAV8-4"), avail)

#--Expanded CDR3 targets and removing microbial probes--
expanded_rows <- grep("Expanded", avail, value = TRUE, ignore.case = TRUE)
expanded_rows <- expanded_rows[!grepl("\\.|strain", expanded_rows, ignore.case = TRUE)]

#--Cardiomyocyte and stress-related genes--
cm_markers_base  <- c("MYH7","MYH6","ACTC1","TNNT2")
cm_markers_extra <- intersect(c("TNNI3","MYL2","MYH7B","TTN","MYL7","RYR2","TNNC1"), avail)
cm_markers <- intersect(unique(c(cm_markers_base, cm_markers_extra)), avail)
stress_markers <- intersect(c("NPPA","NPPB"), avail)

#--CD8 / B / Myeloid / HLA / Apoptosis / Fibrosis--
cd8_genes   <- intersect(c("CD8A","CD8B"), avail)
bcell_genes <- intersect(c("MS4A1","CD79A","CD79B","CD74","CD19","BANK1","IGKC","IGHM"), avail)
myeloid_genes <- intersect(unique(c(macrophage_markers, "ITGAM","S100A8","S100A9")), avail)
hla_genes   <- intersect(c("HLA-A","HLA-B","HLA-C","HLA-E",
                           "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"), avail)
apoptosis_genes <- intersect(c("CASP3","CASP7","BAX","BID","FAS","TNFRSF10A","TNFRSF10B","PMAIP1","BBC3","BCL2L11"), avail)
fibrosis_genes  <- intersect(c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A3","DCN","LUM",
                               "FAP","PDGFRA","PDGFRB","TAGLN","ACTA2","SPARC","FN1"), avail)


REQUIRE_T_FOR_CD4_CD8 <- FALSE   
EXCLUDE_MACRO_FROM_CD4 <- FALSE  
EXCLUDE_MACRO_FROM_CD8 <- FALSE  

#----Core masks----
immune_mask <- present_any("PTPRC")               

t_mask      <- present_any(cd3_genes)


cd4_mask <- present_any("CD4")
if (REQUIRE_T_FOR_CD4_CD8 && length(cd3_genes)) cd4_mask <- cd4_mask & present_any(cd3_genes)
if (EXCLUDE_MACRO_FROM_CD4) cd4_mask <- cd4_mask & (!present_any(macrophage_markers))

#----CD4 CTLs = CD4s with ≥1 cytotoxic gene----
cd4_ctl <- cd4_mask & present_any(cyto_genes)

#----CD4 CTLs with expanded TCR genes----
ctl_tcr <- cd4_ctl & present_any(tcr_panel)

#----CD4 CTLs with Expanded CDR3----
ctl_expanded <- if (length(expanded_rows)) (cd4_ctl & present_any(expanded_rows)) else rep(FALSE, ncol(counts))

#----CD8 cells----
cd8_mask <- present_any(cd8_genes)
if (REQUIRE_T_FOR_CD4_CD8 && length(cd3_genes)) cd8_mask <- cd8_mask & present_any(cd3_genes)
if (EXCLUDE_MACRO_FROM_CD8) cd8_mask <- cd8_mask & (!present_any(macrophage_markers))

#----B cells----
b_mask <- present_any(bcell_genes)

#----Myeloid cells----
myeloid_mask <- present_any(myeloid_genes)

#----HLA+ cells----
hla_mask <- present_any(hla_genes)

#----Apoptosis----
apoptosis_mask <- present_any(apoptosis_genes)

#----Cardiomyocytes----
cm_mask      <- present_any(cm_markers)
stressed_cm  <- cm_mask & present_any(stress_markers)
normal_cm    <- cm_mask & (!stressed_cm)

#----Fibrosis-associated----
fibrosis_mask <- present_any(fibrosis_genes)

#----Densities per mm^2 on ROI-level----
dens <- function(n) n / roi_area_mm2

density_report <- data.frame(
  ROI_area_mm2            = round(roi_area_mm2, 4),
  
  
  Normal_CM_per_mm2       = round(dens(sum(normal_cm)), 3),
  Stressed_CM_per_mm2     = round(dens(sum(stressed_cm)), 3),
  
  
  T_per_mm2               = round(dens(sum(t_mask)), 3),
  CD4_per_mm2             = round(dens(sum(cd4_mask)), 3),
  CD4_CTL_per_mm2         = round(dens(sum(cd4_ctl)), 3),
  CD4_CTL_TCR_per_mm2     = round(dens(sum(ctl_tcr)), 3),
  CD4_CTL_Expanded_per_mm2= round(dens(sum(ctl_expanded)), 3),
  
  
  CD8_per_mm2             = round(dens(sum(cd8_mask)), 3),
  B_per_mm2               = round(dens(sum(b_mask)), 3),
  Myeloid_per_mm2         = round(dens(sum(myeloid_mask)), 3),
  
  
  HLApos_per_mm2          = round(dens(sum(hla_mask)), 3),
  ApoptosisSig_per_mm2    = round(dens(sum(apoptosis_mask)), 3),
  
  
  Fibrosis_like_per_mm2   = round(dens(sum(fibrosis_mask)), 3)
)

print(density_report)




#----Cell Overlay on Tissue----

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(arrow); library(sf); library(ggplot2)
})

#----Load sample----
base_dir <- "Path to Spatial Transcriptomics File Per Patient"  
counts   <- readMM(file.path(base_dir, "matrix.mtx.gz"))
barcodes <- fread(file.path(base_dir, "barcodes.tsv.gz"), header = FALSE)$V1
feats    <- fread(file.path(base_dir, "features.tsv.gz"), header = FALSE)  
rownames(counts) <- trimws(feats$V2)
colnames(counts) <- barcodes

cell_bounds <- read_parquet(file.path(base_dir, "cell_boundaries.parquet"), as_data_frame = TRUE)
stopifnot(all(c("cell_id","vertex_x","vertex_y") %in% names(cell_bounds)))

#----QC,drop zero-transcript & align boundaries----
tot_per_cell <- Matrix::colSums(counts)
keep <- tot_per_cell >= 1
counts <- counts[, keep, drop = FALSE]
cell_bounds <- cell_bounds[cell_bounds$cell_id %in% colnames(counts), ]
stopifnot(nrow(cell_bounds) > 0)


RN <- toupper(trimws(rownames(counts)))
present_any <- function(genes) {
  g <- toupper(trimws(genes))
  idx <- which(RN %in% g)
  if (!length(idx)) return(rep(FALSE, ncol(counts)))
  Matrix::colSums(counts[idx, , drop = FALSE] > 0) > 0
}

#----Gene sets----
cm_markers      <- c("MYH7","MYH6","ACTC1","TNNT2","TNNI3","MYL2")
stress_markers  <- c("NPPA","NPPB")
cd3_genes       <- c("CD3D","CD3E","CD3G","CD247","CD3Z")
macrophage_markers <- c("LYZ","C1QA","C1QB","C1QC","MS4A7","CD68","CSF1R","AIF1","CD163","MRC1","MARCO","CCR2","CX3CR1","APOE","TREM2")
cyto_genes      <- c("GZMB","PRF1","GNLY","NKG7","GZMA","GZMM","GZMH", "ZEB2", "PLEK", "CCL5", "CX3CR1")

#----Expanded CDR3 and removing bacterial probes----
avail <- rownames(counts)
expanded_rows <- grep("Expanded", avail, value = TRUE, ignore.case = TRUE)
expanded_rows <- expanded_rows[!grepl("\\.|strain|Aggregati|Haemophilus|Streptococcus|Staphylococcus|Fusobacterium|Porphyromonas",
                                      expanded_rows, ignore.case = TRUE)]

#----Masks----

#--Cardiomyocytes--
cm_mask      <- present_any(cm_markers)
stressed_cm  <- cm_mask & present_any(stress_markers)
normal_cm    <- cm_mask & (!stressed_cm)

#--T lineage--
t_mask   <- present_any(cd3_genes) & (!present_any(macrophage_markers))

#--CD4 T-cells--
cd4_mask <- present_any("CD4")

#--CD4 CTL--
cd4_ctl  <- cd4_mask & present_any(cyto_genes)

#--CD4 CTL + Expanded CDR3--
ctl_expanded <- if (length(expanded_rows)) (cd4_ctl & present_any(expanded_rows)) else rep(FALSE, ncol(counts))

#----Sanity----
cat("\nMatched markers present on panel:\n")
show_found <- function(label, lst) {
  hit <- intersect(toupper(lst), RN)
  cat(sprintf("  %-12s: %s\n", label, if (length(hit)) paste(hit, collapse=", ") else "NONE"))
}
show_found("CM", cm_markers)
show_found("Stress", stress_markers)
show_found("CD3", cd3_genes)
show_found("Cyto", cyto_genes)
cat(sprintf("  Expanded CDR3 probes detected: %d\n", length(expanded_rows)))

cat("\nCounts by class (cells):\n")
cat(sprintf("  CM_normal           = %d\n", sum(normal_cm)))
cat(sprintf("  CM_stressed         = %d\n", sum(stressed_cm)))
cat(sprintf("  CD4                 = %d\n", sum(cd4_mask)))
cat(sprintf("  CD4_CTL             = %d\n", sum(cd4_ctl)))
cat(sprintf("  CD4_CTL + Expanded  = %d\n", sum(ctl_expanded)))

#----Centroids for plotting----
centroids <- as.data.table(cell_bounds)[, .(
  x = mean(vertex_x, na.rm = TRUE),
  y = mean(vertex_y, na.rm = TRUE)
), by = cell_id]


flag_df <- data.frame(
  cell_id = colnames(counts),
  is_CM_normal     = cm_mask & !stressed_cm,
  is_CM_stressed   = stressed_cm,
  is_CD4           = cd4_mask,
  is_CTL_Expanded  = ctl_expanded,
  stringsAsFactors = FALSE
)
plot_df <- merge(centroids, flag_df, by = "cell_id", all.x = TRUE)
for (nm in c("is_CM_normal","is_CM_stressed","is_CD4","is_CTL_Expanded")) {
  plot_df[[nm]][is.na(plot_df[[nm]])] <- FALSE
}


plot_df$category <- "Other"
plot_df$category[plot_df$is_CM_normal]    <- "CM_normal"
plot_df$category[plot_df$is_CM_stressed]  <- "CM_stressed"
plot_df$category[plot_df$is_CD4]          <- "CD4_T"
plot_df$category[plot_df$is_CTL_Expanded] <- "CD4_CTL_Expanded"
plot_df$category <- factor(plot_df$category,
                           levels = c("CM_normal","CM_stressed","CD4_T","CD4_CTL_Expanded","Other"))

#----Cell Colors----
cols <- c(
  CM_normal        = "#ffb6c1",  
  CM_stressed      = "#ff1493",  
  CD4_T            = "green",  
  CD4_CTL_Expanded = "blue",  
  Other            = "grey70"
)

p <- ggplot() +
  geom_point(
    data = subset(plot_df, category %in% c("CM_normal","CM_stressed","CD4_T","CD4_CTL","Other")),
    aes(x = x, y = y, color = category),
    size = 1.8, alpha = 0.7
  ) +
  geom_point(  
    data = subset(plot_df, category == "CD4_CTL_Expanded"),
    aes(x = x, y = y),
    color = cols["CD4_CTL_Expanded"], size = 2.8, alpha = 0.95
  ) +
  scale_color_manual(values = cols, drop = FALSE) +
  coord_equal() + theme_void() +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(title = "Spatial map")

print(p)

#----saving----
ggsave(file.path(base_dir, "overlay_map_expanded.png"), p, width = 8, height = 7, dpi = 300)






