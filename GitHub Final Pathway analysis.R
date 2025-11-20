
#----Project Direction----
PROJECT_DIR <- "Path to project"
dirs <- file.path(PROJECT_DIR, c("qc_pngs","results/plots","results/tables","logs","rds"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
setwd(PROJECT_DIR)
sink(file.path("logs", paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")), split = TRUE)

## ---- 1) Packages ----
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(Matrix); library(purrr)
  library(SingleCellExperiment); library(scDblFinder); library(glmGamPoi)
  library(DESeq2); library(fgsea); library(msigdbr)
  library(readr); library(forcats); library(stringr)
  library(future)
})
plan(sequential)
options(future.globals.maxSize = 20 * 1024^3)
set.seed(1234)

#----Sample Loading----
samples <- c("CV5","CV12","CV109","CV110","CV39","CV87","CV103","CV106",
             "PD1","PD2","PD3","PD4","PD5","PD6","PD7","PD8")
group_map <- c(
  CV5="NICM_Progressor", CV12="NICM_Progressor", CV109="NICM_Progressor", CV110="NICM_Progressor",
  CV39="NICM_Survivor",  CV87="NICM_Survivor",  CV103="NICM_Survivor",  CV106="NICM_Survivor",
  PD1="PD_HighCTL", PD2="PD_HighCTL", PD3="PD_HighCTL", PD4="PD_HighCTL",
  PD5="PD_LowCTL",  PD6="PD_LowCTL",  PD7="PD_LowCTL",  PD8="PD_LowCTL"
)
rna_base <- "Path to RNA Path"
rna_path <- function(s) file.path(rna_base, paste0("sample_filtered_feature_bc_matrix_", s))
stopifnot(all(dir.exists(vapply(samples, rna_path, character(1)))))

#----QC Check----
mad_hi <- function(x, k=4) { m <- median(x, na.rm=TRUE); m + k*mad(x, na.rm=TRUE) }

plot_qc <- function(seu, s, outdir="qc_pngs"){
  p1 <- VlnPlot(seu, features=c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3, pt.size=0) +
    ggtitle(paste(s, "— violin"))
  p2 <- FeatureScatter(seu, feature1="nCount_RNA", feature2="nFeature_RNA") + ggtitle("nCount vs nFeature")
  p3 <- FeatureScatter(seu, feature1="percent.mt", feature2="nFeature_RNA") + ggtitle("%MT vs nFeature")
  png(file.path(outdir, paste0(s, "_QC.png")), width=1800, height=600, res=150)
  print(p1); print(p2); print(p3)
  dev.off()
}

gate_cd3d_cd4_drop_cd8a <- function(seu) {
  DefaultAssay(seu) <- "SCT"
  X <- GetAssayData(seu, slot="data")
  g <- rownames(X)
  v <- function(gene) if (gene %in% g) X[gene,] else rep(NA_real_, ncol(seu))
  cd3d <- v("CD3D"); cd4 <- v("CD4"); cd8a <- v("CD8A")
  keep <- (cd3d > 1) & (cd4 > 1) & (is.na(cd8a) | cd8a <= 1)
  subset(seu, cells = colnames(seu)[keep])
}

#----Load & QC per sample----
objs <- list(); qc_rows <- list()
for (s in samples) {
  mtx <- Read10X(rna_path(s))
  seu <- CreateSeuratObject(mtx, project=s, min.features=100, min.cells=3)
  seu$sample_id <- s
  seu$group     <- group_map[[s]]
  seu$cohort    <- if (startsWith(s,"CV")) "NICM" else "PD"
  seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern="^MT-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern="^RP[SL]")
  seu[["percent.hb"]]   <- PercentageFeatureSet(seu, pattern="^HB[AB]")
  
  raw <- ncol(seu)
  hard_keep <- with(seu@meta.data,
                    (nFeature_RNA > 500) & (nFeature_RNA < 6500) & (percent.mt < 10))
  seu <- subset(seu, cells = rownames(seu@meta.data)[hard_keep])
  
  nf_hi <- mad_hi(seu$nFeature_RNA, k=4)
  nc_hi <- mad_hi(seu$nCount_RNA,   k=4)
  mt_hi <- min(10, mad_hi(seu$percent.mt, k=3))
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
write_csv(qc_stage1, "results/tables/qc_stage1_counts.csv")

#----Doublets Removal----
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
write_csv(qc_stage2, "results/tables/qc_stage2_counts.csv")

#----Normalization & SCTransform----
sct_list <- lapply(singlets, function(seu){
  SCTransform(seu, assay="RNA", variable.features.n=3000,
              vars.to.regress="percent.mt", method="glmGamPoi", verbose=FALSE)
})
sct_rows <- data.frame(sample_id = names(sct_list),
                       after_SCT = vapply(sct_list, ncol, integer(1)))
qc_stage3 <- left_join(qc_stage2, sct_rows, by="sample_id")
write_csv(qc_stage3, "results/tables/qc_stage3_counts.csv")

#----Filtering out non CD4----
gated_list <- lapply(sct_list, gate_cd3d_cd4_drop_cd8a)
gate_rows <- data.frame(sample_id = names(gated_list),
                        after_gate = vapply(gated_list, ncol, integer(1)))
qc_final <- qc_stage3 |>
  left_join(gate_rows, by="sample_id") |>
  mutate(removed_by_gate = after_SCT - after_gate) |>
  arrange(sample_id)
write_csv(qc_final, "results/tables/qc_stage_all_counts.csv")
print(qc_final)

#--------Pseudobulk & DEGs--------
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(purrr); library(tibble)
  library(DESeq2); library(fgsea); library(msigdbr)
  library(ggplot2); library(forcats); library(readr); library(stringr)
})

message("\n=== Stage 8: Pseudobulk + DE + GSEA ===")

.get_counts_rna <- function(seu){
  if (!"RNA" %in% names(seu@assays)) stop("No RNA assay in object: ", seu$sample_id %||% seu@project.name)
  GetAssayData(seu, assay = "RNA", slot = "counts")
}

#----Pseudobulk per sample----
pb_one <- function(seu){
  stopifnot("sample_id" %in% colnames(seu@meta.data),
            "group"     %in% colnames(seu@meta.data),
            "cohort"    %in% colnames(seu@meta.data))
  cts <- .get_counts_rna(seu)
  if (ncol(cts) == 0) {
    warning("Zero cells after gating for sample_id=", unique(seu$sample_id))
    return(NULL)
  }
  smp <- seu$sample_id
 
  agg <- rowsum(t(as.matrix(cts)), group = smp)  
  list(
    counts = t(agg),  
    meta   = seu@meta.data %>%
      distinct(sample_id, group, cohort) %>%
      arrange(sample_id)
  )
}


make_dds <- function(pb, cohort_keep, grp_ref, grp_alt){
  meta <- pb$meta %>%
    filter(cohort == cohort_keep) %>%
    mutate(group = factor(group, levels = c(grp_ref, grp_alt))) %>%
    arrange(sample_id)
  
  cts  <- pb$counts
  if (is.null(cts) || ncol(cts) == 0) stop("Empty pseudobulk counts matrix")
  
  
  cts <- cts[, colnames(cts) %in% meta$sample_id, drop = FALSE]
  
  
  cat("\n[make_dds] Cohort:", cohort_keep,
      "\n  counts cols:", ncol(cts), " meta rows:", nrow(meta),
      "\n  missing in counts:", paste(setdiff(meta$sample_id, colnames(cts)), collapse=", "),
      "\n  missing in meta:",   paste(setdiff(colnames(cts), meta$sample_id), collapse=", "), "\n")
  
  keep <- intersect(meta$sample_id, colnames(cts))
  if (length(keep) < 2L) stop("Not enough samples after intersection for cohort ", cohort_keep,
                              ". Need ≥2 samples; check gating/QC.")
  
  
  cts  <- cts[, keep, drop = FALSE]
  meta <- meta[match(keep, meta$sample_id), , drop = FALSE]
  rownames(meta) <- meta$sample_id
  stopifnot(identical(colnames(cts), rownames(meta)))
  
 
  cts <- round(as.matrix(cts))
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = meta, design = ~ group)
  
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds
}

#----Run DE with Wald statistic----
run_de_rank <- function(dds, grp_ref, grp_alt){
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", grp_alt, grp_ref))
  tt  <- as.data.frame(res) %>% rownames_to_column("gene")
  stats <- tt$stat; names(stats) <- tt$gene
  stats <- stats[is.finite(stats)]
  stats <- sort(stats, decreasing = TRUE)
  list(res_table = tt, ranks = stats)
}

#----GOBP----
build_msig_all <- function(){
  h_all    <- msigdbr(species="Homo sapiens", category="H") %>%
    select(gs_name, gene_symbol)
  gobp_all <- msigdbr(species="Homo sapiens", category="C5", subcategory="BP") %>%
    select(gs_name, gene_symbol)
  bind_rows(h_all, gobp_all) %>%
    group_by(gs_name) %>%
    summarise(genes = list(unique(gene_symbol)), .groups="drop") %>%
    deframe()
}

#----fgsea----
do_fgsea <- function(ranks, pathways, minSize=10, maxSize=2000, scoreType="std"){
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- ranks[is.finite(ranks)]
  fgseaMultilevel(
    pathways = pathways,
    stats    = ranks,
    minSize  = minSize,
    maxSize  = maxSize,
    scoreType = scoreType
  ) %>%
    arrange(padj) %>%
    mutate(label = sub("^HALLMARK_", "HM:", pathway))
}

#----Preview top enriched or depleted----
plot_lollipop <- function(fg, top_n=15, side=c("enrich","deplete")){
  side <- match.arg(side)
  dat <- fg %>% mutate(dir = ifelse(NES>0, "up", "down"))
  if (side == "enrich") dat <- filter(dat, NES > 0) else dat <- filter(dat, NES < 0)
  dat %>%
    slice_min(padj, n = min(top_n, n())) %>%
    mutate(label = fct_reorder(label, NES)) %>%
    ggplot(aes(x = NES, y = label)) +
    geom_point(aes(size = -log10(padj))) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(x = "NES", y = NULL, size = expression(-log[10](FDR))) +
    theme_classic(base_size = 12)
}

#----NES across contrasts----
plot_concordance <- function(ov){
  ggplot(ov, aes(NES_PD, NES_NICM)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(aes(size = -log10(pmax(FDR_PD, FDR_NICM)))) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "NES (PD HighCTL vs LowCTL)", y = "NES (NICM Progressor vs Survivor)",
         size = expression(-log[10]~FDR)) +
    theme_classic(base_size = 12)
}

#----pseudobulk across all samples----

stopifnot(exists("gated_list"), length(gated_list) > 0)

pb_parts <- map(gated_list, pb_one)
pb_parts <- compact(pb_parts)  

if (length(pb_parts) == 0) stop("All samples empty.")

all_genes   <- Reduce(union, map(pb_parts, ~ rownames(.x$counts)))
all_samples <- unlist(map(pb_parts, ~ colnames(.x$counts)))
M <- matrix(0, nrow = length(all_genes), ncol = length(all_samples),
            dimnames = list(all_genes, all_samples))
meta_rows <- list(); idx <- 1
for (x in pb_parts) {
  g <- rownames(x$counts); s <- colnames(x$counts)
  M[g, idx:(idx + length(s) - 1)] <- x$counts
  meta_rows[[length(meta_rows) + 1]] <- x$meta
  idx <- idx + length(s)
}
pb <- list(counts = M, meta = bind_rows(meta_rows))


write_csv(
  tibble(sample_id = colnames(pb$counts),
         n_genes   = nrow(pb$counts),
         library   = colSums(pb$counts)) %>%
    left_join(pb$meta, by = "sample_id"),
  "results/tables/pseudobulk_library_sizes.csv"
)

#----DESeq2 objects per cohort----

dds_pd <- make_dds(pb, cohort_keep="PD",   grp_ref="PD_LowCTL",     grp_alt="PD_HighCTL")
dds_nc <- make_dds(pb, cohort_keep="NICM", grp_ref="NICM_Survivor", grp_alt="NICM_Progressor")

#----ranks----

pd <- run_de_rank(dds_pd, grp_ref="PD_LowCTL",       grp_alt="PD_HighCTL")
nc <- run_de_rank(dds_nc, grp_ref="NICM_Survivor",   grp_alt="NICM_Progressor")

write_csv(pd$res_table, "results/tables/DE_PD_High_vs_Low.csv")
write_csv(nc$res_table, "results/tables/DE_NICM_Prog_vs_Surv.csv")

#----Gene sets and GSEA----

msig_list <- build_msig_all()

pd_gsea <- do_fgsea(pd$ranks, pathways = msig_list)
nc_gsea <- do_fgsea(nc$ranks, pathways = msig_list)

write_csv(pd_gsea, "results/tables/GSEA_PD_High_vs_Low__ALL_Hallmark_GO-BP.csv")
write_csv(nc_gsea, "results/tables/GSEA_NICM_Prog_vs_Surv__ALL_Hallmark_GO-BP.csv")

#----Plotting----

plot_lollipop <- function(fg, top_n = 15, side = c("enrich","deplete")){
  side <- match.arg(side)
  dat <- fg %>% dplyr::mutate(dir = ifelse(NES > 0, "up", "down"))
  if (side == "enrich") {
    dat <- dplyr::filter(dat, NES > 0)
  } else {
    dat <- dplyr::filter(dat, NES < 0)
  }
  
  take_n <- min(top_n, nrow(dat))
  if (take_n < 1) {
    stop("No pathways to plot for side = '", side, "'.")
  }
  dat %>%
    dplyr::arrange(padj) %>%
    dplyr::slice(1:take_n) %>%  
    dplyr::mutate(label = forcats::fct_reorder(sub("^HALLMARK_", "HM:", pathway), NES)) %>%
    ggplot2::ggplot(ggplot2::aes(x = NES, y = label)) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(padj))) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::labs(x = "NES", y = NULL, size = expression(-log[10](FDR))) +
    ggplot2::theme_classic(base_size = 12)
}

g_pd_up <- plot_lollipop(pd_gsea, top_n = 15, side = "enrich") +
  ggtitle("PD HighCTL vs LowCTL — Top enriched (ALL sets)")
g_nc_up <- plot_lollipop(nc_gsea, top_n = 15, side = "enrich") +
  ggtitle("NICM Progressor vs Survivor — Top enriched (ALL sets)")

ggsave("results/plots/PD_lollipop_up_ALL.png", g_pd_up, width=10, height=10, dpi=300)
ggsave("results/plots/NICM_lollipop_up_ALL.png", g_nc_up, width=10, height=10, dpi=300)



both <- inner_join(
  pd_gsea %>% select(gs_name = pathway, NES_PD = NES, FDR_PD = padj),
  nc_gsea %>% select(gs_name = pathway, NES_NICM = NES, FDR_NICM = padj),
  by = "gs_name"
) %>%
  mutate(rank_PD   = rank(FDR_PD, ties.method="min"),
         rank_NICM = rank(FDR_NICM, ties.method="min"),
         mean_rank = (rank_PD + rank_NICM) / 2) %>%
  arrange(mean_rank)

write_csv(both, "results/tables/GSEA_overlap_PD_vs_NICM__ALL_sets.csv")

p_concord <- plot_concordance(both)
ggsave("results/plots/Concordance_PD_vs_NICM.png", p_concord, width=10, height=5, dpi=300)

cat("\n=== DONE: Stage 8 ===\n",
    "Tables:\n",
    "  - results/tables/DE_PD_High_vs_Low.csv\n",
    "  - results/tables/DE_NICM_Prog_vs_Surv.csv\n",
    "  - results/tables/GSEA_PD_High_vs_Low__ALL_Hallmark_GO-BP.csv\n",
    "  - results/tables/GSEA_NICM_Prog_vs_Surv__ALL_Hallmark_GO-BP.csv\n",
    "  - results/tables/GSEA_overlap_PD_vs_NICM__ALL_sets.csv\n",
    "Figures:\n",
    "  - results/plots/PD_lollipop_up_ALL.png\n",
    "  - results/plots/NICM_lollipop_up_ALL.png\n",
    "  - results/plots/PD_lollipop_down_ALL.png\n",
    "  - results/plots/NICM_lollipop_down_ALL.png\n",
    "  - results/plots/Concordance_PD_vs_NICM.png\n")



#--------PD HighCTL vs LowCTL T-cell Pathway Enrichment--------

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(purrr); library(tibble)
  library(DESeq2); library(fgsea); library(msigdbr)
  library(ggplot2); library(forcats); library(readr); library(stringr)
})

message("\n=== PD HighCTL vs PD LowCTL: start ===")

#----Pseudobulk (genes x samples)----
.get_counts_rna <- function(seu){
  if (!"RNA" %in% names(seu@assays)) stop("No RNA assay in object: ", seu$sample_id %||% seu@project.name)
  GetAssayData(seu, assay = "RNA", slot = "counts")
}
pb_one <- function(seu){
  stopifnot(all(c("sample_id","group","cohort") %in% colnames(seu@meta.data)))
  cts <- .get_counts_rna(seu); smp <- seu$sample_id
  if (ncol(cts) == 0) return(NULL)
  agg <- rowsum(t(as.matrix(cts)), group = smp) 
  list(counts = t(agg),
       meta = seu@meta.data |> distinct(sample_id, group, cohort) |> arrange(sample_id))
}
pb_parts <- gated_list |> map(pb_one) |> compact()
stopifnot(length(pb_parts) > 0)

all_genes   <- Reduce(union, map(pb_parts, ~ rownames(.x$counts)))
all_samples <- unlist(map(pb_parts, ~ colnames(.x$counts)))
M <- matrix(0, nrow=length(all_genes), ncol=length(all_samples),
            dimnames=list(all_genes, all_samples))
meta_rows <- list(); k <- 1
for (x in pb_parts) {
  g <- rownames(x$counts); s <- colnames(x$counts)
  M[g, k:(k+length(s)-1)] <- x$counts
  meta_rows[[length(meta_rows)+1]] <- x$meta
  k <- k + length(s)
}
pb <- list(counts = M, meta = bind_rows(meta_rows))

#----Build DESeq2 for PD----
make_dds_pd <- function(pb){
  meta <- pb$meta %>%
    filter(cohort == "PD") %>%
    mutate(group = factor(group, levels = c("PD_LowCTL","PD_HighCTL"))) %>%
    arrange(sample_id)
  cts  <- pb$counts[, colnames(pb$counts) %in% meta$sample_id, drop = FALSE]
  
  #--sanity--
  cat("PD samples — counts cols:", ncol(cts), " meta rows:", nrow(meta), "\n")
  if (nrow(meta) < 2 || ncol(cts) < 2) stop("Need >=2 PD samples per contrast.")
  
  keep <- intersect(meta$sample_id, colnames(cts))
  cts  <- cts[, keep, drop = FALSE]
  meta <- meta[match(keep, meta$sample_id), , drop = FALSE]
  rownames(meta) <- meta$sample_id
  stopifnot(identical(colnames(cts), rownames(meta)))
  
  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(cts)),
                                colData = meta, design = ~ group)
  dds[rowSums(counts(dds)) >= 10, ]
}
dds_pd <- make_dds_pd(pb)

#----DE (Wald stat)----
run_de_rank <- function(dds){
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds, contrast = c("group","PD_HighCTL","PD_LowCTL"))
  tt  <- as.data.frame(res) |> rownames_to_column("gene")
  stats <- tt$stat; names(stats) <- tt$gene
  stats <- stats[is.finite(stats)]
  stats <- sort(stats, decreasing = TRUE)
  list(res_table = tt, ranks = stats)
}
pd <- run_de_rank(dds_pd)
write_csv(pd$res_table, "results/tables/DE_PD_High_vs_Low.csv")

#----All Gene sets Hallmark----

msig_list <- {
  h_all <- msigdbr(species = "Homo sapiens", category = "H") |>
    dplyr::select(gs_name, gene_symbol)
  
  gobp_all <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") |>
    dplyr::select(gs_name, gene_symbol)
  
  dplyr::bind_rows(h_all, gobp_all) |>
    dplyr::group_by(gs_name) |>
    dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") |>
    tibble::deframe()
}


cat("Gene set count:", length(msig_list), "\n")
cat("Example set size:", length(msig_list[[1]]), "\n")


#----GSEA----
do_fgsea <- function(ranks, pathways, minSize=10, maxSize=2000, scoreType="std"){
  ranks <- ranks[!duplicated(names(ranks))]
  fgseaMultilevel(pathways=pathways, stats=ranks, minSize=minSize, maxSize=maxSize, scoreType=scoreType) |>
    as.data.frame() |>
    mutate(pathway = as.character(pathway), NES = as.numeric(NES), padj = as.numeric(padj)) |>
    arrange(padj)
}
pd_gsea <- do_fgsea(pd$ranks, pathways = msig_list)
write_csv(pd_gsea, "results/tables/GSEA_PD_High_vs_Low__ALL.csv")

#----T-cell–related pathways enriched in PD High----
tcell_pattern <- "(T_CELL|LYMPHOCYTE|LEUKOCYTE|CYTOTOX|KILL|GRANULE|DEGRANUL|TCR|CD3|ACTIVATION|PROLIFER|SELECTION|SIGNALING)"
pd_tcell <- pd_gsea %>%
  filter(NES > 0, str_detect(pathway, regex(tcell_pattern, ignore_case=TRUE))) %>%
  mutate(minuslog10FDR = -log10(pmax(padj, 1e-300)),
         label = str_to_sentence(gsub("_"," ", sub("^HALLMARK_","HM: ", pathway)))) %>%
  arrange(padj)

#----top 12 pathways----
top_n <- min(12, nrow(pd_tcell))

pd_top12 <- pd_tcell %>%
  as_tibble() %>%                        
  dplyr::slice(1:min(n(), top_n)) %>%    
  mutate(
    label = fct_reorder(as.factor(label), as.numeric(minuslog10FDR))
  )
#-----Sanity check for cytotoxicity genes----

table_groups <- colData(dds_pd)$group |> table()
cat("PD group sizes:", paste(names(table_groups), as.integer(table_groups), collapse=" | "), "\n")

ct_markers <- c("GZMB","PRF1","GNLY","NKG7","CX3CR1","IFNG")
tt_pd <- pd$res_table %>% as.data.frame()
rownames(tt_pd) <- tt_pd$gene
tt_pd$gene <- NULL
has_markers <- intersect(ct_markers, rownames(tt_pd))
if (length(has_markers) > 0) {
  dir_ok <- tt_pd[has_markers,"log2FoldChange"] > 0
  cat("Marker direction check (PD High > PD Low):\n")
  print(data.frame(gene=has_markers,
                   log2FC = round(tt_pd[has_markers,"log2FoldChange"],3),
                   ok = dir_ok))
} else {
  warning("None of the cytotoxic markers were found in DE table rownames.")
}


if (top_n > 0) {
  top_term <- pd_top12$pathway[1]
  le <- pd_gsea$leadingEdge[match(top_term, pd_gsea$pathway)][[1]]
  cat("Leading-edge (top term): ", top_term, "\n")
  print(intersect(le, ct_markers))
}

#----Plotting----
plot_pd_high_tcell <- function(df, fdr_cut=0.05,
                               title="PD HighCTL vs LowCTL — T cell pathways enriched in PD High") {
  cutx <- -log10(fdr_cut)
  ggplot(df, aes(x = minuslog10FDR, y = label)) +
    geom_point(size = 5, colour = "#d84a3a") +
    geom_text(aes(label = sprintf("%.1f", NES)), color = "white", size = 2.8, fontface = "bold") +
    geom_vline(xintercept = cutx, linetype = 2) +
    scale_x_continuous(name = expression(-log[10]~"(adjusted p-value)"),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(y = NULL, title = title) +
    theme_classic(base_size = 12)
}
p_pd_high <- plot_pd_high_tcell(pd_top12, fdr_cut = 0.05)
ggsave("results/plots/PD_Tcell_TOP12_PDhigh_only_minuslog10FDR.png",
       p_pd_high, width=7.2, height=5.5, dpi=300)


write_csv(pd_tcell, "results/tables/PD_Tcell_enriched_in_PDHigh_FULL.csv")

message("\n=== PD HighCTL vs PD LowCTL: done ===\n",
        "Plots:\n  results/plots/PD_Tcell_TOP12_PDhigh_only_minuslog10FDR.png\n",
        "Tables:\n  results/tables/GSEA_PD_High_vs_Low__ALL.csv\n",
        "  results/tables/PD_Tcell_enriched_in_PDHigh_FULL.csv\n")


plot_pd_high_tcell <- function(df, fdr_cut = 0.05,
                               title = "PD HighCTL vs LowCTL — T cell pathways enriched in PD High") {
  cutx <- -log10(fdr_cut)
  
  ggplot(df, aes(x = minuslog10FDR, y = label)) +
    
    geom_segment(aes(x = 0, xend = minuslog10FDR, y = label, yend = label),
                 color = "grey80", linewidth = 3) +
    
    geom_point(size = 7, color = "darkred") +
    
    geom_text(aes(label = sprintf("%.1f", NES)), color = "white",
              size = 3, fontface = "bold") +
    
    geom_vline(xintercept = cutx, linetype = 2) +
    
    scale_x_continuous(name = expression(-log[10]~"(adjusted p-value)"),
                       limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(y = NULL, title = title) +
    theme_classic(base_size = 12) +
    theme(panel.grid = element_blank())
}


p_pd_high <- plot_pd_high_tcell(
  pd_top12,
  fdr_cut = 0.05,
  title   = "PD HighCTL vs LowCTL — T cell pathways enriched in PD High"
)


p_pd_high   


ggsave(
  filename = "results/plots/PD_Tcell_TOP12_PDhigh_only_minuslog10FDR.png",
  plot     = p_pd_high,
  width    = 7, height = 4, dpi = 600
)



#--------NICM Progressor vs Survivor Wald-rank--------


#----DESeq2 for NICM----
make_dds_nc <- function(pb){
  meta <- pb$meta %>%
    dplyr::filter(cohort == "NICM") %>%
    dplyr::mutate(group = factor(group, levels = c("NICM_Survivor","NICM_Progressor"))) %>%
    dplyr::arrange(sample_id)
  cts  <- pb$counts[, colnames(pb$counts) %in% meta$sample_id, drop = FALSE]
  
  keep <- intersect(meta$sample_id, colnames(cts))
  cts  <- cts[, keep, drop = FALSE]
  meta <- meta[match(keep, meta$sample_id), , drop = FALSE]
  rownames(meta) <- meta$sample_id
  stopifnot(identical(colnames(cts), rownames(meta)))
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(as.matrix(cts)),
                                        colData = meta, design = ~ group)
  dds[rowSums(DESeq2::counts(dds)) >= 10, ]
}

run_de_rank <- function(dds, contrast_vec){
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = contrast_vec)
  tt  <- as.data.frame(res)
  tt  <- tibble::rownames_to_column(tt, "gene")
  stats <- tt$stat; names(stats) <- tt$gene
  stats <- stats[is.finite(stats)]
  stats <- sort(stats, decreasing = TRUE)
  list(res_table = tt, ranks = stats)
}

dds_nc <- make_dds_nc(pb)

nc      <- run_de_rank(dds_nc, c("group","NICM_Progressor","NICM_Survivor"))

#----GSEA with Wald ranks----
nc_gsea_wald <- do_fgsea(nc$ranks, pathways = msig_list)
readr::write_csv(nc_gsea_wald, "results/tables/GSEA_NICM_Prog_vs_Surv__WALD_ALL.csv")

#----Same 12 PD pathways for Progressors----
keep12 <- pd_top12$pathway

prep_side <- function(gsea_df, label){
  gsea_df %>%
    dplyr::filter(pathway %in% keep12) %>%
    dplyr::transmute(
      pathway,
      label_y = stringr::str_to_sentence(gsub("_"," ", sub("^HALLMARK_","HM: ", pathway))),
      minuslog10FDR = -log10(pmax(padj, .Machine$double.xmin)),
      NES = as.numeric(NES),
      contrast = label
    )
}

pd_side <- prep_side(pd_gsea, "PD HighCTL")
nc_side <- prep_side(nc_gsea, "NICM Progressor")

panel_both <- dplyr::bind_rows(pd_side, nc_side) %>%
  dplyr::group_by(label_y) %>%
  dplyr::mutate(order_key = max(minuslog10FDR[contrast=="PD HighCTL"], na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    label_y  = forcats::fct_reorder(label_y, order_key),
    contrast = factor(contrast, levels = c("PD HighCTL","NICM Progressor")),
    dot_lab  = sprintf("%.1f", NES)  
  )


guides <- panel_both %>%
  dplyr::group_by(label_y) %>%
  dplyr::summarise(xend = max(minuslog10FDR, na.rm = TRUE), .groups = "drop")

#----Plotting----
plot_both <- function(df, guides, fdr_cut = 0.05,
                      title = "T cell pathways — PD HighCTL & NICM Progressor (Wald rank)") {
  cutx <- -log10(fdr_cut)
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = guides,
                          ggplot2::aes(x = 0, xend = xend, y = label_y, yend = label_y),
                          color = "grey80", linewidth = 3) +
    ggplot2::geom_point(data = df,
                        ggplot2::aes(x = minuslog10FDR, y = label_y, color = contrast),
                        size = 6, alpha = 0.95) +
    ggplot2::geom_text(data = df,
                       ggplot2::aes(x = minuslog10FDR, y = label_y, label = dot_lab),
                       color = "white", size = 3, fontface = "bold") +
    ggplot2::geom_vline(xintercept = cutx, linetype = 2) +
    ggplot2::scale_color_manual(values = c("PD HighCTL" = "darkred",     
                                           "NICM Progressor" = "#6A3DAA" 
    ),
    name = "Enriched in") +
    ggplot2::scale_x_continuous(name = expression(-log[10]~"(adjusted p-value)"),
                                limits = c(0, NA),
                                expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(y = NULL, title = title,
                  subtitle = "Numbers in circles = NES; dashed line = FDR 0.05") +
    ggplot2::theme_classic(base_size = 12)
}

p_pd_nc <- plot_both(panel_both, guides, fdr_cut = 0.05)
p_pd_nc   
ggsave("results/plots/PDvsNICM_Tcell_TOP12_WALD.png", p_pd_nc, width = 12, height = 4.0, dpi = 600)
ggsave("results/plots/PDvsNICM_Tcell_TOP12_WALD.svg", p_pd_nc, width = 12, height = 4.0, dpi = 600)

readr::write_csv(panel_both %>% dplyr::arrange(desc(order_key), contrast),
                 "results/tables/PDvsNICM_Tcell_TOP12_WALD_table.csv")


#----Supplementary enrichment GO BP Plots for PD HighCTL vs LowCTL----
suppressPackageStartupMessages({ library(dplyr); library(stringr); library(ggplot2); library(fgsea) })

stopifnot(exists("pd"), exists("msig_list"), exists("pd_gsea"))
stats_pd <- pd$ranks                              
fdr_line <- 0.05


pick_pathway <- function(gsea_df, pattern_regex){
  cand <- gsea_df %>% filter(str_detect(pathway, regex(pattern_regex, ignore_case = TRUE))) %>%
    arrange(padj)
  if (nrow(cand) == 0) stop("No pathway matches regex: ", pattern_regex)
  cand$pathway[1]
}


enrich_plot <- function(stats, pathways, res_tbl, pathway_name,
                        title = NULL, box_pos = c("topright" = TRUE)) {

  if (!pathway_name %in% names(pathways)) stop("Pathway not found in msig list: ", pathway_name)

  
  p <- fgsea::plotEnrichment(pathways[[pathway_name]], stats) +
    labs(title = title %||% stringr::str_to_sentence(gsub("_", " ", pathway_name)),
         x = "Position in the ranked list of genes",
         y = "Running enrichment score") +
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank())

  
  row <- res_tbl[match(pathway_name, res_tbl$pathway), ]
  nes  <- if (is.na(row$NES)) NA_real_ else as.numeric(row$NES)
  padj <- if (is.na(row$padj)) NA_real_ else as.numeric(row$padj)

  
  lab <- sprintf("p.adj = %s\nNES = %s",
                 ifelse(is.finite(padj), format(padj, digits = 3, scientific = TRUE), "NA"),
                 ifelse(is.finite(nes),  format(nes,  digits = 3),                  "NA"))

  
  p <- p + annotate("label", x = Inf, y = Inf, hjust = 1.02, vjust = 1.1,
                    label = lab, size = 3.2, label.size = 0.3, alpha = 0.95)

  p
}

#----By pathway ID----
pw_tcell_activation <- pick_pathway(pd_gsea, "GOBP_T_CELL_ACTIVATION")
pw_cytotox          <- pick_pathway(pd_gsea, "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY")

pw_cell_killing     <- pick_pathway(pd_gsea, "GOBP_CELL_KILLING")

#----make plots----
p_act  <- enrich_plot(stats_pd, msig_list, pd_gsea, pw_tcell_activation,
                      title = "GO: Enrichment of T cell activation")
p_cyto <- enrich_plot(stats_pd, msig_list, pd_gsea, pw_cytotox,
                      title = "GO: Enrichment of leukocyte-mediated cytotoxicity")
p_kill <- enrich_plot(stats_pd, msig_list, pd_gsea, pw_cell_killing,
                      title = "GO: Enrichment of cell killing")


p_act; p_cyto; p_kill


ggsave("results/plots/SUPP_PDHigh_vs_Low_GO_TcellActivation_enrichment.png",
       p_act, width = 4.2, height = 2.8, dpi = 300)
ggsave("results/plots/SUPP_PDHigh_vs_Low_GO_LeukocyteMediatedCytotoxicity_enrichment.png",
       p_cyto, width = 4.2, height = 2.8, dpi = 300)
ggsave("results/plots/SUPP_PDHigh_vs_Low_GO_CellKilling_enrichment.png",
       p_kill, width = 4.2, height = 2.8, dpi = 300)



pick_pathway <- function(gsea_df, pattern_regex){
  cand <- gsea_df %>% filter(str_detect(pathway, regex(pattern_regex, ignore_case = TRUE))) %>%
    arrange(padj)
  if (nrow(cand) == 0) stop("No pathway matches regex: ", pattern_regex)
  cand$pathway[1]
}

enrich_plot <- function(stats, pathways, res_tbl, pathway_name, title = NULL) {
  stopifnot(pathway_name %in% names(pathways))
  p <- fgsea::plotEnrichment(pathways[[pathway_name]], stats) +
    labs(title = title %||% str_to_sentence(gsub("_"," ", pathway_name)),
         x = "Position in the ranked list of genes",
         y = "Running Enrichment Score") +
    theme_classic(base_size = 11) +
    theme(panel.grid = element_blank())
  
  row  <- res_tbl[match(pathway_name, res_tbl$pathway), ]
  nes  <- suppressWarnings(as.numeric(row$NES))
  padj <- suppressWarnings(as.numeric(row$padj))
  lab  <- sprintf("p.adj = %s\nNES = %s",
                  ifelse(is.finite(padj), format(padj, digits = 3, scientific = TRUE), "NA"),
                  ifelse(is.finite(nes),  format(nes,  digits = 3),                  "NA"))
  
  p + annotate("label", x = Inf, y = Inf, hjust = 1.02, vjust = 1.1,
               label = lab, size = 3.2, label.size = 0.3, alpha = 0.95)
}

# ---- choose the six GO:BP pathways (best FDR match for each concept) ----
pw_cell_activation   <- pick_pathway(pd_gsea, "GOBP_CELL_ACTIVATION$")
pw_lymph_activation  <- pick_pathway(pd_gsea, "GOBP_LYMPHOCYTE_ACTIVATION$")
pw_tcell_activation  <- pick_pathway(pd_gsea, "GOBP_T_CELL_ACTIVATION$")
pw_cell_killing      <- pick_pathway(pd_gsea, "GOBP_CELL_KILLING$")
pw_cytotoxicity      <- pick_pathway(pd_gsea, "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY$")
pw_nk_immunity       <- pick_pathway(pd_gsea, "GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY$")

# ---- build plots ----
p_cell_act  <- enrich_plot(stats_pd, msig_list, pd_gsea, pw_cell_activation,
                           "GO: Enrichment of cell activation")
p_lymph_act <- enrich_plot(stats_pd, msig_list, pd_gsea, pw_lymph_activation,
                           "GO: Enrichment of lymphocyte activation")
p_tcell_act <- enrich_plot(stats_pd, msig_list, pd_gsea, pw_tcell_activation,
                           "GO: Enrichment of T cell activation")
p_killing   <- enrich_plot(stats_pd, msig_list, pd_gsea, pw_cell_killing,
                           "GO: Enrichment of cell killing")
p_cytotox   <- enrich_plot(stats_pd, msig_list, pd_gsea, pw_cytotoxicity,
                           "GO: Enrichment of leukocyte-mediated cytotoxicity")
p_nk        <- enrich_plot(stats_pd, msig_list, pd_gsea, pw_nk_immunity,
                           "GO: Enrichment of natural killer cell mediated immunity")

#----preview----
p_cell_act; p_lymph_act; p_tcell_act; p_killing; p_cytotox; p_nk
(p_cell_act | p_lymph_act | p_tcell_act) / (p_killing | p_cytotox | p_nk)

#----save panels----
ggsave("results/plots/SUPP_PDHighLow_GO_CellActivation.png",              p_cell_act,  width=4.2, height=2.8, dpi=300)
ggsave("results/plots/SUPP_PDHighLow_GO_LymphocyteActivation.png",        p_lymph_act, width=4.2, height=2.8, dpi=300)
ggsave("results/plots/SUPP_PDHighLow_GO_TcellActivation.png",             p_tcell_act, width=4.2, height=2.8, dpi=300)
ggsave("results/plots/SUPP_PDHighLow_GO_CellKilling.png",                 p_killing,   width=4.2, height=2.8, dpi=300)
ggsave("results/plots/SUPP_PDHighLow_GO_LeukocyteMediatedCytotoxicity.png", p_cytotox, width=4.2, height=2.8, dpi=300)
ggsave("results/plots/SUPP_PDHighLow_GO_NKcellMediatedImmunity.png",      p_nk,       width=4.2, height=2.8, dpi=300)

