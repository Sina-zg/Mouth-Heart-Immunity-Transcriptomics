# Mouth-Heart-Immunity-Transcriptomics
Reproducible code for the Periodontal Disease→cardiovascular study: scRNA-seq + scTCR-seq (Seurat, scRepertoire, GLIPH2) and spatial transcriptomics (10x). Includes QC, integration, clonotypes, pseudobulk DE, GSEA, and figure scripts.




Volcano Plot
 
This script generates volcano plots based on differential microbial abundance between two groups ( high vs low CD4).
It produces two visualizations: one with labels for the top 10 significant taxa, and one clean version without labels.
 
## What the Script Does
- Reads in taxa codes and corresponding taxonomy.
- Merges taxonomy with differential abundance results.
- Cleans taxonomy strings and extracts Genus and Species.
- Calculates –log10(p-value) and determines significance classification.
- Identifies top 5 upregulated and top 5 downregulated taxa.
- Builds labeled and unlabeled volcano plots.
 
## Required R Packages
- tidyverse
- ggrepel
 
## Required Input Files
### taxa_codes.csv
Must contain:
1. Taxonomy
2. Code
 
Example:
Taxonomy,Code
g__Streptococcus; s__Streptococcus_salivarius,s001
 
### rnaseq_de.csv
Differential analysis results from the RNAseq pipeline.
Required columns:
- Code
- log2FC
- Pvalues
 
Example:
Code,log2FC,Pvalues
s001,1.8,0.0005
 
## Output Files
- volcano_top10.svg
- volcano_top10_unlab.svg
 
These files will be saved to the directory specified in your script.



1. Input Files Required
Place your files in the paths used in the script or update the paths accordingly.
Upstream Processing: All taxonomic abundance tables used in this workflow were generated using the HUMAnN3 microbiome profiling pipeline. This script begins after HUMAnN3 output has been produced.
Mandatory input files:
File	Description
buglist.tsv	HUMAnN output file containing species-level (and genus-level) relative abundances for each sample.
metadata.xlsx	Spreadsheet containing sample information, including Sample IDs and SampleType (e.g., Survivors, Progressors).

2. Software & R Packages
This pipeline uses the following R packages:
●	tidyverse — data cleaning & manipulation

●	vegan — ecological diversity metrics

●	phyloseq — microbiome analysis and object handling

●	ComplexHeatmap & circlize — advanced heatmaps

●	MaAsLin2 — differential abundance testing

●	ggplot2 — all plots and visualizations

●	patchwork — combining multiple plots

To install missing packages:
install.packages(c("tidyverse", "vegan", "phyloseq", "ComplexHeatmap", 
                   "circlize", "viridisLite", "ggplot2", "reshape2", 
                   "dplyr", "scales", "forcats", "patchwork"))

# Install MaAsLin2 (Bioconductor)
BiocManager::install("Maaslin2")


3. Overview of the Pipeline
Step 1: Load and clean HUMAnN taxonomic data
●	Removes unused or failed samples.

●	Extracts species names cleanly from HUMAnN strings (removes t__ tags).

●	Corrects misclassified species (example: Porphyromonas bobii → P. gingivalis).

●	Removes duplicated species and ensures samples are clean.

Step 2: Load and clean metadata
●	Ensures sample IDs match HUMAnN IDs.

●	Cleans sample types (removes “ICM”, empty labels, and excluded samples).

●	Merges metadata with abundance data.


4. PCoA (Principal Coordinates Analysis)
This section:
1.	Converts species table into numeric format

2.	Normalizes abundance to relative abundance

3.	Calculates Bray–Curtis distance

4.	Runs PCoA

5.	Performs PERMANOVA to test significance between groups

6.	Creates a PCoA plot colored by sample group

Output example:
 PCoA_S.svg

5. Heatmap of Top 30 Species
This part uses ComplexHeatmap to generate a high-quality heatmap:
●	Extracts OTU (species) table from phyloseq

●	Applies log1p() scaling for visualization

●	Selects the top 30 most abundant taxa

●	Splits samples by group (e.g., Survivors vs Progressors)

●	Adds color annotations

Output:
 Heatmap_Top30_Taxa.svg

6. Differential Abundance Testing (MaAsLin2)
The script:
1.	Prepares metadata and abundance matrices

2.	Runs MaAsLin2 to test associations between taxa and SampleType

3.	Saves all model output in maaslin2_output/

This allows identification of species or genera significantly associated with different sample groups.

7. Custom Taxon Plots (e.g., Streptococcus sanguinis)
This section produces publication-ready boxplots that include:
●	Jittered sample points

●	Mean ± SE bars

●	Automatically generated p-value, q-value, and significance stars from MaAsLin2

●	Clean minimal theme

Output:
 Boxplot_Ss_custom.svg

8. Genus-Level Analysis
The workflow repeats species-level steps but at the genus level:
●	Cleans genus names

●	Normalizes to relative abundance

●	Creates facet boxplots for selected genera

●	Runs MaAsLin2 for genus-level differential abundance

Output:
 Boxplots_G.svg
 maaslin2_g_output/


9. Output Files Generated
File/Folder	Description
PCoA_S.svg	Species-level PCoA plot
AD_S_SHAN_SIMP.svg	Phyloseq alpha diversity (Observed, Shannon, Simpson)
Heatmap_Top30_Taxa.svg	ComplexHeatmap of top 30 species
maaslin2_output/	MaAsLin2 species-level differential abundance results
Boxplot_Ss_custom.svg	Detailed boxplot for S. sanguinis
Boxplots_G.svg	Genus-level boxplots
maaslin2_g_output/	MaAsLin2 results for genera

How to Run the Script
1.	Open RStudio or run the script in R.

2.	Update the file paths near the top (e.g., "path/to/buglist.tsv").

3.	Ensure input files are present in those paths.

4.	Run the script top to bottom.

5.	All plots and results will appear in your output directories.

----------------------------------------------------------------------

# Complex Heatmap

This script generates a **publication-ready ComplexHeatmap** showing gene expression and TCR clonal features for CD4⁺ T cells from NICM patients, comparing **Progressor** vs **Survivor** groups.

Each column is a single CD4⁺ T cell.  
Each row is a gene of interest.  
Top annotations summarize sample group, patient ID, TCR clone size, and QC metrics.

---

## What this script does

1. **Loads per-patient 10x Genomics RNA data**  
   - Uses `Read10X()` and `CreateSeuratObject()` for each patient in `sample_info`.
   - Annotates each cell with `sample_id` and `group` (Progressor/Survivor).

2. **Performs basic QC and normalization**
   - Computes mitochondrial percentage (`percent.mt`).
   - Keeps cells with:
     - `nFeature_RNA > 500`
     - `nFeature_RNA < 6500`
     - `percent.mt < 10`
   - Normalizes RNA counts with `NormalizeData()`.

3. **Merges all patients into a single Seurat object**
   - Merges count matrices and metadata across samples.
   - Creates `seurat_combined_clean` containing all QC-passing cells.
   - Renormalizes after merging.

4. **Loads and processes TCR data**
   - Reads `filtered_contig_annotations_<sample>.csv` per patient.
   - Keeps only **productive** TCRs (`productive == "true"`).
   - Constructs a unique cell barcode that matches RNA data: `sampleID_barcode`.
   - Defines clones using CDR3 amino acid sequence (`CTaa = cdr3`).
   - Computes clone size (`clone_size` = number of cells per CTaa).
   - Categorizes clones into:
     - `Singleton` (clone_size == 1 or missing)
     - `Small`     (2–10)
     - `Medium`    (11–20)
     - `Large`     (>20)
   - Adds `clone_category` and `tcr_freq` (clone size) to Seurat metadata.

5. **Restricts to CD4⁺ T cells**
   - Keeps cells with `CD3D > 0` **and** `CD4 > 0`.

6. **Selects and scales genes of interest**
   - Uses a predefined list of genes (naïve, activation, cytotoxicity, exhaustion, γδ, CD8, etc.).
   - Keeps only genes present in the RNA assay.
   - Scales these genes using `ScaleData()` and extracts the Z-score matrix from `scale.data`.

7. **Builds the annotation metadata for the heatmap**
   - For each cell (column) the script constructs:
     - `Group` (Progressor / Survivor)
     - `Patient` (sample_id)
     - `CloneSize` (Singleton / Small / Medium / Large)
     - `Mito %` (percent.mt)
     - `Read count` (nCount_RNA)
     - `# of features` (nFeature_RNA)
     - `TCR frequency` (clone size)

8. **Orders cells for visualization**
   - Columns are sorted by:
     1. `Group`
     2. `CloneSize` (Large → Medium → Small → Singleton)
     3. `TCR frequency` (descending)
     4. `Patient`

9. **Builds the ComplexHeatmap**
   - Rows: genes of interest (order fixed, no clustering).
   - Columns: individual cells.
   - Column split by `Group` (Progressor vs Survivor).
   - Top annotation includes:
     - Group
     - Patient
     - CloneSize
     - Mito %
     - Read count
     - # of features
     - TCR frequency
   - Expression values are shown as scaled Z-scores.

10. **Exports the heatmap**
    - High-resolution PNG.
    - Vector PDF.

---

## Required R packages

The script uses:

- `Seurat`
- `dplyr`
- `ComplexHeatmap`
- `circlize`
- `grid`

Make sure these are installed, for example:

```r
install.packages(c("Seurat", "dplyr"))
install.packages("circlize")
install.packages("grid")  # usually comes with base R
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap")

---------------------------------------------------------------------

Pathway Analysis

# Pathway Analysis: Pseudobulk DESeq2 + GSEA

This script takes single-cell RNA-seq data from NICM and periodontitis (PD) patients, extracts CD4⁺ T cells, aggregates them to **pseudobulk** per sample, runs **DESeq2** differential expression, and performs **GSEA** (Hallmark + GO BP) to compare:

- **PD_HighCTL vs PD_LowCTL**
- **NICM_Progressor vs NICM_Survivor**

It then generates tables and plots for main and supplementary pathway figures in the manuscript.

---

## What the script does

1. **Project setup**
   - Sets `PROJECT_DIR` and creates folders: `qc_pngs/`, `results/plots/`, `results/tables/`, `logs/`, `rds/`.
   - Starts a timestamped log file in `logs/`.

2. **Load samples & basic QC**
   - Loads 10x RNA matrices for samples:
     - NICM: `CV5, CV12, CV109, CV110` (Progressors), `CV39, CV87, CV103, CV106` (Survivors)
     - PD: `PD1–PD4` (HighCTL), `PD5–PD8` (LowCTL)
   - Computes `percent.mt`, `percent.ribo`, `percent.hb`.
   - Applies hard + adaptive QC thresholds and saves QC plots to `qc_pngs/`.

3. **Doublet removal & CD4⁺ T-cell gating**
   - Runs `scDblFinder` to keep only singlets.
   - Normalizes with `SCTransform` (glmGamPoi).
   - Gating: keeps CD3D⁺ CD4⁺ CD8A⁻ cells to define CD4⁺ T cells.
   - Writes QC summary tables to `results/tables/qc_stage*_counts.csv`.

4. **Pseudobulk creation**
   - For each sample, sums raw RNA counts across gated CD4⁺ T cells.
   - Builds a gene × sample pseudobulk matrix and metadata (`sample_id`, `group`, `cohort`).
   - Writes library sizes to `results/tables/pseudobulk_library_sizes.csv`.

5. **Differential expression (DESeq2)**
   - Builds separate DESeq2 objects for:
     - PD cohort (`PD_HighCTL` vs `PD_LowCTL`)
     - NICM cohort (`NICM_Progressor` vs `NICM_Survivor`)
   - Uses the Wald statistic as ranking metric.
   - Saves DE tables:
     - `results/tables/DE_PD_High_vs_Low.csv`
     - `results/tables/DE_NICM_Prog_vs_Surv.csv`

6. **GSEA: Hallmark + GO BP**
   - Uses `msigdbr` to collect Hallmark (H) + GO BP (C5) gene sets.
   - Runs `fgseaMultilevel` for PD and NICM contrasts.
   - Saves full GSEA results:
     - `results/tables/GSEA_PD_High_vs_Low__ALL_Hallmark_GO-BP.csv`
     - `results/tables/GSEA_NICM_Prog_vs_Surv__ALL_Hallmark_GO-BP.csv`
   - Produces:
     - Lollipop plots of top enriched pathways per contrast.
     - NES concordance scatter plot (PD vs NICM).

7. **T-cell–focused pathways and overlap**
   - Filters PD GSEA results for T cell / cytotoxic / NK-related terms.
   - Selects top 12 T-cell pathways enriched in PD HighCTL and saves:
     - `results/tables/PD_Tcell_enriched_in_PDHigh_FULL.csv`
     - `results/plots/PD_Tcell_TOP12_PDhigh_only_minuslog10FDR.png`
   - Re-extracts the same 12 pathways for NICM Progressors and builds a 2-condition panel:
     - `results/plots/PDvsNICM_Tcell_TOP12_WALD.png` / `.svg`
     - `results/tables/PDvsNICM_Tcell_TOP12_WALD_table.csv`

8. **Supplementary GO BP enrichment plots**
   - Generates GSEA enrichment curves for key GO BP terms:
     - Cell activation, lymphocyte/T cell activation, cell killing, leukocyte-mediated cytotoxicity, NK cell–mediated immunity.
   - Saves them to `results/plots/SUPP_PDHighLow_GO_*.png` and `results/plots/SUPP_PDHigh_vs_Low_GO_*_enrichment.png`.

---

## Inputs you must set

At the top of the script:

- `PROJECT_DIR <- "Path to project"`  
- `rna_base    <- "Path to RNA Path"`

`rna_base` must point to folders like:

- `sample_filtered_feature_bc_matrix_CV5/`
- `sample_filtered_feature_bc_matrix_PD1/`
- etc.

---

## How to run

1. Install required R packages (`Seurat`, `SingleCellExperiment`, `scDblFinder`, `glmGamPoi`, `DESeq2`, `fgsea`, `msigdbr`, `dplyr`, `ggplot2`, `readr`, `stringr`, `forcats`, `purrr`).
2. Edit `PROJECT_DIR` and `rna_base` to match your directory structure.
3. In R/RStudio:

   ```r
   setwd(PROJECT_DIR)
   source("Pathway_Analysis_Pseudobulk_GSEA.R")
------------------------------------------------------------------------------------------

GLIPH2 Analysis

# GLIPH2 Analysis: Clonotype Preparation, Motif Filtering, and Visualization

This script prepares TCR β–chain clonotype tables for **GLIPH2**, filters GLIPH2 motif outputs with patient- and score-based criteria, and generates publication-ready figures showing motif sharing across **PD HighCTL / PD LowCTL** and **NICM Progressor / Survivor** groups.

---

## Inputs (edit these paths)

- 10x `filtered_contig_annotations` CSV files for each sample:

  ```r
  paths <- list(
    PD_High = c("Sample PD1 filtered_contig_annotations path.csv", ...),
    PD_Low  = c("Sample PD5 filtered_contig_annotations path.csv", ...),
    Prog    = c("Sample NICM1 filtered_contig_annotations path.csv", ...),
    Surv    = c("Sample NICM5 filtered_contig_annotations path.csv", ...)
  )

GLIPH2 result files (one row per motif/pattern):
run1_file <- "Path to gliph2_run1.csv"   # Prog + PD High
run2_file <- "Path to gliph2_run2.csv"   # Surv + PD Low
run3_file <- "Path to gliph2_run3.csv"   # (all samples, for some plots)
res_dir   <- "Path to GLIPH2_results"
fig_dir   <- "Path to figures"


What the script does

Load and standardize contig files

Reads all filtered_contig_annotations_*.csv files for PD_High, PD_Low, Prog, Surv.
Standardizes column names to a 10x-like schema (barcode, chain, cdr3_aa, v_gene, j_gene, productive, reads, umis).
Filters to productive chains only.

Resolve one TRA and one TRB per cell
Within each (sample_id, condition, barcode, chain) group, selects the contig with highest UMI/read support.

Keeps only TRA/TRB and pivots to one row per cell with:
CDR3a, CDR3b, Va, Ja, Vb, Jb (alleles stripped to gene family).

Build clonotype table and GLIPH2 inputs

Motif-level postprocessing of GLIPH2 runs

Loads GLIPH2 outputs (gliph2_run1.csv, gliph2_run2.csv), requiring at least:
pattern (motif), Sample (donor ID), Fisher_Score, Length_Score (for later filtering; names auto-detected)

Classifies each Sample as:
PD (periodontitis, PD1–PD8)
CV (cardiac/NICM, CV samples)

For each run, counts how many distinct donors contribute each motif in PD and CV, and keeps motifs:
Present in ≥ 2 PD donors and ≥ 2 CV donors
With Fisher_Score < 0.05 and Length_Score < 0.05

Writes filtered motif tables to res_dir, e.g.:
run1_motifs_PDHigh>=2_AND_Prog>=2_scores<0.05.csv
run2_motifs_PDLow>=2_AND_Surv>=2_scores<0.05.csv

Motif overlap between groups
Defines two motif sets:
Group1 (Prog + PD High) from run1
Group2 (Surv + PD Low) from run2

Computes:
motifs unique to Group1
motifs unique to Group2
motifs shared between both

Saves lists in res_dir as:
motifs_unique_to_Group1.txt
motifs_unique_to_Group2.txt
motifs_shared_Group1_and_Group2.txt

Frequency-based motif plots (4-layer bar figure)
Converts GLIPH2 run1/run2 outputs into a unified table with:
Subject (donor), condition (Prog / Surv / PD_High / PD_Low), pattern, Frequency

Maps conditions to panels:
Progressor, PD High CTL, PD Low CTL, Survivor

Selects ordered motif patterns:
Left block: motifs enriched in Group1 only
Middle block: shared motifs
Right block: motifs enriched in Group2 only
Builds a 4-row bar plot (pseudo-log2 Y axis) for:
Progressor, PD High CTL, PD Low CTL, Survivor

Saves:
GLIPH2_4layers_FreqLOG2_noTitles.png in fig_dir.
Motif × donor heatmaps
Constructs a donor-level matrix: number of TCRs per (Subject, pattern) for:
Group1 (Progressor + PD High CTL)
Group2 (PD Low CTL + Survivor)
Fills missing combinations with zero and caps fill at 4 TCRs for visualization.
Creates two stacked heatmaps (Group1, Group2) aligned to the same motif order, with vertical dashed lines marking the three motif blocks (Group1-only / shared / Group2-only).
Combines heatmaps + 4-layer bar panels using patchwork and saves:
GLIPH2_heatmaps.png in fig_dir.

-------------------------------------------------------------------------------------
MYH6 / DUF1002 Stimulation – scRNA-seq + scTCR Analysis

This repository contains an end-to-end R pipeline to analyze 10x Genomics 5′ scRNA-seq and scTCR-seq data from CD4⁺ T cells stimulated with the cardiac peptide MYH6, the oral mimic DUF1002, or left unstimulated. The workflow integrates RNA and TCR data across multiple patients, identifies clonotypes expanded upon stimulation, maps them on UMAPs, performs differential expression within expanded clones, and generates publication-grade figures and statistics. 

1. Pipeline Overview
The main script performs the following steps:

Load libraries
Uses Seurat, scRepertoire, data.table, ggplot2, patchwork, and helper utilities for plotting and statistics.

Load and preprocess RNA data
Expects 10x filtered_feature_bc_matrix directories for four patients (CV5, CV12, CV109, CV110), each with three conditions: Unstim, MYH6, and DUF1002.
Paths are specified in the rna_paths list.

For each entry, make_seu():
Reads the 10x matrix with Read10X.
Creates a Seurat object with CreateSeuratObject.
Annotates sample_id, patient ID (e.g., CV5), and Condition (Unstim / MYH6 / DUF1002) from the sample name.

Prefixes cell barcodes with sample_id so they remain unique across samples. 
Per-sample normalization and integration (SCT)
Performs SCTransform on each sample (v2 flavor, 3000 variable features). 

Selects integration features and prepares objects for SCT-based integration.
Uses unstimulated samples (CV5_Unstim, CV12_Unstim, CV109_Unstim, CV110_Unstim) as integration references.
Runs FindIntegrationAnchors and IntegrateData with SCT normalization. 

Dimensional reduction and clustering
Sets the default assay to "integrated".
Runs PCA (40 PCs), UMAP (dims 1–30), neighbor graph construction, and Louvain clustering (resolution 0.6). 

Load and harmonize TCR data
Expects 10x filtered_contig_annotations.csv-style files for each sample, paths defined in contig_paths. 

Reads each CSV (barcode, TCR chains, CDR3s, etc.), adds sample_id, and prefixes barcodes to match Seurat cell names.
Uses scRepertoire combineTCR() to build a combined.TCR object with clonotype information across all samples/conditions. 

Barcode matching and QC
Rebuilds TCR barcodes to exactly match Seurat cell names (sample prefix + core barcode).
Checks which sample_ids overlap between Seurat and combined.TCR.
Summarizes per-sample: number of TCR rows, unique barcodes, CDR3s with non-NA CTaa, number of mapped cells, and % cells with TCRs. 

Clone definition and expansion across conditions
Extracts TRA/TRB amino-acid sequences per row, infers missing chains, and assembles a normalized clone label (Clone). 

Counts clone occurrences per sample and aggregates them per Condition (Unstim, MYH6, DUF1002).
Builds a wide table (per_cond_wide) with counts per clone per condition. 

Identifies expanded clones present in all three conditions but increased in both MYH6 and DUF1002 vs Unstim (shared antigen-responsive clones). 

Alluvial plots of clonal expansion
For the set of expanded clones, builds a long table and orders clones by DUF1002 + MYH6 counts.
Uses ggplot2 (with alluvial geometry) to visualize how clone size changes across Unstim → MYH6 → DUF1002 for each clone, saving high-resolution PNG and PDF outputs. 

UMAP projection of expanded clones (CD4⁺ background)
Infers CD4⁺ T cells either from existing annotations or via CD3D/CD4/CD8A expression. 

Overlays the 27 (or top N) expanded clones on the CD4 UMAP separately for Unstim, MYH6, and DUF1002.
Colors points by expansion bin (singleton / small / medium / large) and saves condition-specific UMAPs as PNG and PDF. 

Cluster-level mapping of expanded-clone cells
Determines which Seurat clusters contain at least one of the expanded-clone cells.
Replots the integrated UMAP, highlighting only those clusters and optionally focusing on a subset of “target” clusters enriched for expanded clones. 

Definition of CD4⁺ cytotoxic clusters
Restricts analysis to CD4⁺ cells.
Computes detection frequency and mean expression for cytotoxic markers GZMB, PRF1, NKG7 by cluster.
Calculates a cytotoxicity score per cluster (Z-score-averaged markers) and flags clusters with high score and high marker detection as CD4 CTL clusters. 

Plots a UMAP of CD4 cells, highlighting CTL-enriched clusters in navy. 

Differential expression in expanded-clone cells
Subsets the Seurat object to cells belonging to expanded clones (e.g., the “27 clones”), creates seu_27.
Ensures RNA assay has a "data" layer (joins layers or runs log-normalization as needed).
Runs FindMarkers (Wilcoxon test) for MYH6 vs Unstim and DUF1002 vs Unstim, focusing on genes expressed in ≥5% of cells with log2FC ≥ 0.1. 

Extracts shared upregulated genes (FDR < 0.05, log2FC > 0 in both comparisons) and previews them in a bar plot ordered by mean log2FC. 

Per-gene summary plots with statistics
For a curated panel of TCR-proximal / activation genes (e.g., JUN, DUSP4, DUSP16, TNFAIP3, BTG1, PDCD4, LDHA, HIF1A, STK17A, IL7R, NFKBIA), the pipeline: 

Aggregates expression per clone and condition (mean log1p(count) per clone).
Computes per-gene paired Wilcoxon tests (MYH6 vs Unstim, DUF1002 vs Unstim) plus rank-biserial effect sizes. 

Generates compact bar + dot plots with overlaid p-values, saving PNG/PDF for each gene. 

Optionally exports a CSV with per-gene statistics (p-values, effect sizes, number of pairs). 

2. Input Data and Directory Structure
Before running the script, you must edit the paths for your data.
2.1 RNA (gene expression) input
The script expects 10x Genomics 5′ filtered_feature_bc_matrix outputs organized as: 

rna_paths <- list(
  CV5_Unstim   = "/path/to/CV5_Unstim/filtered_feature_bc_matrix",
  CV5_MYH6     = "/path/to/CV5_MYH6/filtered_feature_bc_matrix",
  CV5_DUF1002  = "/path/to/CV5_DUF1002/filtered_feature_bc_matrix",
  CV12_Unstim  = "/path/to/CV12_Unstim/filtered_feature_bc_matrix",
  CV12_MYH6    = "/path/to/CV12_MYH6/filtered_feature_bc_matrix",
  CV12_DUF1002 = "/path/to/CV12_DUF1002/filtered_feature_bc_matrix",
  CV109_Unstim = "/path/to/CV109_Unstim/filtered_feature_bc_matrix",
  CV109_MYH6   = "/path/to/CV109_MYH6/filtered_feature_bc_matrix",
  CV109_DUF1002= "/path/to/CV109_DUF1002/filtered_feature_bc_matrix",
  CV110_Unstim = "/path/to/CV110_Unstim/filtered_feature_bc_matrix",
  CV110_MYH6   = "/path/to/CV110_MYH6/filtered_feature_bc_matrix",
  CV110_DUF1002= "/path/to/CV110_DUF1002/filtered_feature_bc_matrix"
)

Each directory should contain the standard 10x files:
barcodes.tsv.gz
features.tsv.gz (or genes.tsv.gz)
matrix.mtx.gz

2.2 TCR (contig) input
The script expects one TCR contig CSV per sample, preferably the 10x filtered_contig_annotations.csv or an equivalent file. Paths are specified in contig_paths: 

contig_paths <- list(
  CV5_Unstim   = "/path/to/CV5_Unstim/filtered_contig_annotations.csv",
  CV5_MYH6     = "/path/to/CV5_MYH6/filtered_contig_annotations.csv",
  CV5_DUF1002  = "/path/to/CV5_DUF1002/filtered_contig_annotations.csv",
  CV12_Unstim  = "/path/to/CV12_Unstim/filtered_contig_annotations.csv",
  CV12_MYH6    = "/path/to/CV12_MYH6/filtered_contig_annotations.csv",
  CV12_DUF1002 = "/path/to/CV12_DUF1002/filtered_contig_annotations.csv",
  CV109_Unstim = "/path/to/CV109_Unstim/filtered_contig_annotations.csv",
  CV109_MYH6   = "/path/to/CV109_MYH6/filtered_contig_annotations.csv",
  CV109_DUF1002= "/path/to/CV109_DUF1002/filtered_contig_annotations.csv",
  CV110_Unstim = "/path/to/CV110_Unstim/filtered_contig_annotations.csv",
  CV110_MYH6   = "/path/to/CV110_MYH6/filtered_contig_annotations.csv",
  CV110_DUF1002= "/path/to/CV110_DUF1002/filtered_contig_annotations.csv"
)

Important: Barcodes in these CSVs must match the 10x RNA barcodes for each sample (the script then adds sample-specific prefixes to match Seurat cells).

3. Software Requirements
R ≥ 4.1 (recommended)
Core packages:
Seurat, SeuratObject
scRepertoire
dplyr, tidyr, data.table
ggplot2, patchwork, scales
ggalluvial, ggnewscale (for some plots; check missing-package errors)
The script begins with:
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

You may need to install any missing packages with install.packages() or BiocManager::install() where appropriate. 

4. How to Run
Clone or copy this repository
git clone <your_repo_url>
cd <your_repo>
Edit paths at the top of the script
Set rna_paths to the actual locations of your 10x RNA matrices.
Set contig_paths to the actual locations of your TCR contig CSVs.
Update any "Path to Save/...png" or "...pdf" paths to desired output directories.
Start R / RStudio and source the script
From an R session:
source("GitHub MYH6:DUF1002 Stim.R")

The script runs sequentially through data loading, integration, TCR mapping, clone expansion analysis, DE tests, and figure generation.

Inspect outputs
Check the console summary tables for:
Per-sample mapping of TCRs to Seurat cells.
Number of expanded clones.
Top shared upregulated genes in MYH6/DUF1002 vs Unstim.

Review PNG/PDF files generated for:
Alluvial plots of clone counts across conditions.
UMAPs showing expanded clones and CTL clusters.
Per-gene bar + dot plots with p-values.

5. Key Outputs
Depending on which sections you run and your chosen output paths, the script will generate:
TCR_alluvial_plot.png / .pdf
– Alluvial plot of clone sizes across Unstim / MYH6 / DUF1002. 

UMAP_CD4_27clones_<Condition>.png / .pdf
– UMAPs of CD4⁺ cells with expanded clones colored by expansion category (singleton / small / medium / large). 

UMAP_selected_clusters_no_annotations.png / .pdf
– Integrated UMAP with clusters harboring expanded clones highlighted. 

UMAP_CD4_cytotoxic_clusters_navy.png / .pdf
– CD4⁺ UMAP highlighting cytotoxic clusters enriched for GZMB/PRF1/NKG7. 

BAR_<GENE>.png / .pdf
– Per-gene bar + dot plots (e.g., BAR_JUN.png, BAR_NFKBIA.png) with p-values for MYH6 vs Unstim and DUF1002 vs Unstim. 

Optional CSV:
A statistics table summarizing p-values and rank-biserial effect sizes for each gene in the panel. 

6. Customization
You can adjust several aspects of the analysis:
Conditions and sample naming
The script infers Condition from the suffix in sample_id (e.g., _Unstim, _MYH6, _DUF1002). If you add new conditions, extend the regex and factor levels accordingly. 

Clonal expansion criteria
By default, “expanded” clones are present in all three conditions and have higher counts in both MYH6 and DUF1002 vs Unstim. Modify the filtering on per_cond_wide to use different thresholds or criteria. 

Number of clones visualized
The script often uses TOP_N or the “top 27” clones for visualization. This can be changed to include more or fewer clones. 

Gene panels
The tcr_direct_genes vector defines the genes for the bar + dot plots; you can add or remove genes as needed. 

Cytotoxic cluster thresholding
You can modify the cytotoxicity score quantile (e.g., 75th percentile) and the minimum detection percentage for GZMB/PRF1/NKG7 used to define enriched CD4 CTL clusters. 

7. Notes and Limitations
The script assumes standard 10x barcoding conventions and that each sample_id in rna_paths has a corresponding TCR file in contig_paths.
Helper objects such as cell_map and helper functions like fix_contig_cols() must be present in the environment if defined elsewhere in your analysis pipeline (e.g., another source file); ensure they are sourced before running sections that depend on them. 

The workflow is tuned for CD4⁺ T-cell biology under MYH6 and DUF1002 peptide stimulation; adapting it to other systems may require updating marker genes, conditions, and QC criteria.

-------------------------------------------------------------------------------------------------------------------

Integrated scRNA-seq + scTCR-seq CD4 T-cell Pipeline

This repository contains an end-to-end R workflow for quality control, integration, annotation, and TCR clonality analysis of 10x Genomics CD4⁺ T cells from:
NICM Progressors vs NICM Survivors

Periodontitis (PD) High-CTL vs PD Low-CTL patients

scRNATCRseq

The pipeline produces:

QC metrics and plots for each sample

An integrated, annotated CD4⁺ T-cell UMAP

Cluster marker tables (top markers per cluster)

Population composition plots (stacked bars) by group

A CD4 CTL “dot plot” comparing Progressor vs PD HighCTL

TCR clonotype calling (TRA+TRB), clone sizes and expansion bins

Clonality barplots and UMAPs colored by clonal expansion

1. Overview of the Workflow
1.1 Sample definitions

The script works with 16 samples:
NICM (CV*): CV5, CV12, CV109, CV110 (Progressors) and CV39, CV87, CV103, CV106 (Survivors)
PD (PD*): PD1, PD2, PD3, PD4 (High CTL) and PD5, PD6, PD7, PD8 (Low CTL)

Each sample is mapped to a clinical group via group_map and to a cohort (NICM vs PD).

scRNATCRseq

2. RNA Pipeline
2.1 Input RNA data

The script expects 10x Genomics filtered_feature_bc_matrix directories for each sample, organized as:
rna_base <- "/Path to Samples/"   # EDIT THIS
rna_path <- function(s) file.path(rna_base, paste0("sample_filtered_feature_bc_matrix_", s))
You must edit rna_base to point to your data. Each directory should contain:

matrix.mtx.gz
barcodes.tsv.gz
features.tsv.gz
 
scRNATCRseq

2.2 Per-sample QC and filtering

For each sample:
Create Seurat object
CreateSeuratObject(mtx, project = s, min.features = 100, min.cells = 3)
Add metadata: sample_id, group, cohort (NICM vs PD). 

Compute QC metrics
percent.mt (mitochondrial), percent.ribo (ribosomal), percent.hb (hemoglobin) using PercentageFeatureSet. 

QC thresholds
Keep cells with:
nFeature_RNA > 500
nFeature_RNA < 6500
percent.mt < 10

Adaptive outlier removal
Use median + k·MAD for:
nFeature_RNA (k = 4)
nCount_RNA (k = 4)
percent.mt capped at 10% (MAD-based but not exceeding 10%)

QC plotting
plot_qc() produces per-sample PNGs (qc_pngs/<sample>_QC.png) with:
Violin plots for nFeature_RNA, nCount_RNA, percent.mt
Scatter plots: nCount vs nFeature, %MT vs nFeature 

QC summary tables
qc_stage1_counts.csv (raw, after hard filter, after adaptive filter, MT cap used)

2.3 Doublet removal
Convert each Seurat object to SingleCellExperiment, run scDblFinder() with dbr = 0.04.
Keep only cells labeled singlet.
Update QC table to qc_stage_all_counts.csv (includes counts after doublet removal and after downstream gating). 

2.4 SCTransform normalization
For each singlet-only object, run:
SCTransform(
  seu,
  assay = "RNA",
  variable.features.n = 3000,
  vars.to.regress = "percent.mt",
  method = "glmGamPoi",
  verbose = FALSE
)
Record remaining cells and append to QC table (after_SCT). 

2.5 CD4⁺ T-cell gating

The function gate_cd3d_cd4_drop_cd8a() gates CD4 T cells at the expression level:
Default assay set to SCT
Keep cells with:
CD3D > 1
CD4 > 1
CD8A ≤ 1 (or missing)

Results:
Per-sample CD4-only objects in gated_list

Final QC table qc_stage_all_counts.csv with "after_gate" and "removed_by_gate" columns
An .rds file containing all QC-clean CD4-only, SCT-normalized objects:
saveRDS(gated_list, "qc_clean_per_sample_SCT_CD4only.rds")
 
---

## 3. Integration, Clustering, and Annotation

### 3.1 Load QC-clean objects


gated_list <- readRDS("Path to qc_clean_per_sample_SCT_CD4only.rds")  # EDIT

3.2 Remove TCR/Ig genes from integration features

For each object:
Default assay: "SCT"
Variable features are filtered to remove TCR/IG genes: ^(TR[ABDG]|IG[HKL])
This keeps integration anchors driven by transcriptional biology rather than clonal identity. 

3.3 SCT integration (RPCA)
SelectIntegrationFeatures(object.list = gated_list, nfeatures = 3000)
PrepSCTIntegration(...)
Run PCA on each object (30 PCs)
Choose the two largest samples as reference (ref_idx)

FindIntegrationAnchors(..., normalization.method="SCT", reduction="rpca")
IntegrateData(..., normalization.method="SCT")

3.4 Dimensionality reduction and clustering

On the integrated object:
RunPCA (30 PCs)
FindNeighbors using "pca" and dims 1–30 → "pca_snn" graph
FindClusters at resolution 0.3
RunUMAP on PCA (dims 1–30)
The integrated assay used for downstream analyses is "SCT"; mitochondrial % is recomputed if missing. 

3.5 Cluster markers

PrepSCTFindMarkers(integrated)
FindAllMarkers with:
only.pos = TRUE
min.pct = 0.25
logfc.threshold = 0.25
Filter out:
TCR/Ig (TR*, IG*)
Mitochondrial (^MT-)
Ribosomal (^RP[SL])

Outputs:
markers_all_clusters_SCT.csv – all filtered markers
markers_top20_per_cluster.csv – top 20 (by log2FC) per cluster
cluster3_top50.csv – detailed markers for cluster 3 (unfiltered logFC threshold, min.pct 0.1) 

3.6 Manual annotation of clusters
Clusters from pca_snn_res.0.3 are mapped to functional labels:
0 → CCR4⁺ memory
1 → Naive/CM
2 → Th1 EM
3 → Treg
4 → CD4 CTL
5 → Activated memory


3.7 Annotated UMAPs

The script generates:

Overall annotated UMAP
UMAP split by cohort (NICM vs PD)
UMAP split by group (NICM_Progressor, NICM_Survivor, PD_HighCTL, PD_LowCTL)
A custom color palette is used for major CD4 subsets (Naive/CM, CCR4+ memory, Th1 EM, Treg, CD4 CTL, Activated memory). 

4. Population Composition Plots
4.1 PD HighCTL vs PD LowCTL

From integrated@meta.data, filter:
cohort == "PD"
group %in% c("PD_HighCTL", "PD_LowCTL")
Reorder annotation and plot stacked bar proportions of CD4 subsets per group:
X-axis: CD4 subset (Naive/CM, Th1 EM, CCR4+ memory, Treg, Activated memory, CD4 CTL)
Fill: PD_HighCTL (red) vs PD_LowCTL (navy)
Y-axis: fraction of cells

Outputs:
PD_High_vs_LowCTL_CD4pop_proportions.png
PD_High_vs_LowCTL_CD4pop_proportions.svg 

4.2 NICM Progressor vs Survivor

Similarly for NICM:
cohort == "NICM"
group %in% c("NICM_Progressor", "NICM_Survivor")
Compute per-cluster proportions and plot as stacked bars:
NICM_Progressor_vs_Survivor_CD4pop_proportions.png/.svg

5. CD4 CTL Signature: Progressor vs PD HighCTL
The pipeline focuses on CD4 CTL cells:
Subset integrated object:
cd4ctl <- subset(integrated, subset = annotation == "CD4 CTL" &
                   group %in% c("NICM_Progressor","PD_HighCTL"))
cd4ctl$plot_group <- recode(cd4ctl$group,
                            "NICM_Progressor"="Progressor",
                            "PD_HighCTL"="PD")
Gene panel (order preserved in plotting):
gene_order <- c("CD28","ZEB2","TOX","LAG3","PDCD1",
                "SELL","CD38","GZMA","IL2RA",
                "CTLA4","CD52","CD69","GZMK","GNLY",
                "GZMH","GZMM","NKG7","EOMES","IFNG","PRF1","GZMB")
Compute:
Average expression per group (AverageExpression)
% of cells expressing each gene per group
Row Z-scores of group averages
Generate a dot plot:
X-axis: Progressor vs PD
Y-axis: genes (reversed order for heatmap-like axis)
Dot size: % cells expressing
Dot color: Z-scored average expression

Outputs:
DotPlot_CD4CTL_Progressor_vs_PD_rowZ.png/.svg
Console print of Pearson r between Progressor and PD mean expressions across this gene set. 

6. TCR Processing and Clonality
6.1 TCR input

The script expects 10x filtered_contig_annotations_*.csv files:
tcr_base <- "Path to Samples TCRs"  # EDIT THIS
tcr_files <- file.path(tcr_base, paste0("filtered_contig_annotations_", samples, ".csv"))
Each file is read and normalized to standard column names (barcode, chain, cdr3_aa, productive, etc.). Only:
is_cell == TRUE
productive == TRUE
Optionally high_confidence == TRUE and full_length == TRUE are retained. 

6.2 One TRA + one TRB per cell
For each sample:
Group by barcode and chain, then keep the top row by umi_count (and reads as tiebreaker).
Spread to wide format with columns TRA_aa, TRB_aa.
Require both chains to be present.

Define:
sample_id = sample name
sample_barcode = sample_id_barcode
CTaa = "<TRA_aa>|<TRB_aa>"
All samples are merged into clono_tbl. 

6.3 Matching TCRs to RNA cells
The script:
Parses RNA cell names to extract 10x barcodes and infer sample_id (from prefixes and/or orig.ident/sample metadata).
Normalizes RNA cell names to "sample_BARCODE-1" format.
Matches clono_tbl$sample_barcode to colnames(integrated) and attaches:
CTaa
TRA_aa / TRB_aa
clone_size (frequency of CTaa across all cells)

expansion_category (categorical bin of clone_size):
Singleton (1)
Small (2–10)
Medium (11–20)
Large (>20)
NA for cells without TCR

6.4 Group and cohort labels
group and cohort are reconstructed if missing.
group_simple maps detailed groups to simplified labels:
NICM_Survivor → Survivor
NICM_Progressor → Progressor
PD_HighCTL → PD HighCTL
PD_LowCTL → PD LowCTL
The script prints diagnostic tables for cohort, group, group_simple, and expansion_category distribution. 

6.5 Clonality composition barplots
Using helper functions mk_df_group and mk_df_cluster:

NICM:
Clonal expansion composition per group (Survivor vs Progressor)
Clonal expansion per cluster (or annotation), with Progressor vs Survivor facets

PD:
Clonal expansion composition per group (PD LowCTL vs PD HighCTL)
Clonal expansion per cluster, faceted by PD group
Expansion categories are colored with a dedicated palette and plotted as proportional stacked bars (0–100%). PNG + SVG are saved via save_both() into out_dir (Path to Save files — edit this). 

6.6 UMAPs colored by clonal expansion
For each group_simple in:
Survivor
Progressor
PD LowCTL
PD HighCTL

the script:
Subsets the integrated object to that group.
Builds a data frame with UMAP embeddings, expansion_category, group_simple.
Removes cells without TCR or with expansion_category == "NA".
Plots UMAP with points colored by expansion category; axis limits fixed to full-cohort ranges.
Each UMAP is saved as PNG/SVG, e.g.:

UMAP_Survivor_expansion.png
UMAP_Progressor_expansion.png
UMAP_PD_LowCTL_expansion.png
UMAP_PD_HighCTL_expansion.png 

7. Requirements
R ≥ 4.1 recommended
CRAN packages:
Seurat, dplyr, ggplot2, Matrix, purrr, glmGamPoi, future, scales, tibble, tidyr, plyr
Bioconductor:
SingleCellExperiment, scDblFinder
Suggested:
svglite (for SVG export)

The script attempts to install missing packages at the top using install.packages() and BiocManager::install().

8. How to Use This Pipeline
Place the R script in your project directory.
Edit paths:
rna_base (RNA matrices)
tcr_base (TCR contigs)
Paths for readRDS("...qc_clean_per_sample_SCT_CD4only.rds") and out_dir for saving plots.
Run the script in R or RStudio.
-----------------------------------------------------------------------------------------------------------
