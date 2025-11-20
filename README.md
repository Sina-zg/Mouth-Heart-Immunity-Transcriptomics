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
