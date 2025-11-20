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
Complex Heatmap – CD4⁺ T-cell Gene × Clonotype Landscape

This script builds a publication-grade ComplexHeatmap of CD4⁺ T cells from NICM patients, integrating gene expression and TCR clonality to contrast Progressor vs Survivor groups. It loads per-patient 10x RNA (filtered_feature_bc_matrix) and TCR (filtered_contig_annotations_*.csv) data, performs basic QC (nFeature/nCount/percent.mt) and normalization, merges samples into a single Seurat object, and restricts the analysis to CD3D⁺CD4⁺ cells. TCR contigs are processed to define CDR3-based clonotypes (CTaa), compute clone size per clonotype, and assign each cell a clone-size category (Singleton, Small, Medium, Large), which is added to the Seurat metadata along with standard QC metrics.

A curated gene panel (naïve, activation, cytotoxic, exhaustion, etc.) is then selected, scaled (Z-scores) and extracted as a cell × gene matrix. Columns (cells) are ordered by group (Progressor vs Survivor), clone-size category, TCR frequency, and patient ID, while rows (genes) retain a fixed biological order. The script assembles a ComplexHeatmap with genes as rows and single cells as columns, split by clinical group and decorated with rich top annotations (group, patient, clone size, mitochondrial %, read depth, feature counts, TCR frequency). The result is exported as both high-resolution PNG and vector PDF, providing a compact, interpretable view of how CD4⁺ T-cell phenotypes and clonal expansion differ between Progressors and Survivors.

---------------------------------------------------------------------

Pathway Analysis – Pseudobulk DESeq2 + GSEA

This script performs pathway-level analysis of CD4⁺ T-cell transcriptional programs in NICM and periodontitis cohorts using a pseudobulk DESeq2 + GSEA framework. Starting from single-cell RNA-seq (10x) for NICM Progressor/Survivor and PD HighCTL/LowCTL patients, it applies QC and CD4⁺ T-cell gating, collapses raw counts to pseudobulk profiles per sample, and fits DESeq2 models to compare PD_HighCTL vs PD_LowCTL and NICM_Progressor vs NICM_Survivor. Wald statistics from these contrasts are then used as ranking metrics for GSEA against Hallmark and GO Biological Process gene sets, generating full enrichment tables, summary plots, and cross-cohort NES concordance.

Downstream, the script focuses on T cell– and cytotoxicity-related pathways, extracting a core panel of immune signatures enriched in PD HighCTL and re-evaluating the same pathways in NICM Progressors to visualize shared and distinct programs (e.g., via two-condition bar/scatter panels). It also produces supplementary enrichment curves for key GO BP terms (T-cell activation, cell killing, leukocyte-mediated cytotoxicity, NK cell immunity) for inclusion in main and supplementary figures. To use it, set the project and RNA base directories, ensure the required R packages (Seurat, scDblFinder, DESeq2, fgsea, msigdbr, etc.) are installed, and run the script end-to-end in R/RStudio; all QC metrics, DE tables, and GSEA results are written to organized results/ subdirectories.

------------------------------------------------------------------------------------------

GLIPH2 Analysis – Motif Filtering and Cross-Cohort Visualization

This script prepares 10x scTCR-seq data for GLIPH2, filters GLIPH2 motif outputs with stringent donor- and score-based criteria, and visualizes convergent TCR motifs across NICM and periodontitis cohorts. It reads filtered_contig_annotations files for PD HighCTL / PD LowCTL and NICM Progressor / Survivor samples, standardizes the TCR tables, and resolves one high-confidence TRA and TRB chain per cell. From these, it builds GLIPH2-ready clonotype inputs (CDR3β with V/J gene calls) and, after GLIPH2 is run externally, re-imports the GLIPH2 CSV outputs for downstream analysis.

GLIPH2 motif tables (e.g., Prog+PD High, Surv+PD Low, and an all-samples run) are then filtered to retain motifs supported by multiple independent donors in each group and passing Fisher_Score and Length_Score thresholds (e.g., < 0.05). The script classifies samples as cardiac (NICM) or periodontal (PD), constructs motif sets for Group1 (Progressor + PD HighCTL) and Group2 (Survivor + PD LowCTL), and derives motifs unique to each group as well as those shared. It summarizes motif usage by donor and condition, generating publication-ready 4-layer bar plots (motif frequencies per condition on a pseudo-log scale) and aligned motif×donor heatmaps for the two groupings.

To use the script, point the input paths to your filtered_contig_annotations files and GLIPH2 result CSVs, set output directories for results and figures, and run the file in R/RStudio from top to bottom. Sample grouping, score cutoffs, and plotting aesthetics are defined in a few central variables and can be easily tuned for alternative cohorts, GLIPH2 runs, or motif-selection criteria.

-------------------------------------------------------------------------------------
MYH6 / DUF1002 Stimulation – scRNA-seq + scTCR-seq Pipeline

This script analyzes 10x scRNA-seq and scTCR-seq data from CD4⁺ T cells stimulated with the cardiac peptide MYH6, the oral mimic DUF1002, or left unstimulated across multiple donors. It reads per-sample 10x filtered_feature_bc_matrix and filtered_contig_annotations files, performs QC and SCT-based integration of all conditions, and builds a unified CD4⁺ T-cell atlas with UMAP and clustering. TCR contigs are processed to assign one TRA and one TRB chain per cell, define amino acid clonotypes (CTaa), and quantify clone sizes across Unstim, MYH6, and DUF1002.

Using these integrated RNA–TCR profiles, the pipeline identifies antigen-responsive clonotypes that are present in all three conditions but show expansion under MYH6 and DUF1002 stimulation. It generates alluvial plots of clonal abundance across conditions, UMAP overlays highlighting expanded clones within CD4 CTL–enriched clusters, and differential expression analyses restricted to cells from expanded clones (MYH6 vs Unstim and DUF1002 vs Unstim). A curated activation/TCR-proximal gene panel is then summarized with compact bar/dot plots and paired statistics, providing a focused transcriptional signature of peptide-responsive CD4 CTLs.

To run the script, set the RNA and TCR path lists to your 10x outputs, choose an output directory, and execute the file in R/RStudio from top to bottom. Sample names, condition labels (Unstim/MYH6/DUF1002), marker panels, DE thresholds, and clone-selection criteria are all defined in a small number of vectors/filters and can be easily adapted to different donors, peptides, or stimulation designs.

-------------------------------------------------------------------------------------------------------------------
scRNA-seq + scTCR-seq CD4 T-cell Pipeline

This script performs an integrated analysis of 10x scRNA-seq and scTCR-seq data from CD4⁺ T cells across NICM Progressor vs Survivor patients and periodontitis HighCTL vs LowCTL groups. It takes 10x filtered_feature_bc_matrix and filtered_contig_annotations files per sample, applies stringent QC (mitochondria, features, counts, doublets), gates on CD4 T cells at the expression level, and uses SCT-based integration to build a unified CD4⁺ T-cell atlas.

After integration, the pipeline clusters and annotates major CD4 subsets (naive/CM, memory, Th1-like, Treg, CD4 CTL, activated memory), generates UMAPs, and quantifies population composition across clinical groups. It then processes TCR contigs to assign one TRA and one TRB chain per cell, constructs combined clonotypes (CTaa), computes clone sizes and expansion categories, and overlays this information onto the integrated object. Downstream outputs include stacked barplots of clonal expansion by group and cluster, UMAPs colored by expansion category, and a focused CD4 CTL “signature” comparison (e.g., Progressor vs PD HighCTL) using a curated gene panel.

To use the script, point the RNA and TCR paths to your per-sample 10x outputs, set the output directory, and run the file in R/RStudio from top to bottom. Sample and group labels, marker panels, QC thresholds, and clonal expansion bins are defined in a few central vectors and can be easily adapted to other cohorts or disease comparisons.

-----------------------------------------------------------------------------------------------------------

Spatial Transcriptomics – Densities and CD4 CTL Overlay

This script analyzes Xenium (or similar) spatial transcriptomics data to (i) quantify cell-type densities per mm² and (ii) generate a spatial map highlighting cardiomyocytes and CD4⁺ cytotoxic T cells with expanded TCRs.

It reads the standard matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz, and cell_boundaries.parquet files for each region of interest, removes zero-transcript cells, and uses the convex hull of cell boundaries to estimate the tissue area in mm². Canonical marker panels (cardiomyocytes, stressed CM, CD4/CD8 T cells, B cells, myeloid cells, HLA, apoptosis, fibrosis, cytotoxic/CD4 CTL, and expanded TCR probes) are then used to classify cells based on ≥1 detected transcript.
For each ROI, the script reports a compact density table (cells per mm²) for normal and stressed cardiomyocytes, T-cell subsets (including CD4 CTL and CD4 CTL with expanded CDR3/TCR genes), B cells, myeloid cells, HLA⁺ cells, apoptosis-signature cells, and fibrosis-like cells. In parallel, it computes cell centroids from the polygon boundaries and produces a spatial overlay plot in which normal and stressed cardiomyocytes, CD4 T cells, and expanded CD4 CTLs are overlaid on the tissue with distinct colors, exporting a high-resolution PNG.

To use the script, set base_dir to the per-patient spatial output directory and run both blocks in R. The gene panels, positivity definition (currently ≥1 molecule), and color scheme can be easily adjusted to match specific biology or figure style requirements.

----------------------------------------------------------------------------------------------------
Xenium Spatial – Progressor vs Survivor Tissue Programs

This script analyzes Xenium In Situ spatial transcriptomics from NICM Progressor and Survivor hearts to map major tissue compartments and compare their abundance and transcriptional states between groups. It loads per-patient Xenium outs directories, applies QC using control probes and gene-count thresholds, normalizes with SCTransform, and builds Harmony-corrected UMAP embeddings. Using curated gene panels (cardiomyocyte, endothelial, immune, fibroblast/myofibroblast, stress, fibrosis), it computes signature scores per cell, assigns each cell to a coarse focus_category (Cardiomyocyte, Stressed CM, Endothelium, Immune, Fibroblast/Myofibroblast), and removes low-confidence “islands” via DBSCAN on the UMAP to retain only well-supported tissue programs.

The cleaned object is then used to generate a split Harmony UMAP (Progressor vs Survivor) colored by focus_category, and a stacked percentage bar plot summarizing the relative proportions of each compartment per clinical group. Finally, the script aggregates SCT-normalized expression by group × focus_category and builds a horizontal ComplexHeatmap of Z-scored expression for a focused gene panel (cardiomyocyte/stress, endothelium, ECM/fibrosis, T-lineage, cytotoxicity, antigen presentation), with rows grouped by biological module and rows annotated by group (Survivor vs Progressor). Outputs include high-resolution UMAPs, composition bar plots, and PNG/PDF heatmaps suitable for main-figure panels.
