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

