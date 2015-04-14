# Illumina_BeadArray
An R script for Illumina BeadArray Data Processing.

The repository provides means for the processing of Illumina BeadArray bead level data to produce heatmaps, volcano plots, Venn diagrams for any specified combination of controls/treatments and p-values/log fold changes.

## Contents
Repository contains each version in a separate folder. Current version is marked accordingly ("-current").
#### Main Script
_code_Analysis_Illumina_BeadArray_from_IDAT-MAIN-v22.R

#### Dependencies Folder 
_code_FUNCTION#01__loadPackages-v1.R
_code_FUNCTION#02_read.norm.filt.data-v2.R
_code_FUNCTION#03_QC-v5.R
_code_FUNCTION#04_design.Matrix.prepare-v1.R
_code_FUNCTION#05_contrasts.Matrix.prepare-v2.R
_code_FUNCTION#06_fit.data-v2.R
_code_FUNCTION#07_produce.volcano.plot-v1.R
_code_FUNCTION#08_make.volcano-v2.R
_code_FUNCTION#09_venn-v2.R
_code_FUNCTION#10_vennDiagram.mod-v1.R
_code_FUNCTION#11_top.tables-v2.R
_code_FUNCTION#12_produce.heatmaps-v2.R
_code_FUNCTION#13_prep.heatmap.sort.data-v1-1.R
_code_FUNCTION#14_plot.heatmap-v1-1.R
_code_FUNCTION#15_heatmap.2.mod-v1-1.R

#### Template for the sample description 
samples.description-TEMPLATE.txt

#### Folder with sample IDAT files
sample IDAT files

## Method description
The script prompts for the folder with dependencies (defaults to current working directory).
Then the script prompts for the folder with bead level data in the form of *.idat files, and a description of the samples as a text file.
Sample IDAT files are provided, along with a template with their description.
  
The data are then converted to ExpressionSetIllumina class, expression values are transformed to log2 scale, and then normalized by quantile method using the beadarray R package.  
The data are further annotated with Gene Symbols, Gene Names, Entrez Gene IDs, and Probe Quality Grades using illuminaMousev2.db annotation data R package.  
After that, the probes with “No match” and “Bad” probe quality grades are removed from the dataset.  
Next, differential expression is assessed by fitting a separate linear model for each gene and subjecting the resultant models to parametric empirical Bayes analysis using limma R package.  
The output is characterized by plotting significance in the form of B-statistic (log-odds) versus log fold changes in expression for each gene.  
Fold change analysis is applied to identify genes with expression ratios above a certain preset value (lfc = 1 by default) between treatments and control at specified significance levels (P-value < 0.01 by default).  
False discovery rate (FDR) is controlled by adjusting p-value using Benjamini-Hochberg algorithm.  
Hierarchical cluster analysis is performed using complete linkage and Euclidean distance as a measure of similarity. 

