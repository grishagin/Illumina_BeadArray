# Illumina_BeadArray
An R script for Illumina BeadArray Data Processing.

The repository provides means for the processing of Illumina BeadArray bead level data to produce heatmaps, volcano plots, Venn diagrams for any specified combination of controls/treatments and p-values/log fold changes.

## Contents
Repository contains each version in a separate folder. Current version is marked accordingly ("-current").
#### Main Script
_code_Analysis_Illumina_BeadArray_from_IDAT-MAIN-v20.R

#### Dependencies 
_code_FUNCTION#1__loadPackages-v1.R  
_code_FUNCTION#2_read.norm.filt.data-v1.R  
_code_FUNCTION#3_QC.CDKN1A.PPIA.HPRT-v3.R  
_code_FUNCTION#4_designMatrix.contMatrix.prep-v1.R  
_code_FUNCTION#5_fit.data-v1.R  
_code_FUNCTION#6_volcano.plot-v1.R  
_code_FUNCTION#7_venn-v1.R  
_code_FUNCTION#8_vennDiagram.mod-v1.R  
_code_FUNCTION#9_top.tables-v1.R  
_code_FUNCTION#10_produce.heatmaps-v1.R  
_code_FUNCTION#11_prep.heatmap.sort.data-v1.R  
_code_FUNCTION#12_plot.heatmap-v1.R  
_code_FUNCTION#13_heatmap.2.mod-v1.R  

#### Template for the sample description 
samples.description.txt

## Method description
The script prompts for the folder with bead level data in the form of *.idat files, and a description of the samples as a text file.  
The data are then converted to ExpressionSetIllumina class, expression values are transformed to log2 scale, and then normalized by quantile method using the beadarray R package.  
The data are further annotated with Gene Symbols, Gene Names, Entrez Gene IDs, and Probe Quality Grades using illuminaMousev2.db annotation data R package.  
After that, the probes with “No match” and “Bad” probe quality grades are removed from the dataset.  
Next, differential expression is assessed by fitting a separate linear model for each gene and subjecting the resultant models to parametric empirical Bayes analysis using limma R package.  
The output is characterized by plotting significance in the form of B-statistic (log-odds) versus log fold changes in expression for each gene.  
Fold change analysis is applied to identify genes with expression ratios above a certain preset value (lfc = 1 by default) between treatments and control at specified significance levels (P-value < 0.01 by default).  
False discovery rate (FDR) is controlled by adjusting p-value using Benjamini-Hochberg algorithm.  
Hierarchical cluster analysis is performed using complete linkage and Euclidean distance as a measure of similarity. 

