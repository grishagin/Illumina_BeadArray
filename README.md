# Illumina_BeadArray
An R script for Illumina BeadArray Data Processing.   
Currently, only MouseRef-8 v2.0 Expression BeadChip is supported.   

The repository provides means for the processing of Illumina BeadArray bead level data to produce heatmaps, volcano plots, Venn diagrams for any specified combination of controls/treatments and p-values/log fold changes.

## Contents

#### _archive
Previous versions.

#### vX-X
Current version of the script and dependencies.
	##### dependencies    
	All R functions.    
	##### _code_Analysis_Illumina_BeadArray_from_IDAT-MAIN-vXX.R       
	Main R script.
	##### _code_Analysis_Illumina_BeadArray_from_IDAT-MAIN-vXX.Rexec    
	"Executable" R script (to be opened with Rscript.exe).     
	##### Parameters_vX-X.xlsx    
	File with Analysis Parameters. Contents are self-explanatory.    
    ##### phenoData_vX-X.xlsx    
	File with Sample Description. Contents are self-explanatory.   

#### sample IDAT files
Folder with sample Illumina BeadArray bead level data (*.idat files).



## Method description   
Prerequisites:
1) R;
2) java.  
The script will check for all the necessary packages, download and install them.  

The script prompts for the following:   
- folder with dependencies (defaults to current working directory);  
- file with analysis parameters;   
- folder with bead level data (*.idat files);   
- file with sample description (phenoData).    
  
Sample IDAT files are provided, along with a template with their description.
  
The data are converted to ExpressionSetIllumina class, expression values are transformed to log2 scale, and then normalized by quantile method using the beadarray R package.  
The data are further annotated with Gene Symbols, Gene Names, Entrez Gene IDs, and Probe Quality Grades using illuminaMousev2.db annotation data R package.  
After that, the probes with “No match” and “Bad” probe quality grades are removed from the dataset.  
Next, differential expression is assessed by fitting a separate linear model for each gene and subjecting the resultant models to parametric empirical Bayes analysis using limma R package.  
The output is characterized by plotting significance in the form of B-statistic (log-odds) versus log fold changes in expression for each gene.  
Fold change analysis is applied to identify genes with expression ratios above the specified threshold value between treatments and control at specified significance levels.  
False discovery rate (FDR) is controlled by adjusting p-value using Benjamini-Hochberg algorithm.  
Hierarchical cluster analysis is performed using complete linkage and Euclidean distance as a measure of similarity. 

