* INTRODUCTION *

This directory (bone_marrow_primary_comp) processes gene expression data from
TARGET-ALL-P2. This data has been subsetted to contain only bone marrow
(primary blood derived cancer) samples and adjusted for replicates. The data
is processed using the preprocessing, normalization, and filtering functions
from TCGAbiolinks. 


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2 
- R packages:
  * TCGAbiolinks (version 2.25.3)
	* tidyverse (version 1.3.2)
	* EDASeq (version 2.32.0)
	* readxl (version 1.4.1)


* PROCESS IN THE SCRIPT *

Process in script TARGET_process_bmp_comp.R:

Step 1: Load data
TARGET-ALL-P2 gene expression data is loaded from TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda

Step 2: Preprocess data
Gene expression data is preprocessed using TCGAanalyze_Preprocessing function
from TCGAbiolinks

Step 3: Normalize data
Preprocessed gene expression data is normalized using TCGAanalyze_Normalization 
function from TCGAbiolinks. This is done using updated geneInfoHT table

Step 4: Filter data
Normalized data is filtered using TCGAanalyze_Filtering function from
TCGAbiolinks

Process in script TARGET_norm_genes_bmp_comp.R:

Step 1: Load data
Gene expression data of TARGET-ALL-P2 bone marrow primary samples and where replicates have
been adjusted for, preprocessed and normalized data are loaded

Step 2: Investigate lost genes following normalization step
Genes that are lost following normalization step are investigated 

Process in script TARGET_norm_genes_plot_bmp_comp.R:

Step 1: Load data
Excel sheet containing information about lost genes following normalization step is loaded

Step 2: Visualize data
Barplot of different gene types of excluded genes is created 


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the script from there. 

TARGET_process_bmp_comp.R, TARGET_norm_genes_bmp_comp.R and TARGET_norm_genes_plot_bmp_comp.R
are the scripts that should be run. 

TARGET_process_bmp_comp.R should be run first in the following way from the terminal:
Rscript TARGET_processing/bone_marrow_primary_comp/TARGET_process_bmp_comp.R

TARGET_norm_genes_bmp_comp.R should be run second in the following way from the terminal:
Rscript TARGET_processing/bone_marrow_primary_comp/TARGET_norm_genes_bmp_comp.R

TARGET_norm_genes_plot_bmp_comp.R should be run third in the following way from the terminal:
Rscript TARGET_processing/bone_marrow_primary_comp/TARGET_norm_genes_plot_bmp_comp.R


* OUTPUT *

Output from script TARGET_process_bmp_comp.R:

./ALL_P2_bmp_prep_data.rda:
- Preprocessed gene expression data from step 2 of * PROCESS IN THE SCRIPT *

./preprocessing_output_AAIC.png:
- Output from preprocessing function from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_norm_data.rda:
- Normalized gene expression data from step 3 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_filt_data.rda:
- Filtered gene expression data from step 4 of * PROCESS IN THE SCRIPT *

Output from script TARGET_norm_genes_bmp_comp.R:

./ALL_P2_bmp_norm_excl_genes.csv:
- Table containing information about lost genes following normalization step
  from step 2 of * PROCESS IN THE SCRIPT *
  This file contains ENSEMBL IDs of excluded genes and if they are found in
  updated geneInfoHT table. 
  These ENSEMBL IDs are subsequently looked up in the biomart database to find
  information about these excluded genes: https://www.ensembl.org/index.html. 
  This information is contained in the excel file ALL_P2_bmp_norm_excl_genes_biomart.xlsx

Output from script TARGET_norm_genes_plot_bmp_comp.R:

./ALL_P2_bmp_excl_genes_norm_plot.pdf
- Barplot showing number of gene types of genes lost following normalization step
  from step 2 of * PROCESS IN THE SCRIPT *



