* INTRODUCTION *

This directory (gene_contributions_pca) performs PCA analysis of 
TARGET-ALL-P2 gene expression data where this data has been subject to
the following workflow. 
1) subsetted to contain bone marrow primary samples
2) adjusted for replicates
3) preprocessed
4) normalized using updated geneInfoHT table
5) filtered 
6) voom transformed
7) batch correction for year of diagnosis variable

This directory contains scripts for performing these analyses.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2
- R packages:
	* factoextra (version 1.0.7)
	* FactoMineR (version 2.7)
	* tidyverse  (version 1.3.2)
	* biomaRt    (version 2.54.0)


* PROCESS IN THE SCRIPT *

Step 1: Load data
TARGET-ALL-P2 batch corrected data and information (metadata) about samples are loaded. 

Step 2: Perform PCA
PCA is performed on batch corrected data to find contributions of features (genes)
to principal components 1 and 2. 


* RUNNING THE SCRIPT * 

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_pca_gene_contrib.R is the script that should be run. 

TARGET_ALL_P2_bmp_pca_gene_contrib.R should be run in the following way from the terminal:
Rscript TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/TARGET_ALL_P2_bmp_pca_gene_contrib.R

This directory also contains the script TARGET_ALL_P2_bmp_pca_gene_contrib_functions.R which
contains functions used to perform MDS analysis and is sourced by script 
TARGET_ALL_P2_bmp_pca_gene_contrib.R


* OUTPUT *

./TARGET_ALL_P2_bmp_pca_contrib_PC1_corr.pdf
- Barplot showing contributions in percentage of features (genes)
  to principal component 1 from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_pca_contrib_PC2_corr.pdf
- Barplot showing contributions in percentage of features (genes)
  to principal component 2 from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_pca_corr_genes_info.csv
- csv file showing contributions in percentages of top 20 features 
  (genes) to principal component 1 and of top 20 features (genes)
  to principal component 2 from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_pca_corr_plot.pdf
- PCA visualization of samples visualized by subtype from 
  step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_pca_corr.rda
- Results of PCA from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_pca_screeplot_corr.pdf
- Scree plot showing percentage of explained variance of
  principal components from step 2 of * PROCESS IN THE SCRIPT *

