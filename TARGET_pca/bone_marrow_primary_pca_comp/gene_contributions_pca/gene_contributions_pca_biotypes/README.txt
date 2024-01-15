* INTRODUCTION *

This directory (gene_contributions_pca_biotypes) visualizes gene contributions to
PCA dimensions where PCA has been performed on 
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
	* tidyverse  (version 1.3.2)
	* viridis    (version 0.6.2)


* PROCESS IN THE SCRIPT *

Step 1: Load data
Top 40 genes contributing to principal components 1 and 2 and results of 
differential expression analysis are loaded.

Step 2: Perform PCA
Top 40 genes contributing to principal components 1 and 2 are compared 
with differential expression analysis results and this is afterwards
visualized in a barplot.


* RUNNING THE SCRIPT * 

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes.R is the script that should be run. 

TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes.R should be run in the following way from the terminal:
Rscript TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/gene_contributions_pca_biotypes/TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes.R


* OUTPUT *

./TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes.csv
- csv file showing the contributions in percentage of the top 40 genes
  contributing to principal components 1 and 2 together with their biotype
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes_info.csv
- csv file containing comparison between PCA gene contributions and results
  from differential expression analysis from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes_plot.pdf
- Barplot showing top 40 genes contributing to principal components 1 and 2
  and results from differential expression analysis from step 2 of * PROCESS IN THE SCRIPT *


