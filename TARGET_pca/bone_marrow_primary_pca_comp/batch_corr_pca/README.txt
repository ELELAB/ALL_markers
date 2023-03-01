* INTRODUCTION *

This directory (batch_corr_pca) performs MDS analysis of TARGET-ALL-P2
gene expression data where this data has been subject to the 
following workflow:
1) subsetted to contain bone marrow primary samples
2) adjusted for replicates
3) preprocessed
4) normalized using updated geneInfoHT table
5) filtered 
6) voom transformed
7) batch corrected for year of diagnosis


* REQUIREMENTS * 

Following R and R packages are needed:
- R 4.2
- R packages:
        * tidyverse (version 1.3.2)
        * patchwork (version 1.1.2)


* PROCESS IN THE SCRIPT * 

Step 1: Load data
Gene expression data that has been batch corrected for year of diagnosis is loaded. 
Table containing information (subtype, year of diagnosis, 
tissue portion, and number isolation from aliquot) about samples is also loaded. 

Step 2: Perform MDS analysis
MDS analysis of batch corrected data is performed coloring for subtype, year of diagnosis, 
tissue portion, and number isolation from aliquot. Plots are also created where the MDS
plot of subtype is side by side with MDS plots of each of the three variables (year of diagnosis, 
tissue portion, and number isolation from aliquot). Plots are saved.


* RUNNING THE SCRIPT * 

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_pca_batch_comp.R is the script that should be run. This script sources the
functions script TARGET_ALL_P2_bmp_pca_functions.R. 

TARGET_ALL_P2_bmp_pca_batch_comp.R should be run in the following way from the terminal:
Rscript TARGET_pca/bone_marrow_primary_pca_comp/batch_corr_pca/TARGET_ALL_P2_bmp_pca_batch_comp.R


* OUTPUT *

./ALL_P2_bmp_subtype_MDS_batch.pdf
- MDS plot coloring for subtypes from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_year_MDS_batch.pdf
- MDS plot coloring for year of diagnosis from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_portion_MDS_batch.pdf
- MDS plot coloring for tissue portion from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_isolation_MDS_batch.pdf
- MDS plot coloring for number isolation from aliquot from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_subtype_year_MDS_batch.pdf
- MDS plots coloring for subtype and year of diagnosis side by side from step
  2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_subtype_portion_MDS_batch.pdf
- MDS plots coloring for subtype and tissue portion side by side from step 2 of
  * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_subtype_isolation_MDS_batch.pdf
- MDS plots coloring for subtype and number isolation from aliquot side by side
  from step 2 of * PROCESS IN THE SCRIPT *
