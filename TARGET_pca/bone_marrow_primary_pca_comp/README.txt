* INTRODUCTION *

This directory (bone_marrow_primary_pca_comp) performs MDS analysis of 
TARGET-ALL-P2 gene expression data where this data has been subject to
two workflows. These two workflows share the same first six step
which also makes up the first workflow:
1) subsetted to contain bone marrow primary samples
2) adjusted for replicates
3) preprocessed
4) normalized using updated geneInfoHT table
5) filtered 
6) voom transformed
In the second workflow, a subsequent step is performed:
7) batch correction for year of diagnosis variable

The MDS analyses belonging to these different workflows are stored in
the two subdirectories:
init_pca (workflow 1) and batch_corr_pca (workflow 2). 
More details about analyses belonging to workflow 1 and 2 can be
found in the README files of the respective directory. 

This directory contains script and resulting output of sample information
needed to perform MDS analysis. 


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2
- R packages:
        * tidyverse (version 1.3.2)
	* patchwork (version 1.1.2)


* PROCESS IN THE SCRIPT *

Step 1: Load data
TARGET-ALL-P2 data subsetted to contain bone marrow primary samples and where replicates
have been adjusted for is loaded. Voom transformed data is also loaded. Clinical data
belonging to TARGET-ALL-P2 project is also loaded. 

Step 2: Create table of sample information
A table containing barcodes, subtype, year of diagnosis, tissue portion, and number
isolation from aliquot is created for the samples. This table is needed for MDS
analyses and is subsequently saved. 


* RUNNING THE SCRIPT * 

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_info_comp.R is the script that should be run. 

TARGET_ALL_P2_bmp_info_comp.R should be run in the following way from the terminal:
Rscript TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_info_comp.R

This directory also contains the script TARGET_ALL_P2_bmp_pca_functions.R which
contains functions used to perform MDS analysis and is sourced by scripts 
performing MDS analysis. 


* OUTPUT *

./TARGET_ALL_P2_bmp_batch_info.csv
- Table containing barcodes, subtype, year of diagnosis, tissue portion, and number
  isolation from aliquot for the samples from step 2 of * PROCESS IN THE SCRIPT *




