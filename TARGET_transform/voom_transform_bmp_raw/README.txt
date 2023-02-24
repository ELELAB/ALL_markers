* INTRODUCTION *

This directory (voom_transform_bmp_raw) contains voom transformed
gene expression data of TARGET-ALL-P2 project. This data has been 
subject to the following workflow prior to voom transformation:
1) subsetted to contain only bone marrow (primary blood derived cancer) samples 


* REQUIREMENTS * 

Following R and R packages are needed:
- R 4.2
- R packages:
  * tidyverse (version 1.3.2)
	* limma (3.54.0)


* PROCESS IN THE SCRIPT *

Step 1: Load data
Raw read counts of data is loaded. 

Step 2: Voom transform data
Raw data is voom transformed. 

Step 3: Save data
Voom transformed data is saved. 


* RUNNING THE SCRIPT * 

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_raw_voom.R is the script that should be run. 

TARGET_ALL_P2_bmp_raw_voom.R should be run in the following way from the terminal:
Rscript TARGET_transform/voom_transform_bmp_raw/TARGET_ALL_P2_bmp_raw_voom.R


* OUTPUT * 

./TARGET_ALL_P2_bmp_raw_voom.rda
- Voom transformed data from step 2 of * PROCESS IN THE SCRIPT * 
