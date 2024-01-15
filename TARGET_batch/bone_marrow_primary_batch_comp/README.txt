* INTRODUCTION *

This directory (TARGET_batch) performs batch correction based
on identified batch factors of TARGET-ALL-P2 gene expression
data. This data has been subsetted to contain only bone marrow
primary samples, adjusted for replicates, preprocessed, 
normalized, filtered, and voom transformed. Batch factors were
identified using multidimensional scaling analysis.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2
- R packages:
	* tidyverse (version 1.3.2)
	* TCGAbiolinks (version 2.25.3) 
	* sva (version 3.46.0)


* PROCESS IN THE SCRIPT *

Step 1: Load data
Voom transformed data TARGET-ALL-P2 project and table containing information about
samples are loaded. 

Step 2: Get information about samples according to batch factor and condition
Number of samples in each batch factor (year of diagnosis and tissue portion) and 
each condition (subtype) are obtained. 

Step 3: Batch correction
Voom transformed data is batch corrected using year of diagnosis as batch factor.

Step 4: Save data
Above data is saved. 


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_batch_comp.R is the script that should be run. This script sources the script:
TARGET_ALL_P2_bmp_batch_functions.R.

TARGET_ALL_P2_bmp_batch_comp.R should be run in the following way from the terminal:
Rscript TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_batch_comp.R


* OUTPUT *

./TARGET_ALL_P2_bmp_samples_portion.csv
- Number of samples in each tissue portion batch factor from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_samples_portion_subtype.csv
- Number of samples in each tissue portion batch factor and subtype from step 2 of * PROCESS IN
  THE SCRIPT *

./TARGET_ALL_P2_bmp_samples_year.csv
- Number of samples in each year of diagnosis batch factor from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_samples_year_subtype.csv
- Number of samples in each year of diagnosis batch factor and subtype from step 2 of * PROCESS IN
  THE SCRIPT *

./TARGET_ALL_P2_bmp_combat_year_plots.pdf
- Plots outputted by COMBAT batch correction method from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_corrected_year.rda
- Corrected gene expression values based on year of diagnosis batch factor from step 3 of * PROCESS
  IN THE SCRIPT *

