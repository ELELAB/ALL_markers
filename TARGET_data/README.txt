* INTRODUCTION *

This directory (TARGET_data) contains gene expression data
of the TARGET-ALL-P2 project. Information about the data 
such as number of genes, number of samples, subtype information, 
and age of patients are saved. This data is subsetted to 
contain only samples obtained from bone marrow (primary blood 
derived cancer). Additionally, meta data and variables about
the data is saved. 


* REQUIREMENTS *

Following R and R packages are needed: 
- R 4.2
- R packages:
	* tidyverse (version 1.3.2)
	* TCGAbiolinks (version 2.25.3) 
	* SummarizedExperiment (version 1.28.0)
	* patchwork (version 1.1.2)


* PROCESS IN THE SCRIPT *

Step 1: Gene expression data of TARGET-ALL-P2 project is loaded
Gene expression data of TARGET-ALL-P2 project is loaded from file ALL-P2_exp_data.rda

Step 2: Get information about data
-Clinical data of TARGET-ALL-P2 project is obtained
-Number of genes, number of samples, and number of samples belonging to different variables
 are obtained: gender, subtype, and vital status
-The ages of TARGET-ALL-P2 patients are obtained
-Number of samples belonging to different combinations between subtype and tissue source
 are obtained

Step 3: Create summary tables of meta data/variables available in TARGET-ALL-P2 project
Two summary tables containing information about meta data/variables are made. For each
variable, the tables include number of NA values, number of unique values, and names
of unique values. One summary table is made of clinical data of the TARGET-ALL-P2 project. 
The other summary table is made of variables available in the SummarizedExperiment object
of the TARGET-ALL-P2 project.

Step 4: Visualize age and sample overview of TARGET-ALL-P2 data
-Age distribution of TARGET-ALL-P2 patients is created
-Number of samples belonging to different combinations of subtype and tissue source is
 created

Step 5: Subset gene expression data to bone marrow primary samples
-The TARGET-ALL-P2 gene expression data is subsetted to contain only samples obtained
 from bone marrow (primary)
-Raw data of the subsetted gene expression data is obtained
-Number of genes, number of samples, and number of samples belonging to different variables
 are obtained: gender, subtype, and vital status
-The ages of TARGET-ALL-P2 patients in subsetted data are obtained

Step 6: Visualize age of TARGET-ALL-P2 patients in bone marrow primary subsetted data
-Age distribution of TARGET-ALL-P2 patients in bone marrow primary subsetted data is 
 created
-Plot of number of samples belonging to different combinations of tissue source and subtype
 is plotted together with age distribution of patients in subsetted data

Step 7: Save data
All of above are saved. 


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the script from there. 

get_TARGET_data.R is the script that should be run. This script sources the script:
get_TARGET_data_functions.R. 

get_TARGET_data.R should be run in the following way from the terminal:
Rscript TARGET_data/get_TARGET_data.R


* OUTPUT *

./TARGET_ALL_P2_clinical_data.rda
- clinical data of TARGET-ALL-P2 project from step 2 of * PROCESS IN
  THE SCRIPT * 

./TARGET_ALL_P2_clinical_data_headers.csv
- headers of clinical data of TARGET-ALL-P2 project from step 3 of 
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_clinical_data_summary.csv
- summary table of meta data/variables available in clinical data of
  TARGET-ALL-P2 project from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_exp_data_headers.csv
- headers of meta data/variables in available in SummarizedExperiment
  object of TARGET-ALL-P2 project from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_exp_data_summary.csv
- summary table of meta data/variable available in SummarizedExperiment
  object of TARGET-ALL-P2 project from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_exp_data_info.csv
- a table containing the number of genes, number of samples, and number of 
  samples of variables: gender, subtype, and vital status in TARGET-ALL-P2
  project from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_exp_data_age_table.csv
- a table containing the ages of TARGET-ALL-P2 patients from step 2
  of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_exp_data_subtype_tissue.csv
- a table containing the number of ALL-P2 samples belonging to combinations
  between subtype (B/T ALL) and tissue source type (blood/bone marrow,
  primary/recurrent cancer) from step 2 of * PROCESS IN THE SCRIPT *

./ALL-P2_bone_marrow_primary_exp_data.rda:
- gene expression data subsetted from TARGET-ALL-P2 expression data to
  contain only samples obtained from bone marrow primary
  from step 5 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_raw_read_counts.rda:
- raw read counts of TARGET-ALL-P2 gene expression data subsetted to
  contain only bone marrow primary samples from step 5 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_exp_data_info.csv
- a table containing the number of genes, number of samples, and number of
  samples of variables: gender, subtype, and vital status in TARGET-ALL-P2
  project subsetted to contain only bone marrow primary samples from step 
  5 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_exp_data_age_table.csv
- a table containing the ages of TARGET-ALL-P2 patients subsetted to contain
  only bone marrow primary samples from step 5 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_exp_data_age_plot.pdf
- a plot showing the age distribution of TARGET-ALL-P2 patients from 
  step 4 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_exp_data_subtype_tissue_plot.pdf
- a plot showing the number of ALL-P2 samples belonging to combinations
  between subtype (B/T ALL) and tissue source type (blood/bone marrow, 
  primary/recurrent cancer) from step 4 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_exp_data_age_plot.df
- a plot showing the age distribution of TARGET-ALL-P2 patients subsetted
  to contain only bone marrow primary samples from step 6 of * PROCESS
  IN THE SCRIPT *

./TARGET_ALL_P2_bmp_exp_data_age_subtype_tissue_plot.df
- a plot containing the combinations of samples by subtype and tissue source
  together with age distribution of patients in bone marrow primary subsetted
  data from step 6 of * PROCESS IN THE SCRIPT * 

