* INTRODUCTION *

This directory, bone_marrow_primary_survival_comp, performs
survival analysis on predicted subtype and cluster markers.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse		(version 1.3.2)
	* TCGAbiolinks		(version 2.25.3)
	* survival		(version 3.5.0)
	* survminer		(version 0.4.9)
	* survMisc		(version 0.5.6)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_survival_subtype_markers.R:

Step 1: Load data
Expression data, clinical data and defined set of subtype-related
gene expression markers are loaded.

Step 2: Wrangle data
Clinical and metadata of samples are wrangled to get survival
data.

Step 3: Survival analysis
Perform survival analysis on predicted subtype-related gene 
expression markers.

Step 4: Save results
Results from survival analysis are saved.

Process in script TARGET_ALL_P2_bmp_survival_cluster_markers.R:

Step 1: Load data
Expression data, clinical data, cluster labels and selected features
from random forest are loaded.

Step 2: Wrangle data
Clinical and metadata of samples are wrangled to get survival
data.

Step 3: Survival analysis
Perform survival analysis on predicted cluster-related gene 
expression markers.

Step 4: Save results
Results from survival analysis are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_survival_subtype_markers.R is the script that should be run first. 

TARGET_ALL_P2_bmp_survival_subtype_markers.R should be run in the following way from the terminal:
Rscript TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_survival_subtype_markers.R

TARGET_ALL_P2_bmp_survival_cluster_markers.R should be run next.

TARGET_ALL_P2_bmp_survival_cluster_markers.R should be run in the following way from the terminal:
Rscript TARGET_survival/bone_marrow_primary_survival_comp/TARGET_ALL_P2_bmp_survival_cluster_markers.R

TARGET_ALL_P2_bmp_survival_subtype_markers.R and TARGET_ALL_P2_bmp_survival_cluster_markers.R source the functions script:
TARGET_ALL_P2_bmp_survival_functions.R


* OUTPUT *

Output from TARGET_ALL_P2_bmp_survival_subtype_markers.R:

./TARGET_ALL_P2_bmp_subtype_markers_cox_ph.rda:
- Results of proportional hazards assumption test from step 3 of *
  PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_subtype_markers_cox_univariate.rda:
- Results of Cox univariate regression analysis from step 3 of *
  PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_subtype_markers_cox_multivariate.rda:
- Results of Cox multivariate regression analysis from step 3 of *
  PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_subtype_markers_km_<ENSEMBL_GENE_ID>.pdf:
- Kaplan-Meier survival plots of subtype-related expression markers
  from step 3 of * PROCESS IN THE SCRIPT *

Output from TARGET_ALL_P2_bmp_survival_cluster_markers.R:

./TARGET_ALL_P2_bmp_cluster_markers_cox_ph.rda:
- Results of proportional hazards assumption test from step 3 of *
  PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_cluster_markers_cox_univariate.rda:
- Results of Cox univariate regression analysis from step 3 of *
  PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_cluster_markers_km_<ENSEMBL_GENE_ID>.pdf:
- Kaplan-Meier survival plots of cluster-related expression markers
  from step 3 of * PROCESS IN THE SCRIPT *

