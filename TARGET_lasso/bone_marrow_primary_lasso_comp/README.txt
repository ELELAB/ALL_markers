* INTRODUCTION *

This directory, bone_marrow_primary_lasso_comp, contains scripts 
to perform regularized logistic regression (LASSO and elastic net) 
on TARGET-ALL-P2 gene expression data. 
This directory also contains the folder elastic_consensus_biotypes
which contains a script that visualizes the results from regularized
logistic regression where biotypes of the predicted
features (genes) from ENSEMBL have been added.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* glmnet    (version 4.1.6)
	* caret     (version 6.0.93)
	* tidyverse (version 1.3.2)
	* biomaRt   (version 2.54.0)
	* viridis   (version 0.6.2)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_lasso.R:

Step 1: Load data
Gene expression data, consensus DEG data and metadata
about samples are loaded. 

Step 2: Logistic regression
Regularized logistic regression is performed on
gene expression data.

Step 3: Correlation
Correlation between mean coefficient values and mean
log2FC values is found.

Step 4: Visualization
Results of logistic regression are visualized.

Step 5: Save results
Above results are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_lasso.R is the script that should be run first. 

TARGET_ALL_P2_bmp_lasso.R should be run in the following way from the terminal:
Rscript TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso.R

TARGET_ALL_P2_bmp_lasso.R sources the functions script: 
TARGET_ALL_P2_bmp_lasso_functions.R


* OUTPUT *

./TARGET_ALL_P2_bmp_lasso_seeds.csv:
- vector of seeds used from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_runs_filt.rds:
- results of 10 LASSO runs using filtered data from step 2 of 
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_runs_filt.rds:
- results of 10 elastic net runs using filtered data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_runs_corr.rds:
- results of 10 LASSO runs using batch corrected data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_runs_corr.rds:
- results of 10 elastic net runs using batch corrected data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_runs_genes_filt.rds:
- minimal subset of genes found in the 10 LASSO runs including 
  their coefficients using filtered data from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_runs_genes_filt.rds:
- minimal subset of genes found in the 10 elastic net runs 
  including their coefficients using filtered data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_runs_genes_corr.rds:
- minimal subset of genes found in the 10 LASSO runs including 
  their coefficients using batch corrected data 
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_runs_genes_corr.rds:
- minimal subset of genes found in the 10 elastic runs including 
  their coefficients using batch corrected data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_consensus_genes_filt.csv:
- consensus minimal subset of genes using filtered data - lasso
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_consensus_genes_filt.csv:
- consensus minimal subset of genes using filtered data - elastic net
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_consensus_genes_corr.csv:
- consensus minimal subset of genes using batch corrected data - lasso
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_consensus_genes_corr.csv:
- consensus minimal subset of genes using batch corrected data - elastic net
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_pred_error_filt.csv:
- table of prediction errors for each of the 10 LASSO runs using filtered data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_pred_error_filt.csv:
- table of prediction errors for each of the 10 elastic net runs using filtered data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_pred_error_corr.csv:
- table of prediction errors for each of the 10 LASSO runs using batch corrected data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_pred_error_corr.csv:
- table of prediction errors for each of the 10 elastic net runs using batch corrected data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_consensus_DEA_corr.csv:
- table containing coefficients of consensus elastic net genes and comparison 
  with consensus DEA from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_pearson.csv:
- results from Pearson correlation test from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_pred_error_plot_filt.pdf:
- plots of prediction errors for each of the 10 LASSO runs using filtered data
  from step 4 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_pred_error_plot_filt.pdf:
- plots of prediction errors for each of the 10 elastic net runs using filtered data
  from step 4 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_lasso_pred_error_plot_corr.pdf:
- plots of prediction errors for each of the 10 LASSO runs using batch corrected data
  from step 4 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_pred_error_plot_corr.pdf:
- plots of prediction errors for each of the 10 elastic net runs using batch corrected data
  from step 4 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_mean_cv_error_plot_corr.pdf:
- plots of mean cros-validation errors for each of the 10 elastic net runs using batch corrected
  data from step 4 of * PROCESS IN THE SCRIPT *

