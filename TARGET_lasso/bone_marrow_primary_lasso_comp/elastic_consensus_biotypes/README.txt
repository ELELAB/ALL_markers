* INTRODUCTION *

This directory, elastic_consensus_biotypes, contains 
a script that visualizes the results from regularized
logistic regression where biotypes of the predicted
features (genes) from ENSEMBL have been added.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse (version 1.3.2)
	* viridis   (version 0.6.2)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_consensus_biotypes.R:

Step 1: Load data
Results from elastic net logistic regression is loaded.

Step 2: Visualization
Results of elastic net logistic regression are visualized.

Step 3: Save results
Above visualization is saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_consensus_biotypes.R is the script that should be run first. 

TARGET_ALL_P2_bmp_consensus_biotypes.R should be run in the following way from the terminal:
Rscript TARGET_lasso/bone_marrow_primary_lasso_comp/elastic_consensus_biotypes/TARGET_ALL_P2_bmp_consensus_biotypes.R


* OUTPUT *

./TARGET_ALL_P2_bmp_elastic_consensus_genes_biotype.pdf:
- Boxplot of elastic net coefficient values of elastic net
  consensus genes from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_elastic_consensus_biotype.csv:
- Table of elastic net consensus genes including biotypes 
  from step 2 of * PROCESS IN THE SCRIPT *

