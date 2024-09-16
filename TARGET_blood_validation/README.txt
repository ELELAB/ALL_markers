* INTRODUCTION *

This directory, TARGET_blood_validation, performs unsupervised
clustering on TARGET blood samples.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* TCGAbiolinks (version 2.25.3)
	* SummarizedExperiment (version 1.28.0)
	* limma (version 3.54.0)
	* gplots (version 3.1.3.1)
	* tidyverse (version 1.3.2)
	* viridis (version 0.6.2)


* PROCESS IN THE SCRIPT *

Process in script TARGET_blood_validation.R:

Step 1: Load data
TARGET ALL P2 expression data is loaded. 

Step 2: Wrangle data
Expression data is subsetted to contain blood samples and
expression data is preprocessed, normalized, filtered and
voom transformed.

Step 3: Clustering
Expression data of candidate genes is subject to clustering
and visualized in heatmap.

Step 4: Save results
Above results are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_blood_validation.R is the script that should be run. 

TARGET_blood_validation.R should be run in the following way from the terminal:
Rscript TARGET_blood_validation/TARGET_blood_validation.R


* OUTPUT *

Output from TARGET_blood_validation.R:

./ALL_P2_blood_primary_exp_data.rda:
- SummarizedExperiment object of expression data of blood samples from step 2 of
  * PROCESS IN THE SCRIPT *

./ALL_P2_blood_primary_prep_data.rda:
- Preprocessed expression data from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_blood_primary_norm_data.rda:
- Normalized expression data from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_blood_primary_filt_data.rda:
- Filtered expression data from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_blood_primary_voom_data.rda:
- Voom transformed data from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_blood_candidate_genes_clustering.pdf:
- Heatmap of clustering of candidate genes from step 3 of * PROCESS IN THE SCRIPT *

