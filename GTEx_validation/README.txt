* INTRODUCTION *

This directory, GTEx_validation, performs unsupervised
clustering on GTEx blood and bone marrow samples.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* TCGAbiolinks (version 2.25.3)
	* limma (version 3.54.0)
	* gplots (version 3.1.3.1)
	* viridis (version 0.6.2)
	* recount (version 1.24.1)


* PROCESS IN THE SCRIPT *

Process in script GTEx_validation.R:

Step 1: Get data
Blood and bone marrow expression data from GTEx are downloaded.

Step 2: Wrangle data
The expression data is scaled, normalized, filtered and boom
transformed. The processed data is subsetted to contain the
candidate genes of interest and unsupervised hierarchical
clustering is performed. 

Step 3: Save results
Above results are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

GTEx_validation.R is the script that should be run. 

GTEx_validation.R should be run in the following way from the terminal:
Rscript GTEx_validation/GTEx_validation.R


* OUTPUT *

Output from GTEx_validation.R:

./GTEx_blood_candidate_genes_clustering.pdf:
- Heatmap of unsupervised hierarchical clustering of blood samples from step 2 of * PROCESS IN THE SCRIPT *

./GTEx_bone_marrow_candidate_genes_clustering.pdf:
- Heatmap of unsupervised hierarchical clustering of bone marrow samples from step 2 of * PROCESS IN THE SCRIPT *

./GTEx_blood_SE.rda:
- Downloaded GTEx blood data in form of SummarizedExperiment from step 1 of * PROCESS IN THE SCRIPT *

./GTEx_bone_marrow_SE.rda:
- Downloaded GTEx bone marrow data in form of SummarizedExperiment from step 1 of * PROCESS IN THE SCRIPT *

./GTEx_blood_scaled_counts.rda:
- Scaled read count of GTEx blood data from step 2 of * PROCESS IN THE SCRIPT *

./GTEx_bone_marrow_scaled_counts.rda:
- Scaled read count of GTEx bone marrow data from step 2 of * PROCESS IN THE SCRIPT *

./GTEx_blood_norm_data.rda:
- Normalized data of blood samples from step 2 of * PROCESS IN THE SCRIPT *

./GTEx_bone_marrow_norm_data.rda:
- Normalized data of bone marrow samples from step 2 of * PROCESS IN THE SCRIPT *

./GTEx_blood_filt_data.rda:
- Filtered data of blood samples from step 2 of * PROCESS IN THE SCRIPT *

./GTEx_bone_marrow_filt_data.rda:
- Filtered data of bone marrow samples from step 2 of * PROCESS IN THE SCRIPT *

./GTEx_blood_voom_data.rda:
- Voom transformed data of blood samples from step 2 of * PROCESS IN THE SCRIPT *

./GTEx_bone_marrow_voom_data.rda:
- Voom transformed data of bone marrow samples from step 2 of * PROCESS IN THE SCRIPT *
