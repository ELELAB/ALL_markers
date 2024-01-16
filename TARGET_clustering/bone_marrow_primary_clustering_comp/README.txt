* INTRODUCTION *

This directory, bone_marrow_primary_clustering_comp, contains 
scripts that perform unsupervised clustering on gene expression
data using the cola framework and analyzes the results.
This directory also contains a subfolder, cola_on_raw_data,
where results are saved from performing unsupervised clustering
on cola-adjusted raw data.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* cola 			 (version 2.4.0)
	* tidyverse		 (version 1.3.2)
	* caret			 (version 6.0.93)
	* viridis		 (version 0.6.2)
	* gridExtra		 (version 2.3)
	* cowplot		 (version 1.1.1)
	* umap			 (version 0.2.9.0)
	* SummarizedExperiment   (version 1.28.0)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_run_cola.R:

Step 1: Load data
Raw and batch corrected gene expression data are loaded.

Step 2: Unsupervised clustering
Unsupervised clustering using the cola framework is
performed.

Step 3: Save results
Above results are saved.


Process in script TARGET_ALL_P2_bmp_analyze_cola.R:

Step 1: Load data
Results from unsupervised clustering are loaded.

Step 2: Visualization
Results of unsupervised clustering are visualized.

Step 3: Save results
Above results are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_run_cola.R is the script that should be run first. 

TARGET_ALL_P2_bmp_run_cola.R should be run in the following way from the terminal:
Rscript TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_run_cola.R

TARGET_ALL_P2_bmp_analyze_cola.R should be run next. 

TARGET_ALL_P2_bmp_analyze_cola.R should be run in the following way from the terminal:
Rscript TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_analyze_cola.R

TARGET_ALL_P2_bmp_analyze_cola.R sources the functions script:
TARGET_ALL_P2_bmp_cola_functions.R


* OUTPUT *

Output from TARGET_ALL_P2_bmp_run_cola.R:

./cola_on_raw_data/TARGET_ALL_P2_bmp_raw_adj_data_cola.rds:
- Cola-adjusted raw data from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_batch_cola_results.rds:
- Results from unsupervised clustering on batch corrected data
  from step 2 of * PROCESS IN THE SCRIPT *

./cola_on_raw_data/TARGET_ALL_P2_bmp_raw_cola_results.rds:
- Results from unsupervised clustering on cola-adjusted raw data
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_batch_cola_report.<file_extension>:
- Results of unsupervised clustering performed on batch corrected
  data in a collective report from step 2 of * PROCESS IN THE SCRIPT *

./cola_on_raw_data/TARGET_ALL_P2_bmp_raw_cola_report.<file_extension>:
- Results of unsupervised clustering performed on cola-adjusted raw
  data in a collective report from step 2 of * PROCESS IN THE SCRIPT *


Output from TARGET_ALL_P2_bmp_analyze_cola.R:

./TARGET_ALL_P2_bmp_cola_samples_clusters.csv:
- Table containing predicted cluster labels of samples together with 
  probability of membership and actual class label from step 2 of
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_cola_batch_confusion_matrices.pdf:
- Confusion matrix plots of batch corrected data using all methods in 
  one pdf from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_cola_raw_confusion_matrices.pdf:
- Confusion matrix plots of raw data using all methods in one pdf
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_cola_batch_confusion_matrices_collective.pdf:
- Plot that contains confusion matrices together in one plot 
  (batch corrected data) from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_cola_raw_confusion_matrices_collective.pdf:
- Plot that contains confusion matrices together in one plot (raw data)
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_cola_stats.pdf:
- Grouped barplot showing statistical measures of clustering methods
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_umap_SD_mclust.pdf:
- UMAP visualization of clusters using SD:mclust k = 4 method from
  step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_umap_MAD_mclust.pdf:
- UMAP visualization of clusters using MAD:mclust k = 4 method from
  step 2 of * PROCESS IN THE SCRIPT *

