* INTRODUCTION *

This directory, bone_marrow_primary_random_forest_comp, performs
variable selection on predicted clusters of TARGET ALL P2 gene
expression data using random forest. The predicted clusters were
found by performing unsupervised clustering using the cola framework.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse		(version 1.3.2)
	* randomForest		(version 4.7.1.1)
	* varSelRF			(version 0.7.8)
	* caret		 	(version 6.0.93)
	* biomaRt			(version 2.54.0)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_clusters_random_forest.R:

Step 1: Load data
Expression data and information about samples including
cluster labels are loaded. 

Step 2: Random forest
Perform variable selection using random forest on
predicted clusters.

Step 3: Save results
Above results are saved.

Process in script TARGET_ALL_P2_bmp_analyze_clusters_random_forest.R:

Step 1: Load data
Random forest results are loaded. Differential expression analysis
results and housekeeping genes are also loaded.

Step 2: Wrangle data
Intersected variables across 10 random forest seed runs are found. 
Random forest selected genes are compared with results from
differential expression analysis and housekeeping analysis. 

Step 3: Save results
Above results are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_clusters_random_forest.R is the script that should be run first. 

TARGET_ALL_P2_bmp_clusters_random_forest.R should be run in the following way from the terminal:
Rscript TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_clusters_random_forest.R

TARGET_ALL_P2_bmp_analyze_clusters_random_forest.R should be run next.

TARGET_ALL_P2_bmp_analyze_clusters_random_forest.R should be run in the following way from the terminal:
Rscript TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_analyze_clusters_random_forest.R

TARGET_ALL_P2_bmp_clusters_random_forest.R and TARGET_ALL_P2_bmp_analyze_clusters_random_forest.R source the functions script:
TARGET_ALL_P2_bmp_random_forest_functions.R


* OUTPUT *

Output from TARGET_ALL_P2_bmp_clusters_random_forest.R:

./TARGET_ALL_P2_bmp_B_clusters_random_forest_results.rda:
- Random forest results performed on predicted clusters within the 
  original B-cell ALL subtype from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_T_clusters_random_forest_results.rda:
- Random forest results performed on predicted clusters within the 
  original T-cell ALL subtype from step 2 of * PROCESS IN THE SCRIPT *

Output from TARGET_ALL_P2_bmp_analyze_clusters_random_forest.R:

./TARGET_ALL_P2_bmp_selected_features_B_clusters_rf.csv:
- Table containing all selected features from all seed runs from random
  forest performed on predicted clusters within the original B-cell ALL
  subgroup from step 2 of * PROCESS IN THE SCRIPT * 

./TARGET_ALL_P2_bmp_selected_features_T_clusters_rf.csv:
- Table containing all selected features from all seed runs from random
  forest performed on predicted clusters within the original T-cell ALL
  subgroup from step 2 of * PROCESS IN THE SCRIPT * 

