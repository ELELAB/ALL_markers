* INTRODUCTION *

This directory, TARGET_known_markers, compares expression 
of predicted markers with known markers. 


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse		(version 1.3.2)
	* reshape2		(version 1.4.4)
	* viridis		(version 0.6.2)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_known_markers.R:

Step 1: Load data
Batch corrected expression data is loaded. 

Step 2: Calculate correlation between each predicted
marker and each known marker.

Step 3: Visualization
Correlations are visualized in a heatmap.

Step 4: Save results
Above results are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_known_markers.R is the script that should be run. 

TARGET_ALL_P2_bmp_known_markers.R should be run in the following way from the terminal:
Rscript TARGET_known_markers/bone_marrow_primary_known_markers_comp/TARGET_ALL_P2_bmp_known_markers.R


* OUTPUT *

Output from TARGET_ALL_P2_bmp_known_markers.R:

./TARGET_ALL_P2_bmp_corr_known_markers_mat.csv:
- csv file of correlations between predicted and known markers from step 2 of
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_corr_known_markers_heatmap.pdf:
- heatmap of correlations between predicted and known markers from step 3 of
  * PROCESS IN THE SCRIPT *
