* INTRODUCTION *

This directory, bone_marrow_primary_housekeeping_comp, contains scripts 
that compare the consensus DEGs found by limma-voom, edgeR and DESeq2 
to lists of housekeeping genes based on literature reviews.
This directory also contains the folder HKG_lists which contains the
lists of housekeeping genes used for comparison.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse (version 1.3.2)
	* viridis   (version 0.6.2)
	* gridExtra (version 2.3)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_housekeeping.R:

Step 1: Load data
Consensus DEG data and housekeeping gene lists are loaded.

Step 2: Data wrangling
Gene lists are extracted from the loaded data.

Step 3: Comparisons
Intersections between the three lists are found.

Step 4: Save results
Gene intersection lists from step 3 are saved.

Process in script TARGET_ALL_P2_bmp_housekeeping_visualization.R:

Step 1: Load data
Consensus DEG data, housekeeping gene lists, filtered and raw
gene expression data, and metadata information about samples
are loaded.

Step 2: Data wrangling
Loaded data is prepared for subsequent plotting.

Step 3: Visualization
Boxplots and density plot of expression data of identified housekeeping DEGs
are created. 

Step 4: Save results
Visualizations are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_housekeeping.R is the script that should be run first. 

TARGET_ALL_P2_bmp_housekeeping.R should be run in the following way from the terminal:
Rscript TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/TARGET_ALL_P2_bmp_housekeeping.R

TARGET_ALL_P2_bmp_housekeeping_visualization.R should be run next.

TARGET_ALL_P2_bmp_housekeeping_visualization.R should be run in the following way
from the terminal:
Rscript TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/TARGET_ALL_P2_bmp_housekeeping_visualization.R

TARGET_ALL_P2_bmp_housekeeping.R and TARGET_ALL_P2_bmp_housekeeping_visualization.R 
source the functions script: TARGET_ALL_P2_bmp_housekeeping_functions.R


* OUTPUT *

Output from script TARGET_ALL_P2_bmp_housekeeping.R:

./eisenberg_vs_homemade.csv:
- list of common genes in the two housekeeping genes lists from step 3 of * PROCESS IN THE SCRIPT *

./DEG_vs_eisenberg.csv:
- list of common genes for Eisenberg HKGs and consensus DEGs from step 3 of * PROCESS IN THE SCRIPT *

./DEG_vs_homemade.csv:
- list of common genes for tissue relevant HKGs and consensus DEGs from step 3 of * PROCESS IN THE SCRIPT *

Output from script TARGET_ALL_P2_bmp_housekeeping_visualization.R:

./TARGET_ALL_P2_bmp_filt_HK_DEG_boxplots.pdf:
- Boxplots of filtered gene expression values of housekeeping DEGs from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_HK_DEG_boxplots.pdf:
- Boxplots of raw gene expression values of housekeeping DEGs from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_hkg_eisenberg_consensus_DEG_density.pdf:
- Density plot of log2FC values of housekeeping DEGs from step 3 of * PROCESS IN THE SCRIPT *

