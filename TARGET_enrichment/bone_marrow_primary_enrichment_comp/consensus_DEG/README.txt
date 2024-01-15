* INTRODUCTION *

This directory, consensus_DEG, performs enrichment analysis 
of identified consensus differentially expressed genes. 


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse (version 1.3.2)
	* enrichR   (version 3.1)
	* patchwork (version 1.1.2)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_DEG_enrichment.R:

Step 1: Load data
Consensus DEG data is loaded.

Step 2: Enrichment analysis
Enrichment analysis is performed on up- and downregulated
consensus DEGs.

Step 3: Visualizations
Enrichment analysis results are visualized.

Step 4: Save results
Enrichment analysis results and visualizations are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_DEG_enrichment.R is the script that should be run. 

TARGET_ALL_P2_bmp_DEG_enrichment.R should be run in the following way from the terminal:
Rscript TARGET_enrichment/bone_marrow_primary_enrichment_comp/consensus_DEG/TARGET_ALL_P2_bmp_DEG_enrichment.R

TARGET_ALL_P2_bmp_DEG_enrichment.R sources the functions script:
TARGET_enrichment/bone_marrow_primary_enrichment_comp/TARGET_ALL_P2_bmp_enrichment_functions.R


* OUTPUT *

./TARGET_ALL_P2_bmp_enrichment_analyses_DEG.rds
- Enrichment analysis results of consensus DEGs for all used databases
  from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_enrichment_analyses_DEG_up.rds
- Enrichment analysis results of upregulated consensus DEGs for all used
  databases from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_enrichment_analyses_DEG_down.rds
- Enrichment analysis results of downregulated consensus DEGs for all used
  databases from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_consensus_DEG_enrich.pdf
- Enrichment analysis visualization of consensus DEGs from step 3 of 
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_consensus_DEG_up_enrich.pdf
- Enrichment analysis visualization of upregulated consensus DEGs from
  step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_consensus_DEG_down_enrich.pdf
- Enrichment analysis visualization of downregulated consensus DEGs from
  step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_consensus_DEG_up_down_enrich.pdf
- Enrichment analysis visualization of up- and downregulated consensus DEGs from
  step 3 of * PROCESS IN THE SCRIPT *
