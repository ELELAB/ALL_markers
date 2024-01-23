* INTRODUCTION *

This directory, bone_marrow_primary_drug_targets_comp, investigates
if any of the predicted subtype-related and cluster-related gene 
expression marers are annotated as drug targets in the Drug Gene 
Interaction Database (DGIdb).


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse	(version 1.3.2)
	* rDGIdb		(version 1.24.0)



* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_drug_targets_markers.R:

Step 1: Load data
Defined set of subtype-related and cluster-related gene 
expression markers are loaded.

Step 2: Drug target investigation
Gene expression markers are queried in DGIdb and results
are extracted.

Step 3: Visualization
Results from drug target investigation are visualized.

Step 4: Save results
Above results are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_drug_targets_markers.R is the script that should be run. 

TARGET_ALL_P2_bmp_drug_targets_markers.R should be run in the following way from the terminal:
Rscript TARGET_drug_targets/bone_marrow_primary_drug_targets_comp/TARGET_ALL_P2_bmp_drug_targets_markers.R


* OUTPUT *

Output from TARGET_ALL_P2_bmp_drug_targets_markers.R:

./TARGET_ALL_P2_bmp_drugs_cluster_genes.csv:
- Results of drug-gene interactions of cluster-related gene expression
  markers from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_drugs_subtype_genes.csv:
- Results of drug-gene interactions of subtype-related gene expression
  markers from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_drugs_subtype_genes_heatmap.pdf:
- Heatmap showing drug-gene interactions of subtype-related gene expression
  markers from step 3 of * PROCESS IN THE SCRIPT *

