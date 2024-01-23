* INTRODUCTION *

This directory, bone_marrow_primary_compare_comp, contains 
scripts that compare ENSEMBL gene IDs and gene symbols 
selected across various methods and further compares these
with genes reported in the Network of Cancer Genes (NCG) 
database. 


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse		(version 1.3.2)
	* UpSetR		(version 1.4.0)
	* ComplexHeatmap	(version 2.14.0)
	* viridis		(version 0.6.2)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_compare.R:

Step 1: Load data
Results and data from DEA, housekeeping analysis, 
elastic net logistic regression, PCA, and NCG are
loaded. 

Step 2: Compare ENSEMBL gene IDs and gene symbols
ENSEMBL gene IDs and gene symbols are compared 
across different combinations of methods.

Step 3: Visualization
Intersections of different combinations of gene
lists are visualized in UpSet plots.

Step 4: Save results
Above results are saved.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the scripts from there. 

TARGET_ALL_P2_bmp_compare.R is the script that should be run first. 

TARGET_ALL_P2_bmp_compare.R should be run in the following way from the terminal:
Rscript TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_compare.R

TARGET_ALL_P2_bmp_compare.R sources the functions script:
TARGET_ALL_P2_bmp_compare_functions.R


* OUTPUT *

Output from TARGET_ALL_P2_bmp_compare.R:

./TARGET_ALL_P2_bmp_compare_ensembl_ids_upset.pdf:
- UpSet plot showing intersections of ENSEMBL gene IDs between consensus DEGs, 
  HKG DEGs, elastic net genes and PCA genes from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_compare_symbols_upset.pdf:
- UpSet plot showing intersections of gene symbols between consensus DEGs, 
  HKG DEGs, elastic net genes, PCA genes and NCG from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_<comparison>.csv:
- Intersections of ENSEMBL gene IDs between a given combination of methods from
  step 2 of * PROCESS IN THE SCRIPT *

./<comparison>_symbols.csv
- Intersections of gene symbols between a given combination of methods from
  step 2 of * PROCESS IN THE SCRIPT * 

