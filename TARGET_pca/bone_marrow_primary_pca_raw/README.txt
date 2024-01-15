* INTRODUCTION *

This directory (bone_marrow_primary_pca_raw) performs MDS analysis of 
raw TARGET-ALL-P2 gene expression data where this data has been subject to
the following:
1) subset to contain bone marrow primary samples
2) voom transformed 

Replicates are also highlighted in these MDS analyses. 


* REQUIREMENTS * 

Following R and R packages are needed:
- R 4.2
- R packages:
        * tidyverse (version 1.3.2)
	* SummarizedExperiment (version 1.28.0)
        * patchwork (version 1.1.2)
	* ggrepel (version 0.9.2)


* PROCESS IN THE SCRIPT * 

Step 1: Load data 
TARGET-ALL-P2 data subsetted to contain only bone marrow primary samples is loaded. 
Voom transformed raw read counts of TARGET-ALL-P2 bone marrow subsetted data is 
also loaded. Clinical data of TARGET-ALL-P2 project is also loaded. Information
about replicates is also loaded. 

Step 2: Create table of sample information
A table containing barcodes, subtype, year of diagnosis, tissue portion, number
isolation from aliquot, and replicates information is created. This table is
needed for MDS analysis and is subsequently saved. 

Step 3: Perform MDS analysis
MDS analyses of voom transformed raw bone marrow primary data is performed 
coloring for subtypes, year of diagnosis, tissue portion, number isolation
from aliquot, and for replicates are performed. Additionally, MDS analysis
of voom transformed raw bone marrow primary data subsetted to only replicate
samples is also performed. Subtype and variable plots are plotted side by
side. Additionally, replicate plots are plotted side by side. All plots
are saved. 


* RUNNING THE SCRIPT * 

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_raw_pca.R is the script that should be run. This script sources the
functions script TARGET_ALL_P2_bmp_raw_pca_functions.R. 

TARGET_ALL_P2_bmp_raw_pca.R should be run in the following way from the terminal:
Rscript TARGET_pca/bone_marrow_primary_pca_raw/TARGET_ALL_P2_bmp_raw_pca.R


* OUTPUT *

./TARGET_ALL_P2_bmp_batch_info.csv
- Information about samples from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_subtype_MDS.pdf
- MDS plot coloring for subtypes from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_year_MDS.pdf
- MDS plot coloring for year of diagnosis from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_portion_MDS.pdf
- MDS plot coloring for portion from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_isolation_MDS.pdf
- MDS plot coloring for number isolation from aliquot from step 3 of 
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_MDS.pdf
- MDS plot coloring for replicate samples from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_rep_MDS.pdf
- MDS plot of only replicate samples from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_subtype_year_MDS.pdf
- MDS plots of subtype and year of diagnosis side by side from step 3 of 
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_subtype_portion_MDS.pdf
- MDS plots of subtype and portion side by side from step 3 of * PROCESS IN
  THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_subtype_isolation_MDS.pdf
- MDS plots of subtype and number isolation from aliquot from step 3 of 
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_raw_rep_side_MDS.pdf
- MDS plots of replicates side by side from step 3 of * PROCESS IN THE 
  SCRIPT *


