* INTRODUCTION *

This directory (limma_voom) contains scripts and results 
from differential expression analysis by limma-voom.
The input gene expression data has been subject to the
following workflow prior to the differential expression
analysis:
1) subsetted to contain only bone marrow (primary blood derived cancer) samples
2) adjusted for replicates
3) preprocessed
4) normalized using updated geneInfoHT table
5) filtered


* REQUIREMENTS * 

Following R and R packages are needed:
- R 4.2.0
- R packages:
        * tidyverse (version 1.3.1)
	* limma (version 3.54.0)
        * biomaRt (version 2.54.0)


* PROCESS IN THE SCRIPT *

Process in the script TARGET_ALL_P2_bmp_dea_limma.R:

Step 1: Load data and load functions
Filtered data and batch factor data are loaded. 
The functions script TARGET_ALL_P2_bmp_dea_limma_functions.R
is loaded.

Step 2: Wrangle data
Batch factor design matrices are created and data is
voom transformed.

Step 3: Differential expression analysis
Differential expression analysis is run for each batch design. 
Subsequently, the ENSEMBL IDs in the results matrix are converted 
to external gene names with the function.
Submatrices are created for significantly expressed genes, 
upregulated genes and downregulated genes for each batch design.

Step 4: Save data
Results are saved for each batch design. 


* RUNNING THE SCRIPT * 

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_dea_limma.R is the script that should be run. 

TARGET_ALL_P2_bmp_dea_limma.R should be run in the following way from the terminal:
Rscript TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/TARGET_ALL_P2_bmp_dea_limma.R

TARGET_ALL_P2_bmp_dea_limma.R sources the functions script: 
TARGET_ALL_P2_bmp_dea_limma_functions.R


* OUTPUT * 

For each batch design, the following output is produced:

./<batch_design>/ALL_P2_bmp_DEA_full_<batch_design>.rds
- DEA results for all genes from step 3 of * PROCESS IN THE SCRIPT *

./<batch_design>/ALL_P2_bmp_DEA_sig_<batch_design>.rds
- DEA results for genes with significant fdr values from step 3 
of * PROCESS IN THE SCRIPT *

./<batch_design>/ALL_P2_bmp_DEG_sig_<batch_design>.csv
- List of ENSEMBL IDs and external gene names for genes with
significant fdr values from step 3 of * PROCESS IN THE SCRIPT *

./<batch_design>/ALL_P2_bmp_DEG_sig_up_<batch_design>.csv
- List of ENSEMBL IDs and external gene names for upregulated genes with 
significant fdr values from step 3 of * PROCESS IN THE SCRIPT *

./<batch_design>/ALL_P2_bmp_DEG_sig_down_<batch_design>.csv
- List of ENSEMBL IDs and external gene names for downregulated genes with 
significant fdr values from step 3 of * PROCESS IN THE SCRIPT *

./<batch_design>/ALL_P2_bmp_voom_<batch_design>.rds
- Voom transformed data from step 2 of * PROCESS IN THE SCRIPT *

./<batch_design>/ALL_P2_bmp_voom_plot_<batch_design>.png
- Plot from voom transformation from step 2 of * PROCESS IN THE SCRIPT *
