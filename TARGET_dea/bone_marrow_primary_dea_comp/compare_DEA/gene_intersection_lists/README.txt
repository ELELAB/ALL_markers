* INTRODUCTION *

This directory (gene_intersection_lists) intersects
identified differentially expressed genes across
three tools: DESeq2, edgeR and limma-voom.  
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
        * tidyverse (version 1.3.2)


* PROCESS IN THE SCRIPT *

Process in the script TARGET_ALL_P2_bmp_dea_gene_intersection_lists.R:

Step 1: Load data
Results of differential expression analyses performed
by the three different tools are loaded.

Step 2: Wrangle data
Intersections of identified differentially expressed
genes comparing different combinations of methods
are found.

Step 3: Save data
Results of gene intersections are saved.


* RUNNING THE SCRIPT * 

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_dea_gene_intersection_lists.R is the script that should be run. 

TARGET_ALL_P2_bmp_dea_gene_intersection_lists.R should be run in the following way from the terminal:
Rscript TARGET_dea/bone_marrow_primary_dea_comp/gene_intersection_lists/TARGET_ALL_P2_bmp_dea_gene_intersection_lists.R


* OUTPUT * 

Gene intersections acorss different combinations of methods
and for each batch factor design are saved in the respective
batch factor design folder:

./<batch_factor>/ALL_P2_bmp_DEG_sig_<comparison>_<batch_factor>.csv
- Intersections of identified differentially expressed genes
  for a given comparison between methods and for a given
  batch factor design from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_numbers_sig.csv
- Table showing numbers of all significant genes for each
  method from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_numbers_sig_down.csv
- Table showing numbers of significantly downregulated
  genes for each method from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_numbers_sig_up.csv
- Table showing numbers of significantly upregulated genes
  for each method from step 2 of * PROCESS IN THE SCRIPT *
