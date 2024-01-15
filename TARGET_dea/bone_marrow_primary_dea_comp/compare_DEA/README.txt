* INTRODUCTION *

This directory (compare_DEA) compares differential expression
analysis results found using three tools: DESeq2, edgeR
and limma-voom.  
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
        * viridis   (version 0.6.2)


* PROCESS IN THE SCRIPT *

Process in the script TARGET_ALL_P2_bmp_dea_compare.R:

Step 1: Load data
Results of differential expression analyses performed
by the three different tools are loaded.

Step 2: Wrangle data
Tables showing the number of identified differentially
expressed genes and their sign of regulation are created.
A consensus table of genes identified as differentially
expressed across all three methods is created.

Step 3: ANOVA
An ANOVA is performed to compare log2FC values for
differentially expressed genes in the consensus list
of differentially expressed genes.

Step 4: Visualize data
log2FC values of consensus differentially expressed genes
are visualized in boxplots.

Step 5: Save data
Results from the above steps are saved. 


* RUNNING THE SCRIPT * 

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_dea_compare.R is the script that should be run. 

TARGET_ALL_P2_bmp_dea_compare.R should be run in the following way from the terminal:
Rscript TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/TARGET_ALL_P2_bmp_dea_compare.R


* OUTPUT * 

./ALL_P2_bmp_DEG_numbers.csv
- Table containing number of identified differentially expressed genes
  from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_up_down_DEG_numbers.csv
- Table containing number of up- and downregulated differentially
  expressed genes from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_DEA_consensus.rda
- Matrix containing results of the consensus differential expression
  analysis from step 2 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_logFC_DEG_common_year_boxplot.pdf
- Boxplot of log2FC values of consensus differentially expressed
  genes across the three differential expression analysis methods
  from step 4 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_anova_logFC.csv
- ANOVA results from step 3 of * PROCESS IN THE SCRIPT *

./ALL_P2_bmp_qq_residuals.pdf
- qq-plot for pooled residuals from step 3 of * PROCESS IN THE
  SCRIPT *
