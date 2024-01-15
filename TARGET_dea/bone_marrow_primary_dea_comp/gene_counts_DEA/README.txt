* INTRODUCTION *

This directory, gene_counts_DEA, makes visualizations of gene counts for the three tools, DESeq2, 
edgeR and limma-voom. Visualizations are placed the folder ./visualizations.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse (version 1.3.2)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_dea_gene_counts.R:

Step 1: Load data
Results data from DESeq2, edgeR and limma-voom is loaded.

Step 2: Data wrangling
Preparation of results from the each of the tools individually
and in combination.

Step 3: Visualizations
Visualizations of gene count.

Step 4: Save results
All visualization plots are saved in the folder ./visualizations.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_dea_gene_counts.R is the script that should be run. 

TARGET_ALL_P2_bmp_dea_gene_counts.R should be run in the following way from the terminal:
Rscript TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/TARGET_ALL_P2_bmp_dea_gene_counts.R



* OUTPUT *

Output from script TARGET_ALL_P2_bmp_dea_gene_counts.R:

./visualizations/ALL_P2_deseq2_counts_log.pdf
- log2FoldChange vs log(baseMean) (gene counts) for DESeq2 results

./visualizations/ALL_P2_deseq2_counts_no_log.pdf
- log2FoldChange vs baseMean for (gene counts)  DESeq2 results

./visualizations/ALL_P2_edgeR_counts_log.pdf
- logFC vs logCPM (gene counts) for edgeR results

./visualizations/ALL_P2_edgeR_counts_no_log.pdf
- logFC vs 10**logCPM (gene counts) for edgeR results

./visualizations/ALL_P2_limma_counts_log.pdf
- logFC vs logCPM (gene counts) for limma-voom results

./visualizations/ALL_P2_limma_counts_no_log.pdf
- logFC vs 10**logCPM (gene counts) for limma-voom results

./visualizations/ALL_P2_comb_deseq2_edgeR.pdf
- edgeR vs DESeq2 gene counts (not filtered for DEGs)

./visualizations/ALL_P2_comb_edgeR_limma.pdf
- edgeR vs limma-voom gene counts (not filtered for DEGs)

./visualizations/ALL_P2_comb_limma_deseq2
- DESeq2 vs limma-voom gene counts (not filtered for DEGs)

./visualizations/ALL_P2_comb_DEG_deseq2_edgeR.pdf
- edgeR vs DESeq2 gene counts (filtered for DEGs)

./visualizations/ALL_P2_comb_DEG_edgeR_limma.pdf
- edgeR vs limma-voom gene counts (filtered for DEGs)

./visualizations/ALL_P2_comb_DEG_limma_deseq.pdf
- DESeq2 vs limma-voom gene counts (filtered for DEGs)



