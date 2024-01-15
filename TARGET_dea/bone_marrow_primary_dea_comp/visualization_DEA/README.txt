* INTRODUCTION *

This directory, visualization_DEA, makes visualizations of differential
expression analysis results found using the three tools, DESeq2, 
edgeR and limma-voom. Visualizations are placed in the respective
subfolder: ./all_dea_pipelines, ./deseq2, ./edgeR and ./limma_voom.


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2.0 
- R packages:
	* tidyverse		(version 1.3.2)
	* UpSetR    		(version 1.4.0)
	* ComplexHeatmap 	(version 2.14.0)
	* TCGAbiolinks		(version 2.25.3)
	* VennDiagram		(version 1.7.3)
	* viridis		(version 0.6.2)
	* biomaRt		(version 2.54.0)


* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_dea_visualization.R:

Step 1: Load data
Results data from DESeq2, edgeR and limma-voom is loaded.

Step 2: Data wrangling
Preparation of results from the each of the tools individually
and in combination.

Step 3: Visualizations
Visualizations of identified differentially expressed genes
from different tools in UpSet plots.

Step 4: Save results
All visualization plots are saved in the respective subfolder:
./all_dea_pipelines, ./deseq2, ./edgeR and ./limma_voom.


* RUNNING THE SCRIPT *

The user should be in the ALL_markers directory and run the script from there. 

TARGET_ALL_P2_bmp_dea_visualization.R is the script that should be run. 

TARGET_ALL_P2_bmp_dea_visualization.R should be run in the following way from the terminal:
Rscript TARGET_dea/bone_marrow_primary_dea_comp/visualization_DEA/TARGET_ALL_P2_bmp_dea_visualization.R

TARGET_ALL_P2_bmp_dea_visualization.R sources the functions script:
TARGET_ALL_P2_bmp_dea_visualization_functions.R



* OUTPUT *

Output from script TARGET_ALL_P2_bmp_dea_visualization.R:

./all_dea_pipelines:
- UpSet plots comparing identified differentially expressed genes
  across different differential expression analysis tools

./deseq2:
- UpSet plots comparing identified differentially expressed genes
  from DESeq2

./edgeR:
- UpSet plots comparing identified differentially expressed genes
  from edgeR

./limma_voom:
- UpSet plots comparing identified differentially expressed genes
  from limma-voom


