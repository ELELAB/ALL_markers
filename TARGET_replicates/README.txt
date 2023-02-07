* INTRODUCTION *

This directory (TARGET_replicates) explores replicate samples present in
the TARGET-ALL-P2 project. Subsequently, TARGET-ALL-P2 data is adjusted
for replicates. 


* REQUIREMENTS *

Following R and R packages are needed:
- R 4.2
- R packages:
	* tidyverse (version 1.3.2)
	* ggrepel (version 0.9.2)
	* SummarizedExperiment (version 1.28.0)
 

* PROCESS IN THE SCRIPT *

Process in script TARGET_ALL_P2_bmp_replicates.R: 

Step 1: Load data
Data of TARGET-ALL-P2 project subsetted to contain bone marrow primary
samples is loaded. 

Step 2: Find replicate samples
Replicate samples, i.e. patients that have more than one sample, are
found by using patient barcodes. Afterwards, read counts of these
replicate samples are found. 

Step 3: Visualize data
-Total sum of read counts in replicate samples are plotted
-Total sum of read counts of all samples are plotted
-Replicate samples are highlighted in plot showing total sum
 of read counts of all samples
-Portion B samples are highlighted in plot showing total sum
 of read counts of all samples
-Portion B replicate samples are highlighted in plot showing
 total sum of read counts of all samples
-Number isolation from aliquot 2 samples are highlighted in
 plot showing total sum of read counts of all samples
-Number isolation from aliquot 2 replicate samples are highlighted
 in plot showing total sum of read counts of all samples

Step 4: Save data
Above data is saved. 


Process in script TARGET_ALL_P2_bmp_replicates_adj.R:

Step 1: Load data
Data of TARGET-ALL-P2 project subsetted to contain bone marrow primary
samples is loaded. Table containing patient barcodee, full barcode, 
and total read count of replicates is loaded. 

Step 2: Adjust data for replicates
Retain only the sample among replicates for a patient with the 
highest total read count. Subset ALL P2 bone marrow primary data to
contain only retained samples among replicates. 

Step 3: Save data
Above data is saved. 


* RUNNING THE SCRIPTS *

The user should be in the moonlight_ALL_paper directory and run the scripts
from there. 

TARGET_ALL_P2_bmp_replicates.R should be run first. 
TARGET_ALL_P2_bmp_replicates_adj.R should be run second. 

The scripts should be run in the following way from the terminal:
Rscript TARGET_replicates/TARGET_ALL_P2_bmp_replicates.R
Rscript TARGET_replicates/TARGET_ALL_P2_bmp_replicates_adj.R


* OUTPUT *

Output from script TARGET_ALL_P2_bmp_replicates.R:

./TARGET_ALL_P2_bmp_all_read_count.csv
- Table containing total read counts and barcode information of all
  samples from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_replicates.csv
- List of patients that are replicates from step 2 of * PROCESS
  IN THE SCRIPT *

./TARGET_ALL_P2_bmp_rep_read_count.csv
- Table containing total read counts and barcode information of 
  replicates from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_rep_read_count.pdf
- Barplot of total read count of replicates from step 3 of 
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_all_read_count.pdf
- Scatter plot of total read count of all samples from step 3 of 
  * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_all_rep_read_count.pdf
- Scatter plot of total read count of all samples, highlighting replicates
  from step 3 of * PROCESS IN THE SCRIPT * 

./TARGET_ALL_P2_bmp_portion_B_read_count.pdf
- Scatter plot of total read count of all samples, highlighting portion B
  samples from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_portion_B_rep_read_count.pdf
- Scatter plot of total read count of all samples, highlighting portion B
  samples and replicates from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_iso_2_read_count.pdf
- Scatter plot of total read count of all samples, highlighting isolation
  2 samples from step 3 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_iso_2_rep_read_count.pdf
- Scatter plot of total read count of all samples, highlighting isolation
  2 samples and replicates from step 3 of * PROCESS IN THE SCRIPT *


Output from script TARGET_ALL_P2_bmp_replicates_adj.R:

./TARGET_ALL_P2_bmp_kept_rep.csv
- List of kept replicate patients from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_discarded_rep.csv
- List of discarded replicate patients from step 2 of * PROCESS IN THE SCRIPT *

./TARGET_ALL_P2_bmp_exp_no_rep.rda
- ALL P2 bone marrow primary data where replicate samples have been filtered
  from step 2 of * PROCESS IN THE SCRIPT * 

