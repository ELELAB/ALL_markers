### This script voom transforms TARGET-ALL-P2 gene expression data 
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples 
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered 



## Load libraries 
library(tidyverse)
library(limma)


## Load data  ------------------------------------------------------------------

# Load in filtered data (filtering done on normalized data)
ALL_P2_bmp_filt <- get(load("TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_filt_data.rda"))


## Wrangle data  ------------------------------------------------------------------

# Voom transform filtered data
ALL_P2_bmp_voom <- voom(ALL_P2_bmp_filt)


## Save data  ------------------------------------------------------------------

save(ALL_P2_bmp_voom, file = "TARGET_transform/voom_transform_bmp_comp/TARGET_ALL_P2_bmp_voom.rda")



