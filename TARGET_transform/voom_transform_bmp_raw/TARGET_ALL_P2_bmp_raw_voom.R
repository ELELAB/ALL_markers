### This script voom transforms TARGET-ALL-P2 gene expression data 
### This data has been 
### subsetted to contain only bone marrow (primary blood derived cancer) samples 


## Load libraries 
library(tidyverse)
library(limma)


## Load data  ------------------------------------------------------------------

# Load in raw read counts of TARGET-ALL-P2 bone marrow (primary) samples
ALL_P2_bone_marrow_primary_raw <- get(load("TARGET_data/ALL_P2_bmp_raw_read_counts.rda"))


## Wrangle data  ------------------------------------------------------------------

# Voom transform filtered data
ALL_P2_bmp_raw_voom <- voom(ALL_P2_bone_marrow_primary_raw)


## Save data  ------------------------------------------------------------------
save(ALL_P2_bmp_raw_voom, file = "TARGET_transform/voom_transform_bmp_raw/TARGET_ALL_P2_bmp_raw_voom.rda")
