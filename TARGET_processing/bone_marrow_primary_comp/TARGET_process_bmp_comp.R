### This script processes TARGET-ALL-P2 gene expression data 
### The gene expression data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates 



# Load libraries
library(TCGAbiolinks)
library(tidyverse)
library(EDASeq)


## Load data ------------------------------------------------------------------

# Load gene expression data of TARGET-ALL-P2 bone marrow (primary) where replicates have been removed
bone_marrow_primary_exp_data <- get(load("TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda"))


## Wrangle data ------------------------------------------------------------------

# Preprocessing of gene expression data of TARGET-ALL-P2 bone marrow (primary)
ALL_P2_bmp_prep_data <- TCGAanalyze_Preprocessing(object = bone_marrow_primary_exp_data, 
                                                  cor.cut = 0.6,
                                                  filename = "TARGET_processing/bone_marrow_primary_comp/preprocessing_output_AAIC.png")

# Normalization of gene expression data of TARGET-ALL-P2 bone marrow (primary)
# This is done using updated geneInfoHT table
ALL_P2_bmp_norm_data <- TCGAanalyze_Normalization(tabDF = ALL_P2_bmp_prep_data, 
                                                  geneInfo = geneInfoHT, 
                                                  method = "gcContent")

# Filtering of gene expression data of TARGET-ALL-P2 bone marrow (primary)
ALL_P2_bmp_filt_data <- TCGAanalyze_Filtering(tabDF = ALL_P2_bmp_norm_data, 
                                              method = "quantile",
                                              qnt.cut = 0.25)


## Save data ------------------------------------------------------------------

# Save preprocessing matrix
save(ALL_P2_bmp_prep_data, 
     file = "TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_prep_data.rda")

# Save normalization matrix
save(ALL_P2_bmp_norm_data, 
     file = "TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_norm_data.rda")

# Save filtering matrix
save(ALL_P2_bmp_filt_data, 
     file = "TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_filt_data.rda")




