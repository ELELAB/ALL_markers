### This script investigates lost genes following normalization step 
### The processing workflow of the gene expression data is:  
### 1) subsetting to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusting for replicates
### 3) preprocessing
### 4) normalizing 
### 5) filtering 


# Load libraries
library(TCGAbiolinks)
library(tidyverse)




## Load data ------------------------------------------------------------------

# Load gene expression data of TARGET-ALL-P2 bone marrow (primary) where replicates have been removed
bone_marrow_primary_exp_data <- get(load("TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda"))

# Load preprocessed data 
ALL_P2_bmp_prep_data <- get(load("TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_prep_data.rda"))

# Load normalized data
ALL_P2_bmp_norm_data <- get(load("TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_norm_data.rda"))


## Wrangle data ------------------------------------------------------------------

# Create tibble containing set difference, i.e. excluded genes, between preprocessing and normalization matrices
excl_genes <- tibble(anti_join(x = as_tibble(rownames(ALL_P2_bmp_prep_data)), 
                               y = as_tibble(rownames(ALL_P2_bmp_norm_data)),
                               by = "value")) %>%
  dplyr::rename("ensembl_id" = "value")

# Double check that excluded genes are not in geneInfoHT table
excl_genes <- excl_genes %>% 
  mutate("in_geneInfoHT" = excl_genes$ensembl_id %in% rownames(geneInfoHT))


## Save data ------------------------------------------------------------------
write_csv(x = excl_genes, file = "TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_norm_excl_genes.csv")








