### This script performs DEA of TARGET-ALL-P2 gene expression data using the limma-voom pipeline
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered


# Load libraries
library(limma)
library(tidyverse)
library(biomaRt)

# Source functions
source("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/TARGET_ALL_P2_bmp_dea_limma_functions.R")


## Load data ---------------------------------------------------------------------------------------------------------

# Load filtered data 
ALL_P2_filt <- get(load("TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_filt_data.rda"))

# Load table containing batch information about samples
ALL_P2_batch_info <- read_csv(file = "TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")


## Wrangle data ---------------------------------------------------------------------------------------------------------

## Vector of batch factor designs
batch_factor_designs <- list.dirs("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom", 
                                  full.names = FALSE, 
                                  recursive = FALSE)


## Create design matrix for each DEA design differing in included batch factors

# Create design matrix for each batch factor design
design_matrix_list <- map(batch_factor_designs, function(batch_factor_designs) {
  
  create_design_matrix(batch_factor = batch_factor_designs, 
                       condition_types = ALL_P2_batch_info$subtype_abbr, 
                       tissue_portion_aliquot = ALL_P2_batch_info$tissue_portion_aliquot, 
                       years_of_diagnosis = as.character(ALL_P2_batch_info$diagnosis_year))
  
})

# Assign names to resulting list of design matrices
names(design_matrix_list) <- batch_factor_designs


## Voom transform data for each DEA design differing in included batch factors

voom_list <- map2(design_matrix_list, seq(batch_factor_designs), function(design_matrix_list, x) {
  
  file_name_voom_plot <- paste0("ALL_P2_bmp_voom_plot_", batch_factor_designs[x], ".png")
  dir_output_path <- paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x])
  
  voom_transform(file_name = file_name_voom_plot, 
                 dir_output = dir_output_path,
                 processed_data = ALL_P2_filt, 
                 design_matrix = design_matrix_list)
  
})


## Do DEA using limma-voom for each DEA design differing in included batch factors

DEA_limma_voom_list <- map2(design_matrix_list, voom_list, function(design_matrix_list, voom_list) {
  
  limma_DEA(design_matrix = design_matrix_list, 
            voom_object = voom_list)
  
})


## Convert ENSEMBL IDs to gene names in DEA matrix

DEA_limma_voom_list_gene <- map(seq(DEA_limma_voom_list), function(x) {
  
  # Convert full DEA matrix to tibble
  DEA_tbl <- DEA_limma_voom_list[[x]] %>%
    as_tibble(rownames = "ENSEMBL_ID")
  
  # Convert ENSEMBL IDs to gene names and add gene names to full DEA matrix
  ensembl_hugo_conversion(ensembl_id = DEA_tbl$ENSEMBL_ID) %>%
    as_tibble() %>%
    dplyr::rename("ENSEMBL_ID" = `Gene stable ID`) %>%
    full_join(x = .,
              y = DEA_tbl, 
              by = "ENSEMBL_ID")
  
})
names(DEA_limma_voom_list_gene) <- batch_factor_designs
  

## Find significant DEGs based on FDR 

DEA_significant_DEG <- map(DEA_limma_voom_list_gene, function(DEA_limma_voom_list_gene) {
  
  DEG_significant(DEA_table = DEA_limma_voom_list_gene,
                  fdr.cut = 0.05)
  
})


## Split significantly DEGs into up- and downregulated DEGs

DEG_sig_up <- map(DEA_significant_DEG, function(DEA_significant_DEG) {
  
  DEA_significant_DEG %>% 
    filter(logFC > 0) %>%
    dplyr::select(ENSEMBL_ID, `Gene name`)
  
})

DEG_sig_down <- map(DEA_significant_DEG, function(DEA_significant_DEG) {
  
  DEA_significant_DEG %>%
    filter(logFC < 0) %>%
    dplyr::select(ENSEMBL_ID, `Gene name`)
  
})


## Save data ---------------------------------------------------------------------------------------------------------

# Save voom transformed and DEA data for all DEA designs

save_data <- map(seq(batch_factor_designs), function(x) {
  
  # Save voom object
  saveRDS(object = voom_list[[x]], 
          file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_voom_", batch_factor_designs[x], ".rds"))
  
  # Save DEA results
  saveRDS(object = DEA_limma_voom_list_gene[[x]], 
          file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_full_", batch_factor_designs[x], ".rds"))
  
  # Save DEA of significant DEGs
  saveRDS(object = DEA_significant_DEG[[x]], 
          file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_sig_", batch_factor_designs[x], ".rds"))
  
  # Save list of significantly DEGs
  write_csv(x = DEA_significant_DEG[[x]] %>% 
              dplyr::select(ENSEMBL_ID, `Gene name`), 
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_", batch_factor_designs[x], ".csv"))
  
  # Save list of significantly upregulated DEGs
  write_csv(x = DEG_sig_up[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_", batch_factor_designs[x], ".csv"))
  
  # Save list of significantly downregulated DEGs
  write_csv(x = DEG_sig_down[[x]],
            file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_", batch_factor_designs[x], ".csv"))
  
})








