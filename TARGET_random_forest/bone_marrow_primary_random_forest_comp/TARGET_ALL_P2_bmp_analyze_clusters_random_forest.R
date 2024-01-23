### This script analyzes selected features from random forest performed on 
### clusters predicted from unsupervised clustering on TARGET ALL P2 gene 
### expression data


# Load libraries
library(tidyverse)
library(biomaRt)


# Source functions script
source("TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_random_forest_functions.R")


### Load data ---------------------------------------------------------------------------------------------------------

# Load results from applying random forest on B clusters
rf_B_cluster_results <- get(load("TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_B_clusters_random_forest_results.rda"))

# Load results from applying random forest on T clusters
rf_T_cluster_results <- get(load("TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_T_clusters_random_forest_results.rda"))

# Load consensus DEA table
DEA_consensus <- get(load("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda"))

# Load list of 103 housekeeping genes
DEG_eisenberg_hkg <- read_csv("TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/DEG_vs_eisenberg.csv")


### Wrangle data ---------------------------------------------------------------------------------------------------------

## Analyze results of random forest applied on B clusters

# Extract selected features from each seed run 
sel_features_B <- rf_B_cluster_results$RFResults$Sel.vars %>% map_dfr(. %>% 
                                                                        as_tibble(), 
                                                                      .id = "seed_run") %>%
  dplyr::rename(ENSEMBL_ID = value) %>%
  dplyr::mutate(seed_run = seed_run %>%
                  as.numeric())

# Convert ENSEMBL IDs to gene names
sel_features_B_gene <- ensembl_hugo_conversion(ensembl_id = sel_features_B$ENSEMBL_ID) %>%
  as_tibble() %>%
  dplyr::rename("ENSEMBL_ID" = `Gene stable ID`,
                "Gene_name" = `Gene name`) %>%
  full_join(x = .,
            y = sel_features_B, 
            by = "ENSEMBL_ID") %>%
  dplyr::arrange(seed_run)

# Compare random forest selected genes with DEA consensus and with housekeeping genes
rf_DEA_hkg_B <- left_join(x = sel_features_B_gene, 
                          y = DEA_consensus, 
                          by = "ENSEMBL_ID") %>% 
  dplyr::select(-c(AveExpr_limma_voom, t_limma_voom, P.Value_limma_voom, unshrunk.logFC_edgeR,
                   logCPM_edgeR, PValue_edgeR, baseMean_DESeq2, lfcSE_DESeq2, stat_DESeq2, pvalue_DESeq2,
                   Gene_name.y)) %>%
  dplyr::mutate(is_hkg = case_when(ENSEMBL_ID %in% DEG_eisenberg_hkg$ENSEMBL_ID ~ "is_hkg",
                                   TRUE ~ "is_not_hkg"))


## Analyze results of random forest applied on T clusters

# Extract selected features from each seed run 
sel_features_T <- rf_T_cluster_results$VarsSelect

# Convert ENSEMBL IDs to gene names
sel_features_T_gene <- ensembl_hugo_conversion(ensembl_id = sel_features_T) %>%
  as_tibble() %>%
  dplyr::rename("ENSEMBL_ID" = `Gene stable ID`,
                "Gene_name" = `Gene name`) 

# Compare random forest selected genes with DEA consensus and with housekeeping genes
rf_DEA_hkg_T <- left_join(x = sel_features_T_gene, 
                          y = DEA_consensus, 
                          by = "ENSEMBL_ID") %>% 
  dplyr::select(-c(AveExpr_limma_voom, t_limma_voom, P.Value_limma_voom, unshrunk.logFC_edgeR,
                   logCPM_edgeR, PValue_edgeR, baseMean_DESeq2, lfcSE_DESeq2, stat_DESeq2, pvalue_DESeq2,
                   Gene_name.y)) %>%
  dplyr::mutate(is_hkg = case_when(ENSEMBL_ID %in% DEG_eisenberg_hkg$ENSEMBL_ID ~ "is_hkg",
                                   TRUE ~ "is_not_hkg"))


### Save data ---------------------------------------------------------------------------------------------------------

## Save results from applying random forest on B clusters

# Save table containing all selected features from all seed runs including their gene names
write_csv(rf_DEA_hkg_B, file = "TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_selected_features_B_clusters_rf.csv")


## Save results from applying random forest on T clusters

# Save table containing all selected features from all seed runs including their gene names
write_csv(rf_DEA_hkg_T, file = "TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_selected_features_T_clusters_rf.csv")






