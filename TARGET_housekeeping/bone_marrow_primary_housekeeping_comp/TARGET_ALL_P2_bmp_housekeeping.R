### This script compares the consensus DEGs found by limma-voom, edgeR and DESeq2
### with two lists of housekeeping genes. The housekeeping gene lists are based on:
### List 1: Eisenberg and Levanon, 2013, doi: 10.1016/j.tig.2013.05.010
### List 2: Homemade list based on Gupta et al., 2022, doi: 10.1007/s11033-022-07337-w,
###                                Villegas-Ruiz et al., 2019, doi: 10.3390/genes10050376,
###                                Caracausi et al., 2017, doi: 10.3892/mmr.2017.6944

# Load libraries
library(tidyverse)
library(viridis)


## Load data ------------------------------------------------------------------

# Load consensus DEG data
dea_consensus <- get(load("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda"))

# Load Eisenberg and Levanon HKGs
HKG_eisenberg_data <- read_tsv("TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/HKG_lists/HKG_Eisenberg_2013.txt", col_names = FALSE)

# Load home made list with bone marrow and ALL HKGs
HKG_homemade_data <- read_csv("TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/HKG_lists/HKG_homemade.csv")


## Wrangle data ----------------------------------------------------------------

# Select gene lists
DEGs <- dea_consensus %>%
  select(ENSEMBL_ID,
         Gene_name)

HKG_eisenberg <- HKG_eisenberg_data %>%
  select(X1) %>%
  dplyr::rename(Gene_name = X1)

HKG_homemade <- HKG_homemade_data %>%
  select(gene) %>%
  dplyr::rename(Gene_name = gene)


## Comparisons -----------------------------------------------------------------

# Overlap between Eisenberg and homemade gene lists
eisenberg_overlap_homemade <- inner_join(HKG_eisenberg, HKG_homemade, by = "Gene_name") %>% # 4
  arrange(Gene_name)

# Overlap between DEGs and Eisenberg
DEG_overlap_eisenberg <- inner_join(DEGs, HKG_eisenberg, by = "Gene_name") %>% # 103
  arrange(Gene_name)

# Overlap betewen DEGs and homemade
DEG_overlap_homemade <- inner_join(DEGs, HKG_homemade, by = "Gene_name") %>% # 0
  arrange(Gene_name)

# Add log2FC values to DEGs that overlap with Eisenberg
DEG_overlap_eisenberg <- DEG_overlap_eisenberg %>%
  inner_join(x = ., y = dea_consensus, by = "ENSEMBL_ID") %>%
  dplyr::select(ENSEMBL_ID, Gene_name.x, logFC_limma_voom, adj.P.Val_limma_voom, logFC_edgeR, 
         FDR_edgeR, log2FoldChange_DESeq2, padj_DESeq2,
         logFC_mean, logFC_sd, FDR_mean, FDR_sd)
  


## Save results ----------------------------------------------------------------

# Save HKG overlap
write_csv(x = eisenberg_overlap_homemade,
          file = "TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/eisenberg_vs_homemade.csv")

# Save DEG overlap with Eisenberg
write_csv(x = DEG_overlap_eisenberg,
          file = "TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/DEG_vs_eisenberg.csv")

# Save DEG overlap with homemade
write_csv(x = DEG_overlap_homemade,
          file = "TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/DEG_vs_homemade.csv")






