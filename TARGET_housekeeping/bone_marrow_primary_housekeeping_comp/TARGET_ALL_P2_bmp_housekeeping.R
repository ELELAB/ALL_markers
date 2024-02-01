### This script compares the consensus DEGs found by limma-voom, edgeR and DESeq2
### with a housekeeping gene list from: Eisenberg and Levanon, 2013, doi: 10.1016/j.tig.2013.05.010

# Load libraries
library(tidyverse)
library(viridis)


## Load data ------------------------------------------------------------------

# Load consensus DEG data
dea_consensus <- get(load("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda"))

# Load Eisenberg and Levanon HKGs
HKG_eisenberg_data <- read_tsv("TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/HKG_lists/HKG_Eisenberg_2013.txt", col_names = FALSE)

# Load batch corrected gene count data (corrected for year of diagnosis batch factor)
ALL_P2_corr <- get(load("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda"))

# Load SummarizedExperiment of expression data
ALL_P2_se <- get(load("TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda"))


## Wrangle data ----------------------------------------------------------------

# Select gene lists
DEGs <- dea_consensus %>%
  dplyr::select(ENSEMBL_ID,
                Gene_name)

HKG_eisenberg <- HKG_eisenberg_data %>%
  dplyr::select(X1) %>%
  dplyr::rename(Gene_name = X1)


## Comparisons -----------------------------------------------------------------

# Overlap between DEGs and Eisenberg
DEG_overlap_eisenberg <- inner_join(DEGs, HKG_eisenberg, by = "Gene_name") %>% # 103
  dplyr::arrange(Gene_name)

# Add log2FC values to DEGs that overlap with Eisenberg
DEG_overlap_eisenberg <- DEG_overlap_eisenberg %>%
  inner_join(x = ., y = dea_consensus, by = "ENSEMBL_ID") %>%
  dplyr::select(ENSEMBL_ID, Gene_name.x, logFC_limma_voom, adj.P.Val_limma_voom, logFC_edgeR, 
                FDR_edgeR, log2FoldChange_DESeq2, padj_DESeq2,
                logFC_mean, logFC_sd, FDR_mean, FDR_sd)

# Intersect housekeeping genes from Eisenberg with total number of genes in our dataset
ensembl_gene <- tibble("ENSEMBL" = ALL_P2_se@rowRanges@elementMetadata@listData[["ensembl_gene_id"]], 
                       "Gene" = ALL_P2_se@rowRanges@elementMetadata@listData[["external_gene_name"]]) %>% 
  dplyr::filter(ENSEMBL %in% rownames(ALL_P2_corr))
intersect_hkg_total <- intersect(HKG_eisenberg$Gene_name, ensembl_gene$Gene) %>%
  as_tibble()


## Save results ----------------------------------------------------------------

# Save DEG overlap with Eisenberg
write_csv(x = DEG_overlap_eisenberg,
          file = "TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/DEG_vs_eisenberg.csv")

# Save intersection between housekeeping genes from Eisenberg and total genes in our dataset
write_csv(x = intersect_hkg_total, 
          file = "TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/eisenberg_vs_dataset_genes.csv")


