### This script visualizes gene contributions to PCA dimensions
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered
### and if batch corrected data is used:
### 6) voom transformed
### 7) batch corrected for year of diagnosis


# Load libraries 
library(tidyverse)
library(viridis)


## Load data ---------------------------------------------------------------------------------------------------------

# Load top 20 genes contributing to PC 1 and 2 with their annotated biotype
gene_contrib <- read_csv("TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/gene_contributions_pca_biotypes/TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes.csv")

# Load consensus DEA table
DEA_consensus <- get(load("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda"))

# Load limma-voom year of diagnosis DEA results
DEA_limma_voom <- readRDS("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/years_batch/ALL_P2_bmp_DEA_sig_years_batch.rds")

# Load edgeR year of diagnosis DEA results
DEA_edger <- readRDS("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/years_batch/ALL_P2_bmp_DEA_sig_years_batch.rds")

# Load DESeq2 year of diagnosis DEA results
DEA_deseq2 <- readRDS("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/years_batch/ALL_P2_bmp_DEA_sig_years_batch.rds")


## Wrangle data ---------------------------------------------------------------------------------------------------------

# Add column combining ENSEMBL gene ID and gene name to gene_contrib table
# For plotting purposes below 
gene_contrib_biotypes <- gene_contrib %>%
  dplyr::rename("gene_name" = `Gene name`) %>%
  mutate(ENSEMBL_gene_name = case_when(is.na(gene_name) ~ ENSEMBL_ID,
                                       TRUE ~ str_c(ENSEMBL_ID,
                                                    "\n",
                                                    gene_name)))

# Compare genes contributing to PC1 and 2 with consensus DEGs
gene_contrib_DEA <- left_join(x = gene_contrib_biotypes, 
                              y = DEA_consensus, 
                              by = "ENSEMBL_ID") %>% 
  dplyr::select(-Gene_name) %>% 
  mutate(sign_DEG = case_when(logFC_mean > 0 ~ "up consensus DEG",
                              logFC_mean < 0 ~ "down consensus DEG",
                              TRUE ~ "non-consensus DEG"))

# Checking if the non-consensus DEGs among genes contributing to PC2 are DEGs by limma-voom, edgeR or DESeq2

# Intersect genes contributing to PC1 and 2 with limma-voom DEGs
gene_contrib_limma_voom <- intersect(gene_contrib_DEA$ENSEMBL_ID, 
                                     DEA_limma_voom$ENSEMBL_ID)

# Intersect genes contributing to PC1 and 2 with edgeR DEGs
gene_contrib_edger <- intersect(gene_contrib_DEA$ENSEMBL_ID, 
                                DEA_edger$ENSEMBL_ID)

# Intersect genes contributing to PC1 and 2 with DESeq2 DEGs
gene_contrib_deseq2 <- intersect(gene_contrib_DEA$ENSEMBL_ID, 
                                 DEA_deseq2$ENSEMBL_ID)

# Add annotation to gene_contrib_DEA table if genes are found as DEGs by the three DEA tools
gene_contrib_DEA <- gene_contrib_DEA %>%
  mutate(DEG_limma_voom = case_when(ENSEMBL_ID %in% gene_contrib_limma_voom ~ "Is DEG limma-voom"),
         DEG_edger = case_when(ENSEMBL_ID %in% gene_contrib_edger ~ "Is DEG edgeR"),
         DEG_deseq2 = case_when(ENSEMBL_ID %in% gene_contrib_deseq2 ~ "Is DEG deseq2")) 


## Visualize data ---------------------------------------------------------------------------------------------------------

# Barplot of top 40 genes contributing to PC dimensions 1 stratified by biotype
gene_contrib_plot <- ggplot(data = gene_contrib_DEA %>%
                              filter(PCA_dim == "PC_1"), 
                            mapping = aes(x = reorder(ENSEMBL_gene_name, 
                                                      -contrib),
                                          y = contrib,
                                          fill = biotype)) +
  geom_col() +
  geom_text(mapping = aes(label = sign_DEG), 
            position = position_dodge(width = 0.9), 
            hjust = 1.05,
            angle = 90,
            size = 8,
            color = "white") +
  scale_fill_manual(values = viridis_pal(option = "D")(gene_contrib_DEA %>% 
                                                         distinct(biotype) %>% 
                                                         nrow)) + 
  labs(x = "",
       y = "Contributions (%)",
       title = "Contributions of top 40 genes to PC 1") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 18, 
                                    color = "black"),
        axis.text.y = element_text(size = 18, 
                                   color = "black"),
        axis.text.x = element_text(size = 18, 
                                   color = "black", 
                                   angle = 90, 
                                   hjust = 0),
        axis.ticks = element_line(size = 2),
        title = element_text(size = 18),
        panel.border = element_rect(size = 2),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.3, "cm")) 



## Save data ---------------------------------------------------------------------------------------------------------

# Save table containing genes contributing to PC dimensions 1 and 2 including comparison with consensus DEA
write_csv(gene_contrib_DEA, file = "TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/gene_contributions_pca_biotypes/TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes_info.csv")

# Save barplot containing top 40 genes contributing to PC 1 
ggsave(filename = "TARGET_ALL_P2_bmp_pca_gene_contrib_biotypes_plot.pdf", 
       plot = gene_contrib_plot, 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/gene_contributions_pca_biotypes/", 
       width = 26, 
       height = 12)










