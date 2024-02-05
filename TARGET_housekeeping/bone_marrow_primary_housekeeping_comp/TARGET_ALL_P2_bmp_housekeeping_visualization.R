### This script visualizes in boxplots gene expression values of 103 housekeeping DEGs which are the
### overlap between housekeeping genes reported by Eisenberg and Levanon, 2013, 
### doi: 10.1016/j.tig.2013.05.010 and our consensus list of DEGs

# Load libraries
library(tidyverse)
library(gridExtra)
library(viridis)
library(SummarizedExperiment)

# Source functions script
source("TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/TARGET_ALL_P2_bmp_housekeeping_functions.R")


## Load data ---------------------------------------------------------------------------------------------------------

# Load consensus DEG data
dea_consensus <- get(load("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda"))

# Load list of 103 consensus housekeeping DEGs
HK_DEG <- read_csv("TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/DEG_vs_eisenberg.csv")

# Load filtered gene expression data 
ALL_P2_bmp_filt_data <- get(load("TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_filt_data.rda"))

# Load raw gene expression data 
ALL_P2_bmp_raw_data <- get(load("TARGET_replicates/TARGET_ALL_P2_bmp_exp_no_rep.rda"))

# Load info of samples
ALL_P2_info <- read_csv("TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")


## Wrangle data ---------------------------------------------------------------------------------------------------------

## Filter expression data to include only prefiltered genes
ALL_P2_bmp_raw_data_filt <- assays(ALL_P2_bmp_raw_data[rownames(ALL_P2_bmp_raw_data) %in% rownames(ALL_P2_bmp_filt_data),])$"HTSeq - Counts"

# Extract barcodes of samples belonging to B-ALL subtype
B_barcodes <- ALL_P2_info %>% 
  dplyr::filter(subtype_abbr == "B-cell ALL") %>% 
  dplyr::select(barcodes)

# Extract barcodes of samples belonging to T-ALL subtype
T_barcodes <- ALL_P2_info %>% 
  dplyr::filter(subtype_abbr == "T-cell ALL") %>% 
  dplyr::select(barcodes)

# Subset filtered gene expression data to contain only the 103 HK DEGs
ALL_P2_bmp_filt_hkg <- ALL_P2_bmp_filt_data[HK_DEG$ENSEMBL_ID,]

# Subset raw gene expression data to 103 HK DEGs
ALL_P2_bmp_raw_data_filt_hkg <- ALL_P2_bmp_raw_data_filt[HK_DEG$ENSEMBL_ID,]

# Prepare filtered gene expression data of each of the 103 HK DEGs to be visualized as boxplot
wrangled_filt_data <- map(1:nrow(ALL_P2_bmp_filt_hkg), function(gene) {
  
  # Wrangle filtered gene expression data
  wrangle_exp_data(exp_data = ALL_P2_bmp_filt_hkg, 
                   HK_data = HK_DEG, 
                   gene_row_num = gene, 
                   B_ALL_barcodes = B_barcodes$barcodes, 
                   T_ALL_barcodes = T_barcodes$barcodes)
  
})

# Prepare raw gene expression data to be visualized as boxplot
wrangled_raw_data <- map(1:nrow(ALL_P2_bmp_raw_data_filt_hkg), function(gene) {
  
  # Wrangle filtered gene expression data
  wrangle_exp_data(exp_data = ALL_P2_bmp_raw_data_filt_hkg, 
                   HK_data = HK_DEG, 
                   gene_row_num = gene, 
                   B_ALL_barcodes = B_barcodes$barcodes, 
                   T_ALL_barcodes = T_barcodes$barcodes)
  
})


## Visualize data ---------------------------------------------------------------------------------------------------------

# For filtered gene expression data
# Create boxplot of each of the 103 HK DEGs stratified on subtype 
boxplot_filt_data <- map(1:nrow(ALL_P2_bmp_filt_hkg), function(gene) {
  
  boxplot_exp(wrangled_exp_data = wrangled_filt_data[[gene]], 
              subtype_col = wrangled_filt_data[[gene]]$subtype, 
              exp_col = wrangled_filt_data[[gene]]$exp_value, 
              ensembl_col = wrangled_filt_data[[gene]]$ENSEMBL_ID, 
              gene_name_col = wrangled_filt_data[[gene]]$Gene_name,
              ylab = "Filtered gene expression values", 
              plot_title = str_c("Filtered gene expression values of ", 
                                 wrangled_filt_data[[gene]]$ENSEMBL_ID %>% 
                                   unique, 
                                 " stratified on ALL subtypes \nExternal gene name: ",
                                 wrangled_filt_data[[gene]]$Gene_name))
  
})

# For raw gene expression data
# Create boxplot of each of the 103 HK DEGs stratified on subtype 
boxplot_raw_data <- map(1:nrow(ALL_P2_bmp_raw_data_filt_hkg), function(gene) {
  
  boxplot_exp(wrangled_exp_data = wrangled_raw_data[[gene]], 
              subtype_col = wrangled_raw_data[[gene]]$subtype, 
              exp_col = wrangled_raw_data[[gene]]$exp_value, 
              ensembl_col = wrangled_raw_data[[gene]]$ENSEMBL_ID, 
              gene_name_col = wrangled_raw_data[[gene]]$Gene_name,
              ylab = "Raw gene expression values", 
              plot_title = str_c("Raw gene expression values of ", 
                                 wrangled_raw_data[[gene]]$ENSEMBL_ID %>% 
                                   unique, 
                                 " stratified on ALL subtypes \nExternal gene name: ",
                                 wrangled_raw_data[[gene]]$Gene_name))
  
})

# Visualize logFC distribution of 103 housekeeping consensus DEGs
hkg_DEG_density <- dea_consensus %>%
  ungroup() %>%
  dplyr::filter(Gene_name %in% HK_DEG$Gene_name) %>%
  dplyr::select(logFC_limma_voom,
                logFC_edgeR,
                log2FoldChange_DESeq2, 
                logFC_mean) %>%
  dplyr::rename(limma_voom = logFC_limma_voom,
                edgeR = logFC_edgeR,
                DESeq2 = log2FoldChange_DESeq2) %>%
  pivot_longer(data = ., 
               cols = everything(), 
               names_to = "DEA_method",
               values_to = "log2FC_values") %>%
  ggplot(data = ., mapping = aes(x = log2FC_values, 
                                 color = DEA_method)) +
  geom_density(lwd = 2) +
  scale_color_manual(values = viridis_pal(option = "D")(4)) + 
  guides(color = guide_legend(override.aes = list(fill = viridis_pal(option = "D")(4)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 18,
                                   color = "black"),
        axis.title.x = element_text(size = 18,
                                    color = "black"),
        axis.text.y = element_text(size = 18,
                                   color = "black"),
        axis.ticks = element_line(size = 2),
        title = element_text(size = 18,
                             color = "black"),
        panel.border = element_rect(size = 2),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.3, "cm")) +
  labs(x = "log2FC",
       y = "Density",
       title = "Density plot of log2FC values from three DEA methods for housekeeping consensus DEGs")


## Save data ---------------------------------------------------------------------------------------------------------

# Save boxplots of filtered gene expression values of 103 HK DEGs in a combined pdf
ggsave(filename = "TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/TARGET_ALL_P2_bmp_filt_HK_DEG_boxplots.pdf", 
       plot = marrangeGrob(boxplot_filt_data, 
                           nrow = 1, 
                           ncol = 1), 
       device = "pdf",
       width = 10,
       height = 8)

# Save boxplots of raw gene expression values of 103 HK DEGs in a combined pdf
ggsave(filename = "TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/TARGET_ALL_P2_bmp_raw_HK_DEG_boxplots.pdf", 
       plot = marrangeGrob(boxplot_raw_data, 
                           nrow = 1, 
                           ncol = 1), 
       device = "pdf",
       width = 10,
       height = 8)

# Save density plot of log2FC values of housekeeping consensus DEGs
ggsave(filename = "TARGET_ALL_P2_bmp_hkg_eisenberg_consensus_DEG_density.pdf",  
       plot = hkg_DEG_density, 
       path = "TARGET_housekeeping/bone_marrow_primary_housekeeping_comp/", 
       width = 20, 
       height = 10)






