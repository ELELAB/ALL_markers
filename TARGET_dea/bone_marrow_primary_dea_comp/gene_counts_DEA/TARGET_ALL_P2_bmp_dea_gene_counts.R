### This script is used to make visualizations of gene counts for the 
### three DEA tools, DESeq2, edgeR and limma-voom. The results from each
### tool is used as data.


# Load libraries
library(tidyverse)


## Load data ---------------------------------------------------------------

# DESeq2 results
res_deseq2 <- readRDS("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/years_batch/ALL_P2_bmp_DEA_full_years_batch.rds")

# edgeR results
res_edger <- readRDS("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/years_batch/ALL_P2_bmp_DEA_full_years_batch.rds")

# limma voom results
res_limma <- readRDS("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/years_batch/ALL_P2_bmp_DEA_full_years_batch.rds")

## Data wrangling

# Add column denoting whether the gene is DEG column to DESeq2 results
res_deseq2 <- res_deseq2 %>%
  dplyr::mutate(DEG = ifelse(padj < 0.05, 1, 0)) %>%
  dplyr::mutate(DEG = as.factor(DEG)) %>%
  dplyr::select(ENSEMBL_ID,`Gene name`, baseMean, log2FoldChange, DEG) %>%
  dplyr::rename(Gene_name = `Gene name`)

# Add column denoting whether the gene is DEG column to edgeR results
res_edger <- res_edger %>% 
  dplyr::mutate(DEG = ifelse(FDR < 0.05, 1, 0)) %>%
  dplyr::mutate(DEG = as.factor(DEG)) %>%
  dplyr::select(ENSEMBL_ID, `Gene name`, logFC, logCPM, DEG) %>%
  dplyr::rename(Gene_name = `Gene name`)

# Add column denoting whether the gene is DEG column to limma results
res_limma <- res_limma %>% 
  dplyr::mutate(DEG = ifelse(adj.P.Val < 0.05, 1, 0)) %>%
  dplyr::mutate(DEG = as.factor(DEG)) %>%
  dplyr::rename(logCPM = AveExpr)  %>% # AveExpr corresponds to logCPM
  dplyr::select(ENSEMBL_ID, `Gene name`, logFC, logCPM, DEG) %>%
  dplyr::rename(Gene_name = `Gene name`)

# res_combined
res_combined <- res_edger %>%
  inner_join(res_limma, by = c("ENSEMBL_ID", "Gene_name"), suffix = c("_edger", "_limma")) %>%
  inner_join(res_deseq2, by = c("ENSEMBL_ID", "Gene_name"), suffix = c("", "_deseq2")) %>%
  dplyr::rename(baseMean_deseq2 = baseMean, log2FoldChange_deseq2 = log2FoldChange, DEG_deseq2 = DEG) %>%
  dplyr::mutate(common_DEG = ifelse(DEG_deseq2 == 1 & DEG_edger == 1 & DEG_limma == 1, 1, 0)) %>%
  dplyr::mutate(common_DEG = as.factor(common_DEG))


## Visualization ---------------------------------------------------------------

# DESeq2 plots
plot_deseq2_counts_log <- ggplot(res_deseq2, aes(log(baseMean), log2FoldChange, color = DEG)) +
  geom_point() +
  geom_hline(yintercept = c(-1,1), color = "darkgray") + 
  ggtitle("DESeq2") + 
  xlab("log(mean of normalized counts)")
  
plot_deseq2_counts_no_log <- ggplot(res_deseq2, aes(baseMean, log2FoldChange, color = DEG)) +
  geom_point() +
  scale_x_continuous(trans = 'log10') +
  geom_hline(yintercept = c(-1,1), color = "darkgray") + 
  ggtitle("DESeq2") + 
  xlab("mean of normalized counts")
  

# edgeR plots
plot_edger_counts_log <- ggplot(res_edger, aes(logCPM, logFC, color = DEG)) +
  geom_point() +
  geom_hline(yintercept = c(-1,1), color = "darkgray") + 
  ggtitle("edgeR")

plot_edger_counts_no_log <- ggplot(res_edger, aes(10**logCPM, logFC, color = DEG)) +
  geom_point() + 
  scale_x_continuous(trans = 'log10') +
  geom_hline(yintercept = c(-1,1), color = "darkgray") + 
  ggtitle("edgeR")


# limma voom plots
plot_limma_counts_log <- ggplot(res_limma, aes(logCPM, logFC, color = DEG)) +
  geom_point() +
  geom_hline(yintercept = c(-1,1), color = "darkgray") + 
  ggtitle("limma")

plot_limma_counts_no_log <- ggplot(res_limma, aes(10**logCPM, logFC, color = DEG)) +
  geom_point() + 
  scale_x_continuous(trans = 'log10') +
  geom_hline(yintercept = c(-1,1), color = "darkgray") + 
  ggtitle("limma")


# Combined plots with common DEGs (common for all three tools)
comb_plot_deseq2_edger <- ggplot(res_combined, aes(baseMean_deseq2, logCPM_edger, color = common_DEG)) +
  geom_point() +
  scale_x_continuous(trans = 'log10')

comb_plot_limma_edger <- ggplot(res_combined, aes(logCPM_limma, logCPM_edger, color = common_DEG)) +
  geom_point()

comb_plot_limma_deseq2 <- ggplot(res_combined, aes(logCPM_limma, baseMean_deseq2, color = common_DEG)) +
  geom_point() + 
  scale_y_continuous(trans = 'log10')

# Combined plots with common DEGs (common for all three tools)
comb_plot_DEG_deseq2_edger <- ggplot(subset(res_combined, common_DEG == 1), aes(baseMean_deseq2, logCPM_edger, color = common_DEG)) +
  geom_point() +
  scale_x_continuous(trans = 'log10')

comb_plot_DEG_limma_edger <- ggplot(subset(res_combined, common_DEG == 1), aes(logCPM_limma, logCPM_edger, color = common_DEG)) +
  geom_point()

comb_plot_DEG_limma_deseq2 <- ggplot(subset(res_combined, common_DEG == 1), aes(logCPM_limma, baseMean_deseq2, color = common_DEG)) +
  geom_point() + 
  scale_y_continuous(trans = 'log10')

## Save results ---------------------------------------------------------------

# Save plots of deseq2 counts (logarithmic and non-logarithmic)
ggsave("ALL_P2_deseq2_counts_log.pdf", 
       plot = plot_deseq2_counts_log, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

ggsave("ALL_P2_deseq2_counts_no_log.pdf", 
       plot = plot_deseq2_counts_no_log, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

# Save plots of edgeR counts (logarithmic and non-logarithmic)
ggsave("ALL_P2_edgeR_counts_log.pdf", 
       plot = plot_edger_counts_log, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

ggsave("ALL_P2_edgeR_counts_no_log.pdf", 
       plot = plot_edger_counts_no_log, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

# Save plots of limma-voom counts (logarithmic and non-logarithmic)
ggsave("ALL_P2_limma_counts_log.pdf", 
       plot = plot_limma_counts_log, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

ggsave("ALL_P2_limma_counts_no_log.pdf", 
       plot = plot_limma_counts_no_log, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

# Save combined plots (not filtered for DEGs)
ggsave("ALL_P2_comb_deseq2_edgeR.pdf", 
       plot = comb_plot_deseq2_edger, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

ggsave("ALL_P2_comb_edgeR_limma.pdf", 
       plot = comb_plot_limma_edger, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

ggsave("ALL_P2_comb_limma_deseq2.pdf", 
       plot = comb_plot_limma_deseq2, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")


# Save combined plots(filtered for DEGs)
ggsave("ALL_P2_comb_DEG_deseq2_edgeR.pdf", 
       plot = comb_plot_DEG_deseq2_edger, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

ggsave("ALL_P2_comb_DEG_edgeR_limma.pdf", 
       plot = comb_plot_DEG_limma_edger, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

ggsave("ALL_P2_comb_DEG_limma_deseq.pdf", 
       plot = comb_plot_DEG_limma_deseq2, device = "pdf", path = "TARGET_dea/bone_marrow_primary_dea_comp/gene_counts_DEA/visualizations/")

