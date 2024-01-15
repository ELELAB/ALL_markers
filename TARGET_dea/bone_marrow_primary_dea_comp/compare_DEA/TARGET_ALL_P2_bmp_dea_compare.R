### This script compares DEA results of TARGET-ALL-P2 gene expression data from three different DEA pipelines
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered


# Load libraries
library(tidyverse)
library(viridis)


## Load data ---------------------------------------------------------------------------------------------------------

# Vector of batch factor designs
batch_factor_designs <- list.dirs("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom", 
                                  full.names = FALSE, 
                                  recursive = FALSE)


# Load DEA of limma-voom pipeline of all batch designs
limma_voom_DEA <- map(seq(batch_factor_designs), function(x) {
  
  # Read DEA of significant DEGs
  readRDS(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_sig_", batch_factor_designs[x], ".rds"))
  
})
names(limma_voom_DEA) <- batch_factor_designs

# Load DEA of edgeR pipeline of all batch designs
edger_DEA <- map(seq(batch_factor_designs), function(x) {
  
  # Read DEA of significant DEGs
  readRDS(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_sig_", batch_factor_designs[x], ".rds"))
  
})
names(edger_DEA) <- batch_factor_designs

# Load DEA of DESeq2 pipeline of all batch designs
deseq2_DEA <- map(seq(batch_factor_designs), function(x) {
  
  # Read DEA of significant DEGs
  readRDS(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEA_sig_", batch_factor_designs[x], ".rds"))
  
})
names(deseq2_DEA) <- batch_factor_designs

# Load list of significantly upregulated DEGs for all DEA designs performed using limma-voom
DEGs_up_limma_voom <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly upregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_up_limma_voom) <- str_c(batch_factor_designs, "_up")

# Load list of significantly downregulated DEGs for all DEA designs performed using limma-voom
DEGs_down_limma_voom <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly downregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/limma_voom/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_down_limma_voom) <- str_c(batch_factor_designs, "_down")

# Load list of significantly upregulated DEGs for all DEA designs performed using edgeR
DEGs_up_edger <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly upregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_up_edger) <- str_c(batch_factor_designs, "_up")

# Load list of significantly downregulated DEGs for all DEA designs performed using edgeR
DEGs_down_edger <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly downregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/edgeR/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_down_edger) <- str_c(batch_factor_designs, "_down")

# Load list of significantly upregulated DEGs for all DEA designs performed using DESeq2
DEGs_up_deseq2 <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly upregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_up_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_up_deseq2) <- str_c(batch_factor_designs, "_up")

# Load list of significantly downregulated DEGs for all DEA designs performed using DESeq2
DEGs_down_deseq2 <- map(seq(batch_factor_designs), function(x) {
  
  # Read list of significantly downregulated DEGs
  read_csv(file = paste0("TARGET_dea/bone_marrow_primary_dea_comp/deseq2/", batch_factor_designs[x], "/ALL_P2_bmp_DEG_sig_down_", batch_factor_designs[x], ".csv"))
  
}) 
names(DEGs_down_deseq2) <- str_c(batch_factor_designs, "_down")



## Wrangle data ---------------------------------------------------------------------------------------------------------

## Create tables of number of DEGs

# Create table comparing number of identified DEGs between the three DEA pipelines
num_DEG_all_DEA <- tibble("Design" = batch_factor_designs,
                          "DEG_limma_voom" = map(limma_voom_DEA, nrow) %>%
                            unlist,
                          "DEG_edgeR" = map(edger_DEA, nrow) %>%
                            unlist,
                          "DEG_DESeq2" = map(deseq2_DEA, nrow) %>%
                            unlist)

# Create table comparing number of identified up- and downregulated DEGs between the three DEA pipelines
num_DEG_up_down_all_DEA <- tibble("Design" = batch_factor_designs,
                                  "DEG_up_limma_voom" = map(DEGs_up_limma_voom, function(x) x[[1]] %>% length) %>%
                                    unlist,
                                  "DEG_down_limma_voom" = map(DEGs_down_limma_voom, function(x) x[[1]] %>% length) %>%
                                    unlist,
                                  "DEG_up_edgeR" = map(DEGs_up_edger, function(x) x[[1]] %>% length) %>%
                                    unlist,
                                  "DEG_down_edgeR" = map(DEGs_down_edger, function(x) x[[1]] %>% length) %>%
                                    unlist,
                                  "DEG_up_DESeq2" = map(DEGs_up_deseq2, function(x) x[[1]] %>% length) %>%
                                    unlist,
                                  "DEG_down_DESeq2" = map(DEGs_down_deseq2, function(x) x[[1]] %>% length) %>%
                                    unlist)

## Create consensus DEA matrix

# Intersect DEGs between the three DEA pipelines for year of diagnosis batch factor design
DEG_common <- Reduce(intersect, list(limma_voom_DEA$years_batch$ENSEMBL_ID,
                                     edger_DEA$years_batch$ENSEMBL_ID,
                                     deseq2_DEA$years_batch$ENSEMBL_ID))

# Subset limma-voom DEA (year of diagnosis) to contain only common DEGs
DEA_common_limma_voom <- limma_voom_DEA$years_batch %>%
  filter(ENSEMBL_ID %in% DEG_common) %>% 
  dplyr::rename_with(.fn = function(x) str_c(x, "_limma_voom"), 
                     .cols = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"))

# Subset edgeR DEA (year of diagnosis) to contain only common DEGs
DEA_common_edger <- edger_DEA$years_batch %>%
  filter(ENSEMBL_ID %in% DEG_common) %>%
  dplyr::rename_with(.fn = function(x) str_c(x, "_edgeR"),
                     .cols = c("logFC", "unshrunk.logFC", "logCPM", "PValue", "FDR"))

# Subset DESeq2 DEA (year of diagnosis) to contain only common DEGs
DEA_common_deseq2 <- deseq2_DEA$years_batch %>%
  filter(ENSEMBL_ID %in% DEG_common) %>%
  dplyr::rename_with(.fn = function(x) str_c(x, "_DESeq2"),
                     .cols = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))

# Create one DEA consensus matrix and add mean logFC
DEA_common_matrix <- full_join(x = DEA_common_limma_voom, 
                               DEA_common_edger, 
                               by = "ENSEMBL_ID") %>%
  full_join(x = ., 
            y = DEA_common_deseq2, 
            by = "ENSEMBL_ID") %>%
  rowwise() %>%
  mutate("logFC_mean" = mean(c(logFC_limma_voom, logFC_edgeR, log2FoldChange_DESeq2)),
         "logFC_sd" = sd(c(logFC_limma_voom, logFC_edgeR, log2FoldChange_DESeq2)),
         "FDR_mean" = mean(c(adj.P.Val_limma_voom, FDR_edgeR, padj_DESeq2)),
         "FDR_sd" = sd(c(adj.P.Val_limma_voom, FDR_edgeR, padj_DESeq2))) %>%
  dplyr::select(-c(`Gene name.y`, `Gene name`)) %>%
  dplyr::rename(Gene_name = `Gene name.x`)

## ANOVA: Comparing logFC values for common DEGs -------------------------------

# Create full continuous variable
y <- c(DEA_common_matrix$logFC_limma_voom, DEA_common_matrix$logFC_edgeR, DEA_common_matrix$log2FoldChange_DESeq2)

# Create factor indications
treatment <- factor(c(rep("limma",length(DEA_common_matrix$logFC_limma_voom)),
                      rep("edgeR", length(DEA_common_matrix$logFC_edgeR)),
                      rep("DESeq2", length(DEA_common_matrix$log2FoldChange_DESeq2))))

## Checking assumptions

# Checking normal distribution of samples

# qq-plot for pooled residuals
residuals  <- lm(y~treatment)$residuals
pdf("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_qq_residuals.pdf")
qq_residuals <- qqnorm(residuals)
qq_residuals <- qqline(residuals)
dev.off()

# qq-plot for each pipeline
qqnorm(DEA_common_matrix$logFC_limma_voom)
qqline(DEA_common_matrix$logFC_limma_voom)

qqnorm(DEA_common_matrix$logFC_edgeR)
qqline(DEA_common_matrix$logFC_edgeR)

qqnorm(DEA_common_matrix$log2FoldChange_DESeq2)
qqline(DEA_common_matrix$log2FoldChange_DESeq2)


# Checking equality of variance of samples
# See logFC_boxplot: They seem to have the same variances

## ANOVA
anova_logFC <- anova(lm(y ~ treatment)) %>%
  as_tibble() %>%
  add_column(" " = c("treatment","residuals"), .before = "Df")


## Visualize data ---------------------------------------------------------------------------------------------------------

# Visualize logFC values of common DEGs 
logFC_boxplot <- DEA_common_matrix %>%
  dplyr::select(c(logFC_limma_voom, 
                  logFC_edgeR, 
                  log2FoldChange_DESeq2,
                  logFC_mean)) %>%
  pivot_longer(data = ., 
               cols = everything(), 
               names_to = "method", 
               values_to = "logFC") %>%
  mutate(method = str_remove(method, 
                             "(.*?)_")) %>%
  ggplot(data = ., 
         mapping = aes(x = method, 
                       y = logFC)) +
  geom_boxplot(color = "black",
               fill = viridis_pal(option = "viridis")(4),
               lwd = 1.5,
               outlier.size = 3) +
  labs(x = "", 
       y = "log2FoldChange",
       title = "log2FoldChange values of common DEGs between three DEA methods") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14, 
                                   color = "Black"),
        axis.text.y = element_text(size = 14,
                                   color = "Black"),
        axis.title.y = element_text(size = 14,
                                    color = "black"),
        title = element_text(size = 14))


## Save data ---------------------------------------------------------------------------------------------------------

# Save table containing number of identified DEGs
write_csv(num_DEG_all_DEA, file = "TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEG_numbers.csv")

# Save table containing number of identified up- and downregulated DEGs
write_csv(num_DEG_up_down_all_DEA, file = "TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_up_down_DEG_numbers.csv")

# Save consensus/common DEA matrix
save(DEA_common_matrix, file = "TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda")

# Save boxplot of logFC values of common DEGs between all three DEA methods
ggsave(filename = "ALL_P2_bmp_logFC_DEG_common_year_boxplot.pdf", 
       plot = logFC_boxplot, 
       path = "TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/", 
       width = 20, 
       height = 8)

# Save ANOVA results
write_csv(anova_logFC, file = "TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_anova_logFC.csv")



