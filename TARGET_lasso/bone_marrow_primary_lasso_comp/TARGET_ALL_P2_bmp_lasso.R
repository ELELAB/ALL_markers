### This script performs LASSO and elastic net logistic regression on TARGET ALL P2 gene expression data
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
library(glmnet)
library(caret)
library(tidyverse)
library(biomaRt)
library(viridis)


# Source functions script
source("TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_functions.R")


### Load data ---------------------------------------------------------------------------------------------------------

# Load filtered gene count data
ALL_P2_filt <- get(load("TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_filt_data.rda"))

# Load batch corrected gene count data (corrected for year of diagnosis batch factor)
ALL_P2_corr <- get(load("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda"))

# Load table containing batch information about samples
ALL_P2_info <- read_csv(file = "TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")

# Load consensus DEA table
DEA_consensus <- get(load("TARGET_dea/bone_marrow_primary_dea_comp/compare_DEA/ALL_P2_bmp_DEA_consensus.rda"))


### Wrangle data ---------------------------------------------------------------------------------------------------------

# Binarize vector containing group of samples (subtype)
group_data <- case_when(ALL_P2_info$subtype_abbr == "B-cell ALL" ~ 0,
                        ALL_P2_info$subtype_abbr == "T-cell ALL" ~ 1)


## Run LASSO 

# Get randomly 10 seeds
#seeds <- sample(1:1000, 10)
seeds <- read_csv(file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_seeds.csv") %>%
  pull

# Empty list to store results of 10 LASSO runs (one run for each seed)
LASSO_res_filt <- list()
LASSO_res_filt_05 <- list()
LASSO_res_corr <- list()
LASSO_res_corr_05 <- list()

# Run LASSO regression 10 times (one run for each seed)
for (idx in 1:length(seeds)) {
  
  # Perform LASSO regression using filtered count data 
  LR_filt <- LASSOFeature(my.seed = seeds[[idx]], 
                          my.data = ALL_P2_filt, 
                          my.group = group_data, 
                          my.LAorEN = 1, 
                          n.cv.folds = 5)
  
  # Perform elastic net regression using filtered count data 
  LR_filt_05 <- LASSOFeature(my.seed = seeds[[idx]], 
                             my.data = ALL_P2_filt, 
                             my.group = group_data, 
                             my.LAorEN = 0.5, 
                             n.cv.folds = 5)
  
  # Perform LASSO regression using batch corrected count data
  LR_corr <- LASSOFeature(my.seed = seeds[[idx]], 
                          my.data = ALL_P2_corr, 
                          my.group = group_data, 
                          my.LAorEN = 1, 
                          n.cv.folds = 5)
  
  # Perform elastic net regression using batch corrected count data (alpha = 0.5)
  LR_corr_05 <- LASSOFeature(my.seed = seeds[[idx]], 
                             my.data = ALL_P2_corr, 
                             my.group = group_data, 
                             my.LAorEN = 0.5, 
                             n.cv.folds = 5)
  
  # Save results of LASSO / elastic net regression on filtered data to list
  LASSO_res_filt[[idx]] <- LR_filt
  LASSO_res_filt_05[[idx]] <- LR_filt_05
  
  # Save results of LASSO / elastic net regression on batch corrected data to list
  LASSO_res_corr[[idx]] <- LR_corr
  LASSO_res_corr_05[[idx]] <- LR_corr_05
  
}
names(LASSO_res_filt) <- names(LASSO_res_filt_05) <- names(LASSO_res_corr) <- names(LASSO_res_corr_05) <- str_c("seed_run_", 
                                                                                                                seq(1:10))


## Extract results from 10 LASSO runs

# For each seed run, get all those genes that do not have a coefficient equal to 0 
# These genes would create the minimal subset of genes of interest
gene_subsets_filt <- get_coef(LASSO_res = LASSO_res_filt)

gene_subsets_filt_05 <- get_coef(LASSO_res = LASSO_res_filt_05)

gene_subsets_corr <- get_coef(LASSO_res = LASSO_res_corr)

gene_subsets_corr_05 <- get_coef(LASSO_res = LASSO_res_corr_05)


# Intersect the 10 minimal subset of genes,
#i.e. one minimal subset of genes from each seed run
gene_select_filt <- intersect_genes(LASSO_genes = gene_subsets_filt)

gene_select_filt_05 <- intersect_genes(LASSO_genes = gene_subsets_filt_05)

gene_select_corr <- intersect_genes(LASSO_genes = gene_subsets_corr)

gene_select_corr_05 <- intersect_genes(LASSO_genes = gene_subsets_corr_05)


# Extract the prediction errors (misclassification error) from the output of the 10 runs 
cv_error_runs_filt <- get_cv_error(LASSO_res = LASSO_res_filt)

cv_error_runs_filt_05 <- get_cv_error(LASSO_res = LASSO_res_filt_05)

cv_error_runs_corr <- get_cv_error(LASSO_res = LASSO_res_corr)

cv_error_runs_corr_05 <- get_cv_error(LASSO_res = LASSO_res_corr_05)


# Compute average coefficient value for the consensus genes across 10 seed runs (batch corrected data)
gene_select_coef_corr_05 <- map(gene_subsets_corr_05, function(x) {
  x %>% 
    as_tibble() %>% 
    filter(ENSEMBL_ID %in% gene_select_corr_05$selected_genes)
}) %>% 
  bind_cols() %>% 
  dplyr::select(1, 
                contains("s1")) %>% 
  mutate(mean = rowMeans(across(where(is.numeric)))) 
names(gene_select_coef_corr_05) <- c("ENSEMBL_ID", 
                                     str_c("seed_run_", 
                                           seq(1:10),
                                           "_s1"),
                                     "mean_s1")

# Add standard deviation of coefficients across 10 seed runs to tibble
gene_select_coef_corr_05 <- gene_select_coef_corr_05 %>% 
  rowwise() %>% 
  mutate(sd_s1 = sd(c(seed_run_1_s1, seed_run_2_s1, seed_run_3_s1, seed_run_4_s1, seed_run_5_s1,
                      seed_run_6_s1, seed_run_7_s1, seed_run_8_s1, seed_run_9_s1, seed_run_10_s1)))
  
# Add external gene names to elastic net consensus genes across 10 seed runs (batch corrected data)
gene_select_corr_name_05 <- ensembl_hugo_conversion(ensembl_id = gene_select_coef_corr_05$ENSEMBL_ID) %>% 
  as_tibble() %>% 
  dplyr::rename(ENSEMBL_ID = `Gene stable ID`,
                gene_name = `Gene name`) %>% 
  full_join(x = .,
            y = gene_select_coef_corr_05,
            by = "ENSEMBL_ID") %>% 
  arrange(desc(abs(mean_s1))) 

# Compare elastic net consensus genes with DEA consensus
gene_select_en_DEA <- inner_join(x = gene_select_corr_name_05, 
                                 y = DEA_consensus, 
                                 by = "ENSEMBL_ID") %>% 
  dplyr::select(-Gene_name)

## Correlation test between mean coefficient values and mean log2FC values

# Check assumption of linear relationship between the two variables
#plot(gene_select_en_DEA$mean_s1, gene_select_en_DEA$logFC_mean)

# Check assumption of normal distribution of each variable
#hist(gene_select_en_DEA$mean_s1)
#hist(gene_select_en_DEA$logFC_mean)

pearson_cc <- cor.test(gene_select_en_DEA$mean_s1, 
                       gene_select_en_DEA$logFC_mean, 
                       method = "pearson")

pearson_cc_tbl <- tibble("corr_coef" = pearson_cc$estimate, 
                         "p-value" = pearson_cc$p.value,
                         "95%_CI" = str_c(pearson_cc$conf.int[1], ",",
                                      pearson_cc$conf.int[2]))


# Filter DEA consensus tibble into upregulated DEGs and downregulated DEGs
DEA_consensus_up <- DEA_consensus %>% 
  filter(logFC_mean > 0)
DEA_consensus_down <- DEA_consensus %>% 
  filter(logFC_mean < 0)


### Visualize data ---------------------------------------------------------------------------------------------------------

## Visualize prediction errors

# Create bar plot of the prediction errors (misclassification error) from each of the 10 runs
cv_plot_filt <- cv_error_plot(seed_runs_col = cv_error_runs_filt$seed_run, 
                              pred_error_col = cv_error_runs_filt$pred_error)

cv_plot_filt_05 <- cv_error_plot(seed_runs_col = cv_error_runs_filt_05$seed_run, 
                                 pred_error_col = cv_error_runs_filt_05$pred_error)

cv_plot_corr <- cv_error_plot(seed_runs_col = cv_error_runs_corr$seed_run, 
                              pred_error_col = cv_error_runs_corr$pred_error)

cv_plot_corr_05 <- cv_error_plot(seed_runs_col = cv_error_runs_corr_05$seed_run, 
                                 pred_error_col = cv_error_runs_corr_05$pred_error)

## Visualize mean cross-validation errors

# Create bar plot of mean cross-validation errors (mean across 5 folds and across all tested lambdas)

mean_cv_error_plot <- map(LASSO_res_corr_05, 
                     function(x) {
                       x[["lasso_fit"]][["cvm"]] }) %>%
  map(., 
      function(x) { 
        mean(x) }) %>% 
  unlist() %>%
  as_tibble() %>%
  dplyr::rename(mean_cv_error_lambdas = value) %>%
  mutate(seed_run = str_c("Run ", seq(1:10)),
         seed_run = factor(seed_run, levels = unique(seed_run))) %>%
  ggplot(data = ., mapping = aes(x = seed_run, y = mean_cv_error_lambdas)) +
  geom_bar(stat = "identity",
           fill = viridis_pal(option = "viridis")(1)) + 
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 0.05, 0.01)) +
  labs(x = "",
       y = "Mean CV error across CV folds and lambda values",
       title = "Mean CV error for each of the 10 seed runs") +
  theme(axis.text.x = element_text(size = 16, 
                                   color = "black"),
        axis.text.y = element_text(size = 16, 
                                   color = "black"),
        axis.title.y = element_text(size = 16, 
                                    color = "black"),
        title = element_text(size = 16,
                             color = "black"))



### Save data ---------------------------------------------------------------------------------------------------------

# Save vector of seeds used
write_csv(as_tibble(seeds), file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_seeds.csv")

# Save results of 10 LASSO runs using filtered data
saveRDS(LASSO_res_filt, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_runs_filt.rds")

# Save results of 10 elastic net runs using filtered data
saveRDS(LASSO_res_filt_05, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_elastic_runs_filt.rds")

# Save results of 10 LASSO runs using batch corrected data
saveRDS(LASSO_res_corr, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_runs_corr.rds")

# Save results of 10 elastic net runs using batch corrected data
saveRDS(LASSO_res_corr_05, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_elastic_runs_corr.rds")

# Save minimal subset of genes found in the 10 LASSO runs including their coefficients using filtered data
saveRDS(gene_subsets_filt, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_runs_genes_filt.rds")

# Save minimal subset of genes found in the 10 elastic net runs including their coefficients using filtered data
saveRDS(gene_subsets_filt_05, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_elastic_runs_genes_filt.rds")

# Save minimal subset of genes found in the 10 LASSO runs including their coefficients using batch corrected data
saveRDS(gene_subsets_corr, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_runs_genes_corr.rds")

# Save minimal subset of genes found in the 10 elastic runs including their coefficients using batch corrected data
saveRDS(gene_subsets_corr_05, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_elastic_runs_genes_corr.rds")

# Save consensus minimal subset of genes using filtered data - lasso
write_csv(gene_select_filt, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_consensus_genes_filt.csv")

# Save consensus minimal subset of genes using filtered data - elastic net
write_csv(gene_select_filt_05, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_elastic_consensus_genes_filt.csv")

# Save consensus minimal subset of genes using batch corrected data - lasso
write_csv(gene_select_corr, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_consensus_genes_corr.csv")

# Save consensus minimal subset of genes using batch corrected data - elastic net
write_csv(gene_select_corr_05, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_elastic_consensus_genes_corr.csv")

# Save table of prediction errors for each of the 10 LASSO runs using filtered data
write_csv(cv_error_runs_filt, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_pred_error_filt.csv")

# Save table of prediction errors for each of the 10 elastic net runs using filtered data
write_csv(cv_error_runs_filt_05, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_elastic_pred_error_filt.csv")

# Save table of prediction errors for each of the 10 LASSO runs using batch corrected data
write_csv(cv_error_runs_corr, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_lasso_pred_error_corr.csv")

# Save table of prediction errors for each of the 10 elastic net runs using batch corrected data
write_csv(cv_error_runs_corr_05, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_elastic_pred_error_corr.csv")

# Save table containing coefficients of consensus elastic net genes and comparison with consensus DEA
write_csv(gene_select_en_DEA, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_elastic_consensus_DEA_corr.csv")

# Save results from Pearson correlation test
write_csv(pearson_cc_tbl, file = "TARGET_lasso/bone_marrow_primary_lasso_comp/TARGET_ALL_P2_bmp_pearson.csv")

# Save plots of prediction errors for each of the 10 LASSO runs using filtered data
ggsave(plot = cv_plot_filt, 
       width = 18, 
       height = 8, 
       filename = "TARGET_ALL_P2_bmp_lasso_pred_error_plot_filt.pdf", 
       path = "TARGET_lasso/bone_marrow_primary_lasso_comp/")

# Save plots of prediction errors for each of the 10 elastic net runs using filtered data
ggsave(plot = cv_plot_filt_05, 
       width = 18, 
       height = 8, 
       filename = "TARGET_ALL_P2_bmp_elastic_pred_error_plot_filt.pdf", 
       path = "TARGET_lasso/bone_marrow_primary_lasso_comp/")

# Save plots of prediction errors for each of the 10 LASSO runs using batch corrected data
ggsave(plot = cv_plot_corr, 
       width = 18, 
       height = 8, 
       filename = "TARGET_ALL_P2_bmp_lasso_pred_error_plot_corr.pdf", 
       path = "TARGET_lasso/bone_marrow_primary_lasso_comp/")

# Save plots of prediction errors for each of the 10 elastic net runs using batch corrected data
ggsave(plot = cv_plot_corr_05, 
       width = 18, 
       height = 8, 
       filename = "TARGET_ALL_P2_bmp_elastic_pred_error_plot_corr.pdf", 
       path = "TARGET_lasso/bone_marrow_primary_lasso_comp/")

# Save plots of mean cros-validation errors for each of the 10 elastic net runs using batch corrected data
ggsave(plot = mean_cv_error_plot,
       width = 10,
       height = 6,
       filename = "TARGET_ALL_P2_bmp_elastic_mean_cv_error_plot_corr.pdf",
       path = "TARGET_lasso/bone_marrow_primary_lasso_comp/")






