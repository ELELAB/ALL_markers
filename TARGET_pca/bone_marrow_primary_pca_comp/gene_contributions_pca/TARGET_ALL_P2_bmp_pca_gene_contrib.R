### This script performs PCA on TARGET ALL P2 gene expression data to investigate gene contributions
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
library(factoextra)
library(FactoMineR)
library(tidyverse)
library(biomaRt)

# Source functions
source("TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/TARGET_ALL_P2_bmp_pca_gene_contrib_functions.R")



## Load data ---------------------------------------------------------------------------------------------------------

# Load batch corrected gene count data (corrected for year of diagnosis batch factor)
ALL_P2_corr <- get(load("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda"))

# Load table containing batch information about samples
ALL_P2_info <- read_csv(file = "TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")


## Wrangle data ---------------------------------------------------------------------------------------------------------

# Binarize vector containing group of samples (subtype)
group_data <- case_when(ALL_P2_info$subtype_abbr == "B-cell ALL" ~ 0,
                        ALL_P2_info$subtype_abbr == "T-cell ALL" ~ 1)

# Perform PCA on batch corrected data
res.pca.corr <- PCA(t(ALL_P2_corr), 
                    graph = FALSE, 
                    ncp = 10, 
                    scale = FALSE)

# Visualize eigenvalues/variances
fviz_screeplot(res.pca.corr, 
               addlabels = TRUE, 
               ylim = c(0, 50), 
               ggtheme = theme_classic(), 
               ncp = 20) 
ggsave(filename = paste0("TARGET_ALL_P2_bmp_pca_screeplot_corr.pdf"), 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/")

# Contributions of variables to PC1
fviz_contrib(res.pca.corr, 
             choice = "var", 
             axes = 1, 
             top = 40, 
             ggtheme = theme_classic()) 
ggsave(filename = paste0("TARGET_ALL_P2_bmp_pca_contrib_PC1_corr.pdf"), 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/")

# Contributions of variables to PC2
fviz_contrib(res.pca.corr, 
             choice = "var", 
             axes = 2, 
             top = 40, 
             ggtheme = theme_classic()) 
ggsave(filename = paste0("TARGET_ALL_P2_bmp_pca_contrib_PC2_corr.pdf"), 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/")

# Create PCA plot, visualize samples by groups
PCA.labels <- "none"
cols <- c("orange", "darkgreen")
fviz_pca_ind(res.pca.corr,
             label = PCA.labels, # hide individual labels
             habillage = as.factor(group_data), # color by groups
             palette = cols,
             addEllipses = TRUE, # concentration ellipses
             repel = TRUE,
             ggtheme = theme_classic(),
             labelsize = 2)
ggsave(filename = paste0("TARGET_ALL_P2_bmp_pca_corr_plot.pdf"), 
       path = "TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/")

# Extract ENSEMBL IDs of top 40 genes contributing to PC1 and PC2 
gene_contrib_axes <- map(c(1:2), function(x) {
  
  facto_summarize(X = res.pca.corr, 
                  element = "var", 
                  result = "contrib", 
                  axes = x) %>%
    as_tibble() %>%
    arrange(desc(contrib)) %>%
    slice_max(n = 40, 
              order_by = contrib) %>% 
    dplyr::rename(ENSEMBL_ID = name) %>% 
    mutate(PCA_dim = rep(str_c("PC_",x), 40))
  
}) %>% 
  bind_rows()

# Convert ENSEMBL IDs to gene names
gene_contrib_axes_names <- ensembl_hugo_conversion(ensembl_id = gene_contrib_axes$ENSEMBL_ID) %>% 
  as_tibble() %>%
  dplyr::rename("ENSEMBL_ID" = `Gene stable ID`) %>%
  full_join(x = gene_contrib_axes,
            y = ., 
            by = "ENSEMBL_ID")



## Save data ---------------------------------------------------------------------------------------------------------

# Save results of PCA
save(res.pca.corr, 
     file = "TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/TARGET_ALL_P2_bmp_pca_corr.rda")

# Save information of top 20 genes contributing to PC1 and PC2
write_csv(gene_contrib_axes_names, file = "TARGET_pca/bone_marrow_primary_pca_comp/gene_contributions_pca/TARGET_ALL_P2_bmp_pca_corr_genes_info.csv")





