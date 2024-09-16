### This script performs validation using unsupervised clustering on
### TARGET blood samples

## Load libraries -------------------------------------------------------------------

library(TCGAbiolinks) #2.25.3
library(SummarizedExperiment) #1.28.0
library(limma) #3.54.0
library(gplots) #3.1.3.1
library(tidyverse) #1.3.2
library(viridis) #0.6.2


## Load data ------------------------------------------------------------------------

# Load ALL expression data as summarized experiment
ALL_P2_exp_data <- get(load("TARGET_data/ALL-P2_exp_data.rda"))


## Wrangle data ------------------------------------------------------------------------

# Subset data to contain only blood samples
blood_primary_exp_data <- subset(ALL_P2_exp_data, 
                                 select = (definition == "Primary Blood Derived Cancer - Peripheral Blood"))

# Create tibble of barcodes and subtypes
barcodes_subtypes <- blood_primary_exp_data@colData@listData %>% 
  as_tibble() %>% 
  dplyr::select(barcode, primary_diagnosis)

# Preprocess data
prep_data <- TCGAanalyze_Preprocessing(object = blood_primary_exp_data, 
                                       cor.cut = 0.6)

# Normalize data
norm_data <- TCGAanalyze_Normalization(tabDF = prep_data, 
                                       geneInfo = geneInfoHT, 
                                       method = "gcContent")

# Filter data
filt_data <- TCGAanalyze_Filtering(tabDF = norm_data, 
                                   method = "quantile",
                                   qnt.cut = 0.25)

# Voom transform filtered data
voom_data <- voom(filt_data)

# Retain only genes of interest which are the candidate genes 
# in the expression data
cand_genes <- c("ENSG00000118523", "ENSG00000128218", "ENSG00000164100", 
                "ENSG00000164330", "ENSG00000199683", "ENSG00000199831", 
                "ENSG00000200087",  "ENSG00000200312", "ENSG00000200959",  
                "ENSG00000201901", "ENSG00000202058", "ENSG00000223806", 
                "ENSG00000227706",  "ENSG00000271394")
cand_genes_exp <- voom_data$E[which(rownames(voom_data$E) %in% cand_genes), ]

# Create functions to use for clustering of expression data of candidate genes 
dist_func <- function(c) { dist(c, method = "euclidean") }
clustering_func <- function(c) { hclust(c, method = "complete") }

# Create vector of colors representing B-ALL samples as magenta and 
# T-ALL samples as cyan
col_vec <- colnames(cand_genes_exp) %>%
  as_tibble() %>% 
  dplyr::rename(barcode = value) %>% 
  left_join(barcodes_subtypes, by = "barcode") %>% 
  dplyr::mutate(color = str_replace(str_replace(primary_diagnosis, 
                                                "Precursor B-cell lymphoblastic leukemia", "#440154FF"), 
                                    "T lymphoblastic leukemia/lymphoma", "#FDE725FF"))

# Perform hiercharcical clustering of expression data of candidate genes and
# visualize in a heatmap
pdf("TARGET_blood_validation/TARGET_ALL_P2_blood_candidate_genes_clustering.pdf", height = 10, width = 10)
heatmap.2(as.matrix(cand_genes_exp), scale = "none", col = "viridis", 
          distfun = dist_func, hclustfun = clustering_func, labCol = FALSE, 
          ColSideColors = col_vec$color, trace = "none", density.info = "none", 
          margins = c(10, 9))
legend("topright", legend = c("B-cell ALL", "T-cell ALL"), 
       col = c("#440154FF", "#FDE725FF"), lty = 1, lwd = 10, border = FALSE, bty = "n", 
       y.intersp = 0.7, x.intersp = 1, cex = 0.9)
dev.off()


## Save data ------------------------------------------------------------------------

# Save SE of blood samples
save(blood_primary_exp_data, file = "TARGET_blood_validation/ALL_P2_blood_primary_exp_data.rda")

# Save preprocessed data
save(prep_data, file = "TARGET_blood_validation/ALL_P2_blood_primary_prep_data.rda")

# Save normalized data
save(norm_data, file = "TARGET_blood_validation/ALL_P2_blood_primary_norm_data.rda")

# Save filtered data
save(filt_data, file = "TARGET_blood_validation/ALL_P2_blood_primary_filt_data.rda")

# Save voom transform filtered data
save(voom_data, file = "TARGET_blood_validation/ALL_P2_blood_primary_voom_data.rda")

