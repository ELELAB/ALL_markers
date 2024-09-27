### This script performs unsupervised hierarchical clustering of GTEx normal
### blood and bone marrow expression data

## Load libraries -------------------------------------------------------------------

library(TCGAbiolinks) #2.25.3
library(recount) #1.24.1
library(limma) #3.54.0
library(gplots) #3.1.3.1
library(viridis) #0.6.2

## Get data ------------------------------------------------------------------------

# Download blood and bone marrow GTEx data
blood_gtex <- TCGAquery_recount2(project = "gtex", tissue = "blood")
bm_gtex <- TCGAquery_recount2(project = "gtex", tissue = "bone_marrow") 


## Wrangle data ------------------------------------------------------------------------

# Scale recount data and get the count matrix
blood_gtex_count_scale <- assays(scale_counts(blood_gtex$gtex_blood, round = TRUE))$counts
bm_gtex_count_scale <- assays(scale_counts(bm_gtex$gtex_bone_marrow, round = TRUE))$counts

# Removing version number in rownames of count matrix
rownames(blood_gtex_count_scale) <- gsub("\\..*", "", rownames(blood_gtex_count_scale))
rownames(bm_gtex_count_scale) <- gsub("\\..*", "", rownames(bm_gtex_count_scale))

# Normalize data
blood_norm_data <- TCGAanalyze_Normalization(tabDF = blood_gtex_count_scale, 
                                             geneInfo = geneInfoHT, 
                                             method = "gcContent")
bm_norm_data <- TCGAanalyze_Normalization(tabDF = bm_gtex_count_scale, 
                                          geneInfo = geneInfoHT, 
                                          method = "gcContent")

# Filter data
blood_filt_data <- TCGAanalyze_Filtering(tabDF = blood_norm_data, 
                                         method = "quantile",
                                         qnt.cut = 0.25)
bm_filt_data <- TCGAanalyze_Filtering(tabDF = bm_norm_data, 
                                      method = "quantile",
                                      qnt.cut = 0.25)

# Voom transform filtered data
blood_voom_data <- voom(blood_filt_data)
bm_voom_data <- voom(bm_filt_data)

# Retain only genes of interest which are the candidate genes 
# in the expression data
cand_genes <- c("ENSG00000118523", "ENSG00000128218", "ENSG00000164100", 
                "ENSG00000164330", "ENSG00000199683", "ENSG00000199831", 
                "ENSG00000200087",  "ENSG00000200312", "ENSG00000200959",  
                "ENSG00000201901", "ENSG00000202058", "ENSG00000223806", 
                "ENSG00000227706",  "ENSG00000271394")
blood_cand_genes_exp <- blood_voom_data$E[which(rownames(blood_voom_data$E) %in% cand_genes), ]
bm_cand_genes_exp <- bm_voom_data$E[which(rownames(bm_voom_data$E) %in% cand_genes), ]

# Create functions to use for clustering of expression data of candidate genes 
dist_func <- function(c) { dist(c, method = "euclidean") }
clustering_func <- function(c) { hclust(c, method = "complete") }

# Perform hiercharcical clustering of expression data of candidate genes and
# visualize in a heatmap

# Blood
pdf("GTEx_validation/GTEx_blood_candidate_genes_clustering.pdf", height = 10, width = 10)
heatmap.2(as.matrix(blood_cand_genes_exp), scale = "none", col = "viridis", 
          distfun = dist_func, hclustfun = clustering_func, labCol = FALSE, 
          trace = "none", density.info = "none", margins = c(10, 9))
dev.off()

# Bone marrow
pdf("GTEx_validation/GTEx_bone_marrow_candidate_genes_clustering.pdf", height = 10, width = 10)
heatmap.2(as.matrix(bm_cand_genes_exp), scale = "none", col = "viridis", 
          distfun = dist_func, hclustfun = clustering_func, labCol = FALSE, 
          trace = "none", density.info = "none", margins = c(10, 9))
dev.off()


## Save data ------------------------------------------------------------------------

# Save GTEx data
save(blood_gtex, file = "GTEx_validation/GTEx_blood_SE.rda")
save(bm_gtex, file = "GTEx_validation/GTEx_bone_marrow_SE.rda")

# Save scaled data
save(blood_gtex_count_scale, file = "GTEx_validation/GTEx_blood_scaled_counts.rda")
save(bm_gtex_count_scale, file = "GTEx_validation/GTEx_bone_marrow_scaled_counts.rda")

# Save normalized data
save(blood_norm_data, file = "GTEx_validation/GTEx_blood_norm_data.rda")
save(bm_norm_data, file = "GTEx_validation/GTEx_bone_marrow_norm_data.rda")

# Save filtered data
save(blood_filt_data, file = "GTEx_validation/GTEx_blood_filt_data.rda")
save(bm_filt_data, file = "GTEx_validation/GTEx_bone_marrow_filt_data.rda")

# Save voom transform filtered data
save(blood_voom_data, file = "GTEx_validation/GTEx_blood_voom_data.rda")
save(bm_voom_data, file = "GTEx_validation/GTEx_bone_marrow_voom_data.rda")
