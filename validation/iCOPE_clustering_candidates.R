### This script performs clustering of gene expression data 



## Load libraries -------------------------------------------------------------------
library(gplots) #3.1.3
library(tidyverse) #2.0.0
library(viridis) #0.6.4


## Load data ------------------------------------------------------------------------

# Load voom transformed filtered data
voom_filt_data <- get(load("transform/iCOPE_ALL_filt_data_voom.rda"))

# Load meta data
sample_meta_data_df <- get(load("get_data/iCOPE_ALL_meta_data.rda"))


## Wrangle and visualize data ---------------------------------------------------------------------

# Retain only genes of interest which are the candidate genes 
# in the expression data
cand_genes <- c("ENSG00000118523", "ENSG00000128218", "ENSG00000164100", 
                "ENSG00000164330", "ENSG00000199683", "ENSG00000199831", 
                "ENSG00000200087",  "ENSG00000200312", "ENSG00000200959",  
                "ENSG00000201901", "ENSG00000202058", "ENSG00000223806", 
                "ENSG00000227706",  "ENSG00000271394")
cand_genes_exp <- voom_filt_data$E[which(rownames(voom_filt_data$E) %in% cand_genes), ]

# Retain only those candidate genes that seem to separate the two subtypes
cand_genes_subset <- c("ENSG00000164100", "ENSG00000227706", "ENSG00000223806",
                       "ENSG00000118523", "ENSG00000128218", "ENSG00000164330")
cand_genes_subset_exp <- voom_filt_data$E[which(rownames(voom_filt_data$E) %in% cand_genes_subset), ]

# Create functions to use for clustering of expression data of candidate genes 
dist_func <- function(c) { dist(c, method = "euclidean") }
clustering_func <- function(c) { hclust(c, method = "complete") }

# Create vector of colors representing B-ALL samples as magenta and 
# T-ALL samples as cyan
col_vec <- str_replace(str_replace(sample_meta_data_df$subtype, "B", "#440154FF"), 
                       "T", "#FDE725FF")

# Perform hiercharcical clustering of expression data of candidate genes and
# visualize in a heatmap
pdf("clustering/iCOPE_candidate_genes_clustering.pdf", height = 10, width = 10)
heatmap.2(as.matrix(cand_genes_exp), scale = "none", col = "viridis", 
          distfun = dist_func, hclustfun = clustering_func, labCol = FALSE, 
          ColSideColors = col_vec, trace = "none", density.info = "none", 
          margins = c(10, 9))
legend("topright", legend = c("B-cell ALL","T-cell ALL"), 
       col = c("#440154FF", "#FDE725FF"), lty = 1, lwd = 10, border = FALSE, bty = "n", 
       y.intersp = 0.7, x.intersp = 1, cex = 0.9)
dev.off()

# Perform hiercharcical clustering of expression data of a subset of the 
# candidate genes and visualize in a heatmap
pdf("clustering/iCOPE_candidate_genes_subset_clustering.pdf", height = 8, width = 8)
heatmap.2(as.matrix(cand_genes_subset_exp), scale = "none", col = "viridis", 
          distfun = dist_func, hclustfun = clustering_func, labCol = FALSE, 
          ColSideColors = col_vec, trace = "none", density.info = "none", 
          margins = c(15, 15))
legend("topright", legend = c("B-cell ALL","T-cell ALL"), 
       col = c("#440154FF", "#FDE725FF"), lty = 1, lwd = 10, border = FALSE, bty = "n", 
       y.intersp = 0.7, x.intersp = 1, cex = 0.9)
dev.off()



