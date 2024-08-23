### This script investigates correlation between predicted subtype-specific
### expression markers and known B- and T-cell markers


# Load libraries
library(tidyverse)
library(viridis)
library(reshape2)


### Load data ---------------------------------------------------------------------------------------------------------

# Load batch corrected gene count data (corrected for year of diagnosis batch factor)
ALL_P2_corr <- get(load("TARGET_batch/bone_marrow_primary_batch_comp/TARGET_ALL_P2_bmp_corrected_year.rda"))


### Wrangle data ---------------------------------------------------------------------------------------------------------

# Subset expression data to contain predicted subtype-specific expression markers
cand_genes <- c("ENSG00000118523", "ENSG00000128218", "ENSG00000164100", 
                "ENSG00000164330", "ENSG00000199683", "ENSG00000199831", 
                "ENSG00000200087",  "ENSG00000200312", "ENSG00000200959",  
                "ENSG00000201901", "ENSG00000202058", "ENSG00000223806", 
                "ENSG00000227706",  "ENSG00000271394")
cand_genes_exp <- corrected_counts_year[which(rownames(corrected_counts_year) %in% cand_genes), ]

# Subset expression data to contain known B/T-cell markers
known_markers <- c("ENSG00000177455", "ENSG00000012124", "ENSG00000272398",
                   "ENSG00000174059", "ENSG00000004468", "ENSG00000105369",
                   "ENSG00000158477", "ENSG00000116824", "ENSG00000010610",
                   "ENSG00000110448", "ENSG00000173762")
known_markers_exp <- corrected_counts_year[which(rownames(corrected_counts_year) %in% known_markers), ]

# Combine the two expression matrices
comb_markers_exp <- rbind(cand_genes_exp, known_markers_exp)

# Calculate correlation matrix between all expression markers (both predicted and known)
cor_matrix <- cor(t(comb_markers_exp), method = "spearman")

# Subset rows of correlation matrix to include only predicted markers
cor_matrix <- cor_matrix[which(rownames(cor_matrix) %in% cand_genes), ]

# Subset columns of correlation matrix to include only known markers
cor_matrix <- cor_matrix[, which(colnames(cor_matrix) %in% known_markers)]

# Reformat the correlation matrix
matrix_melt <- melt(cor_matrix)
colnames(matrix_melt) <- c("Predicted_markers", "Known_markers", "Correlation")


### Visualize data ---------------------------------------------------------------------------------------------------------

# Create heatmap of correlation values
cor_heatmap <- ggplot(matrix_melt, mapping = aes(x = Known_markers, y = Predicted_markers, fill = Correlation)) +
  geom_tile() +
  scale_fill_viridis(option = "D") +
  labs(title = "Heatmap of correlation values between markers") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 90, hjust = -1),
        axis.title.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black"),
        title = element_text(size = 8, color = "black"),
        panel.border = element_blank())


### Save data ---------------------------------------------------------------------------------------------------------

# Save correlation matrix
write_csv(matrix_melt, file = "TARGET_known_markers/bone_marrow_primary_known_markers_comp/TARGET_ALL_P2_bmp_corr_known_markers_mat.csv")

# Save correlation heatmap
ggsave(filename = "TARGET_ALL_P2_bmp_corr_known_markers_heatmap.pdf", 
       plot = cor_heatmap, 
       path = "TARGET_known_markers/bone_marrow_primary_known_markers_comp/", 
       width = 5, 
       height = 4)

