### This script plots gene types of lost genes following normalization step 
### The processing workflow of the gene expression data is:  
### 1) subsetting to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusting for replicates
### 3) preprocessing
### 4) normalizing 
### 5) filtering 


# Load libraries
library(tidyverse)
library(readxl)




## Load data ------------------------------------------------------------------

# Load in excel sheet containing information from Biomart about the excluded genes
genes_biomart <- read_xlsx(path = "TARGET_processing/bone_marrow_primary_comp/ALL_P2_bmp_norm_excl_genes_biomart.xlsx", 
                           col_names = TRUE, 
                           na = "NA")


## Visualize data ------------------------------------------------------------------

# Create barplot of different gene types of excluded genes
gene_types <- genes_biomart %>%
  dplyr::count(gene_type_archive_biomart) %>%
  ggplot(data = ., 
         mapping = aes(x = reorder(gene_type_archive_biomart, -n), 
                       y = n)) +
  geom_col(fill = "darkgreen") +
  labs(x = "Gene type", 
       y = "Number of genes",
       title = "Number of excluded gene types following normalization") +
  scale_y_continuous(breaks = seq(from = 0, to = 110, by = 10)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1, 
                                   colour = "black",
                                   size = 16),
        axis.text.y = element_text(colour = "black", 
                                   size = 16),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.ticks.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 16),
        panel.grid.major.y = element_line(linetype = 1))


## Save data ------------------------------------------------------------------
ggsave(filename = "ALL_P2_bmp_excl_genes_norm_plot.pdf",
       plot = gene_types, 
       path = "TARGET_processing/bone_marrow_primary_comp/",
       width = 12, 
       height = 8)






