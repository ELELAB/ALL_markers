### This script visualizes ENSEMBL gene IDs selected from elastic net logistic regression 
### performed on TARGET ALL P2 gene expression data
### This data has been 
### 1) subsetted to contain only bone marrow (primary blood derived cancer) samples
### 2) adjusted for replicates
### 3) preprocessed
### 4) normalized using updated geneInfoHT table
### 5) filtered
### 6) voom transformed
### 7) batch corrected for year of diagnosis


# Load libraries
library(tidyverse)
library(viridis)


## Load data ---------------------------------------------------------------------------------------------------------

# Load table containing selected genes from elastic net logistic regression 
# together with DEA and biotype information
elastic_consensus_biotypes <- read_csv("TARGET_lasso/bone_marrow_primary_lasso_comp/elastic_consensus_biotypes/TARGET_ALL_P2_bmp_elastic_consensus_biotype.csv")


## Visualize data ---------------------------------------------------------------------------------------------------------

## Boxplot of elastic net coefficient values of elastic net consensus genes
## stratified by sign of DEG and biotype
consensus_en_boxplot <- elastic_consensus_biotypes %>% 
  dplyr::select(ENSEMBL_ID, 
                gene_name,
                biotype,
                contains("seed"),
                mean_s1,
                logFC_mean) %>% 
  pivot_longer(data = ., 
               cols = contains("seed"), 
               names_to = "seed_run", 
               values_to = "coefficient") %>% 
  mutate(sign_DEG = case_when(logFC_mean > 0 ~ "up DEG",
                              logFC_mean < 0 ~ "down DEG"),
         ENSEMBL_gene_name = case_when(is.na(gene_name) ~ ENSEMBL_ID,
                                       TRUE ~ str_c(ENSEMBL_ID,
                                                    "\n",
                                                    gene_name))) %>% 
  arrange(mean_s1) %>% 
  ggplot(data = ., 
         mapping = aes(x = reorder(ENSEMBL_gene_name, 
                                   mean_s1),
                       y = coefficient,
                       color = sign_DEG,
                       fill = biotype)) +
  geom_boxplot() +
  geom_hline(yintercept = 0,
             linetype = "dashed") +
  scale_fill_manual(values = viridis_pal(alpha = 1,
                                         option = "H")(elastic_consensus_biotypes %>% 
                                                         distinct(biotype) %>% 
                                                         nrow)) +
  scale_color_manual(values = c("darkgreen",
                                "darkorange")) +
  labs(x = "",
       y = "Coefficient value",
       title = "Boxplot of elastic net coefficient values across 10 seed runs for 31 elastic net consensus genes") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 18, 
                                    color = "black"),
        axis.text.y = element_text(size = 18, 
                                   color = "black"),
        axis.text.x = element_text(size = 18, 
                                   color = "black", 
                                   angle = 90, 
                                   hjust = 0),
        axis.ticks = element_line(size = 2),
        title = element_text(size = 18),
        panel.border = element_rect(size = 2),
        legend.text = element_text(size = 18),
        legend.key.size = unit(1.3, "cm")) 




## Save data ---------------------------------------------------------------------------------------------------------

# Save table of elastic net consensus genes including biotypes
write_csv(elastic_consensus_biotypes, 
          file = "TARGET_lasso/bone_marrow_primary_lasso_comp/elastic_consensus_biotypes/TARGET_ALL_P2_bmp_elastic_consensus_biotype.csv")

# Save biotypes barplot
ggsave(filename = "TARGET_ALL_P2_bmp_elastic_consensus_genes_biotype.pdf",
       plot = consensus_en_boxplot, 
       path = "TARGET_lasso/bone_marrow_primary_lasso_comp/elastic_consensus_biotypes/",
       width = 24, 
       height = 12)
  
  
  
  





  
  