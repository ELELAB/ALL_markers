## This script finds investigates drug-gene interactions
## among the predicted gene expression markers between
## ALL subtypes (B- and T-ALL) and between predicted 
## clusters of T-ALL



## Load libraries ------------------------------------
library(tidyverse)
library(rDGIdb)


## Load data ------------------------------------

# Load defined small gene set of subtype-related markers
subtype_genes_info <- read_csv("TARGET_compare_genes/bone_marrow_primary_compare_comp/TARGET_ALL_P2_bmp_consensus_DEGs_elastic_pca.csv")

# Load predicted markers between predicted T-ALL clusters
cluster_genes_info <- read_csv("TARGET_random_forest/bone_marrow_primary_random_forest_comp/TARGET_ALL_P2_bmp_selected_features_T_clusters_rf.csv")


## Wrangle data ----------------------------------

# Extract subtype genes from full table
subtype_genes <- subtype_genes_info$Gene_name 

# Extract cluster genes from full table
cluster_genes <- cluster_genes_info$Gene_name.x


## Analyze data ----------------------------------

# Query DGIdb for biomarker-drug interactions
# using only cancer-specific source databases
drug_gene_interact_subtype <- queryDGIdb(subtype_genes, 
                                         sourceDatabases = c("CGI", "CIViC", "COSMIC",
                                                             "CancerCommons", 
                                                             "ClearityFoundationBiomarkers",
                                                             "ClearityFoundationClinicalTrial",
                                                             "DoCM", "JAX-CKB", "MyCancerGenome",
                                                             "MyCancerGenomeClinicalTrial",
                                                             "NCI", "OncoKB", "TALC"))
drug_gene_interact_cluster <- queryDGIdb(cluster_genes, 
                                         sourceDatabases = c("CGI", "CIViC", "COSMIC",
                                                             "CancerCommons", 
                                                             "ClearityFoundationBiomarkers",
                                                             "ClearityFoundationClinicalTrial",
                                                             "DoCM", "JAX-CKB", "MyCancerGenome",
                                                             "MyCancerGenomeClinicalTrial",
                                                             "NCI", "OncoKB", "TALC"))

# Extract table showing number of drugs each driver gene interacts with
drug_gene_count_subtype <- drug_gene_interact_subtype@byGene %>% 
  as_tibble() 
drug_gene_count_cluster <- drug_gene_interact_cluster@byGene %>% 
  as_tibble() 

# Extract table showing which drugs interact with which driver genes
# Add annotation if gene and drug interact (yes = interaction)
drug_gene_names_subtype <- drug_gene_interact_subtype@detailedResults %>% 
  as_tibble() %>% 
  dplyr::mutate(Interaction = "Yes")
drug_gene_names_cluster <- drug_gene_interact_cluster@detailedResults %>% 
  as_tibble() %>% 
  dplyr::mutate(Interaction = "Yes")


## Visualize data ----------------------------------

# Visualize drug-gene interactions in heatmap
drug_heatmap <- drug_gene_names_subtype %>%
  ggplot(., mapping = aes(x = Drug, y = Gene,
                          fill = Interaction)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("#21908CFF")) +
  scale_x_discrete(position = "top", guide = guide_axis(angle = 90)) +
  scale_y_discrete(position = "right") +
  labs(x = "", y = "") +
  theme_classic() +
  theme(line = element_blank(),
        axis.text.x = element_text(size = 12,
                                   color = "black"),
        axis.text.y = element_text(size = 12,
                                   color = "black"),
        title = element_text(size = 18,
                             color = "black"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, "cm"))



## Save data ------------------------------------

# Save results of drug-gene interactions of subtype genes
write_csv(drug_gene_names_subtype, file = "TARGET_drug_targets/bone_marrow_primary_drug_targets_comp/TARGET_ALL_P2_bmp_drugs_subtype_genes.csv")

# Save results of drug-gene interactions of cluster genes
write_csv(drug_gene_names_cluster, file = "TARGET_drug_targets/bone_marrow_primary_drug_targets_comp/TARGET_ALL_P2_bmp_drugs_cluster_genes.csv")

# Save heatmap showing drug-gene interactions 
ggsave(filename = "TARGET_ALL_P2_bmp_drugs_subtype_genes_heatmap.pdf",
       plot = drug_heatmap,
       path = "TARGET_drug_targets/bone_marrow_primary_drug_targets_comp/",
       width = 8,
       height = 2.8)





