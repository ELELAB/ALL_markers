### This script performs PCA of filtered 
### gene expression data focusing on the candidate genes



## Load libraries -------------------------------------------------------------------
library(stats) #4.2.0
library(patchwork) #1.1.3
library(factoextra) #1.0.7
library(FactoMineR) #2.6
library(tidyverse) #2.0.0
library(viridis) #0.6.4

## Source functions
source("pca/iCOPE_pca_functions.R")


## Load data ------------------------------------------------------------------------

# Load voom transformed filtered data
voom_filt_data <- get(load("transform/iCOPE_ALL_filt_data_voom.rda"))

# Load meta data
sample_meta_data_df <- get(load("get_data/iCOPE_ALL_meta_data.rda"))


## Wrangle data ---------------------------------------------------------------------

# MDSPlot requires a vector of group IDs that is used for coloring and labeling. 
# The length of this vector should be equal to the number of columns in the exp dataframe.

# Get vector of subtypes
subtypes <- sample_meta_data_df$subtype

# Define candidate genes
cand_genes <- c("ENSG00000118523", "ENSG00000128218", "ENSG00000164100", 
                "ENSG00000164330", "ENSG00000199683", "ENSG00000199831", 
                "ENSG00000200087",  "ENSG00000200312", "ENSG00000200959",  
                "ENSG00000201901", "ENSG00000202058", "ENSG00000223806", 
                "ENSG00000227706",  "ENSG00000271394")

# Subset voom transformed filtered data to include only candidate genes
voom_filt_cand <- voom_filt_data$E[cand_genes,]

# Subset voom transformed filtered data to exclude candidate genes, i.e.
# so that the matrix includes all genes except the 14 candidate genes
voom_filt_wo_cand <- voom_filt_data$E[!(rownames(voom_filt_data$E) %in% cand_genes), ]

# Perform PCA of expression data using only candidate genes
pca_cand <- PCA(t(voom_filt_cand), 
                graph = FALSE, 
                ncp = 10, 
                scale.unit = FALSE)

# Extract contributions to PC1 and PC2 of the 14 genes from PCA performed
# on 14 candidate genes
gene_contrib_axes <- map(c(1:2), function(x) {
  
  facto_summarize(X = pca_cand, 
                  element = "var", 
                  result = "contrib", 
                  axes = x) %>%
    as_tibble() %>%
    dplyr::arrange(desc(contrib)) %>%
    slice_max(n = 14, 
              order_by = contrib) %>% 
    dplyr::rename(ENSEMBL_ID = name) %>% 
    dplyr::mutate(PCA_dim = rep(str_c("PC_",x), 14))
  
}) %>% 
  bind_rows() 

# Add external gene names to gene contribution tibble
translation_rules <- tribble(
  ~ENSEMBL_ID,      ~gene_name,
  "ENSG00000118523", "CCN2",
  "ENSG00000128218", "VPREB3",
  "ENSG00000164100", "NDST3",
  "ENSG00000164330", "EBF1",
  "ENSG00000199683", "RN7SKP185",
  "ENSG00000199831", "RN7SKP291",
  "ENSG00000200087", "SNORA73B",
  "ENSG00000200312", "RN7SKP255",
  "ENSG00000200959", "SNORA74A",
  "ENSG00000201901", "RN7SKP48",
  "ENSG00000202058", "RN7SKP80",
  "ENSG00000223806", "LINC00114",
  "ENSG00000227706", "novel gene",
  "ENSG00000271394", "7SK")
gene_contrib_axes_gene <- gene_contrib_axes %>%
  left_join(translation_rules, by = c("ENSEMBL_ID" = "ENSEMBL_ID")) %>%
  dplyr::mutate(ENSEMBL_gene = dplyr::coalesce(paste(ENSEMBL_ID, "\n", gene_name), 
                                               ENSEMBL_ID)) 

# Perform PCA of expression data using all genes except candidate genes
pca_wo_cand <- PCA(t(voom_filt_wo_cand), 
                   graph = FALSE, 
                   ncp = 10, 
                   scale.unit = FALSE)

# Perform PCA of expression data of all genes
pca_all <- PCA(t(voom_filt_data$E), 
               graph = FALSE, 
               ncp = 10, 
               scale.unit = FALSE)


## Visualize data ---------------------------------------------------------------------

## Visualize MDS

# MDS of expression data using only candidate genes
MDS_cand <- MDSPlot(my.data = voom_filt_cand,
                    my.group = subtypes, 
                    my.labels = sample_meta_data_df$ids, 
                    my.cols = c("orange", "darkgreen"))
ggsave(filename = "iCOPE_ALL_filt_subtype_MDS_candidate_genes.pdf", 
       plot = MDS_cand, 
       path = "pca/", 
       width = 12, 
       height = 8)

# MDS of expression data excluding candidate genes, i.e. 
# so that the matrix includes all genes except the 14 candidate genes
MDS_wo_cand <- MDSPlot(my.data = voom_filt_wo_cand,
                       my.group = subtypes, 
                       my.labels = sample_meta_data_df$ids, 
                       my.cols = c("orange", "darkgreen"))
ggsave(filename = "iCOPE_ALL_filt_subtype_MDS_wo_candidate_genes.pdf", 
       plot = MDS_wo_cand, 
       path = "pca/", 
       width = 12, 
       height = 8)

# MDS side by side
MDS_with_wo_cand <- MDS_cand + MDS_wo_cand
ggsave(filename = "iCOPE_ALL_filt_subtype_MDS_with_wo_candidate_genes.pdf", 
       plot = MDS_with_wo_cand, 
       path = "pca/", 
       width = 12, 
       height = 8)

## Visualize PCAs
PCA_labels <- "none"
cols <- c("#440154FF", "#FDE725FF")

# Visualize PCA of expression data of only candidate genes
fviz_pca_ind(pca_cand,
             label = PCA_labels, # hide individual labels
             habillage = as.factor(subtypes), # color by groups
             palette = cols,
             addEllipses = TRUE, # concentration ellipses
             repel = TRUE,
             ggtheme = theme_classic(),
             labelsize = 2)
ggsave(filename = "iCOPE_ALL_filt_subtype_PCA_candidate_genes.pdf",
       path = "pca/", width = 5, height = 3)

# Visualize gene contributions to PC1 and PC2 as a stacked barplot
# from PCA performed on candidate genes
gene_contrib_barplot <- ggplot(data = gene_contrib_axes_gene,
                               mapping = aes(fill = factor(PCA_dim, 
                                                           levels = c("PC_2", "PC_1")), 
                                             y = contrib, 
                                             x = reorder(ENSEMBL_gene, -contrib))) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#20A387FF", "#440154FF")) +
  labs(x = "", y = "Contributions (%)",
       title = "Contributions of 14 markers to PC1 and PC2") +
  theme(axis.text.x = element_text(color = "black",
                                   angle = 90),
        legend.title= element_blank())
ggsave(gene_contrib_barplot, 
       filename = "iCOPE_ALL_filt_PCA_contrib_PC1_PC2_ccandidate_genes.pdf",
       path = "pca/", width = 5, height = 3)

# Visualize eigenvalues/variances of dimensions from PCA performed
# on candidate genes
fviz_screeplot(pca_cand, 
               addlabels = TRUE, 
               ylim = c(0, 50), 
               ggtheme = theme_classic(), 
               ncp = 20,
               barfill = "#20A387FF", barcolor = "#20A387FF") 
ggsave(filename = paste0("iCOPE_ALL_filt_PCA_screeplot_candidate_genes.pdf"), 
       path = "pca/", width = 5, height = 3)

# Feature contributions of variables to PC1 from PCA performed
# on candidate genes
fviz_contrib(pca_cand, 
             choice = "var", 
             axes = 1, 
             top = 14, 
             ggtheme = theme_classic(),
             fill = "#20A387FF", color = "#20A387FF") 
ggsave(filename = paste0("iCOPE_ALL_filt_PCA_contrib_PC1_candidate_genes.pdf"), 
       path = "pca/", width = 5, height = 3)

# Feature contributions of variables to PC2 from PCA performed
# on candidate genes
fviz_contrib(pca_cand, 
             choice = "var", 
             axes = 2, 
             top = 14, 
             ggtheme = theme_classic(),
             fill = "#20A387FF", color = "#20A387FF") 
ggsave(filename = paste0("iCOPE_ALL_filt_PCA_contrib_PC2_candidate_genes.pdf"), 
       path = "pca/", width = 5, height = 3)

# Visualize PCA of expression data of all genes except candidate genes
fviz_pca_ind(pca_wo_cand,
             label = PCA_labels, # hide individual labels
             habillage = as.factor(subtypes), # color by groups
             palette = cols,
             addEllipses = TRUE, # concentration ellipses
             repel = TRUE,
             ggtheme = theme_classic(),
             labelsize = 2)
ggsave(filename = "iCOPE_ALL_filt_subtype_PCA_wo_candidate_genes.pdf",
       path = "pca/", width = 5, height = 3)

# Visualize PCA of expression data of all genes
fviz_pca_ind(pca_all,
             label = PCA_labels, # hide individual labels
             habillage = as.factor(subtypes), # color by groups
             palette = cols,
             addEllipses = TRUE, # concentration ellipses
             repel = TRUE,
             ggtheme = theme_classic(),
             labelsize = 2)
ggsave(filename = "iCOPE_ALL_filt_subtype_PCA_all_genes.pdf",
       path = "pca/", width = 5, height = 3)

# Feature contributions of top 50 variables to PC1 from PCA performed
# on all genes
fviz_contrib(pca_all, 
             choice = "var", 
             axes = 1, 
             top = 50, 
             ggtheme = theme_classic(),
             fill = "#20A387FF", color = "#20A387FF") +
  theme(axis.text.x = element_text(angle = 90, color = "black"),
        axis.text.y = element_text(color = "black"))
ggsave(filename = paste0("iCOPE_ALL_filt_PCA_contrib_PC1_all_genes.pdf"), 
       path = "pca/", width = 10, height = 5)

# Feature contributions of top 50 variables to PC2 from PCA performed
# on all genes
fviz_contrib(pca_all, 
             choice = "var", 
             axes = 2, 
             top = 50, 
             ggtheme = theme_classic(),
             fill = "#20A387FF", color = "#20A387FF") +
  theme(axis.text.x = element_text(angle = 90, color = "black"),
        axis.text.y = element_text(color = "black"))
ggsave(filename = paste0("iCOPE_ALL_filt_PCA_contrib_PC2_all_genes.pdf"), 
       path = "pca/", width = 10, height = 5)



