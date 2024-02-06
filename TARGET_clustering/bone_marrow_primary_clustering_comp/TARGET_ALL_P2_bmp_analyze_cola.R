### This script analyzes unsupervised clustering results done on TARGET ALL P2 gene expression data using the cola framework


# Load libraries
library(cola)
library(tidyverse)
library(caret)
library(viridis)
library(gridExtra)
library(cowplot)
library(umap)


# Source functions
source("TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_cola_functions.R")


### Load data ---------------------------------------------------------------------------------------------------------

# Load results of unsupervised clustering performed on batch corrected data using cola 
cola_results_batch <- readRDS("TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_batch_cola_results.rds")

# Load results of unsupervised clustering performed on raw data using cola 
cola_results_raw <- readRDS("TARGET_clustering/bone_marrow_primary_clustering_comp/cola_on_raw_data/TARGET_ALL_P2_bmp_raw_cola_results.rds")

# Load table containing information of samples
ALL_P2_info <- read_csv(file = "TARGET_pca/bone_marrow_primary_pca_comp/TARGET_ALL_P2_bmp_batch_info.csv")


### Wrangle and visualize data ---------------------------------------------------------------------------------------------------------

## Compare predicted two clusters (k=2) with actual labels of B/T ALL (on batch corrected data and on raw data) using all clustering methods

## Batch corrected data

# Extract names of all clustering methods
all_methods_batch <- rownames(suggest_best_k(cola_results_batch)) %>%
  sort()

## Compare predicted two clusters (k=2) with actual labels of B/T ALL on batch corrected data using all clustering methods
confusion_mat_plot_batch_list <- map(all_methods_batch,  function(x) {
  
  # Visualize confusion matrix comparing predicted and actual class labels when using each of the different clustering methods on batch corrected data
  cm_plot_batch <- confusion_mat_plot(cola_results_list = cola_results_batch, 
                                      method_name = x, 
                                      data_info = ALL_P2_info, 
                                      actual_subtype_col = "subtype_abbr", 
                                      name_subtype_1 = "B-cell ALL", 
                                      name_subtype_2 = "T-cell ALL",
                                      title_plot = x)
  
})

# Assign name of method to each element in list
names(confusion_mat_plot_batch_list) <- all_methods_batch


## Plot confusion matrices together in one plot

# Extract legend from first plot (since all individual plots have identical legend)
legend_batch <- get_legend(confusion_mat_plot_batch_list[[1]])

# Remove legends and axis labels from all individual plots and plot all plots together
plot_wo_legend_batch <- map(confusion_mat_plot_batch_list, 
                            function(x) {x + 
                                theme(legend.position = "none")}) %>%
  plot_grid(plotlist = ., 
            ncol = 5, 
            nrow = 4) 

# Add legend to combined plot
plot_comb_batch <- plot_grid(plot_wo_legend_batch, 
                             legend_batch, 
                             rel_widths = c(1, 
                                            .3)) 


## Raw data

# Extract names of all clustering methods
all_methods_raw <- rownames(suggest_best_k(cola_results_raw)) %>%
  sort()

## Compare predicted two clusters (k=2) with actual labels of B/T ALL on batch corrected data using all clustering methods
confusion_mat_plot_raw_list <- map(all_methods_raw,  function(x) {
  
  # Visualize confusion matrix comparing predicted and actual class labels when using each of the different clustering methods on batch corrected data
  cm_plot_raw <- confusion_mat_plot(cola_results_list = cola_results_raw, 
                                    method_name = x, 
                                    data_info = ALL_P2_info, 
                                    actual_subtype_col = "subtype_abbr", 
                                    name_subtype_1 = "B-cell ALL", 
                                    name_subtype_2 = "T-cell ALL",
                                    title_plot = x)
  
})

# Assign name of method to each element in list
names(confusion_mat_plot_raw_list) <- all_methods_raw


## Plot confusion matrices together in one plot

# Extract legend from first plot (since all individual plots have identical legend)
legend_raw <- get_legend(confusion_mat_plot_raw_list[[1]])

# Remove legends from all individual plots and plot all plots together
plot_wo_legend_raw <- map(confusion_mat_plot_raw_list, 
                          function(x) {x + 
                              theme(legend.position = "none")} ) %>%
  plot_grid(plotlist = ., 
            ncol = 5, 
            nrow = 4) 

# Add legend to combined plot
plot_comb_raw <- plot_grid(plot_wo_legend_raw, 
                           legend_raw, 
                           rel_widths = c(1, 
                                          .3))


## Visualize statistical metrics of 10 methods that could 100% correctly 
## cluster the B- and T-cell ALL samples into their own separate clusters 
## (when using batch corrected data)

# Define the 10 methods
methods_correct <- c("ATC:kmeans", "ATC:mclust", "CV:kmeans", "MAD:hclust", "MAD:kmeans", 
                     "MAD:mclust", "SD:hclust", "SD:kmeans", "SD:mclust", "SD:pam")

# Obtain statistic metrics for each of the 10 methods
methods_stats <- map(methods_correct, function(x) {
  
  # Extract results for method 
  res <- cola_results_batch[[x]] 
  
  # Get statistic metrics for each method for its optimal k 
  get_stats(res, 
            k = suggest_best_k(res)) %>%
    as_tibble() %>%
    mutate(method_name = x)
  
}) %>% # Combine statistic metrics for all methods into one table
  bind_rows(.) %>% # Remove methods with optimal k = 2
  dplyr::filter(k != 2) %>% # Combine method name with optimal k into one string for x axis
  dplyr::mutate(x_axis_names = str_c(method_name, 
                                     ", k=", 
                                     k)) %>% # Convert data to long format needed for plotting
  pivot_longer(data = .,
               cols = -c(method_name, 
                         x_axis_names),
               names_to = "metric",
               values_to = "value") 

# Visualize the mean silhouette, concordance, and 1-PAC scores
# These scores are used to determine the best performing method
stats_plot <- ggplot(data = methods_stats %>%
                       dplyr::filter(metric != "k") %>%
                       dplyr::filter(metric != "area_increased"), 
                     mapping = aes(x = x_axis_names, 
                                   y = value, 
                                   fill = metric)) + 
  geom_bar(position = "dodge", 
           stat = "identity") +
  scale_fill_manual(values = viridis_pal(option = "D")(methods_stats %>%
                                                         dplyr::filter(metric != "k") %>%
                                                         distinct(metric) %>%
                                                         pull %>%
                                                         length)) +
  labs(title = "Statistical measures of clustering methods",
       x = "Clustering methods and selected k",
       y = "Statistical measure") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12,
                                   color = "black"),
        axis.title.x = element_text(size = 12,
                                    color = "black"),
        axis.text.y = element_text(size = 12,
                                   color = "black"),
        axis.title.y = element_text(size = 12,
                                    color = "black"),
        axis.ticks = element_blank(),
        title = element_text(size = 12,
                             color = "black"),
        panel.border = element_blank())


## Get table of predicted clusters of samples using SD:mclust k = 4 method 

# Get predicted cluster of each sample and add actual class labels
samples_clusters <- get_classes(cola_results_batch["SD:mclust"], 
                                k = 4) %>%
  bind_cols(get_membership(cola_results_batch["SD:mclust"], 
                           k = 4)) %>%
  rownames_to_column(var = "barcodes") %>%
  full_join(x = .,
            y = ALL_P2_info,
            by = "barcodes") %>%
  dplyr::select(barcodes, class, entropy, silhouette, p1, p2, p3, p4, subtype_abbr)


## UMAP visualization of clusters using SD:mclust k = 4 method

pdf(file = "TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_umap_SD_mclust.pdf",
    width = 10,
    height = 8)

# Perform UMAP
umap_SD_mclust <- dimension_reduction(object = cola_results_batch["SD:mclust"],
                                      k = 4,
                                      method = "UMAP")

dev.off()

## UMAP visualization of clusters using MAD:mclust k = 4 method

pdf(file = "TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_umap_MAD_mclust.pdf",
    width = 10,
    height = 8)

# Perform UMAP
umap_MAD_mclust <- dimension_reduction(object = cola_results_batch["MAD:mclust"],
                                       k = 4,
                                       method = "UMAP")

dev.off()


### Save data ---------------------------------------------------------------------------------------------------------

# Save table containing predicted cluster labels of samples together with probability of membership and actual class label
#write_csv(samples_clusters, 
#          file = "TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_cola_samples_clusters.csv")

# Save confusion matrix plots of batch corrected data using all methods in one pdf
ggsave(filename = "TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_cola_batch_confusion_matrices.pdf", 
       plot = marrangeGrob(confusion_mat_plot_batch_list, 
                           nrow = 1, 
                           ncol = 1), 
       device = "pdf",
       width = 10,
       height = 8)

# Save confusion matrix plots of raw data using all methods in one pdf
ggsave(filename = "TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_cola_raw_confusion_matrices.pdf", 
       plot = marrangeGrob(confusion_mat_plot_raw_list, 
                           nrow = 1, 
                           ncol = 1), 
       device = "pdf",
       width = 10,
       height = 8)

# Save plot that contains confusion matrices together in one plot (batch corrected data)
ggsave(filename = "TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_cola_batch_confusion_matrices_collective.pdf",
       plot = plot_comb_batch,
       device = "pdf",
       width = 10,
       height = 8)

# Save plot that contains confusion matrices together in one plot (raw data)
ggsave(filename = "TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_cola_raw_confusion_matrices_collective.pdf",
       plot = plot_comb_raw,
       device = "pdf",
       width = 10,
       height = 8)

# Save grouped barplot showing statistical measures of clustering methods
ggsave(filename = "TARGET_clustering/bone_marrow_primary_clustering_comp/TARGET_ALL_P2_bmp_cola_stats.pdf",
       plot = stats_plot,
       device = "pdf",
       width = 10,
       height = 8)











