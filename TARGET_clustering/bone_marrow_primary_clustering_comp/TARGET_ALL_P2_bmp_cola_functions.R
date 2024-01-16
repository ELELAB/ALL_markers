### This script contains functions used for cola clustering analysis of TARGET-ALL-P2 gene expression data 


### ---------------------------------------------------------------------------------------------------------

# This function, confusion_mat_plot, compares predicted class labels of samples when two clusters
# are used with the actual class labels through a confusion matrix and visualizes this confusion
# matrix. 
# This function takes seven arguments: 
# cola_results_list:            a list containing results from the cola analysis
#                               where each element in the list is the output of a 
#                               clustering method
# method_name:                  name of the clustering method
# data_info:                    a tibble containing information/metadata of the dataset 
#                               on which clustering was performed. It should include
#                               an identier of the samples (here, barcodes) and the actual
#                               class labels of the samples
# actual_subtype_col:           a character string indicating the name of the column in
#                               data_info which contains the actual class labels of the samples
# name_subtype_1:               a character string indicating the first name of the class
# name_subtype_2:               a character string indicating the second name of the class
# title_plot:                   a character string containing the title of the plot
# The function returns a ggplot visualizing the confusion matrix. 
confusion_mat_plot <- function(cola_results_list, 
                               method_name,
                               data_info,
                               actual_subtype_col,
                               name_subtype_1,
                               name_subtype_2,
                               title_plot) {
  
  # Extract results for specific method
  res <- cola_results_list[method_name]
  
  # Compare predicted classes with actual classes
  compare_classes <- get_classes(res, 
                                 k = 2) %>% 
    rownames_to_column(., 
                       var = "barcodes") %>% 
    as_tibble() %>%
    full_join(x = ., 
              y = data_info, 
              by = "barcodes") %>% 
    dplyr::select(barcodes, 
                  class, 
                  entropy, 
                  silhouette, 
                  actual_subtype_col) %>%
    dplyr::rename(predicted_class = class) %>%
    mutate(actual_class = case_when(data_info[actual_subtype_col] == name_subtype_1 ~ 1,
                                    data_info[actual_subtype_col] == name_subtype_2 ~ 2) %>%
             as.factor(),
           predicted_class = as.factor(predicted_class))
  
  # Create confusion matrix by comparing predicted and actual class labels
  confusion_mat <- confusionMatrix(data = compare_classes$predicted_class,
                                   reference = compare_classes$actual_class)$table %>%
    as_tibble() %>%
    dplyr::rename(Predicted_class = Prediction,
                  Actual_class = Reference,
                  Number_of_samples = n)
  
  # Visualize confusion matrix 
  confusion_mat_visualize <- ggplot(data = confusion_mat,
                                    mapping = aes(x = Actual_class,
                                                  y = Predicted_class)) +
    geom_tile(mapping = aes(fill = Number_of_samples)) +
    geom_text(mapping = aes(label = sprintf("%1.0f",
                                            Number_of_samples)),
              size = 5) +
    scale_fill_viridis(limits = c(0,
                                  250)) +
    labs(title = title_plot) +
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
          panel.border = element_blank(),
          panel.grid = element_blank())
  
  # Return plot of confusion matrix
  return(confusion_mat_visualize)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

# This function, kruskall_wallis_test, performs the Kruskall-Wallis test on gene 
# expression values for each gene in the input data matrix against the predicted
# class labels. This function is used in the get_signatures function in the Cola
# framework to investigate gene signatures. 
# This function takes two arguments: 
# mat:              A data matrix containing with genes in rows and samples
#                   in columns
# class:            A vector containing the predicted class labels for each 
#                   sample. The order must match the order of samples in the
#                   data matrix. 
# The function returns a vector of FDR values. 
kruskall_wallis_test <- function(mat, class) {
  
  # Convert the class vector to factor
  class <- as.numeric(factor(class))
  n_class <- length(unique(class))
  
  # Transpose the data matrix and convert to tibble
  mat <- t(mat) %>%
    as_tibble()
  
  # For each gene (now, each column), perform the Kruskall-Wallis test
  # against the predicted class labels
  kw_test <- map(colnames(mat), function(gene) {
    gene_data <- mat %>%
      dplyr::select(all_of(gene)) %>% 
      pull
    
    gene_pvalue <- tibble(gene_name = gene, 
                          pvalue_kw = kruskal.test(gene_data ~ class)$p.value)
    
  })
  
  # Combine results from all genes into one table
  kw_test <- bind_rows(kw_test)
  
  # Apply multiple testing correction by calculating FDR values
  FDR <- p.adjust(kw_test$pvalue_kw, 
                  method = "BH")
  
  return(FDR)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, ensembl_hugo_conversion, converts ENSEMBL IDs to gene names.
### The function takes one input:
### ensembl_id:     A list of ENSEMBL IDs to be converted into gene names.
### The function returns the ENSEMBL IDs and their corresponding gene name.
### The function is based on code from the CAMPP pipeline.
ensembl_hugo_conversion <- function(ensembl_id) {
  
  # Connect to biomart database and dataset hosted by ENSEMBL (release 106)
  ensembl_106 <- useEnsembl(biomart = 'genes', 
                            dataset = 'hsapiens_gene_ensembl',
                            version = 106)
  
  # Convert ENSEMBL IDs to HUGO gene names
  gene_names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      mart = ensembl_106,
                      values = ensembl_id,
                      filters = 'ensembl_gene_id',
                      uniqueRows = TRUE,
                      bmHeader = T)
  
  return(gene_names)
  
}

### ---------------------------------------------------------------------------------------------------------

