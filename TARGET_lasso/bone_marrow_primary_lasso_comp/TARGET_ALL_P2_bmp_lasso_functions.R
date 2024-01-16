### This script contains functions used for LASSO regression of TARGET-ALL-P2 gene expression data 

### ---------------------------------------------------------------------------------------------------------

# This function, LASSOFeature, performs LASSO regression for gene count data. 
# This function takes six arguments: 
# my.seed:            an integer to denote a seed
# my.data:            a data matrix containing gene counts with genes in rows 
#                     and samples in columns
# my.group:           vector of integers specifying the group of each sample
# my.LAorEN:          alpha parameter, an integer between 0 and 1 where
#                     1 denotes LASSO regression and 0 is ridge regression and
#                     a number between 0 and 1 denotes elastic net regression.
# my.validation:      Boolean, if TRUE validation will be performed on test data.
#                     If FALSE, validation will not be performed on test data.
# my.multinorm:       Boolean, if TRUE multinomial regression will be performed.
#                     If FALSE, binomial regression will be performed.
# The function returns a list of three elements:
# lasso.fit:          The fit returned when fitting the model.
# coefficients:       The coefficients for each of the features input to the model.
# pred_error:         The prediction error when predicting the test data. 
# Code in function is modified from code associated with papers 
# 10.1371/journal.pcbi.1007665 and 10.1002/1878-0261.12850
LASSOFeature <- function(my.seed, 
                         my.data, 
                         my.group, 
                         my.LAorEN,
                         n.cv.folds) {
  
  # Divide my.data into a validation dataset and a training dataset
  
  ll <- list()
  llev <- levels(as.factor(my.group))
  
  # Generate a list of length number of groups 
  # The first element in the list is the indexes of the samples belonging to group 1
  # The second element in the list is the indexes of the samples belonging to group 2
  for (idx in 1:length(llev)) {
    
    pos <- which(my.group == as.character(llev[idx]))
    ll[[idx]] <- pos
    
  }
  
  # Generate randomly a vector of 1/4 of the number of samples 
  # These samples will be used as validation data set
  # It will contain samples from both groups
  my.samp <- unlist(lapply(ll, 
                           function(x) sample(x, 
                                              ceiling((length(x)/4)))))
  
  # Generate the validation data set
  my.data.val <- my.data[,my.samp]
  my.group.val <- as.integer(my.group[my.samp])
  
  # Generate the training data set, i.e. the original data minus the validation data
  my.data <- my.data[,-my.samp]
  my.group <- as.integer(my.group[-my.samp]) 
  
  # Set seed
  set.seed(my.seed)
  
  # Fit binomial logistic regression model with LASSO regularization and do 10-fold cross-validation 
  # Use misclassification error as the criterion for 10-fold CV
  my.fit <- cv.glmnet(x = t(my.data), 
                      y = my.group, 
                      family = "binomial", 
                      type.measure = "class", 
                      nfolds = n.cv.folds, 
                      alpha = my.LAorEN)
  
  # Obtain model coefficients when using lambda = lambda.min
  my.coef <- coef(my.fit, 
                  s = my.fit$lambda.min)
  
  # Force my.coef object to be a matrix
  # my.ma is a vector with length number of genes in dataset and contains the
  # coefficient for each gene as well as the intercept
  my.ma <- as(my.coef, 
              "matrix")
  
  # Use trained model to predict validation dataset (test data) 
  # Compare predicted class for each sample with actual group labels stored in my.group.val
  # Find those samples where there is NOT correspondence between predicted class and actual class
  # Take the mean of comparisons to obtain misclassification error and convert to %
  pred_error <- as.numeric(mean(predict(my.fit, 
                                        t(my.data.val), 
                                        s = my.fit$lambda.min, 
                                        type="class") != my.group.val)) * 100
  
  # Return the minimal subset of genes and the misclassification error
  return(list("lasso_fit" = my.fit,
              "coefficients" = my.ma, 
              "pred_error" = pred_error))
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

# This function, get_coef, extracts those genes that do not have a coefficient
# equal to 0 from each LASSO seed run. 
# This function takes one argument: 
# LASSO_res:          The results of each LASSO seed run stored in a list
#                     where each element in the list is a LASSO seed run
# The function returns a list where each element of the list corresponds to a
# LASSO seed run. Each element of the list contains a list of genes that do not
# have a coefficient equal to 0 and the value of the coefficient. 
get_coef <- function(LASSO_res) { 
  
  # For each seed run, extract those genes that do not have a coefficient equal to 0
  gene_subsets <- map(LASSO_res, 
                      function(x) { 
                        genes = x[["coefficients"]] %>% 
                          as_tibble(rownames = "ENSEMBL_ID") %>% 
                          filter(s1 != 0) 
                      })
  
  return(gene_subsets)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

# This function, intersect_genes, intersects the genes that do not have a coefficient
# equal to 0 across the seed runs. 
# This function takes one argument: 
# LASSO_genes:        A list where each element in the list contains the genes from 
#                     each seed run that do not have a coefficient equal to 0. 
# The function returns a vector of those genes that are the intersection of the seed runs, 
# i.e. between those genes from each seed run that do not have a coefficient equal to 0. 
intersect_genes <- function(LASSO_genes) {
  
  # Intersect the minimal subset of genes from each LASSO seed run
  # and convert to tibble and remove the intercept 
  genes <- Reduce(intersect, 
                  map(LASSO_genes, 
                      function(x) { 
                        x[["ENSEMBL_ID"]] })) %>% 
    as_tibble() %>% 
    filter(!value == "(Intercept)") %>% 
    dplyr::rename("selected_genes" = value)
  
  return(genes)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

# This function, get_cv_error, extracts the prediction errors from each LASSO seed run. 
# This function takes one argument: 
# LASSO_res:        A list where each element in the list contains the results from LASSO
#                   regression. 
# The function returns a tibble with the prediction error from each seed run.  
get_cv_error <- function(LASSO_res) {
  
  # Extract the prediction error from each seed run and convert to tibble
  cv_error_runs <- map(LASSO_res, 
                       function(x) {
                         x[["pred_error"]] }) %>% 
    unlist() %>% 
    as_tibble(rownames = "seed_run") %>% 
    dplyr::rename("pred_error" = value) %>% 
    mutate(seed_run = str_c("Run ", seq(1:10)),
           seed_run = factor(seed_run, levels = unique(seed_run)))
  
  return(cv_error_runs)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

# This function, cv_error_plot, plots the prediction errors from each LASSO seed run. 
# This function takes two arguments: 
# seed_runs_col:        A vector indicating the name of the seed runs.
# pred_error_col:       A vector of the prediction errors from each seed run.
# The function returns a ggplot bar plot of the prediction errors for each seed run. 
cv_error_plot <- function(seed_runs_col,
                          pred_error_col) {
  
  # Barplot of prediction errors for each seed run
  ggplot(mapping = aes(x = seed_runs_col, 
                       y = pred_error_col)) +
    geom_bar(stat = "identity",
             fill = "darkgreen") + 
    theme_minimal() +
    labs(x = "",
         y = "Prediction error (misclassified samples) in %",
         title = "Prediction error in % for each of the 10 seed runs") +
    theme(axis.text.x = element_text(size = 16, 
                                     color = "black"),
          axis.text.y = element_text(size = 16, 
                                     color = "black"),
          axis.title.y = element_text(size = 16, 
                                      color = "black"),
          title = element_text(size = 16,
                               color = "black"))
  
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


### ---------------------------------------------------------------------------------------------------------

### This function, venn_diagram, creates a Venn diagram of intersections between
### different sets. This function takes five inputs:
### sets_list:        A list containing those sets that are to be compared 
### names_sets:       A character vector containing names of each element in the list, sets_list
### title_plot:       The title to be included in the Venn diagram
### dir_output:       The path to the directory where the plot should be saved
### file_name:        File name of plot that will be saved; best if saved as .png
### The function returns a Venn diagram.  
venn_diagram <- function(sets_list,
                         names_sets,
                         title_plot,
                         dir_output,
                         file_name) {
  
  # Assign names to list of sets
  names(sets_list) <- names_sets
  
  # Define colors to use in Venn diagram
  n_col <- length(sets_list)
  palette_col <- viridis_pal(alpha = 0.8,
                             option = "viridis")(n_col)
  
  #Create Venn diagram and save it as pdf in plots directory of TCGA project
  venn.diagram(sets_list, 
               filename = paste0(dir_output, 
                                 "/", 
                                 file_name), 
               category.names = names_sets, 
               output = TRUE, 
               main = title_plot, 
               main.cex = 0.5, 
               main.fontface = "bold", 
               main.fontfamily = "sans", 
               lwd = 1, 
               lty = "blank", 
               fill = palette_col, 
               cex = 0.5, 
               fontface = "bold", 
               fontfamily = "sans",
               cat.cex = 0.5, 
               cat.fontface = "bold", 
               cat.fontfamily = "sans")
  
}

### ---------------------------------------------------------------------------------------------------------









