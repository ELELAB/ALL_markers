## This script contains functions to obtain information about TARGET-ALL-P2 data

### This function, get_meta_data, creates a table containing number of samples of 
### each unique case of user-specified variables stored in a SE object. 
### This function takes two inputs:
### SE_exp_data:    SummarizedExperiment object of data containing the variables
### var_vector:     A string vector containing names of those variables the 
###                 number of samples should be found. 
### The function returns a tibble containing the number of samples for each 
### unique case of a veriable. 
get_meta_data <- function(SE_exp_data, var_vector) {
  
  # Convert column data in SummarizedExperiment object to tibble
  col_data <- as_tibble(SE_exp_data@colData@listData)
  
  # Split vector of variables  
  var_vector_split <- str_split(var_vector, ",") %>% 
    all_of()
  
  # Count number of samples for each variable
  # Combine lists containing number of samples for each variable into one tibble
  num_samples <- map(var_vector_split, function(x) { 
    
    col_data %>%
      dplyr::count(across(x)) %>% 
      pivot_wider(names_from = x, values_from = n)
    
  }) %>% 
    bind_cols()
  
  # Add number of genes and number of total samples 
  num_samples <- num_samples %>% mutate("genes" = nrow(SE_exp_data),
                                        "total_samples" = ncol(SE_exp_data))
  
  return(num_samples)
  
}




