### This script contains functions used for survival analysis 


### ---------------------------------------------------------------------------------------------------------

### This function, cox_reg_ph, tests the proportional hazards assumption of
### Cox regression. 
### This function takes three inputs:
### exp_data:       Matrix of expression values with genes in rows and 
###                 samples in columns.
### gene:           Character vector of the gene name
### clinical_data:  A tibble containing the clinical data of the patient 
###                 barcodes. The patient barcodes must be in a column 
###                 called submitter_id. The tibble must contain the
###                 following columns: time_years with the survival time, 
###                 vital_status_binary with the vital status encoded as 0/1,
###                 age_years with the age of the patients in years, and 
###                 gender_binary indicating the gender of the patients
###                 encoded as 0/1. 
### The function returns a tibble of the results of the proportional hazards
### assumption test. 
cox_reg_ph <- function(exp_data, gene, clinical_data) {
  
  ## Fit a Cox regression model for each gene
  
  # Get a tibble of the expression values of respective gene and 
  # barcodes of samples
  exp_gene <- exp_data[gene, ] %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "barcode") %>% 
    as_tibble() %>% 
    dplyr::rename("exp_value" = ".")
  
  # Merge expression values of respective gene with tibble containing
  # clinical data of tumor samples
  clinical_exp_gene <- clinical_data %>% 
    full_join(x = ., y = exp_gene, by = "barcode") %>%
    dplyr::filter(!is.na(time_years),
                  !is.na(vital_status_binary),
                  !is.na(age_years),
                  !is.na(gender_binary))
  
  # Fit a univariate Cox regression model with expression values of respective gene
  # as explanatory variable
  cox_fit_uni <- coxph(Surv(time_years, vital_status_binary)~exp_value, data = clinical_exp_gene)
  
  # Fit a multivariate Cox regression model with expression values, 
  # age, and gender as explanatory variables 
  cox_fit_multi <- coxph(Surv(time_years, vital_status_binary)~exp_value+age_years+gender_binary, data = clinical_exp_gene)
  
  # Check the proportional hazards assumption of the univariate analysis
  ph_check_uni <- cox.zph(cox_fit_uni)
  ph_check_uni <- ph_check_uni$table %>% 
    as_tibble(rownames = "test") %>% 
    dplyr::mutate("test_type" = "univariate")
  
  # Check the proportional hazards assumption of the multivariate analysis
  ph_check_multi <- cox.zph(cox_fit_multi)
  ph_check_multi <- ph_check_multi$table %>% 
    as_tibble(rownames = "test") %>% 
    dplyr::mutate("test_type" = "multivariate")
  
  # Bind results from univariate and multivariate analyses together
  ph_check_results <- ph_check_uni %>% bind_rows(ph_check_multi) %>% 
    dplyr::mutate("gene" = gene)
  
  
  return(ph_check_results)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, cox_reg_univariate, performs univariate Cox regression.
### This function takes three inputs:
### exp_data:       Matrix of expression values with genes in rows and 
###                 samples in columns.
### gene:           Character vector of the gene name
### clinical_data:  A tibble containing the clinical data of the patient 
###                 barcodes. The patient barcodes must be in a column 
###                 called submitter_id. The tibble must contain the
###                 following columns: time_years with the survival time and 
###                 vital_status_binary with the vital status encoded as 0/1
### The function returns a tibble of the results of univariate 
### Cox regression analysis. 
cox_reg_univariate <- function(exp_data, gene, clinical_data) {
  
  ## Fit a univariate Cox regression model for a gene
  
  # Get a tibble of the expression values of respective gene and 
  # barcodes of samples
  exp_gene <- exp_data[gene, ] %>% 
    BiocGenerics::as.data.frame() %>% 
    rownames_to_column(var = "barcode") %>% 
    as_tibble() %>% 
    dplyr::rename("exp_value" = ".")
  
  # Merge expression values of respective gene with tibble containing
  # clinical data of tumor samples
  clinical_exp_gene <- clinical_data %>% 
    full_join(x = ., y = exp_gene, by = "barcode") %>% 
    dplyr::filter(!is.na(time_years),
                  !is.na(vital_status_binary),
                  !is.na(age_years),
                  !is.na(gender_binary))
  
  # Fit a univariate Cox regression model with expression values of respective gene
  # as explanatory variable
  cox_fit_uni <- coxph(Surv(time_years, vital_status_binary)~exp_value, data = clinical_exp_gene)
  
  # Extract the coefficients and p value
  sum_cox_uni <- summary(cox_fit_uni)$coefficients
  sum_cox_uni <- sum_cox_uni %>% as_tibble(rownames = "variable") %>% 
    dplyr::mutate("cox_model_type" = "univariate",
                  "gene" = gene)
  sum_cox_uni_95CI <- exp(confint(cox_fit_uni)) %>% 
    as_tibble(rownames = "variable") %>% 
    full_join(sum_cox_uni, by = "variable") %>% 
    dplyr::relocate(., c(`2.5 %`, `97.5 %`), .after = `Pr(>|z|)`)
  #cox_p_value_uni <- sum_cox_uni$`Pr(>|z|)`
  
  return(sum_cox_uni_95CI)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, cox_reg_multivariate, performs multivariate Cox regression.
### This function takes three inputs:
### exp_data:       Matrix of expression values with genes in rows and 
###                 samples in columns.
### gene:           Character vector of the gene name
### clinical_data:  A tibble containing the clinical data of the patient 
###                 barcodes. The patient barcodes must be in a column 
###                 called submitter_id. The tibble must contain the
###                 following columns: time_years with the survival time, 
###                 vital_status_binary with the vital status encoded as 0/1,
###                 age_years with the age of the patients in years, and 
###                 gender_binary indicating the gender of the patients 
###                 encoded as 0/1. 
### The function returns a tibble of the results of multivariate 
### Cox regression analysis. 
cox_reg_multivariate <- function(exp_data, gene, clinical_data) {
  
  ## Fit a multivariate Cox regression model for a gene
  
  # Get a tibble of the expression values of respective gene and 
  # barcodes of samples
  exp_gene <- exp_data[gene, ] %>% 
    BiocGenerics::as.data.frame() %>% 
    rownames_to_column(var = "submitter_id") %>% 
    as_tibble() %>% 
    dplyr::rename("exp_value" = ".")
  
  # Merge expression values of respective gene with tibble containing
  # clinical data of tumor samples
  clinical_exp_gene <- clinical_data %>% 
    full_join(x = ., y = exp_gene, by = "submitter_id") %>% 
    dplyr::filter(!is.na(time_years),
                  !is.na(vital_status_binary),
                  !is.na(age_years),
                  !is.na(gender_binary))
  
  # Fit a multivariate Cox regression model with expression values,
  # age, and gender as explanatory variables 
  cox_fit_multi <- coxph(Surv(time_years, vital_status_binary)~exp_value+age_years+gender_binary, data = clinical_exp_gene)
  
  # Extract the coefficients, hazard ratio, p value and 95% CI of the hazard ratio
  sum_cox_multi <- summary(cox_fit_multi)$coefficients
  sum_cox_multi <- sum_cox_multi %>% as_tibble(rownames = "variable") %>% 
    dplyr::mutate("cox_model_type" = "multivariate",
                  "gene" = gene)
  sum_cox_multi_95CI <- exp(confint(cox_fit_multi)) %>% 
    as_tibble(rownames = "variable") %>% 
    full_join(sum_cox_multi, by = "variable") %>% 
    dplyr::relocate(., c(`2.5 %`, `97.5 %`), .after = `Pr(>|z|)`)
  #cox_p_value_multi <- sum_cox_multi$`Pr(>|z|)`
  
  return(sum_cox_multi_95CI)
  
}

### ---------------------------------------------------------------------------------------------------------

### ---------------------------------------------------------------------------------------------------------

### This function, km_survival, performs Kaplan-Meier survival analysis. 
### This function takes three inputs:
### exp_data:       Matrix of expression values with genes in rows and 
###                 samples in columns.
### gene:           Character vector of the gene name
### clinical_data:  A tibble containing the clinical data of the patient 
###                 barcodes. The patient barcodes must be in a column 
###                 called submitter_id. The tibble must contain the
###                 following columns: time_years with the survival time, 
###                 vital_status_binary with the vital status encoded as 0/1,
###                 age_years with the age of the patients in years, and 
###                 gender indicating the gender of the patients. 
### The function returns the fitted KM survival object. 
km_survival <- function(exp_data, gene, clinical_data) {
  
  # Get a tibble of the expression values of respective gene and 
  # barcodes of samples
  exp_gene <- exp_data[gene, ] %>% 
    BiocGenerics::as.data.frame() %>% 
    rownames_to_column(var = "submitter_id") %>% 
    as_tibble() %>% 
    dplyr::rename("exp_value" = ".")
  
  # Merge expression values of respective gene with tibble containing
  # clinical data of tumor samples
  clinical_exp_gene <- clinical_data %>% 
    full_join(x = ., y = exp_gene, by = "submitter_id") %>% 
    dplyr::filter(!is.na(time_years),
                  !is.na(vital_status_binary),
                  !is.na(age_years),
                  !is.na(gender_binary))
  
  # Get median expression level of gene
  median_exp <- median(clinical_exp_gene$exp_value)
  
  # Divide samples into two expression groups:
  # one group with expression values > median expression
  # one group with expression values < median expression
  clinical_exp_gene <- clinical_exp_gene %>% 
    dplyr::mutate(exp_group = case_when(exp_value > median_exp ~ "exp_high",
                                        exp_value < median_exp ~ "exp_low"))
  
  ## Do Kaplan-Meier survival analysis 
  
  # Fit a survival curve
  surv_object <- survfit(Surv(time_years, vital_status_binary)~exp_group, data = clinical_exp_gene)
  
  # Test for significant difference in survival between the two groups
  #surv_diff <- survdiff(Surv(time_years, vital_status_binary)~exp_group, data = clinical_exp_gene)
  
  # Visualize Kaplan-Meier survival analysis
  km_plot <- ggsurvplot(surv_object, data = clinical_exp_gene, 
                        risk.table = TRUE, pval = TRUE, conf.int = TRUE) +
    labs(x = "Time [years]", title = paste0("Survival curve of ", gene)) 
  
  return(km_plot)
  
}

### ---------------------------------------------------------------------------------------------------------
