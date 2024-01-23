## This script contains functions used for random forest


### ---------------------------------------------------------------------------------------------------------

# This function, ForestFeatures, performs random forest on gene expression data.
# This function takes seven arguments: 
# seed:                an integer containing a random seed number. Default is 123.
# data:                a matrix of (transformed and normalized) feature counts from
#                      "seq", "array", "ms" or "other" technology (with feature IDs as row names
#                      and sample IDs as columns). If available, batch corrected data should be used.
# group:               a factor specifying group for each sample (e.g. could be
#                      represented by a column from a metadata file)
# validation:          a boolean indicating if validation will be performed
#                      on test data. If TRUE validation will be performed on test data. 
#                      If FALSE validation will not be performed on test data. Default is FALSE.
# test.train.ratio:    a floating point number between 0 and 1 representing
#                      the ratio of samples to keep as validation dataset. For example, a
#                      test.train.ratio = 0.25 splits 25 percent of the data into a validation dataset,
#                      meaning 75 percent of the data will be kept as the training dataset.
# num.trees.init:      an integer specifying number of trees to use for the first
#                      forest in the feature selection process. Default is 5000.
# num.trees.iterat:    an integer specifying number of trees to use for
#                      all additional forests in the feature selection process. Default is 2000.
# The function returns a list of four elements:
#                      1) a varSelRF object containing results of variable
#                      selection using random forest, 2) a randomForest object containing random forest model
#                      fitted to data, 3) a factor containing predictions of test data using fitted random
#                      forest model and 4) a confusionMatrix object containing confusion matrix of test data
#                      using fitted random forest model.
#                      If validation = FALSE, the last three elements in output will be NA.
# Code in function is copied from CAMPP2 available at https://github.com/ELELAB/CAMPP2
# which is the second version of CAMPP published at doi: 10.1002/1878-0261.12850
ForestFeatures <- function(seed = 123,
                           data,
                           group,
                           validation = FALSE,
                           test.train.ratio,
                           num.trees.init = 5000,
                           num.trees.iterat = 2000) {
  
  group <- as.factor(group)
  
  if (validation == TRUE) {
    
    ll <- list()
    llev <- levels(as.factor(group))
    
    # Generate a list of length number of groups
    # The first element in the list is the indexes of the samples belonging to group 1
    # The n-th element in the list is the indexes of the samples belonging to group n
    # where n is the number of groups
    for (idx in 1:length(llev)) {
      
      pos <- which(group == as.character(llev[idx]))
      ll[[idx]] <- pos
      
    }
    
    # Split randomly a subset of the samples into a validation and training set
    # A ratio of test.train.ratio is used for the splitting
    # These samples will be used as validation data set
    # It will contain samples from both groups
    set.seed(seed)
    samp <- unlist(lapply(ll,
                          function(x) sample(x,
                                             ceiling((length(x) * test.train.ratio)))))
    # Generate the validation data set
    data.val <- t(data[,samp])
    group.val <- group[samp]
    
    # Generate the training data set, i.e. the original data minus the validation data
    data.train <- t(data[,-samp])
    group.train <- group[-samp]
    
    # Carry out variable selection using random forest on training data and out-of-bag errors
    set.seed(seed)
    RFvars <- varSelRF(data.train,
                       group.train,
                       ntree = num.trees.init,
                       ntreeIterat = num.trees.iterat,
                       vars.drop.frac = 0.2)
    
    # Fit random forest algorithm on training data
    set.seed(seed)
    RFforest <- randomForest(data.train,
                             group.train,
                             ntree = num.trees.init,
                             ntreeIterat = num.trees.iterat,
                             vars.drop.frac = 0.2)
    
    # Use random forest to predict test data
    pred <- predict(RFforest,
                    newdata = data.val)
    
    # Create confusion matrix using predicted group classes and actual group classes on test data
    RFpred <- confusionMatrix(data = pred,
                              reference = group.val)
    
    # Combine results of variable selection, predictions, confusion matrix and random forest
    # into one list
    res <- list(RF.var.sel = RFvars,
                RF.model = RFforest,
                RF.pred = pred,
                RF.confusion.mat = RFpred)
    
  } else {
    
    # Transpose data
    data <- t(data)
    
    # Carry out variable selection using random forest and out-of-bag errors
    set.seed(seed)
    RFvars <- varSelRF(data,
                       group,
                       ntree = num.trees.init,
                       ntreeIterat = num.trees.iterat,
                       vars.drop.frac = 0.2)
    
    # Return results from variable selection using random forest
    res <- list(RFvars,
                NA,
                NA,
                NA)
  }
  
  return(res)
  
}

### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

# This function, RFApply, applies random forest on gene expression data.
# This function takes six arguments: 
# data:                 a matrix of (transformed and normalized) feature counts from
#                       "seq", "array", "ms" or "other" technology (with feature IDs as row names
#                       and sample IDs as columns). If available, batch corrected data should be used.
#                       group a factor specifying group for each sample (e.g. could be
#                       represented by a column from a metadata file).
# validation:           a boolean indicating if validation will be performed
#                       on test data. If TRUE validation will be performed on test data. 
#                       If FALSE, validation will not be performed on test data. Default is FALSE.
# test.train.ratio:     a floating point number between 0 and 1 representing
#                       the ratio of samples to keep as validation dataset. For example, a
#                       test.train.ratio = 0.25 splits 25 percent of the data into 
#                       a validation dataset, meaning 75 percent of the data will be 
#                       kept as the training dataset.
# num.trees.init:       an integer specifying number of trees to use for the first
#                       forest in the feature selection process. Default is 5000.
# num.trees.iterat:     an integer specifying number of trees to use for
#                       all additional forests in the feature selection process. Default is 2000.
# The function returns a list of seven elements:
#                       1) a list where each element is a character string
#                       containing selected variables from the random forest feature selection process from
#                       each seed run, 2) a list where each element is out-of-bag errors (of type numeric)
#                       from all random forest models used in the feature selection process for each seed run,
#                       3) a list where each element is the OOB error (of type numeric) from the best random
#                       forest model (i.e. the final selected model) of the feature selection process for each
#                       seed run, 4) a matrix containing average initial importance (before any variable deletion)
#                       for each feature across seed runs, 5) a list where each element is the OOB error (of type
#                       numeric) from fitted random forest model for each seed run (if a random forest was not
#                       fitted then this will be NA), 6) a matrix containing accuracy and 95 percent confidence
#                       interval for predictions of test data using fitted random forest model where each row
#                       in the matrix is a seed run, and 7) a vector of integers containing seeds used in the
#                       procedure.
# Code in function is copied from CAMPP2 available at https://github.com/ELELAB/CAMPP2
# which is the second version of CAMPP published at doi: 10.1002/1878-0261.12850
RFApply <- function(data,
                    group,
                    validation = FALSE,
                    test.train.ratio,
                    num.trees.init = 5000,
                    num.trees.iterat = 2000) {
  
  # Get randomly 10 seeds
  seeds <- sample(1:1000, 10)
  
  # Assign empty lists to store results
  RFvars <- list() # store selected variables from feature selection process
  RFvars.oob <- list() # store out-of-bag errors from feature selection process
  RFsel.vars.oob <- list() # store out-of-bag error from the best random forest model (i.e. the final selected model)
  RFimp <- list() # store initial importance of variables before any variable deletion from feature selection process
  RFoob.fit <- list() # store OOB error from fitted random forest model
  RFacc <- list() # store accuracy and 95% confidence interval for predictions of test data using fitted random forest model
  
  # For each seed
  for (idx in 1:length(seeds)) {
    
    # Run random forest
    RF <- ForestFeatures(seed = seeds[[idx]],
                         data = data,
                         group = group,
                         validation = validation,
                         test.train.ratio = test.train.ratio,
                         num.trees.init = num.trees.init,
                         num.trees.iterat = num.trees.iterat)
    
    ## Extract different results from the random forest variable selection process
    
    # Obtain selected variables found from running feature selection using random forest
    RFvars[[idx]] <- RF[[1]]$selected.vars
    
    # Obtain OOB error from all random forest models used in feature selection process
    RFvars.oob[[idx]] <- RF[[1]]$selec.history$OOB
    
    # Obtain OOB error from the best random forest model (i.e. the final selected model)
    tmp <- which(RF[[1]]$selec.history$Number.Variables ==
                   RF[[1]]$best.model.nvars) # Find index of the RF that is the final selected model
    RFsel.vars.oob[[idx]] <- RF[[1]]$selec.history$OOB[tmp] # Get OOB error of final selected RF model
    
    # Obtain initial importance of features before any variable deletion from feature selection process
    RFimp[[idx]] <- RF[[1]]$initialImportances
    
    
    ## Extract results from the process of fitting a random forest model to the data
    
    # If a random forest model was fitted
    if (class(RF[[2]]) == "randomForest") {
      
      # Obtain OOB error from fitted random forest model
      RFoob.fit[[idx]] <- tail(RF[[2]]$err.rate, n = 1)[,1]
      
    } else {
      
      RFoob.fit[[idx]] <- NA
      
    }
    
    # If predictions were performed on test data
    if (class(RF[[4]]) == "confusionMatrix") {
      
      # Obtain accuracy and 95% confidence interval for predictions of test data using fitted random forest model
      RFacc[[idx]] <- RF[[4]]$overall[c(1, 3, 4)]
      
    } else {
      
      RFacc[[idx]] <- NA
      
    }
    
    # Print message stating how many seed runs are completed
    print(paste0(idx, " out of ", length(seeds), " runs completed."))
    
  }
  
  ## Calculate average importance for each feature across seed runs
  
  # Initialize all.df to be importance of all features from first seed run
  all.df <- RFimp[[1]]
  
  # For each additional seed run
  for (idx in 2:length(RFimp)) {
    
    # Sum the importance of each feature across the seed runs
    all.df <- all.df + RFimp[[idx]]
    
  }
  
  # Divide the summed importance of each feature with the number of seed runs to get average
  RFimp <- all.df / length(RFimp)
  
  # Combine all results into one list
  RFres <- list(RFvars,
                RFvars.oob,
                RFsel.vars.oob,
                RFimp,
                RFoob.fit,
                do.call(rbind,
                        RFacc),
                seeds)
  
  # Assign names to elements in list
  names(RFres) <- c("Sel.vars",
                    "Var.sel.oob",
                    "Sel.rf.oob",
                    "Mean.importance.features",
                    "oob.rf.model",
                    "accuracy.rf.model",
                    "seeds")
  
  return(RFres)
  
}

### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

# This function, RunRF, applies random forest on gene expression data.
# This function takes six arguments: 
# data:                 a matrix of (transformed and normalized) feature counts from
#                       "seq", "array", "ms" or "other" technology (with feature IDs as row names
#                       and sample IDs as columns). If available, batch corrected data should be used.
#                       group a factor specifying group for each sample (e.g. could be
#                       represented by a column from a metadata file).
# group:                a factor specifying group for each sample (e.g. could be
#                       represented by a column from a metadata file)
# split.size:           an integer specifying the minimum number of samples that
#                       the groups must contain in order to carry out random forest classification
#                       and validation
# test.train.ratio:     a floating point number between 0 and 1 representing
#                       the ratio of samples to keep as validation dataset. For example, a
#                       test.train.ratio = 0.25 splits 25 percent of the data into 
#                       a validation dataset, meaning 75 percent of the data will be 
#                       kept as the training dataset.
# num.trees.init:       an integer specifying number of trees to use for the first
#                       forest in the feature selection process. Default is 5000.
# num.trees.iterat:     an integer specifying number of trees to use for
#                       all additional forests in the feature selection process. Default is 2000.
# The function returns a list of two elements:
#                       1) The first element is a list containing the output from the RFApply function
#                       2) The second element is a character vector containing the intersection of
#                       selected variables (e.g., genes) from the feature selection process
# Code in function is copied from CAMPP2 available at https://github.com/ELELAB/CAMPP2
# which is the second version of CAMPP published at doi: 10.1002/1878-0261.12850
RunRF <- function(data,
                  group,
                  split.size,
                  test.train.ratio,
                  num.trees.init = 5000,
                  num.trees.iterat = 2000) {
  
  ## Run random forest
  
  # Length of each group
  group.size <- as.numeric(table(group))
  test.train <- unique(group.size >= split.size)
  
  if (FALSE %in% test.train) {
    
    # Run random forest variable selection
    print("Validation is not going to be done due to the low number (<split.size) of samples in at least one of the sample groups. Only variable selection will be performed.")
    RF.results <- RFApply(data = data, group = group, validation = FALSE, test.train.ratio = test.train.ratio,
                          num.trees.init = num.trees.init, num.trees.iterat = num.trees.iterat)
    
  } else {
    
    # Run random forest classification and variable selection
    RF.results <- RFApply(data = data, group = group, validation = TRUE, test.train.ratio = test.train.ratio,
                          num.trees.init = num.trees.init, num.trees.iterat = num.trees.iterat)
    
  }
  
  # Intersect features selected in each seed run
  intersect.sel.vars <- Reduce(intersect, RF.results$Sel.vars)
  
  return(list("RFResults" = RF.results, "VarsSelect" = intersect.sel.vars))
  
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








