### This script contains functions used for the batch effect correction 
### The batch effect correction is done by COMBAT

### ---------------------------------------------------------------------------------------------------------

### This function, TCGAbatch_Correction_updated, does batch correction and is updated
### from the original function in TCGAbiolinks, TCGAbatch_Correction. The update is the
### inclusion of a design matrix in the call to ComBat when UnpublishedData = TRUE. 
TCGAbatch_Correction_updated <- function(tabDF, batch.factor = NULL, adjustment = NULL, ClinicalDF = data.frame(), 
                                         UnpublishedData = FALSE, AnnotationDF = data.frame(), 
                                         design.matrix) 
{
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("sva is needed. Please install it.", call. = FALSE)
  }
  if (UnpublishedData == TRUE) {
    batch.factor <- as.factor(AnnotationDF$Batch)
    batch_corr <- sva::ComBat(dat = tabDF, batch = batch.factor, mod = design.matrix,
                              par.prior = TRUE, prior.plots = TRUE)
  }
  if (UnpublishedData == FALSE) {
    if (length(batch.factor) == 0 & length(adjustment) == 
        0) 
      message("batch correction will be skipped")
    else if (batch.factor %in% adjustment) {
      stop(paste0("Cannot adjust and correct for the same factor|"))
    }
    my_IDs <- get_IDs(tabDF)
    if (length(batch.factor) > 0 || length(adjustment) > 
        0) 
      if ((nrow(ClinicalDF) > 0 & batch.factor == "Year") || 
          ("Year" %in% adjustment == TRUE & nrow(ClinicalDF) > 
           0)) {
        names(ClinicalDF)[names(ClinicalDF) == "bcr_patient_barcode"] <- "patient"
        ClinicalDF$age_at_diag_year <- floor(ClinicalDF$age_at_diagnosis/365)
        ClinicalDF$diag_year <- ClinicalDF$age_at_diag_year + 
          ClinicalDF$year_of_birth
        diag_yearDF <- ClinicalDF[, c("patient", "diag_year")]
        Year <- merge(my_IDs, diag_yearDF, by = "patient")
        Year <- Year$diag_year
        Year <- as.factor(Year)
      }
    else if (nrow(ClinicalDF) == 0 & batch.factor == 
             "Year") {
      stop("Cannot extract Year data. Clinical data was not provided")
    }
    Plate <- as.factor(my_IDs$plate)
    Condition <- as.factor(my_IDs$condition)
    TSS <- as.factor(my_IDs$tss)
    Portion <- as.factor(my_IDs$portion)
    Sequencing.Center <- as.factor(my_IDs$center)
    design.matrix <- model.matrix(~Condition)
    design.mod.combat <- model.matrix(~Condition)
    options <- c("Plate", "TSS", "Year", "Portion", "Sequencing Center")
    if (length(batch.factor) > 1) 
      stop("Combat can only correct for one batch variable. Provide one batch factor")
    if (batch.factor %in% options == FALSE) 
      stop(paste0(o, " is not a valid batch correction factor"))
    for (o in adjustment) {
      if (o %in% options == FALSE) 
        stop(paste0(o, " is not a valid adjustment factor"))
    }
    adjustment.data <- c()
    for (a in adjustment) {
      if (a == "Sequencing Center") 
        a <- Sequencing.Center
      adjustment.data <- cbind(eval(parse(text = a)), 
                               adjustment.data)
    }
    if (batch.factor == "Sequencing Center") 
      batch.factor <- Sequencing.Center
    batchCombat <- eval(parse(text = batch.factor))
    if (length(adjustment) > 0) {
      adjustment.formula <- paste(adjustment, collapse = "+")
      adjustment.formula <- paste0("+", adjustment.formula)
      adjustment.formula <- paste0("~Condition", adjustment.formula)
      print(adjustment.formula)
      model <- data.frame(batchCombat, row.names = colnames(tabDF))
      design.mod.combat <- model.matrix(eval(parse(text = adjustment.formula)), 
                                        data = model)
    }
    print(unique(batchCombat))
    batch_corr <- sva::ComBat(dat = tabDF, batch = batchCombat, 
                              mod = design.mod.combat, par.prior = TRUE, prior.plots = TRUE)
  }
  return(batch_corr)
}

### ---------------------------------------------------------------------------------------------------------


### ---------------------------------------------------------------------------------------------------------

### This function, batch correction, corrects gene expression data for batch effects using COMBAT. 
### The function takes seven inputs:
### batchTable: A tibble containing information about samples including their batch.
###             The tibble must contain the batch for each sample as a column and 
###             the barcodes of the samples as a column. 
### barcodesCol: The column in batchTable indicating the barcodes of each sample. 
### batchFactor: The column in batchTable indicating the batch of each sample. 
### filePath: The path and filename to where plots that COMBAT outputs should be saved
### plotWidth: The width of the plots that COMBAT outputs
### plotHeight: The height of the plots that COMBAT outputs
### expData: A gene expression matrix containing genes as rows and
###          samples as columns. 
batch_correction <- function(batchTable, barcodesCol, batchFactor, subtype_group, filePath, plotWidth, plotHeight, expData) {
  
  # Create data frame containing sample barcodes and corresponding batch information 
  annotation_counts <- batchTable %>%
    select({{barcodesCol}}, 
           {{batchFactor}}) %>%
    dplyr::rename(Samples = {{barcodesCol}}, 
                  Batch = {{batchFactor}}) %>%
    mutate(row_names = Samples) %>%
    column_to_rownames(var = "row_names")
  
  # Create design matrix
  subtype_design <- model.matrix(~subtype_group)
  
  # pdf to save plot that combat outputs
  pdf(file = filePath, 
      width = plotWidth,
      height = plotHeight)
  
  # Correct gene counts according to batch factor using COMBAT
  corrected_counts <- TCGAbatch_Correction_updated(tabDF = expData, 
                                                   UnpublishedData = TRUE, 
                                                   AnnotationDF = annotation_counts, 
                                                   design.matrix = subtype_design)
  
  # End with dev.off()
  dev.off()
  
  return(corrected_counts)
  
}

### ---------------------------------------------------------------------------------------------------------


