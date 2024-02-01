#quantileNormalization function
.quantileNormalization <-
  function(wd, distribution) {
    n <- nrow(wd)
    m <- ncol(wd)
    if(!missing(distribution))
      if(m != length(distribution))
        stop("The reference distribution has length
                 different from the col dimension of the data matrix.") else
                   distribution  <-  sort(distribution)
    
    o <- matrix(0, n, m)
    for(i in 1:n)
      o[i,] <- order(wd[i,])
    
    j <- 1
    tmp <- rep(0, n)
    
    while(j <= m) {
      if(missing(distribution)) {
        for(i in 1:n)
          tmp[i] <- wd[i,o[i,j]]
        value <- mean(tmp)
      } else value  <- distribution[j]
      for(i in 1:n)
        wd[i,o[i,j]] <- value
      j <- j + 1
    }
    return(wd)
  }

#TCGAanalyze_Normalization_revised function
TCGAanalyze_Normalization_revised <- function(
    tabDF,
    geneInfo,
    method = "geneLength"
) {
  
  #check_package("EDASeq")
  
  # Check if we have a SE, we need a gene expression matrix
  if (is(tabDF, "SummarizedExperiment")) tabDF <- assay(tabDF)
  
  geneInfo <- geneInfo[!is.na(geneInfo[, 1]), ]
  geneInfo <- as.data.frame(geneInfo)
  geneInfo$geneLength <- as.numeric(as.character(geneInfo$geneLength))
  geneInfo$gcContent <- as.numeric(as.character(geneInfo$gcContent))
  
  
  if(any(grepl("\\|",rownames(tabDF)))){
    tmp <- as.character(rownames(tabDF))
    geneNames <- stringr::str_split(tmp, "\\|",simplify = T)
    tmp <- which(geneNames[, 1] == "?")
    geneNames[tmp, 1] <- geneNames[tmp, 2]
    tmp <- table(geneNames[, 1])
    tmp <- which(geneNames[, 1] %in% names(tmp[which(tmp > 1)]))
    geneNames[tmp, 1] <- paste(geneNames[tmp, 1], geneNames[tmp, 2], sep = ".")
    tmp <- table(geneNames[, 1])
    rownames(tabDF) <- geneNames[, 1]
  } else if(any(grepl("ENSG",rownames(tabDF)))){
    rownames(tabDF) <- gsub("\\.[0-9]*","",rownames(tabDF))
  }
  
  
  if (method == "gcContent") {
    rawCounts <- tabDF
    commonGenes <-  intersect(rownames(geneInfo), rownames(rawCounts))
    
    geneInfo <- geneInfo[commonGenes, ]
    rawCounts <- rawCounts[commonGenes, ]
    
    timeEstimated <- format(ncol(tabDF) * nrow(tabDF) / 80000, digits = 2)
    message(
      messageEstimation <- paste(
        "I Need about ",
        timeEstimated,
        "seconds for this Complete Normalization Upper Quantile",
        " [Processing 80k elements /s]  "
      )
    )
    
    ffData  <- as.data.frame(geneInfo)
    rawCounts <- floor(rawCounts)
    
    message("Step 1 of 4: newSeqExpressionSet ...")
    tmp <- EDASeq::newSeqExpressionSet(as.matrix(rawCounts), featureData = ffData)
    
    #fData(tmp)[, "gcContent"] <- as.numeric(geneInfo[, "gcContent"])
    
    message("Step 2 of 4: withinLaneNormalization ...")
    tmp <- EDASeq::withinLaneNormalization(
      tmp, "gcContent",
      which = "upper",
      offset = TRUE
    )
    
    #Filtering step added for Inf values, not present in TCGAbiolinks TCGAanalyze_Normalization function
    message("Step 3 of 4: betweenLaneNormalization ...")
    if(any(is.na(EDASeq::normCounts(tmp))) || any(is.infinite(EDASeq::normCounts(tmp)))) {
      tmp <- tmp[rowSums(is.na(EDASeq::normCounts(tmp))) == 0,]
      tmp <- tmp[rowSums(is.infinite(EDASeq::normCounts(tmp))) == 0,]
    }

    tmp <- EDASeq::betweenLaneNormalization(
      tmp,
      which = "upper",
      offset = TRUE
    )
    normCounts <-  log(rawCounts[rownames(tmp),] + .1) + EDASeq::offst(tmp)
    normCounts <-  floor(exp(normCounts) - .1)
    
    message("Step 4 of 4: .quantileNormalization ...")
    tmp <- t(.quantileNormalization(t(normCounts)))
    tabDF_norm <- floor(tmp)
  }
  
  if (method == "geneLength") {
    tabDF <- tabDF[!duplicated(rownames(tabDF)), !duplicated(colnames(tabDF))]
    tabDF <- tabDF[rownames(tabDF) %in% rownames(geneInfo), ]
    #tabDF <- tabDF[rowMeans(tabDF) > 1,]
    tabDF <- tabDF[which(rowSums(tabDF == 0) < ncol(tabDF)),]
    tabDF <- as.matrix(tabDF)
    
    geneInfo <-  geneInfo[rownames(geneInfo) %in% rownames(tabDF),]
    geneInfo <- geneInfo[!duplicated(rownames(geneInfo)),]
    toKeep <- which(geneInfo[, "geneLength"] != 0)
    geneInfo <- geneInfo[toKeep,]
    tabDF <- tabDF[toKeep,]
    geneInfo <- as.data.frame(geneInfo)
    tabDF <- round(tabDF)
    commonGenes <- intersect(rownames(tabDF), rownames(geneInfo))
    
    tabDF <- tabDF[commonGenes, ]
    geneInfo <- geneInfo[commonGenes, ]
    
    timeEstimated <-  format(ncol(tabDF) * nrow(tabDF) / 80000, digits = 2)
    message(
      messageEstimation <- paste(
        "I Need about ",
        timeEstimated,
        "seconds for this Complete Normalization Upper Quantile",
        " [Processing 80k elements /s]  "
      )
    )
    
    message("Step 1 of 4: newSeqExpressionSet ...")
    
    tabDF_norm <- EDASeq::newSeqExpressionSet(
      tabDF,
      featureData = geneInfo
    )
    message("Step 2 of 4: withinLaneNormalization ...")
    
    tabDF_norm <-  EDASeq::withinLaneNormalization(
      tabDF_norm,
      "geneLength",
      which = "upper",
      offset = FALSE
    )
    
    message("Step 3 of 4: betweenLaneNormalization ...")
    if(any(is.na(EDASeq::normCounts(tabDF_norm)))) {
      tabDF_norm <- tabDF_norm[rowSums(is.na(EDASeq::normCounts(tabDF_norm))) == 0,]
    }
    
    tabDF_norm <- EDASeq::betweenLaneNormalization(
      tabDF_norm,
      which = "upper",
      offset = FALSE
    )
    
    message("Step 4 of 4: exprs ...")
    tabDF_norm <- EDASeq::counts(tabDF_norm)
  }

  #Filtering step for genes with only NA values removed from TCGAbiolinks TCGAanalyze_Normalization function


  return(tabDF_norm)
}

