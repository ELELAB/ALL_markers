## Expression correlation

This directory contains script and results from the expression correlation analysis on 
our available T-ALL/B-ALL pediatric cohort, performed according to what described in 
Supplementary Text S3. It does not contain raw or input expression data as those
cannot be released.

## Requirements

The analysis has been performed on R 4.4.0 using the following R packages:

  - CEMiTool 1.28.0
  - enrichr 3.2
  - annotable 0.2.0
  - tidyverse 2.0.0

## Performing the analysis

The analysis has been performed running the included R script, `run.R`:

```
Rscript run.R
```

## Results

This directory contains results of the analysis, with modules_ALL_B.csv and modules_ALL_T.csv
containing the expression correlation modules identified by CEMiTool and the different csv,
pdf and png file containing the results of the enrichment analysis for each gene module

