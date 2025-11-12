# Usefull links and session information{#sec-end}



## Data sets

We refer here the data sources used in the book:

### E. Coli LFQ spike-in data set

Original study: Shen X, Shen S, Li J, Hu Q, Nie L, Tu C, et al. (2018)
Ionstar enables high-precision, low-missing-data proteomics
quantification in large bio- logical cohorts. Proc.  Natl.  Acad.
Sci. U.S.A. 115, E4767–E4776

Reanalysis study: Sticker A, Goeminne L, Martens L, Clement L. Robust
Summarization and Inference in Proteome-wide Label-free
Quantification. Mol Cell Proteomics. 2020;19(7):1209-1219.

Link to data: https://github.com/statOmics/MSqRobSumPaper/raw/refs/heads/master/spikein/data/maxquant/peptides.zip

Link to data in archive: **TODO** add link to Zenodo

Used in chapters \@ref(sec-basics), and \@ref(sec-benchmarking).

### TMT spike-in data set

Original study: Huang T, Choi M, Tzouros M, Golling S, Pandya NJ,
Banfai B, et al. MSstatsTMT: Statistical Detection of Differentially
Abundant Proteins in Experiments with Isobaric Labeling and Multiple
Mixtures. Mol Cell Proteomics. 2020;19(10):1706-1723.

Reanalysis study: Vandenbulcke S, Vanderaa C, Crook O, Martens L,
Clement L. Msqrob2TMT: Robust linear mixed models for inferring
differential abundant proteins in labeled experiments with arbitrarily
complex design. Mol Cell Proteomics. 2025;24(7):101002.
  
Data source: MassIVE repository (RMSV000000265)

Link to data from archive: https://zenodo.org/records/14767905

Used in chapter \@ref(sec-advanced).

### Francisella data set



### Heart data set


### Mouse diet data set


## License

This book is distributed under a
[Artistic-2.0](https://opensource.org/license/artistic-2-0) license.

**TODO**: add logo and explicit how people can access and modify this work

## Citation

Please cite this book as:

**TODO**: add citation once published

Please cite the `msqrob2` package as: 

> Goeminne L, Gevaert K, Clement L (2016). "Peptide-level Robust 
  Ridge Regression Improves Estimation, Sensitivity, and Specificity
  in Data-dependent Quantitative Label-free Shotgun Proteomics."
  _Molecular & Cellular Proteomics_, *15*(2), 657-668.
  doi:[10.1074/mcp.m115.055897](https://doi.org/10.1074/mcp.m115.055897).

If you opt for a summarised based workflow, you can also cite: 

> Sticker A, Goeminne L, Martens L, Clement L (2020). "Robust 
  Summarization and Inference in Proteome-wide Label-free
  Quantification." _Molecular & Cellular Proteomics_, *19*(7),
  1209-1219.
  doi:[10.1074/mcp.ra119.001624](https://doi.org/10.1074/mcp.ra119.001624).

If you use TMT-based workflows, please cite

> Vandenbulcke S, Vanderaa C, Crook O, Martens L, Clement L.
msqrob2TMT: robust linear mixed models for inferring differential
abundant proteins in labelled experiments with arbitrarily complex
design. bioRxiv. Published online March 29, 2024:2024.03.29.587218.
doi:10.1101/2024.03.29.587218

**TODO**: populate the citation section and 

## Session Info

The following packages have been used to generate this document.




``` r
sessionInfo()
```

```
## R version 4.5.2 (2025-10-31)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 24.04.3 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Europe/Brussels
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] msqrob2book_0.0.99          msqrob2_1.17.1             
##  [3] QFeatures_1.19.3            MultiAssayExperiment_1.35.1
##  [5] SummarizedExperiment_1.39.0 Biobase_2.69.0             
##  [7] GenomicRanges_1.61.0        GenomeInfoDb_1.45.3        
##  [9] IRanges_2.43.0              S4Vectors_0.47.0           
## [11] BiocGenerics_0.55.0         generics_0.1.3             
## [13] MatrixGenerics_1.21.0       matrixStats_1.5.0          
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.2.1        dplyr_1.1.4             fastmap_1.2.0          
##  [4] lazyeval_0.2.2          promises_1.3.2          digest_0.6.37          
##  [7] mime_0.13               lifecycle_1.0.4         cluster_2.1.8.1        
## [10] ProtGenerics_1.41.0     ellipsis_0.3.2          statmod_1.5.0          
## [13] magrittr_2.0.3          compiler_4.5.2          rlang_1.1.6            
## [16] tools_4.5.2             igraph_2.1.4            knitr_1.50             
## [19] S4Arrays_1.9.0          htmlwidgets_1.6.4       pkgbuild_1.4.7         
## [22] DelayedArray_0.35.1     plyr_1.8.9              pkgload_1.4.0          
## [25] abind_1.4-8             BiocParallel_1.44.0     miniUI_0.1.2           
## [28] withr_3.0.2             purrr_1.0.4             desc_1.4.3             
## [31] grid_4.5.2              urlchecker_1.0.1        profvis_0.4.0          
## [34] xtable_1.8-4            MASS_7.3-65             cli_3.6.5              
## [37] crayon_1.5.3            reformulas_0.4.1        remotes_2.5.0          
## [40] rstudioapi_0.17.1       httr_1.4.7              reshape2_1.4.4         
## [43] sessioninfo_1.2.3       BiocBaseUtils_1.11.0    minqa_1.2.8            
## [46] cachem_1.1.0            stringr_1.5.1           splines_4.5.2          
## [49] parallel_4.5.2          AnnotationFilter_1.33.0 XVector_0.49.0         
## [52] vctrs_0.6.5             devtools_2.4.5          boot_1.3-31            
## [55] Matrix_1.7-4            jsonlite_2.0.0          clue_0.3-66            
## [58] limma_3.65.0            tidyr_1.3.1             glue_1.8.0             
## [61] nloptr_2.2.1            codetools_0.2-20        stringi_1.8.7          
## [64] later_1.4.2             UCSC.utils_1.5.0        lme4_1.1-37            
## [67] tibble_3.2.1            pillar_1.10.2           htmltools_0.5.8.1      
## [70] R6_2.6.1                Rdpack_2.6.4            rprojroot_2.0.4        
## [73] evaluate_1.0.3          shiny_1.10.0            lattice_0.22-7         
## [76] rbibutils_2.3           memoise_2.0.1           httpuv_1.6.16          
## [79] Rcpp_1.0.14             SparseArray_1.9.0       nlme_3.1-168           
## [82] xfun_0.52               fs_1.6.6                MsCoreUtils_1.21.0     
## [85] usethis_3.1.0           pkgconfig_2.0.3
```

