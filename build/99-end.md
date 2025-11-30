# Useful links and session information{#sec-end}



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

**TODO**

### Heart data set

**TODO**

### Mouse diet data set

**TODO**

## License

<a rel="license"
href="http://creativecommons.org/licenses/by-sa/4.0/"><img
alt="Creative Commons Licence" style="border-width:0"
src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br
/>This material is licensed under a <a rel="license"
href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons
Attribution-ShareAlike 4.0 International License</a>. You are free to
**share** (copy and redistribute the material in any medium or format)
and **adapt** (remix, transform, and build upon the material) for any
purpose, even commercially, as long as you give appropriate credit and
distribute your contributions under the same license as the original.

## Citation

Please cite this book as:

**TODO**: add citation once published

Please cite the `msqrob2` package as: 

> Goeminne L, Gevaert K, Clement L (2016). "Peptide-level Robust 
  Ridge Regression Improves Estimation, Sensitivity, and Specificity
  in Data-dependent Quantitative Label-free Shotgun Proteomics."
  _Molecular & Cellular Proteomics_, *15*(2), 657-668.
  doi:[10.1074/mcp.m115.055897](https://doi.org/10.1074/mcp.m115.055897).

If you opt for a summarisation-based workflow, you can also cite: 

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
## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] bookdown_0.45               tidyr_1.3.1                
##  [3] scater_1.38.0               scuttle_1.19.0             
##  [5] SingleCellExperiment_1.31.0 patchwork_1.3.2            
##  [7] MsDataHub_1.10.0            impute_1.84.0              
##  [9] ggrepel_0.9.6               ggplot2_4.0.0              
## [11] ExploreModelMatrix_1.22.0   dplyr_1.1.4                
## [13] ComplexHeatmap_2.26.0       BiocFileCache_3.0.0        
## [15] dbplyr_2.5.0                BiocParallel_1.44.0        
## [17] msqrob2book_0.0.99          msqrob2_1.17.1             
## [19] QFeatures_1.19.3            MultiAssayExperiment_1.35.1
## [21] SummarizedExperiment_1.39.0 Biobase_2.69.0             
## [23] GenomicRanges_1.61.0        GenomeInfoDb_1.45.3        
## [25] IRanges_2.43.0              S4Vectors_0.47.0           
## [27] BiocGenerics_0.55.0         generics_0.1.3             
## [29] MatrixGenerics_1.21.0       matrixStats_1.5.0          
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.5.2           later_1.4.2             filelock_1.0.3         
##   [4] tibble_3.2.1            lifecycle_1.0.4         httr2_1.1.2            
##   [7] Rdpack_2.6.4            doParallel_1.0.17       rprojroot_2.0.4        
##  [10] lattice_0.22-7          MASS_7.3-65             magrittr_2.0.3         
##  [13] limma_3.65.0            yaml_2.3.10             remotes_2.5.0          
##  [16] httpuv_1.6.16           sessioninfo_1.2.3       pkgbuild_1.4.7         
##  [19] cowplot_1.1.3           MsCoreUtils_1.21.0      DBI_1.2.3              
##  [22] minqa_1.2.8             RColorBrewer_1.1-3      abind_1.4-8            
##  [25] pkgload_1.4.0           purrr_1.0.4             AnnotationFilter_1.33.0
##  [28] rappdirs_0.3.3          circlize_0.4.16         irlba_2.3.5.1          
##  [31] codetools_0.2-20        DelayedArray_0.35.1     DT_0.33                
##  [34] tidyselect_1.2.1        shape_1.4.6.1           UCSC.utils_1.5.0       
##  [37] farver_2.1.2            lme4_1.1-37             ScaledMatrix_1.17.0    
##  [40] viridis_0.6.5           jsonlite_2.0.0          GetoptLong_1.0.5       
##  [43] BiocNeighbors_2.3.0     ellipsis_0.3.2          iterators_1.0.14       
##  [46] foreach_1.5.2           tools_4.5.2             Rcpp_1.0.14            
##  [49] glue_1.8.0              gridExtra_2.3           SparseArray_1.9.0      
##  [52] BiocBaseUtils_1.11.0    xfun_0.52               usethis_3.1.0          
##  [55] shinydashboard_0.7.3    withr_3.0.2             BiocManager_1.30.25    
##  [58] fastmap_1.2.0           boot_1.3-31             shinyjs_2.1.0          
##  [61] digest_0.6.37           rsvd_1.0.5              R6_2.6.1               
##  [64] mime_0.13               colorspace_2.1-1        RSQLite_2.3.11         
##  [67] httr_1.4.7              htmlwidgets_1.6.4       S4Arrays_1.9.0         
##  [70] pkgconfig_2.0.3         gtable_0.3.6            blob_1.2.4             
##  [73] S7_0.2.0                XVector_0.49.0          htmltools_0.5.8.1      
##  [76] profvis_0.4.0           ProtGenerics_1.41.0     rintrojs_0.3.4         
##  [79] clue_0.3-66             scales_1.4.0            png_0.1-8              
##  [82] reformulas_0.4.1        knitr_1.50              rstudioapi_0.17.1      
##  [85] reshape2_1.4.4          rjson_0.2.23            nlme_3.1-168           
##  [88] curl_6.2.2              nloptr_2.2.1            cachem_1.1.0           
##  [91] GlobalOptions_0.1.2     stringr_1.5.1           BiocVersion_3.22.0     
##  [94] parallel_4.5.2          miniUI_0.1.2            vipor_0.4.7            
##  [97] AnnotationDbi_1.71.0    desc_1.4.3              pillar_1.10.2          
## [100] vctrs_0.6.5             urlchecker_1.0.1        promises_1.3.2         
## [103] BiocSingular_1.25.0     beachmat_2.25.0         xtable_1.8-4           
## [106] cluster_2.1.8.1         beeswarm_0.4.0          evaluate_1.0.3         
## [109] cli_3.6.5               compiler_4.5.2          rlang_1.1.6            
## [112] crayon_1.5.3            plyr_1.8.9              fs_1.6.6               
## [115] ggbeeswarm_0.7.2        stringi_1.8.7           viridisLite_0.4.2      
## [118] Biostrings_2.77.0       lazyeval_0.2.2          devtools_2.4.5         
## [121] Matrix_1.7-4            ExperimentHub_2.99.0    bit64_4.6.0-1          
## [124] KEGGREST_1.49.0         statmod_1.5.0           shiny_1.10.0           
## [127] AnnotationHub_3.99.0    rbibutils_2.3           igraph_2.1.4           
## [130] memoise_2.0.1           bit_4.6.0
```
