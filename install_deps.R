
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## Install dependencies required within the chapters
BiocManager::install(
    c(
        "BiocParallel",
        "BiocFileCache",
        "dplyr",
        "ExploreModelMatrix",
        "ggplot2",
        "ggrepel",
        "patchwork",
        "scater",
    ),
    ask = FALSE, udpate = TRUE
)

## Install specific package versions
BiocManager::install("cvanderaa/QFeatures", ref = "uniquePrecId")
BiocManager::install("statOmics/msqrob2")

## Install dependencies required to render the book
BiocManager::install(
    c(
        "lgatto/msmbstyle",
        "bookdown"
    ),
    ask = FALSE, udpate = TRUE
)