---
title: "Statistical analysis of mass spectrometry-based proteomics data"
subtitle: "A dive into the msqrob2 universe"
author: "Christophe Vanderaa, Stijn Vandenbulcke, Lieven Clement"
date: "2025-11-16"
output:
  msmbstyle::msmb_html_book:
    highlight: tango
    toc: TRUE
    toc_depth: 1
    split_by: chapter
    margin_references: TRUE
    css: style.css
bibliography: [refs.bib]
link-citations: yes
---

# Preamble



This book provides comprehensive hands-on tutorials on how to apply the
`msqrob2` software for the statistical analysis of mass spectrometry
(MS)-based proteomics data. It includes the latest improvements of the
software that enable statistical modelling for a wide panel of use
cases. The book will first introduce the general concepts behind
the statistical analysis of proteomics data and `msqrob2`. Further
chapters will demonstrate the application of `msqrob2` for assessing
different biological questions starting from datasets with different
experimental designs, acquisition strategies, instruments, and search
engines. The book aims to help proteomics researchers tailoring their
statistical analysis workflow to their specific datasets and research
questions.

## Why msqrob2?

MS-based proteomics experiments often imposes a complex correlation
structure among observations. Addressing this correlation is key for
correct statistical inference and reliable biomarker discovery. This
`msqrob2` book provides a set of (mixed) model-based workflows
dedicated to differential abundance analysis for label-free as well
as labeled MS-based proteomics data. The key features of `msqrob2`
workflows are:

1. Modularity: all core functions rely on the `QFeatures` class, a
   standardised data structure, meaning that output of a function can
   be fed as input to any other function. Hence, different functions
   are assembled as modular blocks into a complete data analysis
   workflows that can be easily adapted to the peculiarities of any
   MS-based proteomics data set. Therefore, the approach extends well
   beyond the use case presented in this chapter
2. Flexibility: the `msqrob2` modelling approach relies on the
   `lme4::lmer()` model specification syntax, meaning that any linear
   model can be specified. For fixed effects, this includes modelling
   categorical and numerical variables, as well as their interaction.
   Moreover, `msqrob2` can model both sample-specific and
   feature-specific (e.g. peptide or protein) covariates, which
   unlocks the inference to experiments with arbitrarily complex
   designs as well as to correct explicitly for feature-specific
   properties.
3. Performance: thanks to the inclusion of robust ridge regression, we
   demonstrated improved performance of `msqrob2` workflows upon the
   competing software [@Goeminne2016-tr;@Sticker2020-rl;@Vandenbulcke2025-sj].

## Outline

The book is divided in three parts.

### Concepts

This parts introduces the user to the concepts for . It explains the
theory behind proteomics data analysis and provides extensive
description of the code. While this part is conceptual, the concepts
are illustrated using a real spike-in study.

- Chapter \@ref(sec-basics) introduces the basic concepts for MS-based 
  proteomics analysis. We **recommend** user to first read this 
  chapter before reading any other chapter. 
- Chapter \@ref(sec-advanced) builds upon the previous chapter and 
  introduces more more concepts that will be used in later chapters
  that involve more complex designs and analyses.

### Benchmarking

This part illustrates how to benchmark data analysis workflows and
demonstrates how the guidelines presented in this book were derived.
The main sections of the chapters in this part are intended for
advanced users with R programming skills. However, the conclusions in
each chapter are more accessible and also intended for entry-level
users that want to understand how to apply the guidelines and
recommendations to their analyses.

- Chapter \@ref(sec-benchmarking) explains how to conduct a 
  benchmarking experiment to assess different data analysis workflows.
  As an example, the chapter compares the performance when starting
  from the different MaxQuant input files: the evidence file, the
  peptides file, and the protein-group file.
- Chapter **TODO** summarise Stijn's chapter.

### Use cases{#sec-overview_use_cases}

This part contains a set of chapter that illustrate the data analysis
for a range of experimental designs and technological setups.

- Chapter \@sec-francisella analyses the francisella use case: a 
  MaxQuant LFQ DDA dataset with technical replication.
- Chapter \@sec-heart_chapter analyses the heart use case: a MaxQuant
  LFQ DDA dataset with a more complex design.
- Chapter \@sec-mouse_diet analyses the mouse diet use case: a Skyline
  TMT DDA dataset.

**TODO** make table with dataset descriptions to let users easily find
their relevant use case.

## Targeted audience and assumed background

The course material is targeted to either proteomics practitioners or
data analysts/bioinformaticians that would like to learn how to
analyse proteomics data. 

A working knowledge of R (R syntax, commonly used functions, basic
data structures such as data frames, vectors, matrices, ... and their
manipulation) is required. Familiarity with MS or proteomics in
general is recommended, this would allow for a better understanding of
the modelling assumptions taken throughout this book. Familiarity with
other Bioconductor omics data classes and the tidyverse syntax is
useful.

We **highly recommend** reading the [quantitative proteomics
chapter](https://rformassspectrometry.github.io/book/sec-quant.html)
of the R for mass spectrometry book.

## Setup

To install all the necessary package, please use the latest release of
R and execute:


``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
   "BiocParallel",
   "BiocFileCache",
   "ComplexHeatmap",
   "dplyr",
   "ExploreModelMatrix",
   "ggpattern",
   "ggplot2",
   "ggrepel",
   "impute",
   "MsDataHub",
   "patchwork",
   "scater",
   "tidyr",
   "bookdown"
))
```

All software versions used to generate this document are recorded at
the end of the book in chapter \@ref(sec-end).

## Citation

If you need to cite this book, please use the following reference:

**TODO** add citation when available on Zenodo (and paper)

## Acknowledgments

We thank the [R for Mass
Spectrometry](https://www.rformassspectrometry.org/) initiative: 

- For developing the `QFeatures` package on which `msqrob2` depends
- For openly sharing their
  [book](https://rformassspectrometry.github.io/book), which we used
  as a template.

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
