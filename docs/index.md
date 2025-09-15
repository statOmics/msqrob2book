---
title: "Statistical analysis for proteomics data"
subtitle: "A dive into the msqrob2 universe"
author: "Christophe Vanderaa, Stijn Vandenbulcke, Lieven Clement"
date: "2025-09-05"
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



The aim of this book is to provide thorough documentation of the
`msqrob2` software for the statistical analysis of mass spectrometry
(MS)-based proteomics. It includes the latest improvements of the
software that enable statistical modelling for a wide panel of use
cases. The book will first introduce the general concepts behind
`msqrob2` and statistical analysis of proteomics data. Further
chapters will demonstrate the application of `msqrob2` for datasets
that were obtained using different acquisition strategies, different
instruments, and different search engines to answers different
biological questions. We hope that this book will help proteomics
researchers tailoring their statistical analysis workflow to their
specific dataset and research question.


## Targeted audience and assumed background

The course material is targeted to either proteomics practitioners or
data analysts/bioinformaticians that would like to learn how to
analyse proteomics data. Familiarity with MS or proteomics in general
is recommended, this would allow for a better understanding of the
modelling assumptions taken throughout this book.

A working knowledge of R (R syntax, commonly used functions, basic
data structures such as data frames, vectors, matrices, ... and their
manipulation) is required. Familiarity with other Bioconductor omics
data classes and the tidyverse syntax is useful, and we highly
recommend reading [Chapter
5](https://rformassspectrometry.github.io/book/sec-quant.html) of the
R for mass spectrometry book.

## Setup

To install all the necessary package, please use the latest release of
R and execute:

**TODO** install dependencies
**TODO** provide a Docker image


``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```

## Citation

If you need to cite this book, please use the following reference:

**TODO** add citation when available on Zenodo (and paper)

## Acknowledgments

We thank the [R for Mass
Spectrometry](https://www.rformassspectrometry.org/) initiative for
openly sharing their
[book](https://rformassspectrometry.github.io/book), which we used as
a template when writing this book.

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
