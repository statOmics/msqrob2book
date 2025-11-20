# msqrob2book

The repository provides the source code for generating the book: 

*"Statistical analysis of mass spectrometry-based proteomics data".*

The book is focused around the msqrob2 software.

<div>
<img src="figs/msqrob2.png" alt="msqrob2 sticker" width="300px" />
<p style="font-size: 10px; color: grey;">Sticker license</p>
</div>

**TODO** add sticker license 

## Compiling the book

Each chapter of the book is compiled/knit separately. The resulting
`.md` file (one per chapter) is saved in the `build` subdirectory.
From the `.md` files, the book is rendered using `bookdown` into html
files in the `docs` directory.

To compile the book, follow the instructions and the run the code in 
R:

1. Install all dependencies.

```r
installed.packages("BiocManager")
BiocManager::install(version = "devel")
remotes::install_deps(dependencies = TRUE, repos = BiocManager::repositories())
BiocManager::install("lgatto/msmbstyle")
```

2. Compile the chapters you want (obviously rename `CHAPTER` with the
   chapter name you want to compile). This will create intermediate 
   markdown (`.md`) files in the `build/` folder. 

```r
knitr::knit("CHAPTER.Rmd"", output = "build/CHAPTER.md")
``` 

4. Copy additional files (only needed if changed) to the `build/` 
   directory: `refs.bib`, `_bookdown.yml`, and the full `figure/` 
   directory (contains the figures generated within the code chunks by
   `knitr`).

3. Assemble the book as HTML.

```
bookdown::render_book("build/index.md", output_dir = "docs")
```

This process is automatised using `make` which you can run in a 
(`bash`) shell. Credits to
https://github.com/rstudio/bookdown/issues/1262:

```
make all
```

The book can also be compiled through [GitHub
Actions](https://github.com/statOmics/msqrob2book/actions/workflows/publish_book.yaml)

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
