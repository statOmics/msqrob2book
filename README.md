# msqrob2book

The repository provides the source code for generating the book: 

"Statistical analysis of mass spectrometry-based proteomics data".

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

2. Compile each chapter.

```r
knitr::knit("CHAPTER.Rmd"", output = "build/CHAPTER.md", quiet = FALSE)
``` 

3. Assemble the book as HTML.

```
bookdown::render_book("build/index.md", output_dir = "docs")
```

This process is automatised using `make` in a (`bash`) shell. Credits
to https://github.com/rstudio/bookdown/issues/1262:

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

> Vandenbulcke S, Vanderaa C, Crook O, Martens L, Clement L.
msqrob2TMT: robust linear mixed models for inferring differential
abundant proteins in labelled experiments with arbitrarily complex
design. bioRxiv. Published online March 29, 2024:2024.03.29.587218.
doi:10.1101/2024.03.29.587218

**TODO**: populate the citation section and 
