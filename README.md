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

This process is automatised using `make` in a (`bash`) shell:

```
make all
```

The book can also be compiled through [GitHub
Actions](https://github.com/statOmics/msqrob2book/actions/workflows/publish_book.yaml)
