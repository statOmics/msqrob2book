# --- Directories ---
SRC_DIR := .
BUILD_DIR := build
DOCS_DIR := docs

# --- Files ---
RMDS := $(wildcard $(SRC_DIR)/*.Rmd)
MDS := $(patsubst $(SRC_DIR)/%.Rmd, $(BUILD_DIR)/%.md, $(RMDS))
HTML := $(DOCS_DIR)/index.html

# --- Default target ---
all: $(HTML)

# --- Rule: Knit Rmd -> build/md ---
$(BUILD_DIR)/%.md: $(SRC_DIR)/%.Rmd | $(BUILD_DIR)
	Rscript -e "knitr::knit('$<', output='$@')"

# --- Ensure build directory exists ---
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# --- Rule: Render bookdown HTML from built md files ---
$(HTML): $(MDS) _bookdown.yml refs.bib | $(DOCS_DIR)
	cp refs.bib $(BUILD_DIR)/refs.bib
	cp _bookdown.yml $(BUILD_DIR)/_bookdown.yml
	cp figure/* $(BUILD_DIR)/figure
	Rscript -e "setwd('$(BUILD_DIR)');bookdown::render_book('index.md', output_dir = '../$(DOCS_DIR)')"

# --- Ensure docs directory exists ---
$(DOCS_DIR):
	mkdir -p $(DOCS_DIR)

# --- Clean up unecessary files ---
clean:
	build/msqrob2book.rds