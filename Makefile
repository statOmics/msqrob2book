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
	Rscript -e "knitr::knit('$<', output='$@', quiet=FALSE)"

# --- Ensure build directory exists ---
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# --- Rule: Render bookdown HTML from built md files ---
$(HTML): $(MDS) _bookdown.yml | $(DOCS_DIR)
	Rscript -e "setwd('$(BUILD_DIR)');bookdown::render_book('index.md', output_dir = '../$(DOCS_DIR)')"

# --- Ensure docs directory exists ---
$(DOCS_DIR):
	mkdir -p $(DOCS_DIR)

# --- Clean up generated files ---
clean:
