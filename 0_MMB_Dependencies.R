# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# Define required CRAN packages
cran_packages <- c(
  "readr", "dplyr", "tidyr", "tibble", "ggsurvfit", "survival", "Boruta",
  "DataExplorer", "broom", "knitr", "ggfortify", "BiocManager"
)

# Install missing CRAN packages
missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

# Load CRAN packages
invisible(lapply(cran_packages, require, character.only = TRUE))

# Install and load Bioconductor package
if (!"limma" %in% rownames(installed.packages())) {
  BiocManager::install("limma")
}
library(limma)
