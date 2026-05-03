# ============================================================
# 00_install_packages.R
# Install all required packages for the MR pipeline
# Run this script once before the analysis
# ============================================================

# CRAN packages
cran_pkgs <- c(
  "data.table",   # fast file I/O
  "dplyr",        # data manipulation
  "tidyr",        # reshaping
  "ggplot2",      # base plotting
  "ggrepel",      # non-overlapping labels
  "scales",       # axis formatting
  "patchwork",    # multi-panel figures
  "viridis",      # colour scales
  "RColorBrewer", # colour palettes
  "stringr",      # string manipulation
  "forcats",      # factor handling
  "broom",        # tidy model outputs
  "jsonlite",     # JSON I/O
  "R.utils"       # gzip support
)

install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("qvalue"), ask = FALSE)

# TwoSampleMR — from MRCIEU GitHub
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/MRInstruments")  # optional: access IEU OpenGWAS

cat("\n✅  All packages installed successfully.\n")
