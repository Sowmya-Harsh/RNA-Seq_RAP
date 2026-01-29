library(rix)

latest_date <- "2026-01-29"

rix(
  date = latest_date,
  r_pkgs = c("BiocManager", "dplyr", "ggplot2", "DESeq2", "pasilla", "knitr", 
             "testthat", "apeglm", "KernSmooth", "tibble", "quarto",
            "pheatmap", "RColorBrewer", "gridExtra", "tidyr"),
  system_pkgs = c("pandoc"),
  ide = "none",
  project_path = ".",
  overwrite = TRUE
)
