library(rix)
dates <- available_dates()
latest_date <- tail(dates, 1)

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