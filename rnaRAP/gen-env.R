
library(rix)

rix(
  date = "2025-10-14",
  r_pkgs = c(
    "DESeq2",
    "pasilla",
    "devtools",
    "quarto"
  ),
  system_pkgs = NULL,
  git_pkgs = NULL,
  ide = "none",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)
