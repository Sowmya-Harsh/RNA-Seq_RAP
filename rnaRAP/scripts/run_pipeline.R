#!/usr/bin/env Rscript

# scripts/run_pipeline.R

library(DESeq2)
library(rnaRAP)

message("1) Loading example data...")
d <- load_example_data()
cts <- d$cts
coldata <- d$coldata

message("2) Running DESeq2...")
res <- run_deseq2_condition(cts, coldata)

message("3) Saving results and MA plot...")
dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, file = "data_processed/deseq2_results.rds")

plot_ma_simple(res, out_path = "reports/ma_plot.png")

message("Done.")

