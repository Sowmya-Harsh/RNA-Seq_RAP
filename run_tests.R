#!/usr/bin/env Rscript
# run_tests.R: Test rnaRAP core functions (Chapter 5)

library(testthat)
library(pasilla)
source("analysis_core.R")

cat(" Testing rnaRAP...\n\n")

test_that("load_example_data works", {
  data <- load_example_data()
  expect_true(is.list(data))
  expect_true("cts" %in% names(data))
  expect_true("coldata" %in% names(data))
  expect_true(ncol(data$cts) > 0)
})

test_that("top_genes_table works", {
  mock_res <- data.frame(
    row.names = paste0("gene_", 1:20),
    baseMean = runif(20, 50, 500),
    log2FoldChange = rnorm(20),
    padj = runif(20)
  )
  mock_res$padj[1:5] <- 0.01
  result <- top_genes_table(mock_res, n = 5, outdir = tempdir())
  expect_equal(nrow(result), 5)
  expect_true(all(result$padj < 0.05))
})

cat(" ALL TESTS PASS!\n")
cat("Coverage: load_example_data() + top_genes_table()\n")
