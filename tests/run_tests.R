#!/usr/bin/env Rscript
# Turn-key test runner for the CNA module.
#
# Usage (from the package root):
#   JRS=/Applications/jamovi.app/Contents/Frameworks/R.framework/Versions/4.5-arm64/Resources/bin/Rscript
#   $JRS tests/run_tests.R
#
# Resolves the build library (so CNA/Nestimate/cograph load) and runs the whole
# tests/testthat/ suite, reporting a non-zero exit on any failure.

# --- locate this script's directory so paths work regardless of cwd ----------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(file_arg))
  normalizePath(file.path(dirname(file_arg), "..")) else normalizePath(".")

# --- put the freshly-built package library first ------------------------------
build_lib <- list.files(file.path(root, "build"),
                        pattern = "^R[0-9].*", full.names = TRUE)
build_lib <- build_lib[dir.exists(build_lib)]
if (length(build_lib)) .libPaths(c(build_lib[[1]], .libPaths()))

if (!requireNamespace("CNA", quietly = TRUE))
  stop("CNA not found on .libPaths(); run jmvtools::install('.') first so ",
       "build/R*/CNA exists.")

library(testthat)
# The Monte Carlo equivalence sweeps are marked skip_on_cran(). This is the
# local full run, so they run unless CNA_QUICK=true asks for the quick suite.
if (!nzchar(Sys.getenv("NOT_CRAN")) && !identical(Sys.getenv("CNA_QUICK"), "true"))
  Sys.setenv(NOT_CRAN = "true")
options(CNA_TEST_ROOT = root,            # tests resolve data/ relative to this
        testthat.progress.max_fails = Inf)  # report every failure, not the first 10

reporter <- testthat::SummaryReporter$new()
res <- testthat::test_dir(file.path(root, "tests", "testthat"),
                          reporter = reporter, stop_on_failure = FALSE)
df <- as.data.frame(res)
fail <- sum(df$failed) + sum(df$error)
cat(sprintf("\n==== CNA tests: %d passed, %d failed, %d skipped across %d test blocks ====\n",
            sum(df$passed), fail, sum(df$skipped), nrow(df)))
if (fail > 0) quit(status = 1)
