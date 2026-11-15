# testthat entrypoint for the CNA jamovi module.
#
# These tests need the Nestimate + cograph engine, which lives in jamovi's
# bundled R, not system R. Run them with jamovi's Rscript and the build lib
# on the path (see tests/run_tests.R for a turn-key runner):
#
#   JRS=/Applications/jamovi.app/Contents/Frameworks/R.framework/Versions/<ver>/Resources/bin/Rscript
#   $JRS tests/run_tests.R
#
# `tests/run_tests.R` prepends build/R*/  to .libPaths so loadNamespace("CNA"),
# Nestimate and cograph all resolve.

library(testthat)
library(CNA)

test_check("CNA")
