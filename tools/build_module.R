#!/usr/bin/env Rscript
# tools/build_module.R
#
# Builds CNA_<version>.jmo WITHOUT installing it into jamovi, optionally
# bundling engine packages from local source checkouts.
#
# jamovi's compiler (jmc --build) installs the DESCRIPTION Remotes (Nestimate,
# cograph) from GitHub HEAD into build/R*/ and zips that library into the .jmo.
# When the engine you want to ship exists only locally (for example a Nestimate
# release that is not pushed yet), pass its source with --engine: it is
# installed into the build library and replaces the GitHub copy in the .jmo.
# The engines are pure R, so any R 4.5 build can install them.
#
# Usage, from the repository root, with jamovi's Rscript:
#   <jamovi Rscript> tools/build_module.R
#   <jamovi Rscript> tools/build_module.R --engine ~/Documents/Github/Nestimate
#
# Afterwards install by sideloading the .jmo in jamovi (Modules > Sideload), or
# unzip it into the jamovi modules directory (see README).

stopifnot("run from the repository root" = file.exists("DESCRIPTION"))

args <- commandArgs(trailingOnly = TRUE)
engines <- args[which(args == "--engine") + 1L]
engines <- normalizePath(engines[!is.na(engines)], mustWork = TRUE)
r_for_install <- Sys.getenv("CNA_INSTALL_R", Sys.which("R"))
jamovi_home <- getOption("jamovi_home", "/Applications/jamovi.app")

version <- read.dcf("DESCRIPTION")[, "Version"][[1]]
# absolute: the archive is rebuilt from inside a temporary directory
jmo <- file.path(normalizePath("."), sprintf("CNA_%s.jmo", version))

run <- function(cmd, args, what) {
  status <- system2(cmd, args)
  if (!identical(status, 0L))
    stop(sprintf("%s failed (exit status %s)", what, status), call. = FALSE)
  invisible(status)
}

# 1. compile the module and its dependency library (no install into jamovi)
node <- get("node", envir = asNamespace("jmvtools"))()
run(node, c(jmvtools:::jmcPath(), "--build", shQuote("."), "--home", shQuote(jamovi_home),
            "--rpath", shQuote(R.home(component = "bin"))), "jamovi compiler")
lib <- list.files("build", pattern = "^R[0-9]", full.names = TRUE)[[1L]]

# 2. replace engine packages with local builds and repack the .jmo
if (length(engines)) {
  pkgs <- vapply(engines, function(d) read.dcf(file.path(d, "DESCRIPTION"))[, "Package"][[1]],
                 character(1))
  invisible(Map(function(dir, pkg) {
    message(sprintf("Bundling %s from %s", pkg, dir))
    run(r_for_install, c("CMD", "INSTALL", "--no-multiarch", "-l", shQuote(lib), shQuote(dir)),
        paste("installing", pkg))
  }, engines, pkgs))

  staging <- tempfile("cna-jmo-")
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  utils::unzip(jmo, exdir = staging)
  invisible(lapply(pkgs, function(pkg) {
    target <- file.path(staging, "CNA", "R", pkg)
    unlink(target, recursive = TRUE)
    file.copy(file.path(lib, pkg), file.path(staging, "CNA", "R"), recursive = TRUE)
  }))
  unlink(jmo)
  # zip from inside the staging directory so the archive root is CNA/
  run("sh", c("-c", shQuote(sprintf("cd %s && zip -r9Xq %s CNA", shQuote(staging), shQuote(jmo)))),
      "zip")
}

# 3. report what the .jmo carries
check <- tempfile("cna-jmo-check-")
utils::unzip(jmo, exdir = check)
bundled <- list.files(file.path(check, "CNA", "R"))
versions <- vapply(bundled, function(p)
  read.dcf(file.path(check, "CNA", "R", p, "DESCRIPTION"))[, "Version"][[1]], character(1))
unlink(check, recursive = TRUE)
cat(sprintf("\nBuilt %s (%.1f MB) with:\n", basename(jmo), file.size(jmo) / 1e6))
cat(sprintf("  %-12s %s\n", names(versions), versions), sep = "")
