#!/usr/bin/env Rscript
# Compatibility entry point: the comprehensive revision supersedes the partial plot pass.
a <- commandArgs(FALSE)
here <- dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", a[grepl("^--file=", a)][1]), fixed=TRUE)))
source(file.path(here, "07_complete_reference_outputs.R"))
