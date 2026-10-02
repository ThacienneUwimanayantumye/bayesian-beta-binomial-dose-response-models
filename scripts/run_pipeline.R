#!/usr/bin/env Rscript
#' Reproduce the thesis hierarchical beta-binomial fit.
#' Same model file, same priors, same seed and MCMC settings.
#'
#'   Rscript scripts/run_pipeline.R
#'   Rscript scripts/run_pipeline.R --prepare-only
#'   Rscript scripts/run_pipeline.R --refit

args <- commandArgs(trailingOnly = TRUE)
refit <- "--refit" %in% args
prepare_only <- "--prepare-only" %in% args

root <- tryCatch(here::here(), error = function(e) normalizePath("."))
if (!file.exists(file.path(root, "R", "config.R"))) {
  root <- normalizePath(file.path(".."))
}
setwd(root)

source(file.path("R", "config.R"))
source(file.path("R", "prepare_data.R"))
source(file.path("R", "fit.R"))
source(file.path("R", "plots.R"))

prepared <- prepare_bugs_data()
message(
  "Host status ", HOST_STATUS, ": ", prepared$bugs$N_total,
  " dose groups, ", prepared$bugs$K, " strains (",
  paste(prepared$serovars, collapse = ", "), ")."
)
print(prepared$raw_subset[, c("t", "S", "log10dose", "N", "Y")])

if (prepare_only) {
  message("Stopped after data prep (--prepare-only). OpenBUGS was not called.")
  quit(save = "no", status = 0)
}

results <- load_or_fit(prepared, refit = refit)
table <- export_posterior_summary(results)
print(table)

if (file.exists(fit_path()) || exists("results")) {
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("gridExtra", quietly = TRUE) &&
      requireNamespace("MASS", quietly = TRUE) &&
      requireNamespace("tidyr", quietly = TRUE)) {
    write_all_figures(results, prepared)
    message("Figures written to ", path_output())
  } else {
    message("Skip figures: install ggplot2, gridExtra, MASS, tidyr, dplyr to plot.")
  }
}

message("Done. Posterior summary: ", file.path(path_derived(), "posterior_summary.csv"))
