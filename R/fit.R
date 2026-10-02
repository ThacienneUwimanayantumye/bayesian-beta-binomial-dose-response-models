#' Fit the hierarchical beta-binomial model with the thesis MCMC settings.
#' Requires OpenBUGS and R2OpenBUGS. debug=FALSE so a script can run headless;
#' that flag does not change the posterior.

fit_hbb <- function(prepared,
                    model_file = path_bugs_model(),
                    n.chains = N_CHAINS,
                    n.iter = N_ITER,
                    n.burnin = N_BURNIN,
                    n.thin = N_THIN,
                    debug = FALSE,
                    seed = RNG_SEED) {
  if (!requireNamespace("R2OpenBUGS", quietly = TRUE)) {
    stop(
      "R2OpenBUGS is not installed. Install OpenBUGS, then ",
      "install.packages('R2OpenBUGS')."
    )
  }
  if (!file.exists(model_file)) {
    stop("BUGS model not found: ", model_file)
  }
  set.seed(seed)
  R2OpenBUGS::bugs(
    data = prepared$bugs,
    inits = bugs_inits(prepared$bugs$K),
    parameters.to.save = MONITORED,
    model.file = model_file,
    n.chains = n.chains,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.thin = n.thin,
    DIC = TRUE,
    debug = debug
  )
}

fit_path <- function() file.path(path_derived(), "bugs_fit.rds")

load_or_fit <- function(prepared, refit = FALSE) {
  cached <- fit_path()
  dir.create(path_derived(), recursive = TRUE, showWarnings = FALSE)
  if (!refit && file.exists(cached)) {
    message("Using cached fit: ", cached)
    return(readRDS(cached))
  }
  message("Fitting OpenBUGS model (seed ", RNG_SEED, ")...")
  results <- fit_hbb(prepared)
  saveRDS(results, cached)
  results
}

export_posterior_summary <- function(results, path = file.path(path_derived(), "posterior_summary.csv")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  table <- signif(results$summary, 2)
  utils::write.csv(table, path)
  table
}
