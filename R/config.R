#' Paths and MCMC settings from the thesis notebook.
#' Do not change n.chains / n.iter / n.burnin / n.thin / seed: those define the fit.

project_root <- function() {
  if (requireNamespace("here", quietly = TRUE)) {
    tryCatch(here::here(), error = function(e) getwd())
  } else {
    getwd()
  }
}

path_raw <- function() file.path(project_root(), "data", "raw")
path_derived <- function() file.path(project_root(), "data", "derived")
path_bugs_model <- function() {
  file.path(project_root(), "inst", "bugs", "hierarchical_dose_response_model_sigmapriors.txt")
}
path_figures <- function() file.path(project_root(), "docs", "figures")
path_output <- function() file.path(project_root(), "output")

HOST_STATUS <- "Normal"
RNG_SEED <- 123L

# Original bugs() call in Hierarchical_Beta_Binomial_model.Rmd
N_CHAINS <- 3L
N_ITER <- 10000L
N_BURNIN <- 1000L
N_THIN <- 2L

MONITORED <- c(
  "w_0", "z_0", "sigma_w", "sigma_z",
  "w", "z",
  "alpha0", "beta0", "alpha", "beta",
  "znew", "wnew", "alphanew", "betanew"
)
