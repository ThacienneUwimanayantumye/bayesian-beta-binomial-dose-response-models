#' Outbreak rows used in the thesis: host status "Normal" only.
#' Strain order is first-appearance in that subset (not alphabetical).

load_salmonella <- function(path = file.path(path_raw(), "salmonella.csv")) {
  read.csv(path, header = TRUE)
}

prepare_bugs_data <- function(raw = load_salmonella(), host_status = HOST_STATUS) {
  data_in <- subset(raw, S == host_status)
  serovars <- unique(data_in$t)
  strain <- as.integer(factor(data_in$t, levels = serovars, labels = seq_along(serovars)))
  list(
    raw_subset = data_in,
    serovars = serovars,
    bugs = list(
      y = data_in$Y,
      d = 10^data_in$log10dose,
      n = data_in$N,
      strain = strain,
      N_total = length(data_in$Y),
      K = length(serovars)
    )
  )
}

bugs_inits <- function(K) {
  function() {
    list(
      w_0 = 0,
      z_0 = 0,
      alpha0 = 0,
      beta0 = 0,
      sigma_w = 1,
      sigma_z = 1,
      w = rep(0, K),
      z = rep(0, K),
      alpha = rep(0, K),
      beta = rep(0, K),
      wnew = 0,
      znew = 0,
      alphanew = 0,
      betanew = 0
    )
  }
}

infection_prob <- function(alpha, beta, dose) {
  1 - exp(
    lgamma(alpha + beta) + lgamma(beta + dose) -
      lgamma(beta) - lgamma(alpha + beta + dose)
  )
}

curve_quantile <- function(alpha, beta, dose_levels, p = 0.5) {
  vapply(dose_levels, function(dose) {
    stats::quantile(infection_prob(alpha, beta, dose), probs = p, names = FALSE)
  }, numeric(1))
}

curve_median <- function(alpha, beta, dose_levels) {
  vapply(dose_levels, function(dose) {
    stats::median(infection_prob(alpha, beta, dose))
  }, numeric(1))
}

ed50_dose <- function(median_probs, doses, target_prob = 0.5) {
  doses[[which.min(abs(median_probs - target_prob))]]
}
