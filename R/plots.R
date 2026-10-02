#' Thesis figures. Geometry matches Hierarchical_Beta_Binomial_model.Rmd.
#' Writes to output/ so docs/figures/ (the original PDFs) stay untouched.

ensure_output <- function() {
  dir.create(path_output(), recursive = TRUE, showWarnings = FALSE)
  path_output()
}

dose_grid <- function(d, n = 100) {
  max_dose <- max(d)
  log10_dose <- seq(0, log10(max_dose), length.out = n)
  list(log10_dose = log10_dose, dose_levels = 10^log10_dose)
}

plot_traces <- function(results, path = file.path(ensure_output(), "HBB_trace.pdf")) {
  posterior_samples <- results$sims.list
  trace_data <- list(
    w_0 = posterior_samples$w_0,
    z_0 = posterior_samples$z_0,
    sigma_w = posterior_samples$sigma_w,
    sigma_z = posterior_samples$sigma_z,
    w = posterior_samples$w,
    z = posterior_samples$z
  )

  trace_data_long <- do.call(rbind, lapply(names(trace_data), function(param) {
    samples <- trace_data[[param]]
    if (is.matrix(samples)) {
      num_chains <- N_CHAINS
      num_iterations <- nrow(samples) / num_chains
      if (num_iterations %% 1 != 0) stop("Number of iterations is not an integer.")
      samples_long <- as.data.frame(samples)
      colnames(samples_long) <- paste0(param, "_", seq_len(ncol(samples_long)))
      samples_long <- tidyr::pivot_longer(
        samples_long,
        cols = names(samples_long),
        names_to = "parameter",
        values_to = "value"
      )
      samples_long$chain <- rep(1:num_chains, each = num_iterations * ncol(samples))
      samples_long$iteration <- rep(1:num_iterations, num_chains * ncol(samples))
    } else {
      num_chains <- N_CHAINS
      num_iterations <- length(samples) / num_chains
      if (num_iterations %% 1 != 0) stop("Number of iterations is not an integer.")
      samples_long <- data.frame(
        parameter = param,
        value = samples,
        chain = rep(1:num_chains, each = num_iterations),
        iteration = rep(1:num_iterations, num_chains)
      )
    }
    samples_long
  }))

  p <- ggplot2::ggplot(trace_data_long, ggplot2::aes(x = iteration, y = value, color = factor(chain))) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      title = "Trace Plots for Posterior Samples (by Chain)",
      x = "Iteration", y = "Parameter Value", color = "Chain"
    ) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(path, p, width = 11, height = 11)
  p
}

create_contour_plot <- function(alpha, beta, strain_id) {
  density_est <- MASS::kde2d(alpha, beta, n = 100)
  density_df <- with(density_est, expand.grid(x = x, y = y))
  density_df$z <- as.vector(density_est$z)
  ggplot2::ggplot(data = density_df, ggplot2::aes(x = x, y = y, z = z)) +
    ggplot2::geom_contour_filled(ggplot2::aes(fill = ggplot2::after_stat(level))) +
    ggplot2::labs(
      title = strain_id,
      x = expression(log10(alpha)),
      y = expression(log10(beta)),
      fill = "Density"
    ) +
    ggplot2::theme_minimal()
}

plot_contours <- function(results, serovars, path = file.path(ensure_output(), "HBB_cont.pdf")) {
  a_0 <- results$sims.list$alpha0
  b_0 <- results$sims.list$beta0
  a <- results$sims.list$alpha
  b <- results$sims.list$beta
  a_new <- results$sims.list$alphanew
  b_new <- results$sims.list$betanew
  alpha_matrix <- cbind(log10(as.matrix(a_0)), log10(a), log10(a_new))
  beta_matrix <- cbind(log10(as.matrix(b_0)), log10(b), log10(b_new))
  serovars_new <- c("Overall_strain", serovars, "New_strain")
  K <- ncol(alpha_matrix)
  plot_list <- vector("list", K)
  for (k in 1:K) {
    plot_list[[k]] <- create_contour_plot(alpha_matrix[, k], beta_matrix[, k], serovars_new[k])
  }
  plot_grobs <- lapply(plot_list, ggplot2::ggplotGrob)
  main_title <- grid::textGrob(
    "Contour plot of log10(alpha) vs log10(beta)",
    gp = grid::gpar(fontsize = 16, fontface = "bold")
  )
  p <- gridExtra::grid.arrange(grobs = plot_grobs, ncol = 2, top = main_title)
  ggplot2::ggsave(path, p, width = 10, height = 11)
  p
}

create_histogram_with_ci <- function(u, strain_id) {
  ci_lower <- stats::quantile(u, probs = 0.025)
  ci_upper <- stats::quantile(u, probs = 0.975)
  u_df <- data.frame(u = u)
  strain_title <- paste(strain_id, "(95% CI: [", signif(ci_lower, 3), ", ", signif(ci_upper, 3), "])")
  ggplot2::ggplot(u_df, ggplot2::aes(x = u)) +
    ggplot2::geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    ggplot2::geom_vline(xintercept = ci_lower, color = "red", linetype = "dashed", size = 1) +
    ggplot2::geom_vline(xintercept = ci_upper, color = "red", linetype = "dashed", size = 1) +
    ggplot2::labs(title = strain_title, x = "u", y = "Frequency") +
    ggplot2::theme_minimal()
}

plot_u_histograms <- function(results, serovars, path = file.path(ensure_output(), "HBB_U_with_CI_in_Title.pdf")) {
  w_matrix <- cbind(results$sims.list$w_0, results$sims.list$w, results$sims.list$wnew)
  u_matrix <- 1 / (1 + exp(-w_matrix))
  serovars_new <- c("Overall_strain", serovars, "New_strain")
  K <- ncol(w_matrix)
  histogram_list <- vector("list", K)
  for (k in 1:K) {
    histogram_list[[k]] <- create_histogram_with_ci(u_matrix[, k], serovars_new[k])
  }
  p <- do.call(gridExtra::grid.arrange, c(histogram_list, ncol = 2, top = "Histograms of u by Strains with 95% Credible Intervals"))
  ggplot2::ggsave(path, p, width = 11, height = 11)
  p
}

plot_infection_curves <- function(results, prepared, path = file.path(ensure_output(), "HBB_Pinf_with_ED50_segments.pdf")) {
  bugs <- prepared$bugs
  serovars <- prepared$serovars
  grid <- dose_grid(bugs$d)
  log10_dose <- grid$log10_dose
  dose_levels <- grid$dose_levels
  a_0 <- results$sims.list$alpha0
  b_0 <- results$sims.list$beta0
  a <- results$sims.list$alpha
  b <- results$sims.list$beta
  a_new <- results$sims.list$alphanew
  b_new <- results$sims.list$betanew
  K <- length(serovars)

  overall_prob_median <- curve_median(a_0, b_0, dose_levels)
  overall_prob_lower <- curve_quantile(a_0, b_0, dose_levels, 0.025)
  overall_prob_upper <- curve_quantile(a_0, b_0, dose_levels, 0.975)
  new_prob_median <- curve_median(a_new, b_new, dose_levels)
  new_prob_lower <- curve_quantile(a_new, b_new, dose_levels, 0.025)
  new_prob_upper <- curve_quantile(a_new, b_new, dose_levels, 0.975)

  predicted_medians <- sapply(1:K, function(k) curve_median(a[, k], b[, k], dose_levels))
  predicted_lowers <- sapply(1:K, function(k) curve_quantile(a[, k], b[, k], dose_levels, 0.025))
  predicted_uppers <- sapply(1:K, function(k) curve_quantile(a[, k], b[, k], dose_levels, 0.975))

  plot_data <- data.frame(
    log10_dose = rep(log10_dose, K),
    dose = rep(dose_levels, K),
    strain = factor(rep(serovars, each = length(log10_dose))),
    median_prob = as.vector(predicted_medians),
    lower_95CI = as.vector(predicted_lowers),
    upper_95CI = as.vector(predicted_uppers)
  )
  overall_data_expanded <- data.frame(
    log10_dose = log10_dose,
    dose = dose_levels,
    strain = factor(rep("Overall_strain", length(log10_dose))),
    median_prob = overall_prob_median,
    lower_95CI = overall_prob_lower,
    upper_95CI = overall_prob_upper
  )
  new_strain_data_expanded <- data.frame(
    log10_dose = log10_dose,
    dose = dose_levels,
    strain = factor(rep("New Strain", length(log10_dose))),
    median_prob = new_prob_median,
    lower_95CI = new_prob_lower,
    upper_95CI = new_prob_upper
  )
  plot_data_all <- rbind(plot_data, overall_data_expanded, new_strain_data_expanded)

  plot_list <- lapply(1:K, function(g) {
    group_indices <- which(bugs$strain == g)
    data.frame(
      log10_dose = log10(bugs$d[group_indices]),
      dose = bugs$d[group_indices],
      strain = as.factor(serovars[g]),
      observed_prob = bugs$y[group_indices] / bugs$n[group_indices],
      n = bugs$n[group_indices]
    )
  })
  observed_data <- do.call(rbind, plot_list)
  overall_observed_data <- data.frame(
    log10_dose = log10(bugs$d),
    dose = bugs$d,
    strain = factor(rep("Overall_strain", length(bugs$d))),
    observed_prob = bugs$y / bugs$n,
    n = bugs$n
  )
  new_strain_data <- data.frame(
    log10_dose = log10(bugs$d),
    dose = bugs$d,
    strain = factor(rep("New Strain", length(bugs$d))),
    observed_prob = NA_real_,
    n = NA_real_
  )
  observed_data_all <- rbind(observed_data, overall_observed_data, new_strain_data)

  ed50_doses <- sapply(1:K, function(k) ed50_dose(predicted_medians[, k], dose_levels))
  overall_ed50_dose <- ed50_dose(overall_prob_median, dose_levels)
  new_strain_ed50_dose <- ed50_dose(new_prob_median, dose_levels)
  Ed50_doses <- c(ed50_doses, overall_ed50_dose, new_strain_ed50_dose)
  all_strains <- c(serovars, "Overall_strain", "New Strain")
  ed50_labels <- data.frame(
    strain = factor(all_strains),
    ed50 = log10(Ed50_doses),
    label = signif(log10(Ed50_doses), 4),
    y_pos = -0.25
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = plot_data_all,
      ggplot2::aes(x = log10_dose, ymin = lower_95CI, ymax = upper_95CI, fill = strain),
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = plot_data_all,
      ggplot2::aes(x = log10_dose, y = median_prob, color = strain),
      lwd = 0.5
    ) +
    ggplot2::geom_point(
      data = observed_data_all,
      ggplot2::aes(x = log10_dose, y = observed_prob, shape = "Observed Data"),
      size = 2, fill = "black", na.rm = TRUE
    ) +
    ggplot2::geom_segment(
      data = ed50_labels,
      ggplot2::aes(x = ed50, xend = ed50, y = 0, yend = 0.5, color = strain),
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      data = ed50_labels,
      ggplot2::aes(x = 0, xend = ed50, y = 0.5, yend = 0.5, color = strain),
      linetype = "dashed"
    ) +
    ggplot2::geom_text(
      data = ed50_labels,
      ggplot2::aes(x = ed50, y = y_pos, label = label, color = strain),
      angle = 45, hjust = -0.1, vjust = 0.5, size = 3
    ) +
    ggplot2::labs(
      x = "log10 Dose",
      y = "Probability of Infection",
      title = "Beta-Binomial Probability of Infection by Dose Level and Strain"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(name = "Strain", values = c(grDevices::rainbow(K), "purple", "blue")) +
    ggplot2::scale_fill_manual(name = "Strain", values = c(grDevices::rainbow(K), "purple", "blue")) +
    ggplot2::scale_shape_manual(name = "Data Type", values = c("Observed Data" = 21)) +
    ggplot2::facet_wrap(~strain, scales = "free_y") +
    ggplot2::theme(legend.position = "right")
  ggplot2::ggsave(path, p, width = 10, height = 6)
  invisible(list(plot = p, ed50 = ed50_labels))
}

write_all_figures <- function(results, prepared) {
  plot_traces(results)
  plot_contours(results, prepared$serovars)
  plot_u_histograms(results, prepared$serovars)
  plot_infection_curves(results, prepared)
}
