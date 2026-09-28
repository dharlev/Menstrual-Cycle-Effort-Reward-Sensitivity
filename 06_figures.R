# Analytical figures from trial data, model estimates and posterior draws
posterior_plot_data <- function(fit, sessions) {
  x <- fit$draws(variables = c("effSens_sess", "rewSens_sess", "delta_eff", "delta_rew"), format = "matrix")
  rows <- list()
  for (v in c("effSens", "rewSens")) for (ph in PHASES) {
    ids <- sessions$session_idx[as.character(sessions$phase) == ph]
    values <- rowMeans(x[, paste0(v, "_sess[", ids, "]"), drop = FALSE])
    rows[[paste(v, ph)]] <- data.frame(parameter = v, phase = ph,
      mean = mean(values), lower = unname(stats::quantile(values, .025)),
      upper = unname(stats::quantile(values, .975)))
  }
  list(means = dplyr::bind_rows(rows), contrasts = summarise_draws(x, c("delta_eff", "delta_rew")))
}

save_figure <- function(plot, output_dir, name, width = 10, height = 7) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(output_dir, paste0(name, ".pdf")), plot, width = width, height = height)
  ggplot2::ggsave(file.path(output_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 300)
}

make_figures <- function(inputs, secondary, fit, prepared, behaviour, recovery, output_dir) {
  palette <- c(Follicular = "#4FA8D5", Luteal = "#F28E72")
  phase_labels <- c(Follicular = "Late-follicular", Luteal = "Mid-luteal")
  theme <- ggplot2::theme_classic(base_size = 11) + ggplot2::theme(legend.position = "bottom")
  paired_plot <- function(d, variable, label) {
    ggplot2::ggplot(d, ggplot2::aes(x = phase, y = .data[[variable]], group = participant_id)) +
      ggplot2::geom_line(colour = "grey80", linewidth = .4) +
      ggplot2::geom_point(ggplot2::aes(colour = phase), alpha = .65) +
      ggplot2::stat_summary(ggplot2::aes(group = phase), fun.data = function(y) {
        se <- stats::sd(y) / sqrt(length(y)); m <- mean(y); h <- stats::qt(.975, length(y) - 1) * se
        data.frame(y = m, ymin = m - h, ymax = m + h)
      }, geom = "pointrange", colour = "black", linewidth = .5) +
      ggplot2::scale_colour_manual(values = palette) + ggplot2::scale_x_discrete(labels = phase_labels) +
      ggplot2::labs(x = NULL, y = label, colour = NULL) + theme
  }
  scatter <- function(d, x, y, xlab, ylab) {
    d <- d[is.finite(d[[x]]) & is.finite(d[[y]]), ]
    ggplot2::ggplot(d, ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = phase, fill = phase)) +
      ggplot2::geom_point(alpha = .6) + ggplot2::geom_smooth(method = "lm", formula = y ~ x, linewidth = .7) +
      ggplot2::scale_colour_manual(values = palette) + ggplot2::scale_fill_manual(values = palette) +
      ggplot2::labs(x = xlab, y = ylab, colour = NULL, fill = NULL) + theme
  }
  curves <- inputs$ep_phase |>
    dplyr::group_by(participant_id, phase, force) |>
    dplyr::summarise(mean_rating = mean(scale_rating), .groups = "drop")
  curve_means <- curves |> dplyr::group_by(phase, force) |>
    dplyr::summarise(mean_rating = mean(mean_rating), .groups = "drop")
  a <- ggplot2::ggplot(curves, ggplot2::aes(force, mean_rating, colour = phase,
                                         group = interaction(participant_id, phase))) +
    ggplot2::geom_line(alpha = .1) +
    ggplot2::geom_line(data = curve_means, ggplot2::aes(group = phase), linewidth = 1.2) +
    ggplot2::geom_point(data = curve_means, ggplot2::aes(group = phase), size = 2) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(x = "Force (% MVC)", y = "Effort rating", colour = NULL) + theme
  b <- paired_plot(secondary$slopes, "effort_slope", "Effort differentiation slope")
  mod_slopes <- ep_session_slopes(inputs$ep_mod)
  c <- scatter(dplyr::inner_join(mod_slopes, dplyr::distinct(secondary$ep_hrv, participant_id, phase, hrv_mvc_z),
                                 by = c("participant_id", "phase")),
                "hrv_mvc_z", "effort_slope", "MVC log RMSSD (z)", "Effort differentiation slope")
  d <- scatter(dplyr::inner_join(mod_slopes, dplyr::distinct(secondary$ep_ema, participant_id, phase, valence_sd_z),
                                 by = c("participant_id", "phase")),
                "valence_sd_z", "effort_slope", "Valence variability (z)", "Effort differentiation slope")
  save_figure(patchwork::wrap_plots(a, b, c, d, ncol = 2) + patchwork::plot_annotation(tag_levels = "A"),
               output_dir, "Figure_3_effort_perception", 11, 8)
  choices <- prepared$trials |>
    dplyr::mutate(hard_reward = ifelse(effort1 > effort2, reward1, reward2),
      hard_effort = pmax(effort1, effort2),
      choose_hard = as.integer((effort1 > effort2 & option_chosen == 0) | (effort2 > effort1 & option_chosen == 1))) |>
    dplyr::group_by(participant_id, phase, hard_reward, hard_effort) |>
    dplyr::summarise(p_hard = mean(choose_hard), .groups = "drop") |>
    dplyr::group_by(phase, hard_reward, hard_effort) |>
    dplyr::summarise(p_hard = mean(p_hard), n_total = dplyr::n(), .groups = "drop") |>
    tidyr::pivot_wider(names_from = phase, values_from = c(p_hard, n_total)) |>
    dplyr::mutate(difference = 100 * (p_hard_Luteal - p_hard_Follicular),
                  n_total = dplyr::coalesce(n_total_Follicular, 0L) + dplyr::coalesce(n_total_Luteal, 0L))
  write_table(choices, dirname(output_dir), "decision_space_values")
  surface_path <- file.path(output_dir, "Figure_2A_decision_space.png")
  draw_choice_diff_3d(dplyr::mutate(choices, diff_pp = difference), surface_path,
                      "Location of phase differences in decision space")
  a <- patchwork::wrap_elements(full = grid::rasterGrob(png::readPNG(surface_path), interpolate = TRUE))
  pp <- posterior_plot_data(fit, prepared$sessions)
  write_table(pp$means, dirname(output_dir), "posterior_phase_means")
  pp$means$label <- paste(ifelse(pp$means$parameter == "effSens", "Effort", "Reward"),
                          ifelse(pp$means$phase == "Follicular", "(LF)", "(ML)"))
  b <- ggplot2::ggplot(pp$means, ggplot2::aes(mean, label, colour = phase)) +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = label)) +
    ggplot2::geom_point() + ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(x = "Sensitivity: posterior mean and 95% CrI", y = NULL, colour = NULL) + theme
  pp$contrasts$label <- ifelse(pp$contrasts$parameter == "delta_eff", "Effort", "Reward")
  c <- ggplot2::ggplot(pp$contrasts, ggplot2::aes(mean, label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_segment(ggplot2::aes(x = lower, xend = upper, yend = label), colour = "#8E6BBE") +
    ggplot2::geom_point(colour = "#8E6BBE") + ggplot2::labs(x = "Log sensitivity: ML - LF (95% CrI)", y = NULL) + theme
  dmplot <- secondary$dm_ema |> dplyr::mutate(effSens_z = z_score(effSens))
  d <- scatter(dmplot, "arousal_mean_z", "effSens_z", "Mean arousal (z)", "Effort sensitivity (z)")
  e <- scatter(dmplot, "valence_mean_z", "effSens_z", "Mean valence (z)", "Effort sensitivity (z)")
  save_figure(patchwork::wrap_plots(a, b, c, d, e, ncol = 2) + patchwork::plot_annotation(tag_levels = "A"),
               output_dir, "Figure_2_decision_making", 12, 11)
  panels <- list(paired_plot(behaviour, "hard_pct", "Hard choices (%)"),
                 paired_plot(behaviour, "success_pct", "Execution success (%)"))
  if (!is.null(recovery)) for (v in c("effSens", "rewSens")) {
    panels[[length(panels) + 1L]] <- ggplot2::ggplot(recovery,
      ggplot2::aes(.data[[paste0("true_", v)]], .data[[v]])) +
      ggplot2::geom_point(alpha = .6) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::labs(x = paste("Generating", v), y = paste("Recovered", v)) + theme
  }
  save_figure(patchwork::wrap_plots(panels, ncol = 2) + patchwork::plot_annotation(tag_levels = "A"),
               output_dir, "Figure_S2_behaviour_recovery", 10, if (is.null(recovery)) 4 else 8)
  ep_ami <- dplyr::left_join(mod_slopes,
    inputs$ep_mod |> dplyr::group_by(participant_id, phase) |>
      dplyr::summarise(AMI_Score = dplyr::first(AMI_Score), .groups = "drop"),
    by = c("participant_id", "phase"))
  pooled <- function(d, x, y, label) {
    d <- d[is.finite(d[[x]]) & is.finite(d[[y]]), ]
    ggplot2::ggplot(d, ggplot2::aes(.data[[x]], .data[[y]])) + ggplot2::geom_point(alpha = .55) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
      ggplot2::labs(x = "AMI total score", y = label) + theme
  }
  save_figure(patchwork::wrap_plots(pooled(ep_ami, "AMI_Score", "effort_slope", "Effort differentiation slope"),
    pooled(secondary$dm, "AMI_Score", "effSens", "Effort sensitivity")) + patchwork::plot_annotation(tag_levels = "A"),
    output_dir, "Figure_S3_apathy", 10, 4)
  ppc <- as.data.frame(fit$draws(variables = c("ybar_rep", "ybar_obs"), format = "matrix"))
  save_figure(ggplot2::ggplot(ppc, ggplot2::aes(ybar_rep)) + ggplot2::geom_histogram(bins = 40, fill = "grey70") +
    ggplot2::geom_vline(xintercept = ppc$ybar_obs[1], colour = "#B83B41") +
    ggplot2::labs(x = "Replicated proportion choosing option 2", y = "Posterior draws") + theme,
    output_dir, "DM_posterior_predictive_check", 6, 4)
}

draw_choice_diff_3d <- function(dat, out_path, title_txt) {
  ord <- dat[order(dat$hard_reward, dat$hard_effort), ]
  fit <- loess(
    diff_pp ~ hard_reward * hard_effort,
    data = ord,
    weights = n_total,
    span = 0.60,
    degree = 2,
    control = loess.control(surface = "direct")
  )
  xg <- sort(unique(ord$hard_reward))
  yg <- sort(unique(ord$hard_effort))
  grid_df <- expand.grid(hard_reward = xg, hard_effort = yg)
  grid_df$z <- as.vector(predict(fit, newdata = grid_df))
  grid_df$z[!is.finite(grid_df$z)] <- mean(ord$diff_pp, na.rm = TRUE)
  rng <- max(abs(c(grid_df$z, ord$diff_pp)), na.rm = TRUE)
  zlim <- c(-rng, rng)
  zmat <- matrix(grid_df$z, nrow = length(xg), ncol = length(yg))
  surf_pal <- colorRampPalette(c("#313695", "#4575B4", "#74ADD1", "#FEE08B", "#F46D43", "#A50026"))
  zbreaks <- seq(zlim[1], zlim[2], length.out = 101)
  facet_cols <- matrix(
    surf_pal(100)[cut(zmat[-1, -1], breaks = zbreaks, include.lowest = TRUE)],
    nrow = nrow(zmat) - 1,
    ncol = ncol(zmat) - 1
  )

  png(out_path, width = 1500, height = 1200, res = 220)
  layout(matrix(c(1, 2), nrow = 2), heights = c(1.0, 9.0))
  par(mar = c(0.4, 3.0, 0.0, 3.0))
  plot.new()
  xleft <- seq(0.10, 0.90, length.out = 100)
  xright <- seq(0.10 + 0.80 / 100, 0.90 + 0.80 / 100, length.out = 100)
  rect(xleft, 0.22, xright, 0.50, col = surf_pal(100), border = NA, xpd = NA)
  legend_ticks <- pretty(zlim, n = 5)
  tick_x <- 0.10 + 0.80 * ((legend_ticks - min(zlim)) / diff(range(zlim)))
  segments(tick_x, 0.52, tick_x, 0.58, xpd = NA)
  text(tick_x, 0.62, labels = round(legend_ticks, 1), cex = 0.88, xpd = NA)
  text(0.50, 0.82, "Mid-luteal - late-follicular difference in hard choice (%)", cex = 0.92, font = 2, xpd = NA)

  par(mar = c(0.2, 0.2, 0.2, 0.2))
  zfloor <- matrix(min(zlim), nrow = length(xg), ncol = length(yg))
  persp(
    x = xg, y = yg, z = zfloor, theta = 50, phi = 25, expand = 0.78, ticktype = "detailed",
    col = matrix(adjustcolor("grey78", 0.72), nrow = length(xg) - 1, ncol = length(yg) - 1),
    shade = 0, border = adjustcolor("grey58", 0.40), xlab = "", ylab = "", zlab = "",
    main = "", zlim = zlim, box = FALSE, axes = FALSE
  )
  par(new = TRUE)
  persp(
    x = xg, y = yg, z = zmat, theta = 50, phi = 25, expand = 0.78, ticktype = "detailed",
    col = facet_cols, shade = 0.10, border = adjustcolor("grey30", 0.55),
    xlab = "Reward (points)", ylab = "Effort (%MVC)", zlab = "Mid-luteal - late-follicular (%)",
    main = "", zlim = zlim, box = TRUE, axes = TRUE
  )
  title(title_txt, line = -0.45, cex.main = 1.10, font.main = 2)
  dev.off()
}
