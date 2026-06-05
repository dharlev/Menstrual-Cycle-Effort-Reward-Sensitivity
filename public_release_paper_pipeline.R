#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggdist)
  library(patchwork)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(signal)
  library(pracma)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

project_dir <- "."
raw_dir <- file.path(project_dir, "data", "raw")
processed_dir <- file.path(project_dir, "data", "processed")
figures_dir <- file.path(project_dir, "outputs", "figures")
models_dir <- file.path(project_dir, "outputs", "models")
tables_dir <- file.path(project_dir, "outputs", "tables")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

phase_levels <- c("Follicular", "Luteal")
phase_labels <- c(Follicular = "Late-follicular", Luteal = "Mid-luteal")

col_fol <- "#4FA8D5"
col_fol_soft <- "#D5EAF6"
col_lut <- "#F28E72"
col_lut_soft <- "#F8D5C8"
col_contrast <- "#8E6BBE"

ep_outlier_sd <- 2.5
min_rr_baseline <- 30L
min_rr_active <- 10L
n_boot <- 4000L
seed_base <- 20260605L

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

read_csv_clean <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

require_columns <- function(df, cols, name) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("%s is missing columns: %s", name, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
}

z_score <- function(x) {
  as.numeric(scale(x))
}

boot_mean_draws <- function(x, n = n_boot, seed = seed_base) {
  x <- x[is.finite(x)]
  set.seed(seed)
  replicate(n, mean(sample(x, length(x), replace = TRUE), na.rm = TRUE))
}

summary_to_draws <- function(mean, lower, upper, n = n_boot, seed = seed_base) {
  sd_approx <- (upper - lower) / (2 * qnorm(0.975))
  set.seed(seed)
  rnorm(n, mean = mean, sd = sd_approx)
}

tidy_lmm <- function(model, label) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(model = label)
}

theme_pub <- function(base_size = 11.5) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank(),
      legend.title = element_text(face = "bold")
    )
}

tag_panel <- function(p, tag) {
  p + labs(tag = tag) +
    theme_pub() +
    theme(plot.tag = element_text(face = "bold", size = 15))
}

save_plot_pair <- function(plot_obj, stem, width, height) {
  ggsave(file.path(figures_dir, paste0(stem, ".png")), plot_obj, width = width, height = height, dpi = 320)
  ggsave(file.path(figures_dir, paste0(stem, ".pdf")), plot_obj, width = width, height = height)
}

# -----------------------------------------------------------------------------
# ECG / HRV
# -----------------------------------------------------------------------------

bandpass_ecg <- function(ecg, fs, low = 0.5, high = 40, order = 2) {
  wn <- c(low, high) / (fs / 2)
  filt <- signal::butter(order, wn, type = "pass")
  as.numeric(signal::filtfilt(filt, ecg))
}

detect_r_peaks <- function(ecg_filt, fs, min_distance_sec = 0.30, min_prominence = NULL) {
  if (is.null(min_prominence)) {
    min_prominence <- 0.35 * stats::sd(ecg_filt, na.rm = TRUE)
  }

  peak_mat <- pracma::findpeaks(
    ecg_filt,
    minpeakdistance = round(min_distance_sec * fs),
    minpeakheight = stats::median(ecg_filt, na.rm = TRUE) + min_prominence
  )

  if (is.null(peak_mat)) {
    return(integer(0))
  }

  as.integer(peak_mat[, 2])
}

clean_rr_intervals <- function(rr_sec) {
  rr_sec <- rr_sec[is.finite(rr_sec)]
  rr_sec <- rr_sec[rr_sec >= 0.30 & rr_sec <= 1.50]
  if (length(rr_sec) < 3) {
    return(numeric(0))
  }
  rr_med <- stats::median(rr_sec, na.rm = TRUE)
  rr_sec[rr_sec >= 0.80 * rr_med & rr_sec <= 1.20 * rr_med]
}

compute_rmssd <- function(rr_sec) {
  if (length(rr_sec) < 2) {
    return(NA_real_)
  }
  sqrt(mean(diff(rr_sec)^2, na.rm = TRUE))
}

process_ecg_segment <- function(path, sampling_rate) {
  ecg <- read_csv_clean(path)
  signal_col <- names(ecg)[1]
  ecg_vec <- as.numeric(ecg[[signal_col]])
  ecg_filt <- bandpass_ecg(ecg_vec, fs = sampling_rate)
  r_idx <- detect_r_peaks(ecg_filt, fs = sampling_rate)
  rr_raw <- diff(r_idx) / sampling_rate
  rr_clean <- clean_rr_intervals(rr_raw)

  tibble(
    r_peak_n = length(r_idx),
    rr_clean_n = length(rr_clean),
    rmssd = compute_rmssd(rr_clean)
  )
}

build_hrv_table <- function(ecg_index) {
  require_columns(
    ecg_index,
    c("subject", "session_id", "phase", "sampling_rate", "baseline_file", "active_file"),
    "ecg_index"
  )

  purrr::pmap_dfr(
    ecg_index,
    function(subject, session_id, phase, sampling_rate, baseline_file, active_file, ...) {
      baseline <- process_ecg_segment(file.path(raw_dir, baseline_file), sampling_rate)
      active <- process_ecg_segment(file.path(raw_dir, active_file), sampling_rate)

      tibble(
        subject = subject,
        session_id = session_id,
        phase = phase,
        hrv_baseline = baseline$rmssd,
        hrv_active = active$rmssd,
        rr_clean_n_baseline = baseline$rr_clean_n,
        rr_clean_n_active = active$rr_clean_n,
        ecg_ok_baseline = baseline$r_peak_n > 1 && baseline$rr_clean_n >= min_rr_baseline,
        ecg_ok_active = active$r_peak_n > 1 && baseline$rr_clean_n >= 1 && active$rr_clean_n >= min_rr_active
      )
    }
  )
}

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------

phase_schedule <- read_csv_clean(file.path(raw_dir, "phase_schedule.csv"))
er_trials_raw <- read_csv_clean(file.path(raw_dir, "er_trials.csv"))
dm_trials_raw <- read_csv_clean(file.path(raw_dir, "dm_trials.csv"))
clinical_raw <- read_csv_clean(file.path(raw_dir, "clinical_subject_phase.csv"))
ema_raw <- read_csv_clean(file.path(raw_dir, "ema_subject_phase.csv"))

require_columns(phase_schedule, c("subject", "session_id", "phase", "phase_verified", "phase_target"), "phase_schedule")
require_columns(er_trials_raw, c("subject", "session_id", "phase", "trial_num", "force", "scale_rating", "is_practice", "vas_rt", "surpassed_goal"), "er_trials")
require_columns(dm_trials_raw, c("subject", "session_id", "phase", "trial_num", "effort1", "effort2", "reward1", "reward2", "option_chosen", "is_practice", "surpassed_goal", "choice_rt"), "dm_trials")
require_columns(clinical_raw, c("subject", "session_id", "phase", "DRSP_Score", "AMI_Score"), "clinical_subject_phase")
require_columns(ema_raw, c("subject", "session_id", "phase", "arousal_mean_z", "arousal_sd_z", "valence_mean_z", "valence_sd_z"), "ema_subject_phase")

hrv_path <- file.path(raw_dir, "hrv_subject_phase.csv")
ecg_index_path <- file.path(raw_dir, "ecg_index.csv")

if (file.exists(hrv_path)) {
  hrv_raw <- read_csv_clean(hrv_path)
  require_columns(
    hrv_raw,
    c("subject", "session_id", "phase", "hrv_baseline", "hrv_active", "rr_clean_n_baseline", "rr_clean_n_active", "ecg_ok_baseline", "ecg_ok_active"),
    "hrv_subject_phase"
  )
} else if (file.exists(ecg_index_path)) {
  ecg_index <- read_csv_clean(ecg_index_path)
  hrv_raw <- build_hrv_table(ecg_index)
  write_csv(hrv_raw, file.path(processed_dir, "hrv_subject_phase_from_ecg.csv"))
} else {
  stop("Provide either data/raw/hrv_subject_phase.csv or data/raw/ecg_index.csv.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Cohorts
# -----------------------------------------------------------------------------

strict_phase_subjects <- phase_schedule %>%
  filter(phase_verified == TRUE, phase %in% phase_levels, phase_target %in% phase_levels) %>%
  distinct(subject, phase) %>%
  count(subject, name = "n_phase") %>%
  filter(n_phase == 2) %>%
  pull(subject)

# -----------------------------------------------------------------------------
# EP
# -----------------------------------------------------------------------------

er_trials <- er_trials_raw %>%
  filter(subject %in% strict_phase_subjects, is_practice == 0) %>%
  mutate(
    subject = factor(subject),
    phase = factor(phase, levels = phase_levels),
    force_c = force - mean(force, na.rm = TRUE),
    trial_num_c = trial_num - mean(trial_num, na.rm = TRUE)
  ) %>%
  left_join(clinical_raw, by = c("subject", "session_id", "phase")) %>%
  left_join(ema_raw, by = c("subject", "session_id", "phase"))

ep_qc <- er_trials %>%
  group_by(subject, phase) %>%
  summarise(force_rating_cor = suppressWarnings(cor(force, scale_rating, use = "complete.obs")), .groups = "drop") %>%
  group_by(subject) %>%
  summarise(mean_force_rating_cor = mean(force_rating_cor, na.rm = TRUE), .groups = "drop") %>%
  mutate(cor_z = z_score(mean_force_rating_cor))

ep_excluded_subjects <- ep_qc %>%
  filter(cor_z <= -ep_outlier_sd) %>%
  pull(subject)

er_model_data <- er_trials %>%
  filter(!subject %in% ep_excluded_subjects)

ep_subject_phase <- er_model_data %>%
  group_by(subject, phase) %>%
  summarise(
    effort_slope = coef(lm(scale_rating ~ force))[2],
    success_pct = 100 * mean(surpassed_goal, na.rm = TRUE),
    mean_rt = mean(vas_rt, na.rm = TRUE),
    .groups = "drop"
  )

ep_curve_subject <- er_model_data %>%
  group_by(subject, phase, force) %>%
  summarise(mean_rating = mean(scale_rating, na.rm = TRUE), .groups = "drop")

ep_curve_phase <- ep_curve_subject %>%
  group_by(phase, force) %>%
  summarise(
    mean_rating = mean(mean_rating, na.rm = TRUE),
    se = sd(mean_rating, na.rm = TRUE) / sqrt(sum(is.finite(mean_rating))),
    .groups = "drop"
  )

ep_primary_test <- ep_subject_phase %>%
  select(subject, phase, effort_slope) %>%
  pivot_wider(names_from = phase, values_from = effort_slope) %>%
  summarise(test = list(t.test(Follicular, Luteal, paired = TRUE))) %>%
  pull(test) %>%
  .[[1]]

fit_ep_moderator <- function(moderator, label) {
  dat <- er_model_data %>%
    filter(is.finite(.data[[moderator]])) %>%
    mutate(mod_z = z_score(.data[[moderator]]))

  fit <- lmer(
    scale_rating ~ force_c * phase * mod_z + trial_num_c * phase + (1 | subject),
    data = dat,
    REML = TRUE
  )

  tidy_lmm(fit, label)
}

ep_moderators <- bind_rows(
  fit_ep_moderator("DRSP_Score", "EP_DRSP"),
  fit_ep_moderator("arousal_mean_z", "EP_arousal_mean"),
  fit_ep_moderator("arousal_sd_z", "EP_arousal_sd"),
  fit_ep_moderator("valence_mean_z", "EP_valence_mean"),
  fit_ep_moderator("valence_sd_z", "EP_valence_sd")
)

hrv_clean <- hrv_raw %>%
  filter(subject %in% levels(er_model_data$subject)) %>%
  mutate(
    phase = factor(phase, levels = phase_levels),
    hrv_baseline_ok = ecg_ok_baseline == TRUE & rr_clean_n_baseline >= min_rr_baseline,
    hrv_active_ok = ecg_ok_active == TRUE & rr_clean_n_active >= min_rr_active,
    hrv_baseline_z = z_score(log(hrv_baseline)),
    hrv_active_z = z_score(log(hrv_active))
  )

ep_hrv_trials <- er_model_data %>%
  inner_join(
    hrv_clean %>%
      filter(hrv_baseline_ok, hrv_active_ok) %>%
      select(subject, session_id, phase, hrv_baseline_z, hrv_active_z),
    by = c("subject", "session_id", "phase")
  )

fit_ep_hrv <- function(moderator, label) {
  fit <- lmer(
    scale_rating ~ force_c * phase * .data[[moderator]] + trial_num_c * phase + (1 | subject),
    data = ep_hrv_trials,
    REML = TRUE
  )
  tidy_lmm(fit, label)
}

ep_hrv_models <- bind_rows(
  fit_ep_hrv("hrv_baseline_z", "EP_hrv_baseline"),
  fit_ep_hrv("hrv_active_z", "EP_hrv_active")
)

# -----------------------------------------------------------------------------
# DM
# -----------------------------------------------------------------------------

dm_trials <- dm_trials_raw %>%
  filter(subject %in% strict_phase_subjects, is_practice == 0) %>%
  mutate(
    subject = factor(subject),
    phase = factor(phase, levels = phase_levels)
  ) %>%
  left_join(clinical_raw, by = c("subject", "session_id", "phase")) %>%
  left_join(ema_raw, by = c("subject", "session_id", "phase"))

catch_trials <- dm_trials %>%
  filter((reward1 > reward2 & effort1 < effort2) | (reward2 > reward1 & effort2 < effort1)) %>%
  mutate(correct = case_when(
    reward1 > reward2 & effort1 < effort2 & option_chosen == 0 ~ 1L,
    reward2 > reward1 & effort2 < effort1 & option_chosen == 1 ~ 1L,
    TRUE ~ 0L
  )) %>%
  group_by(subject, phase) %>%
  summarise(n_errors = sum(correct == 0, na.rm = TRUE), .groups = "drop")

dm_excluded_subjects <- catch_trials %>%
  filter(n_errors >= 2) %>%
  distinct(subject) %>%
  pull(subject)

dm_trials_strict <- dm_trials %>%
  filter(!subject %in% dm_excluded_subjects)

dm_choice_space <- dm_trials_strict %>%
  mutate(
    hard_reward = ifelse(effort1 > effort2, reward1, reward2),
    hard_effort = ifelse(effort1 > effort2, effort1, effort2),
    choose_hard = ifelse(
      (effort1 > effort2 & option_chosen == 0) | (effort2 > effort1 & option_chosen == 1),
      1, 0
    )
  ) %>%
  group_by(subject, phase, hard_reward, hard_effort) %>%
  summarise(p_choose_hard = mean(choose_hard, na.rm = TRUE), .groups = "drop") %>%
  group_by(phase, hard_reward, hard_effort) %>%
  summarise(p_choose_hard = mean(p_choose_hard, na.rm = TRUE), n_total = n(), .groups = "drop") %>%
  pivot_wider(names_from = phase, values_from = c(p_choose_hard, n_total)) %>%
  mutate(
    diff_pp = 100 * (p_choose_hard_Luteal - p_choose_hard_Follicular),
    n_total = coalesce(n_total_Follicular, 0) + coalesce(n_total_Luteal, 0)
  ) %>%
  filter(is.finite(diff_pp))

dm_session_params <- read_csv_clean(file.path(raw_dir, "dm_session_parameters.csv")) %>%
  filter(subject %in% strict_phase_subjects, !subject %in% dm_excluded_subjects) %>%
  mutate(
    subject = factor(subject),
    phase = factor(phase, levels = phase_levels),
    effSens_z = z_score(effSens),
    rewSens_z = z_score(rewSens)
  ) %>%
  left_join(clinical_raw, by = c("subject", "session_id", "phase")) %>%
  left_join(ema_raw, by = c("subject", "session_id", "phase")) %>%
  left_join(
    hrv_clean %>% select(subject, session_id, phase, hrv_baseline_z, hrv_active_z),
    by = c("subject", "session_id", "phase")
  )

require_columns(dm_session_params, c("effSens", "rewSens"), "dm_session_parameters")

dm_phase_summary <- read_csv_clean(file.path(raw_dir, "dm_phase_population_effect.csv"))
require_columns(dm_phase_summary, c("term", "estimate", "ci_low", "ci_high"), "dm_phase_population_effect")

fit_dm_moderator <- function(outcome, moderator, label) {
  dat <- dm_session_params %>%
    filter(is.finite(.data[[outcome]]), is.finite(.data[[moderator]])) %>%
    mutate(
      outcome_z = z_score(.data[[outcome]]),
      mod_z = z_score(.data[[moderator]])
    )

  fit <- lmer(
    outcome_z ~ phase * mod_z + (1 | subject),
    data = dat,
    REML = TRUE
  )

  tidy_lmm(fit, label)
}

dm_moderators <- bind_rows(
  fit_dm_moderator("effSens", "DRSP_Score", "DM_eff_DRSP"),
  fit_dm_moderator("effSens", "arousal_mean_z", "DM_eff_arousal_mean"),
  fit_dm_moderator("effSens", "arousal_sd_z", "DM_eff_arousal_sd"),
  fit_dm_moderator("effSens", "valence_mean_z", "DM_eff_valence_mean"),
  fit_dm_moderator("effSens", "valence_sd_z", "DM_eff_valence_sd"),
  fit_dm_moderator("effSens", "hrv_baseline_z", "DM_eff_hrv_baseline"),
  fit_dm_moderator("effSens", "hrv_active_z", "DM_eff_hrv_active"),
  fit_dm_moderator("rewSens", "DRSP_Score", "DM_rew_DRSP"),
  fit_dm_moderator("rewSens", "arousal_mean_z", "DM_rew_arousal_mean"),
  fit_dm_moderator("rewSens", "arousal_sd_z", "DM_rew_arousal_sd"),
  fit_dm_moderator("rewSens", "valence_mean_z", "DM_rew_valence_mean"),
  fit_dm_moderator("rewSens", "valence_sd_z", "DM_rew_valence_sd"),
  fit_dm_moderator("rewSens", "hrv_baseline_z", "DM_rew_hrv_baseline"),
  fit_dm_moderator("rewSens", "hrv_active_z", "DM_rew_hrv_active")
)

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------

write_csv(er_model_data, file.path(processed_dir, "er_model_data.csv"))
write_csv(ep_subject_phase, file.path(processed_dir, "ep_subject_phase.csv"))
write_csv(ep_moderators, file.path(models_dir, "ep_moderators.csv"))
write_csv(ep_hrv_models, file.path(models_dir, "ep_hrv_models.csv"))
write_csv(hrv_clean, file.path(processed_dir, "hrv_subject_phase_clean.csv"))

write_csv(dm_trials_strict, file.path(processed_dir, "dm_trials_strict.csv"))
write_csv(dm_session_params, file.path(processed_dir, "dm_session_params.csv"))
write_csv(dm_choice_space, file.path(processed_dir, "dm_choice_space.csv"))
write_csv(dm_moderators, file.path(models_dir, "dm_moderators.csv"))

analysis_counts <- tibble(
  metric = c(
    "strict_phase_subjects",
    "ep_subjects_after_qc",
    "ep_hrv_subjects_after_qc",
    "dm_subjects_after_catch_qc"
  ),
  value = c(
    length(strict_phase_subjects),
    n_distinct(er_model_data$subject),
    n_distinct(ep_hrv_trials$subject),
    n_distinct(dm_trials_strict$subject)
  )
)

write_csv(analysis_counts, file.path(tables_dir, "analysis_counts.csv"))

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------

ep_curve_plot <- ggplot(ep_curve_subject, aes(x = force, y = mean_rating, group = interaction(subject, phase), color = phase)) +
  geom_line(alpha = 0.08, linewidth = 0.40, show.legend = FALSE) +
  geom_ribbon(
    data = ep_curve_phase,
    aes(x = force, ymin = mean_rating - 1.96 * se, ymax = mean_rating + 1.96 * se, fill = phase, group = phase),
    inherit.aes = FALSE,
    alpha = 0.18,
    colour = NA
  ) +
  geom_line(data = ep_curve_phase, aes(group = phase), inherit.aes = FALSE, linewidth = 1.35) +
  geom_point(data = ep_curve_phase, size = 2.1, inherit.aes = FALSE) +
  scale_color_manual(values = c(Follicular = col_fol, Luteal = col_lut), labels = phase_labels) +
  scale_fill_manual(values = c(Follicular = col_fol_soft, Luteal = col_lut_soft), labels = phase_labels) +
  scale_x_continuous(breaks = sort(unique(ep_curve_phase$force))) +
  labs(x = "Required force (%MVC)", y = "Subjective effort rating") +
  theme_pub()

ep_slope_plot <- ggplot(ep_subject_phase, aes(x = phase, y = effort_slope, group = subject)) +
  geom_line(color = "grey82", linewidth = 0.7, alpha = 0.9) +
  geom_point(aes(color = phase), size = 2.2, alpha = 0.85) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 1.4) +
  stat_summary(aes(group = 1), fun.data = mean_cl_normal, geom = "errorbar", width = 0.08, color = "black", linewidth = 0.9) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 21, fill = "white", color = "black", size = 3.4, stroke = 1.1) +
  scale_color_manual(values = c(Follicular = col_fol, Luteal = col_lut), labels = phase_labels) +
  scale_x_discrete(labels = phase_labels) +
  labs(x = NULL, y = "Effort differentiation slope") +
  theme_pub()

dm_eff_fol_draw <- boot_mean_draws(dm_session_params$effSens[dm_session_params$phase == "Follicular"], seed = seed_base)
dm_eff_lut_draw <- boot_mean_draws(dm_session_params$effSens[dm_session_params$phase == "Luteal"], seed = seed_base + 1L)
dm_rew_fol_draw <- boot_mean_draws(dm_session_params$rewSens[dm_session_params$phase == "Follicular"], seed = seed_base + 2L)
dm_rew_lut_draw <- boot_mean_draws(dm_session_params$rewSens[dm_session_params$phase == "Luteal"], seed = seed_base + 3L)

dm_eff_delta_draw <- summary_to_draws(
  mean = dm_phase_summary$estimate[dm_phase_summary$term == "effort_sensitivity_phase_difference"],
  lower = dm_phase_summary$ci_low[dm_phase_summary$term == "effort_sensitivity_phase_difference"],
  upper = dm_phase_summary$ci_high[dm_phase_summary$term == "effort_sensitivity_phase_difference"],
  seed = seed_base + 4L
)

dm_rew_delta_draw <- summary_to_draws(
  mean = dm_phase_summary$estimate[dm_phase_summary$term == "reward_sensitivity_phase_difference"],
  lower = dm_phase_summary$ci_low[dm_phase_summary$term == "reward_sensitivity_phase_difference"],
  upper = dm_phase_summary$ci_high[dm_phase_summary$term == "reward_sensitivity_phase_difference"],
  seed = seed_base + 5L
)

dm_phase_draws <- bind_rows(
  tibble(term = "Late-follicular mean effort sensitivity", family = "follicular", value = dm_eff_fol_draw),
  tibble(term = "Mid-luteal mean effort sensitivity", family = "luteal", value = dm_eff_lut_draw),
  tibble(term = "Late-follicular mean reward sensitivity", family = "follicular", value = dm_rew_fol_draw),
  tibble(term = "Mid-luteal mean reward sensitivity", family = "luteal", value = dm_rew_lut_draw),
  tibble(term = "Mid-luteal - late-follicular effort sensitivity", family = "contrast", value = dm_eff_delta_draw),
  tibble(term = "Mid-luteal - late-follicular reward sensitivity", family = "contrast", value = dm_rew_delta_draw)
)

dm_posterior_plot <- ggplot(dm_phase_draws, aes(x = value, y = term, fill = family)) +
  ggdist::stat_gradientinterval(.width = c(.9, .5), slab_size = 1) +
  scale_fill_manual(values = c(follicular = col_fol, luteal = col_lut, contrast = col_contrast)) +
  labs(x = "Posterior parameter estimate", y = NULL) +
  theme_pub()

save_plot_pair(tag_panel(ep_curve_plot, "A") | tag_panel(ep_slope_plot, "B"), "ep_core_panels", 11, 4.8)
save_plot_pair(dm_posterior_plot, "dm_posterior_panels", 8.5, 5.2)

message("Pipeline completed.")
