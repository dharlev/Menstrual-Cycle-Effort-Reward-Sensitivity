# Configuration and shared analysis functions
INPUT_DIR <- "data"
OUTPUT_DIR <- "outputs"
TRIALS_PER_SESSION <- 44L
PHASES <- c("Follicular", "Luteal")
EMA_VARIABLES <- c("arousal_mean", "arousal_sd", "valence_mean", "valence_sd")
QUESTIONNAIRES <- c("PHQ9_Score", "GAD7_Score", "AMI_Score", "CIRENS_Score")
PRIMARY_STAN_FILE <- "dm_hierarchical_phase.stan"
VISIT_STAN_FILE <- "dm_hierarchical_phase_visit.stan"
PRIMARY_SAMPLING <- list(
  chains = 4L, parallel_chains = 4L,
  iter_warmup = 1000L, iter_sampling = 1000L,
  seed = 123L, adapt_delta = 0.95, max_treedepth = 12L
)
VISIT_SAMPLING <- list(
  chains = 4L, parallel_chains = 4L,
  iter_warmup = 1000L, iter_sampling = 1000L,
  seed = 20260910L, adapt_delta = 0.99, max_treedepth = 12L
)

require_condition <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

require_columns <- function(data, columns) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) {
    stop("Input is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  }
}

read_analysis_table <- function(input_dir, filename, columns,
                                key = c("participant_id", "phase")) {
  path <- file.path(input_dir, filename)
  require_condition(file.exists(path), paste("Missing input:", path))
  x <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE,
                       col_types = readr::cols(participant_id = readr::col_character()))
  require_columns(x, columns)
  require_condition(nrow(readr::problems(x)) == 0L, paste("CSV parsing error:", filename))
  x <- x[, columns, drop = FALSE]
  require_condition(nrow(x) > 0L, paste("Empty input:", filename))
  require_condition(!anyNA(x$participant_id) && all(nzchar(trimws(x$participant_id))),
                    paste("Blank participant ID:", filename))
  if ("phase" %in% names(x)) {
    x$phase <- dplyr::recode(x$phase, LF = "Follicular", ML = "Luteal")
    require_condition(!anyNA(x$phase) && all(x$phase %in% PHASES),
                      paste("Invalid phase:", filename))
    x$phase <- factor(x$phase, levels = PHASES)
  }
  if (length(key)) require_condition(!anyDuplicated(x[key]), paste("Duplicate key:", filename))
  x
}

z_score <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  require_condition(is.finite(s) && s > 0, "Cannot standardise a constant/empty variable.")
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

write_table <- function(x, output_dir, name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(as.data.frame(x), file.path(output_dir, paste0(name, ".csv")))
}

paired_summary <- function(data, variable, label = variable) {
  w <- data |>
    dplyr::select(participant_id, phase, dplyr::all_of(variable)) |>
    tidyr::pivot_wider(names_from = phase, values_from = dplyr::all_of(variable))
  require_columns(w, PHASES)
  w <- w[is.finite(w$Follicular) & is.finite(w$Luteal), ]
  require_condition(nrow(w) >= 2L, paste("Insufficient pairs:", label))
  test <- stats::t.test(w$Follicular, w$Luteal, paired = TRUE)
  overall <- (w$Follicular + w$Luteal) / 2
  data.frame(measure = label, n = nrow(w), overall_mean = mean(overall),
    overall_sd = stats::sd(overall), LF_mean = mean(w$Follicular),
    LF_sd = stats::sd(w$Follicular), LF_min = min(w$Follicular), LF_max = max(w$Follicular),
    ML_mean = mean(w$Luteal), ML_sd = stats::sd(w$Luteal),
    ML_min = min(w$Luteal), ML_max = max(w$Luteal),
    LF_minus_ML = unname(test$estimate), CI_low = test$conf.int[1],
    CI_high = test$conf.int[2], t = unname(test$statistic),
    df = unname(test$parameter), raw_p = test$p.value)
}

fit_mixed_model <- function(formula, data, name, output_dir, reml = TRUE,
                             sum_coding = FALSE, bobyqa = TRUE) {
  data$participant_id <- factor(data$participant_id)
  data$phase <- factor(data$phase, levels = PHASES)
  contrasts(data$phase) <- if (sum_coding) stats::contr.sum(2) else stats::contr.treatment(2)
  # Named treatment contrasts keep LF/ML coefficient names readable.
  if (!sum_coding) contrasts(data$phase) <- stats::contr.treatment(PHASES)
  control <- if (bobyqa) lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)) else lme4::lmerControl()
  fit <- lmerTest::lmer(formula, data = data, REML = reml, control = control,
                        na.action = stats::na.fail)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, file.path(output_dir, paste0(name, ".rds")))
  coefficients <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) |>
    dplyr::mutate(model = name, n_participants = dplyr::n_distinct(data$participant_id),
                  n_sessions = dplyr::n_distinct(paste(data$participant_id, data$phase)),
                  n_rows = nrow(data), REML = reml, singular = lme4::isSingular(fit))
  write_table(coefficients, output_dir, paste0(name, "_coefficients"))
  writeLines(capture.output(summary(fit)), file.path(output_dir, paste0(name, "_summary.txt")))
  list(fit = fit, coefficients = coefficients, data = data)
}

read_secondary_inputs <- function(input_dir) {
  ep_columns <- c("participant_id", "phase", "trial_num", "force", "scale_rating")
  ep_phase <- read_analysis_table(input_dir, "ep_phase_trials.csv",
    c(ep_columns, "vas_rt", "surpassed_goal"), c("participant_id", "phase", "trial_num"))
  ep_mod <- read_analysis_table(input_dir, "ep_moderator_trials.csv",
    c(ep_columns, "AMI_Score"), c("participant_id", "phase", "trial_num"))
  for (x in list(ep_phase, ep_mod)) {
    require_condition(all(is.finite(x$force)) && all(is.finite(x$scale_rating)) &&
      all(is.finite(x$trial_num)), "EP force, rating and trial number must be finite.")
    require_condition(all(x$force %in% c(10, 30, 50, 70)) &&
      all(x$scale_rating >= 0 & x$scale_rating <= 100), "Invalid EP force or rating scale.")
    counts <- dplyr::count(x, participant_id, phase)
    require_condition(all(counts$n == 16L), "EP inputs require 16 test trials per session.")
  }
  paired_counts <- dplyr::count(dplyr::distinct(ep_phase, participant_id, phase), participant_id)
  require_condition(all(paired_counts$n == 2L), "EP phase input requires two phases per participant.")
  ep_hrv <- read_analysis_table(input_dir, "ep_hrv_sessions.csv",
    c("participant_id", "phase", "hrv_baseline", "hrv_mvc", "rr_n_baseline", "rr_n_mvc",
      "baseline_ok", "mvc_ok"))
  ep_reference <- read_analysis_table(input_dir, "ep_ema_reference.csv", "participant_id", "participant_id")
  dm_covariates <- read_analysis_table(input_dir, "dm_moderator_sessions.csv",
    c("participant_id", "phase", "effSens", "rewSens", "hrv_baseline", "hrv_mvc", "AMI_Score"))
  questionnaires <- read_analysis_table(input_dir, "questionnaires.csv",
    c("participant_id", "phase", QUESTIONNAIRES))
  ema <- read_analysis_table(input_dir, "ema_observations.csv",
    c("participant_id", "observation_id", "day_before_LF", "day_before_ML",
      "state_energy", "state_fatigue", "state_mood", "state_anxiety"),
    c("participant_id", "observation_id"))
  items <- unlist(ema[c("state_energy", "state_fatigue", "state_mood", "state_anxiety")])
  require_condition(all(is.na(items) | (is.finite(items) & items >= 0 & items <= 1)),
                    "EMA items must use the 0-1 scale, with NA for missing values.")
  list(ep_phase = ep_phase, ep_mod = ep_mod, ep_hrv = ep_hrv,
       ep_reference = ep_reference, dm_covariates = dm_covariates,
       questionnaires = questionnaires, ema = ema)
}

check_packages <- function() {
  packages <- c("cmdstanr", "posterior", "dplyr", "tidyr", "readr", "lme4",
                "lmerTest", "broom.mixed", "ggplot2", "patchwork", "png")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  require_condition(!length(missing), paste("Install required packages:", paste(missing, collapse = ", ")))
  invisible(packages)
}
