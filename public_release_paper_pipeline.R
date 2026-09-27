#!/usr/bin/env Rscript
# Menstrual-cycle phase, effort perception, and effort-reward decisions
#
# Run: Rscript public_release_paper_pipeline.R
# Input: de-identified analysis tables in data/raw/ (schemas in README.md).
# Output: fitted models, numerical tables, and figures in outputs/.
# Sections 1-7 specify the decision model; sections 8-13 prepare covariates,
# fit the remaining analyses, and produce tables and figures.

# ============================================================================
# 1. Analysis configuration
# ============================================================================

OUTPUT_DIR <- "outputs"
TRIALS_PER_SESSION <- 44L

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

# ============================================================================
# 2. Input checks, session indexing, and offer normalisation
# ============================================================================
# Each trial must have the same pair of offers across sessions.
# Reward is divided by the larger reward WITHIN each offered pair.
# Effort is divided by 100; it is squared later in the Stan value function.

prepare_choice_data <- function(input_file, include_visit = TRUE) {
  x <- readr::read_csv(
    input_file, show_col_types = FALSE, progress = FALSE,
    col_types = readr::cols(participant_id = readr::col_character(),
                           phase = readr::col_character())
  )
  require_condition(nrow(readr::problems(x)) == 0L,
                    "Input contains CSV parsing errors; check column types.")
  required <- c("participant_id", "phase", "trial_num", "option_chosen",
                "reward1", "reward2", "effort1", "effort2")
  if (include_visit) {
    required <- c(required, "visit", "is_force_exerted", "surpassed_goal")
  }
  require_columns(x, required)
  x <- x[, required]
  x$phase <- dplyr::recode(x$phase, LF = "Follicular", ML = "Luteal")
  choice_columns <- c("participant_id", "phase", "trial_num", "option_chosen",
                      "reward1", "reward2", "effort1", "effort2")
  require_condition(nrow(x) > 0L && !anyNA(x[choice_columns]),
                    "Choice input is empty or has missing required values.")
  require_condition(all(nzchar(trimws(x$participant_id))), "Participant IDs cannot be blank.")
  require_condition(all(x$phase %in% c("Follicular", "Luteal")),
                    "phase must be LF/Follicular or ML/Luteal.")
  numeric_columns <- c("trial_num", "option_chosen", "reward1", "reward2", "effort1", "effort2")
  require_condition(all(vapply(x[numeric_columns], is.numeric, logical(1))),
                    "Trial numbers, choices, rewards, and efforts must be numeric.")
  require_condition(all(is.finite(as.matrix(x[numeric_columns]))),
                    "Numeric choice input must be finite.")
  require_condition(all(x$option_chosen %in% c(0, 1)), "Choices must be coded 0 or 1.")
  require_condition(all(x$trial_num %in% seq_len(TRIALS_PER_SESSION)),
                    "Trial numbers must be integers from 1 through 44.")
  require_condition(all(x$reward1 >= 0 & x$reward2 >= 0), "Rewards cannot be negative.")
  require_condition(all(x$effort1 >= 0 & x$effort1 <= 100 &
                        x$effort2 >= 0 & x$effort2 <= 100),
                    "Effort must be recorded as percent MVC, from 0 to 100.")
  require_condition(!anyDuplicated(x[c("participant_id", "phase", "trial_num")]),
                    "Participant/phase/trial rows must be unique.")

  x$phase <- factor(x$phase, levels = c("Follicular", "Luteal"))
  counts <- dplyr::count(x, participant_id, phase, name = "n_trials")
  require_condition(all(counts$n_trials == TRIALS_PER_SESSION) &&
                      all(table(counts$participant_id) == 2L),
                    "Each participant must have two complete phase sessions.")

  subjects <- x |>
    dplyr::distinct(participant_id) |>
    dplyr::arrange(participant_id) |>
    dplyr::mutate(subj_idx = dplyr::row_number())
  sessions <- x |>
    dplyr::distinct(participant_id, phase) |>
    dplyr::left_join(subjects, by = "participant_id") |>
    dplyr::arrange(participant_id, phase) |>
    dplyr::mutate(session_idx = dplyr::row_number(),
                  phase01 = as.integer(phase == "Luteal"))

  template <- x |>
    dplyr::distinct(trial_num, reward1, reward2, effort1, effort2) |>
    dplyr::arrange(trial_num)
  require_condition(nrow(template) == TRIALS_PER_SESSION,
                    "Offer pairs must be identical across sessions at each trial number.")
  pair_max <- pmax(template$reward1, template$reward2)
  denominator <- ifelse(pair_max == 0, 1, pair_max)
  template$rew1 <- template$reward1 / denominator
  template$rew2 <- template$reward2 / denominator
  template$eff1 <- template$effort1 / 100
  template$eff2 <- template$effort2 / 100
  indexed <- dplyr::left_join(
    x, dplyr::select(sessions, participant_id, phase, session_idx),
    by = c("participant_id", "phase")
  )
  y <- matrix(NA_integer_, nrow = nrow(sessions), ncol = TRIALS_PER_SESSION)
  y[cbind(indexed$session_idx, indexed$trial_num)] <- as.integer(indexed$option_chosen)
  require_condition(!anyNA(y), "Choice matrix is incomplete.")

  stan_data <- list(
    N_subj = nrow(subjects), N_sess = nrow(sessions), T = TRIALS_PER_SESSION,
    subj = as.integer(sessions$subj_idx), phase = sessions$phase01, y = y,
    rew1 = as.array(template$rew1), rew2 = as.array(template$rew2),
    eff1 = as.array(template$eff1), eff2 = as.array(template$eff2)
  )
  if (include_visit) {
    require_condition(!anyNA(x$visit) && all(x$visit %in% c(1, 2)),
                      "Chronological visit must be 1 or 2 on every trial.")
    visits <- dplyr::distinct(x, participant_id, phase, visit)
    require_condition(nrow(visits) == nrow(sessions) &&
                        !anyDuplicated(visits[c("participant_id", "visit")]),
                      "Each phase session must map to exactly one physical visit.")
    sessions <- dplyr::left_join(sessions, visits, by = c("participant_id", "phase")) |>
      dplyr::arrange(session_idx)
  }
  list(trials = x, sessions = sessions, template = template, stan_data = stan_data)
}

# ============================================================================
# 3. Primary hierarchical Bayesian choice model
# ============================================================================
# For each option: V = reward sensitivity * normalised reward
#                    - effort sensitivity * squared normalised effort.
# P(choose option 2) = inverse-logit(V2 - V1).
# Log sensitivities have population intercepts and ML-minus-LF phase effects.
# Four correlated participant effects represent reward/effort intercepts and
# reward/effort phase slopes. Priors are specified in the Stan model below.
# LF = 0 and ML = 1. Every test choice enters the likelihood, including catches.

PRIMARY_STAN_CODE <- r"(
data {
  int<lower=1> N_subj;
  int<lower=1> N_sess;
  int<lower=1> T;

  array[N_sess] int<lower=1, upper=N_subj> subj;
  array[N_sess] int<lower=0, upper=1> phase;

  array[N_sess, T] int<lower=0, upper=1> y;

  array[T] real rew1;
  array[T] real rew2;
  array[T] real eff1;
  array[T] real eff2;
}

parameters {
  vector[2] mu_base;                 // baseline means: [reward, effort]
  vector[2] mu_phase;                // phase means:    [reward, effort]

  vector<lower=0>[4] sigma_subj;     // SDs for 4 subject-level effects
  cholesky_factor_corr[4] L_subj;    // correlations between them
  matrix[4, N_subj] z_subj;          // subject latent z-scores
}

transformed parameters {
  matrix[4, N_subj] subj_eff;

  vector[N_sess] rewSens;
  vector[N_sess] effSens;

  subj_eff = diag_pre_multiply(sigma_subj, L_subj) * z_subj;

  for (n in 1:N_sess) {
    int s = subj[n];
    real ph = phase[n];

    real theta_rew =
      mu_base[1] +
      subj_eff[1, s] +
      (mu_phase[1] + subj_eff[3, s]) * ph;

    real theta_eff =
      mu_base[2] +
      subj_eff[2, s] +
      (mu_phase[2] + subj_eff[4, s]) * ph;

    rewSens[n] = exp(theta_rew);
    effSens[n] = exp(theta_eff);
  }
}

model {
  mu_base ~ normal(0, 1);
  mu_phase ~ normal(0, 0.5);

  sigma_subj ~ exponential(1);
  L_subj ~ lkj_corr_cholesky(2.0);

  to_vector(z_subj) ~ std_normal();

  for (n in 1:N_sess) {
    for (t in 1:T) {
      real V1;
      real V2;

      V1 = rewSens[n] * rew1[t] - effSens[n] * square(eff1[t]);
      V2 = rewSens[n] * rew2[t] - effSens[n] * square(eff2[t]);

      y[n, t] ~ bernoulli_logit(V2 - V1);
    }
  }
}

generated quantities {
  real delta_rew = mu_phase[1];
  real delta_eff = mu_phase[2];

  real delta_rew_exp = exp(mu_phase[1]);
  real delta_eff_exp = exp(mu_phase[2]);

  corr_matrix[4] Omega_subj;
  vector[N_sess] rewSens_sess;
  vector[N_sess] effSens_sess;

  real ybar_rep = 0;
  real ybar_obs = 0;
  int total = 0;

  Omega_subj = multiply_lower_tri_self_transpose(L_subj);

  for (n in 1:N_sess) {
    rewSens_sess[n] = rewSens[n];
    effSens_sess[n] = effSens[n];
  }

  for (n in 1:N_sess) {
    for (t in 1:T) {
      real V1;
      real V2;
      real p;
      int y_rep;

      V1 = rewSens[n] * rew1[t] - effSens[n] * square(eff1[t]);
      V2 = rewSens[n] * rew2[t] - effSens[n] * square(eff2[t]);

      p = inv_logit(V2 - V1);
      y_rep = bernoulli_rng(p);

      ybar_rep += y_rep;
      ybar_obs += y[n, t];
      total += 1;
    }
  }

  ybar_rep /= total;
  ybar_obs /= total;
}

)"

# ============================================================================
# 4. Chronological-visit extension of the hierarchical model
# ============================================================================
# The likelihood, phase effects, participant effects, and their priors are the
# same as above. A population visit-2 coefficient is added to each log
# sensitivity, with a Normal(0, 0.5) prior. visit2 = physical visit minus one.

VISIT_STAN_CODE <- r"(
data {
  int<lower=1> N_subj;
  int<lower=1> N_sess;
  int<lower=1> T;

  array[N_sess] int<lower=1, upper=N_subj> subj;
  array[N_sess] int<lower=0, upper=1> phase;

  array[N_sess, T] int<lower=0, upper=1> y;
  array[N_sess] int<lower=0,upper=1> visit2;

  array[T] real rew1;
  array[T] real rew2;
  array[T] real eff1;
  array[T] real eff2;
}

parameters {
  vector[2] mu_base;                 // baseline means: [reward, effort]
  vector[2] mu_visit;
  vector[2] mu_phase;                // phase means:    [reward, effort]

  vector<lower=0>[4] sigma_subj;     // SDs for 4 subject-level effects
  cholesky_factor_corr[4] L_subj;    // correlations between them
  matrix[4, N_subj] z_subj;          // subject latent z-scores
}

transformed parameters {
  matrix[4, N_subj] subj_eff;

  vector[N_sess] rewSens;
  vector[N_sess] effSens;

  subj_eff = diag_pre_multiply(sigma_subj, L_subj) * z_subj;

  for (n in 1:N_sess) {
    int s = subj[n];
    real ph = phase[n];

    real theta_rew =
      mu_base[1] +
      subj_eff[1, s] +
      (mu_phase[1] + subj_eff[3, s]) * ph + mu_visit[1]*visit2[n];

    real theta_eff =
      mu_base[2] +
      subj_eff[2, s] +
      (mu_phase[2] + subj_eff[4, s]) * ph + mu_visit[2]*visit2[n];

    rewSens[n] = exp(theta_rew);
    effSens[n] = exp(theta_eff);
  }
}

model {
  mu_base ~ normal(0, 1);
  mu_phase ~ normal(0, 0.5);
  mu_visit ~ normal(0,0.5);

  sigma_subj ~ exponential(1);
  L_subj ~ lkj_corr_cholesky(2.0);

  to_vector(z_subj) ~ std_normal();

  for (n in 1:N_sess) {
    for (t in 1:T) {
      real V1;
      real V2;

      V1 = rewSens[n] * rew1[t] - effSens[n] * square(eff1[t]);
      V2 = rewSens[n] * rew2[t] - effSens[n] * square(eff2[t]);

      y[n, t] ~ bernoulli_logit(V2 - V1);
    }
  }
}

generated quantities {
  real delta_rew = mu_phase[1];
  real delta_eff = mu_phase[2];

  real delta_rew_exp = exp(mu_phase[1]);
  real delta_eff_exp = exp(mu_phase[2]);

  corr_matrix[4] Omega_subj;
  vector[N_sess] rewSens_sess;
  vector[N_sess] effSens_sess;

  real ybar_rep = 0;
  real ybar_obs = 0;
  int total = 0;

  Omega_subj = multiply_lower_tri_self_transpose(L_subj);

  for (n in 1:N_sess) {
    rewSens_sess[n] = rewSens[n];
    effSens_sess[n] = effSens[n];
  }

  for (n in 1:N_sess) {
    for (t in 1:T) {
      real V1;
      real V2;
      real p;
      int y_rep;

      V1 = rewSens[n] * rew1[t] - effSens[n] * square(eff1[t]);
      V2 = rewSens[n] * rew2[t] - effSens[n] * square(eff2[t]);

      p = inv_logit(V2 - V1);
      y_rep = bernoulli_rng(p);

      ybar_rep += y_rep;
      ybar_obs += y[n, t];
      total += 1;
    }
  }

  ybar_rep /= total;
  ybar_obs /= total;
}

)"

# ============================================================================
# 5. Model fitting and diagnostics
# ============================================================================
# Both model specifications are embedded above. The script writes them to the
# output directory, compiles them, and samples the posterior from the choices.
# Posterior summaries and sampler diagnostics are saved for assessment.
# Diagnostics never change the observations, seed, or model automatically.

fit_choice_model <- function(stan_code, stan_data, sampling, output_dir, model_name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  model_file <- file.path(output_dir, paste0(model_name, ".stan"))
  writeLines(stan_code, model_file, useBytes = TRUE)
  cmdstanr::write_stan_json(stan_data, file.path(output_dir, "stan_data.json"))
  model <- cmdstanr::cmdstan_model(model_file)
  fit <- do.call(model$sample, c(
    list(data = stan_data, output_dir = output_dir, output_basename = model_name,
         refresh = 500L), sampling
  ))
  fit$save_object(file.path(output_dir, "fit.rds"))
  readr::write_csv(fit$summary(), file.path(output_dir, "posterior_diagnostics.csv"))
  diagnostics <- fit$diagnostic_summary()
  writeLines(capture.output(diagnostics), file.path(output_dir, "sampler_diagnostics.txt"))
  readr::write_csv(as.data.frame(sampling), file.path(output_dir, "sampling_settings.csv"))
  invisible(fit)
}

# ============================================================================
# 6. Posterior contrasts and participant/session estimates
# ============================================================================
# Population contrasts are on the log-sensitivity scale, ML minus LF.
# Intervals are equal-tailed 95% credible intervals; prob_positive is the
# proportion of actual posterior draws above zero. Session sensitivities are
# posterior means on the positive sensitivity scale, not log sensitivities.

summarise_draws <- function(draws, parameters) {
  do.call(rbind, lapply(parameters, function(parameter) {
    values <- draws[, parameter]
    data.frame(
      parameter = parameter,
      mean = mean(values),
      lower = unname(stats::quantile(values, 0.025)),
      upper = unname(stats::quantile(values, 0.975)),
      prob_positive = mean(values > 0)
    )
  }))
}

export_choice_results <- function(fit, sessions, output_dir, include_visit = FALSE) {
  parameters <- c("delta_eff", "delta_rew")
  if (include_visit) parameters <- c(parameters, "mu_visit[1]", "mu_visit[2]")
  draws <- fit$draws(variables = c("delta_eff", "delta_rew",
                                  if (include_visit) "mu_visit",
                                  "effSens_sess", "rewSens_sess"), format = "matrix")
  contrasts <- summarise_draws(draws, parameters)
  contrasts$contrast <- c("phase_effort_ML_minus_LF", "phase_reward_ML_minus_LF",
                         if (include_visit) c("visit_reward_2_minus_1", "visit_effort_2_minus_1"))
  estimates <- sessions
  estimates$effSens <- vapply(sessions$session_idx, function(i) {
    mean(draws[, paste0("effSens_sess[", i, "]")])
  }, numeric(1))
  estimates$rewSens <- vapply(sessions$session_idx, function(i) {
    mean(draws[, paste0("rewSens_sess[", i, "]")])
  }, numeric(1))
  readr::write_csv(contrasts, file.path(output_dir, "population_contrasts.csv"))
  readr::write_csv(estimates, file.path(output_dir, "session_parameters.csv"))
  readr::write_csv(as.data.frame(draws[, parameters, drop = FALSE]),
                   file.path(output_dir, "contrast_draws.csv"))
  print(contrasts, digits = 6)
  invisible(contrasts)
}

# ============================================================================
# 7. Behavioural chronological-visit models
# ============================================================================
# Each participant/session contributes:
#   success_pct = 100 * successful executed trials / recorded executed trials
#   hard_pct    = 100 * choices of the higher-effort option / all test choices
# Missing or unrecognised execution flags are not counted as executed. Their
# choices remain in the choice model. Success must be known on executed trials.
# Equal-effort offers, if present, do not count as higher-effort choices.
# Both outcomes use outcome ~ phase + visit2 + (1 | participant).
# Models use REML, Satterthwaite degrees of freedom and two-sided t inference;
# confidence intervals are 95% Wald intervals from broom.mixed::tidy.
# Coefficients for these outcomes are percentage-point differences.

parse_binary_flag <- function(values) {
  values <- tolower(as.character(values))
  ifelse(values %in% c("1", "true"), 1L,
         ifelse(values %in% c("0", "false"), 0L, NA_integer_))
}

prepare_behaviour_data <- function(trials) {
  x <- trials
  x$executed <- parse_binary_flag(x$is_force_exerted)
  x$success <- parse_binary_flag(x$surpassed_goal)
  require_condition(!anyNA(x$success[x$executed %in% 1L]),
                    "Success must be recorded on every executed trial.")
  x$hard <- ifelse(x$option_chosen == 0, x$effort1 > x$effort2, x$effort2 > x$effort1)
  b <- x |>
    dplyr::group_by(participant_id, phase, visit) |>
    dplyr::summarise(
      n_choices = dplyr::n(),
      n_executed = sum(executed %in% 1L),
      n_execution_unknown = sum(is.na(executed)),
      n_successful = sum(success[executed %in% 1L]),
      success_pct = 100 * mean(success[executed %in% 1L]),
      hard_pct = 100 * mean(hard), .groups = "drop"
    ) |>
    dplyr::mutate(participant_id = factor(participant_id),
                  phase = factor(phase, levels = c("Follicular", "Luteal")),
                  visit2 = as.integer(visit) - 1L)
  require_condition(all(b$n_executed > 0L) && all(is.finite(b$success_pct)),
                    "Every session needs at least one recorded executed trial for success analysis.")
  if (any(b$n_execution_unknown > 0L)) {
    message(sum(b$n_execution_unknown),
            " execution flags are missing/unrecognised; denominators are saved per session.")
  }
  b
}

fit_behaviour_models <- function(data, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(data, file.path(output_dir, "behaviour_session_summaries.csv"))
  results <- lapply(c("success_pct", "hard_pct"), function(outcome) {
    formula <- stats::as.formula(paste(outcome, "~ phase + visit2 + (1 | participant_id)"))
    fit <- lmerTest::lmer(
      formula, data = data, REML = TRUE,
      contrasts = list(phase = stats::contr.treatment(c("Follicular", "Luteal"), base = 1))
    )
    saveRDS(fit, file.path(output_dir, paste0(outcome, "_fit.rds")))
    writeLines(capture.output(summary(fit)),
               file.path(output_dir, paste0(outcome, "_model_summary.txt")))
    broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE,
                       conf.level = 0.95, conf.method = "Wald") |>
      dplyr::mutate(outcome = outcome)
  }) |>
    dplyr::bind_rows()
  readr::write_csv(results, file.path(output_dir, "behaviour_fixed_effects.csv"))
  print(dplyr::filter(results, term == "visit2"))
  invisible(results)
}

# ============================================================================
# 8. Analysis tables and common summaries
# ============================================================================
# LF is the reference phase. Input tables contain the sessions included in
# each analysis, with visit matching and task quality control completed.
# EP phase comparisons use paired sessions; moderator models use their
# available-case trial/session tables. Missing covariates are not imputed.

PHASES <- c("Follicular", "Luteal")
EMA_VARIABLES <- c("arousal_mean", "arousal_sd", "valence_mean", "valence_sd")
QUESTIONNAIRES <- c("PHQ9_Score", "GAD7_Score", "AMI_Score", "CIRENS_Score")

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

# ============================================================================
# 9. EMA composites and HRV covariates
# ============================================================================
# EMA means and sample SDs use individual prompts on days -3, -2 and -1.
# A mean requires one valid observation; an SD requires two. No daily averaging
# is performed. Each composite requires both of its component items.
# DM EMA scores are standardised over all available session windows; EP scores
# are standardised within the EP reference sample, before trial-level joins.
# HRV inputs are validated RMSSD values in milliseconds. EP requires both
# recordings to pass QC (>=30 resting and >=10 MVC RR intervals) in both phases.
# Log RMSSD is standardised over eligible sessions before joining EP trials.

aggregate_ema <- function(observations, reference_ids = NULL) {
  x <- observations |>
    dplyr::mutate(arousal = state_energy + (1 - state_fatigue),
                  valence = state_mood + (1 - state_anxiety))
  if (!is.null(reference_ids)) x <- dplyr::filter(x, participant_id %in% reference_ids)
  x <- dplyr::bind_rows(
    dplyr::filter(x, day_before_LF %in% -3:-1) |> dplyr::mutate(phase = "Follicular"),
    dplyr::filter(x, day_before_ML %in% -3:-1) |> dplyr::mutate(phase = "Luteal"))
  x |>
    dplyr::group_by(participant_id, phase) |>
    dplyr::summarise(n_prompts = dplyr::n(), n_arousal = sum(is.finite(arousal)),
      n_valence = sum(is.finite(valence)), arousal_mean = mean(arousal, na.rm = TRUE),
      arousal_sd = stats::sd(arousal, na.rm = TRUE), valence_mean = mean(valence, na.rm = TRUE),
      valence_sd = stats::sd(valence, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(phase = factor(phase, levels = PHASES),
      dplyr::across(dplyr::all_of(EMA_VARIABLES), z_score, .names = "{.col}_z"))
}

prepare_ep_hrv <- function(hrv, participant_ids) {
  x <- hrv |>
    dplyr::filter(participant_id %in% participant_ids,
      baseline_ok == TRUE, mvc_ok == TRUE, rr_n_baseline >= 30L, rr_n_mvc >= 10L,
      is.finite(hrv_baseline), hrv_baseline > 0, is.finite(hrv_mvc), hrv_mvc > 0) |>
    dplyr::group_by(participant_id) |>
    dplyr::filter(dplyr::n_distinct(phase) == 2L) |>
    dplyr::ungroup() |>
    dplyr::mutate(hrv_baseline_z = z_score(log(hrv_baseline)), hrv_mvc_z = z_score(log(hrv_mvc)))
  x
}

# ============================================================================
# 10. Effort perception, moderator models, and questionnaires
# ============================================================================
# EP: an OLS rating~force slope per session, followed by a two-sided paired
# LF-minus-ML t-test. Trial-level moderator models use a participant intercept
# and REML. Rating, force and trial number are standardised over the full EP
# moderator trial table before selecting rows with available covariates.
# HRV models include trial order and its phase interaction; EMA models do not.
# DM: session posterior means are entered into ML participant-intercept models.
# HRV models use standardised outcomes and log-HRV; EMA uses raw sensitivities
# and sum-coded phase (LF=+1, ML=-1). Computational uncertainty is not propagated
# into these second-stage regressions. AMI models pool phases.
# DM HRV/AMI models use the session estimates in dm_moderator_sessions.csv;
# DM EMA and symptom models use estimates extracted from the primary fit.

ep_session_slopes <- function(trials) {
  trials |>
    dplyr::group_by(participant_id, phase) |>
    dplyr::summarise(effort_slope = unname(stats::coef(stats::lm(scale_rating ~ force))[2]),
                      .groups = "drop")
}

run_secondary_models <- function(inputs, dm_parameters, output_dir) {
  model_dir <- file.path(output_dir, "models")
  table_dir <- file.path(output_dir, "tables")
  data_dir <- file.path(output_dir, "processed")
  models <- list()
  ep <- inputs$ep_mod |>
    dplyr::mutate(rating_z = z_score(scale_rating), force_z = z_score(force), trial_z = z_score(trial_num))
  slopes <- ep_session_slopes(inputs$ep_phase)
  write_table(slopes, data_dir, "ep_session_slopes")
  write_table(paired_summary(slopes, "effort_slope"), table_dir, "ep_phase_comparison")
  hrv <- prepare_ep_hrv(inputs$ep_hrv, unique(ep$participant_id))
  ema_ep <- aggregate_ema(inputs$ema, inputs$ep_reference$participant_id)
  ema_dm <- aggregate_ema(inputs$ema)
  ep_hrv <- dplyr::inner_join(ep, hrv, by = c("participant_id", "phase"), relationship = "many-to-one")
  ep_ema <- dplyr::inner_join(ep, ema_ep, by = c("participant_id", "phase"), relationship = "many-to-one")
  for (v in c("hrv_baseline_z", "hrv_mvc_z")) {
    name <- paste0("EP_", v)
    f <- stats::as.formula(paste("rating_z ~ force_z * phase *", v, "+ trial_z * phase + (1 | participant_id)"))
    models[[name]] <- fit_mixed_model(f, ep_hrv, name, model_dir)
  }
  for (v in paste0(EMA_VARIABLES, "_z")) {
    name <- paste0("EP_", v)
    d <- dplyr::filter(ep_ema, is.finite(.data[[v]]))
    f <- stats::as.formula(paste("rating_z ~ force_z * phase *", v, "+ (1 | participant_id)"))
    models[[name]] <- fit_mixed_model(f, d, name, model_dir)
  }
  d <- ep |> dplyr::filter(is.finite(AMI_Score)) |> dplyr::mutate(ami_z = z_score(AMI_Score))
  models$EP_AMI <- fit_mixed_model(rating_z ~ force_z * ami_z + trial_z + (1 | participant_id),
                                  d, "EP_AMI", model_dir)
  covariates <- inputs$dm_covariates
  require_condition(nrow(dplyr::anti_join(covariates, dm_parameters,
    by = c("participant_id", "phase"))) == 0L, "DM covariates contain sessions outside the fitted sample.")
  dm <- covariates |>
    dplyr::mutate(hrv_reactivity = log(hrv_mvc / hrv_baseline))
  dm_ema <- dplyr::left_join(dm_parameters, ema_dm, by = c("participant_id", "phase"),
                            relationship = "one-to-one")
  for (y in c("effSens", "rewSens")) {
    for (v in c("hrv_baseline", "hrv_mvc", "hrv_reactivity")) {
      d <- dm |> dplyr::filter(is.finite(.data[[v]]), is.finite(.data[[y]]))
      if (v != "hrv_reactivity") d <- dplyr::filter(d, .data[[v]] > 0)
      d <- d |> dplyr::mutate(outcome_z = z_score(.data[[y]]),
        pred_z = z_score(if (v == "hrv_reactivity") .data[[v]] else log(.data[[v]])))
      name <- paste("DM", y, v, sep = "_")
      models[[name]] <- fit_mixed_model(outcome_z ~ phase * pred_z + (1 | participant_id),
                                       d, name, model_dir, reml = FALSE)
    }
    for (v in paste0(EMA_VARIABLES, "_z")) {
      d <- dm_ema |> dplyr::filter(is.finite(.data[[v]]), is.finite(.data[[y]])) |>
        dplyr::mutate(outcome = .data[[y]], pred_z = .data[[v]])
      name <- paste("DM", y, v, sep = "_")
      models[[name]] <- fit_mixed_model(outcome ~ phase * pred_z + (1 | participant_id),
                                       d, name, model_dir, reml = FALSE, sum_coding = TRUE, bobyqa = FALSE)
    }
    d <- dm |> dplyr::filter(is.finite(AMI_Score), is.finite(.data[[y]])) |>
      dplyr::mutate(ami_z = z_score(AMI_Score), outcome_z = z_score(.data[[y]]))
    name <- paste0("DM_", y, "_AMI")
    models[[name]] <- fit_mixed_model(outcome_z ~ ami_z + (1 | participant_id),
                                     d, name, model_dir, reml = FALSE)
  }
  # Questionnaire totals are matched to the same physical visit as each DM session.
  q <- dplyr::semi_join(inputs$questionnaires, dm_parameters, by = c("participant_id", "phase"))
  q_table <- dplyr::bind_rows(lapply(QUESTIONNAIRES, function(v) paired_summary(q, v)))
  q_table$FDR_p <- stats::p.adjust(q_table$raw_p, method = "BH")
  write_table(q_table, table_dir, "table1_questionnaires")
  symptoms <- dplyr::inner_join(dm_parameters, q, by = c("participant_id", "phase"), relationship = "one-to-one")
  for (v in c("PHQ9_Score", "GAD7_Score")) {
    d <- symptoms |> dplyr::filter(is.finite(.data[[v]])) |> dplyr::mutate(score_z = z_score(.data[[v]]))
    name <- paste0("DM_effSens_", v)
    models[[name]] <- fit_mixed_model(effSens ~ phase * score_z + (1 | participant_id),
                                     d, name, model_dir, reml = FALSE, sum_coding = TRUE, bobyqa = FALSE)
  }
  write_table(dplyr::bind_rows(lapply(models, `[[`, "coefficients")), table_dir, "all_mixed_model_coefficients")
  write_table(ema_ep, data_dir, "ep_ema_session_summary")
  write_table(ema_dm, data_dir, "dm_ema_session_summary")
  write_table(hrv, data_dir, "ep_hrv_sessions")
  write_table(dm_ema, data_dir, "dm_sessions_with_ema")
  write_table(ep, data_dir, "ep_moderator_trials")
  # Simple slopes come from the fitted interaction models, with their covariance.
  simple_slopes <- list()
  for (name in names(models)) {
    if (grepl("AMI$", name)) next
    m <- models[[name]]$fit
    terms <- names(lme4::fixef(m))
    sum_coded <- any(grepl("phase1:", terms))
    target <- if (startsWith(name, "EP_")) grep("^force_z:phaseLuteal:", terms, value = TRUE) else
      grep("^phase(Luteal|1):", terms, value = TRUE)
    require_condition(length(target) == 1L, paste("Ambiguous interaction:", name))
    main <- sub("phase(Luteal|1):", "", target)
    for (phase in c("Follicular", "Luteal", "ML_minus_LF")) {
      L <- stats::setNames(rep(0, length(terms)), terms)
      if (phase != "ML_minus_LF") L[main] <- 1
      L[target] <- if (sum_coded) switch(phase, Follicular = 1, Luteal = -1, ML_minus_LF = -2) else
        switch(phase, Follicular = 0, Luteal = 1, ML_minus_LF = 1)
      simple_slopes[[paste(name, phase)]] <- cbind(data.frame(model = name, contrast = phase),
                                                 lmerTest::contest1D(m, L, confint = TRUE))
    }
  }
  write_table(dplyr::bind_rows(simple_slopes), table_dir, "moderator_simple_slopes")
  list(models = models, slopes = slopes, ep_trials = ep, ep_hrv = ep_hrv,
       ep_ema = ep_ema, dm = dm, dm_ema = dm_ema, questionnaires = q_table)
}

# ============================================================================
# 11. Posterior predictive checks and parameter recovery
# ============================================================================
# Recovery simulates one full choice dataset from a randomly selected posterior
# draw (seed 456), then fits the same model (seed 789, 800 warmup + 800 samples).
# Recovery points correspond to participant-session parameters.

simulate_recovery <- function(fit, prepared) {
  draws <- fit$draws(variables = c("effSens_sess", "rewSens_sess"), format = "matrix")
  set.seed(456)
  index <- sample.int(nrow(draws), 1L)
  n <- prepared$stan_data$N_sess
  eff <- draws[index, paste0("effSens_sess[", seq_len(n), "]")]
  rew <- draws[index, paste0("rewSens_sess[", seq_len(n), "]")]
  dat <- prepared$stan_data
  y <- matrix(0L, n, dat$T)
  for (i in seq_len(n)) for (t in seq_len(dat$T)) {
    v1 <- rew[i] * dat$rew1[t] - eff[i] * dat$eff1[t]^2
    v2 <- rew[i] * dat$rew2[t] - eff[i] * dat$eff2[t]^2
    y[i, t] <- stats::rbinom(1L, 1L, stats::plogis(v2 - v1))
  }
  dat$y <- y
  list(data = dat, true_effSens = unname(eff), true_rewSens = unname(rew), draw_index = index)
}

run_recovery <- function(fit, prepared, output_dir) {
  sim <- simulate_recovery(fit, prepared)
  settings <- PRIMARY_SAMPLING
  settings$iter_warmup <- 800L
  settings$iter_sampling <- 800L
  settings$seed <- 789L
  refit <- fit_choice_model(PRIMARY_STAN_CODE, sim$data, settings, output_dir, "dm_recovery")
  export_choice_results(refit, prepared$sessions, output_dir)
  r <- readr::read_csv(file.path(output_dir, "session_parameters.csv"), show_col_types = FALSE)
  r$true_effSens <- sim$true_effSens
  r$true_rewSens <- sim$true_rewSens
  summary <- dplyr::bind_rows(lapply(c("effSens", "rewSens"), function(v) data.frame(
    parameter = v, r = stats::cor(r[[v]], r[[paste0("true_", v)]]),
    RMSE = sqrt(mean((r[[v]] - r[[paste0("true_", v)]])^2)), simulated_draw = sim$draw_index)))
  write_table(r, output_dir, "recovery_sessions")
  write_table(summary, output_dir, "recovery_summary")
  r
}

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

# ============================================================================
# 12. Figures
# ============================================================================
# Posterior intervals use posterior draws. Phase means average session
# sensitivities within each draw. Descriptive scatterplot lines use OLS;
# inferential coefficients and simple slopes are saved from the mixed models.

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
    dplyr::summarise(p_hard = mean(p_hard), .groups = "drop") |>
    tidyr::pivot_wider(names_from = phase, values_from = p_hard) |>
    dplyr::mutate(difference = 100 * (Luteal - Follicular))
  write_table(choices, dirname(output_dir), "decision_space_values")
  a <- ggplot2::ggplot(choices, ggplot2::aes(hard_reward, hard_effort, fill = difference)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B83B41") +
    ggplot2::labs(x = "Hard-option reward", y = "Hard-option effort (% MVC)", fill = "ML - LF\npercentage points") + theme
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

# ============================================================================
# 13. Run the pipeline
# ============================================================================

run_paper_pipeline <- function(input_dir = file.path("data", "raw"), output_dir = OUTPUT_DIR,
                                include_visit = TRUE, include_recovery = TRUE) {
  packages <- c("cmdstanr", "posterior", "dplyr", "tidyr", "readr", "lme4", "lmerTest",
                "broom.mixed", "ggplot2", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  require_condition(!length(missing), paste("Install required packages:", paste(missing, collapse = ", ")))
  inputs <- read_secondary_inputs(input_dir)
  prepared <- prepare_choice_data(file.path(input_dir, "dm_trials.csv"), include_visit = TRUE)
  behaviour <- prepare_behaviour_data(prepared$trials)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_table(prepared$sessions, file.path(output_dir, "processed"), "dm_session_index")
  write_table(prepared$template, file.path(output_dir, "processed"), "dm_offer_template")
  write_table(data.frame(analysis = c("DM", "EP_phase", "EP_moderators"),
    participants = c(prepared$stan_data$N_subj, dplyr::n_distinct(inputs$ep_phase$participant_id),
                      dplyr::n_distinct(inputs$ep_mod$participant_id)),
    sessions = c(prepared$stan_data$N_sess, nrow(dplyr::distinct(inputs$ep_phase, participant_id, phase)),
                  nrow(dplyr::distinct(inputs$ep_mod, participant_id, phase))),
    trials = c(nrow(prepared$trials), nrow(inputs$ep_phase), nrow(inputs$ep_mod))),
    file.path(output_dir, "tables"), "analysis_counts")
  primary_dir <- file.path(output_dir, "dm_primary")
  primary <- fit_choice_model(PRIMARY_STAN_CODE, prepared$stan_data, PRIMARY_SAMPLING, primary_dir, "dm_primary")
  export_choice_results(primary, prepared$sessions, primary_dir)
  parameters <- readr::read_csv(file.path(primary_dir, "session_parameters.csv"), show_col_types = FALSE)
  secondary <- run_secondary_models(inputs, parameters, output_dir)
  write_table(behaviour, file.path(output_dir, "processed"), "dm_behaviour")
  write_table(dplyr::bind_rows(paired_summary(behaviour, "success_pct"), paired_summary(behaviour, "hard_pct")),
               file.path(output_dir, "tables"), "dm_behaviour_phase_comparisons")
  if (include_visit) {
    visit_data <- prepared$stan_data
    visit_data$visit2 <- as.integer(prepared$sessions$visit) - 1L
    visit_dir <- file.path(output_dir, "dm_visit")
    visit <- fit_choice_model(VISIT_STAN_CODE, visit_data, VISIT_SAMPLING, visit_dir, "dm_phase_visit")
    export_choice_results(visit, prepared$sessions, visit_dir, include_visit = TRUE)
    fit_behaviour_models(behaviour, file.path(output_dir, "behaviour_visit"))
  }
  recovery <- if (include_recovery) run_recovery(primary, prepared, file.path(output_dir, "dm_recovery")) else NULL
  make_figures(inputs, secondary, primary, prepared, behaviour, recovery, file.path(output_dir, "figures"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  write_table(data.frame(component = c("CmdStan", packages), version = c(as.character(cmdstanr::cmdstan_version()),
    vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)))), output_dir, "software_versions")
  message("Pipeline complete: ", normalizePath(output_dir, winslash = "/"))
  invisible(list(primary = primary, secondary = secondary))
}

if (sys.nframe() == 0L) run_paper_pipeline()
