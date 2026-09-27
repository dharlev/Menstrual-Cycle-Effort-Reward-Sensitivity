
# ============================================================================
# 1. Analysis configuration
# ============================================================================

INPUT_FILE <- file.path("data", "raw", "dm_trials.csv")
OUTPUT_DIR <- file.path("outputs", "dm_analysis")
RUN_VISIT_ANALYSES <- TRUE
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

check_packages <- function(include_visit) {
  packages <- c("cmdstanr", "posterior", "dplyr", "readr")
  if (include_visit) packages <- c(packages, "lme4", "lmerTest", "broom.mixed")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(packages)
}

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
# These checks concern data structure and coding, not the direction, size,
# or statistical significance of any result. No rows are selected by outcome.
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
# 8. Run the analysis in order and save outputs
# ============================================================================
# Output folders:
#   primary/     Stan source, fitted model, actual posterior draws, population
#                contrasts, session estimates, and sampling diagnostics.
#   visit/       Corresponding outputs with chronological visit adjustment.
#   behaviour/   Session denominators, fitted mixed models, fixed-effect tables.
# Top-level files record sample counts, the offer template, session indexing,
# and software versions. There are no expected-result targets in the analysis.

run_analysis <- function(input_file = INPUT_FILE, output_dir = OUTPUT_DIR,
                         include_visit = RUN_VISIT_ANALYSES) {
  packages <- check_packages(include_visit)
  prepared <- prepare_choice_data(input_file, include_visit)
  # Validate behavioural input before starting the more expensive sampling.
  if (include_visit) behaviour <- prepare_behaviour_data(prepared$trials)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  versions <- data.frame(
    component = c("R", "CmdStan", packages),
    version = c(as.character(getRversion()), as.character(cmdstanr::cmdstan_version()),
                vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)))
  )
  readr::write_csv(versions, file.path(output_dir, "software_versions.csv"))
  readr::write_csv(data.frame(
    participants = prepared$stan_data$N_subj, sessions = prepared$stan_data$N_sess,
    choices = nrow(prepared$trials), choices_per_session = TRIALS_PER_SESSION
  ), file.path(output_dir, "sample_counts.csv"))
  readr::write_csv(prepared$sessions, file.path(output_dir, "session_index.csv"))
  readr::write_csv(prepared$template, file.path(output_dir, "offer_template.csv"))

  primary_dir <- file.path(output_dir, "primary")
  primary_fit <- fit_choice_model(PRIMARY_STAN_CODE, prepared$stan_data,
                                   PRIMARY_SAMPLING, primary_dir, "dm_primary")
  primary_results <- export_choice_results(primary_fit, prepared$sessions, primary_dir)

  if (include_visit) {
    visit_data <- prepared$stan_data
    visit_data$visit2 <- as.integer(prepared$sessions$visit) - 1L
    visit_dir <- file.path(output_dir, "visit")
    visit_fit <- fit_choice_model(VISIT_STAN_CODE, visit_data, VISIT_SAMPLING,
                                  visit_dir, "dm_phase_visit")
    export_choice_results(visit_fit, prepared$sessions, visit_dir, include_visit = TRUE)
    fit_behaviour_models(behaviour, file.path(output_dir, "behaviour"))
  }
  message("Analysis complete. Outputs: ", output_dir)
  invisible(primary_results)
}

if (sys.nframe() == 0L) run_analysis()
