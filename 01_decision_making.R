# Primary decision-making model and chronological-visit sensitivity
# All 44 test choices, including the two catch trials, enter the likelihood.
# Reward is divided by the maximum within each offer pair; effort is divided by 100.
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
  require_condition(nrow(subjects) == 44L, "The primary DM input must contain the 44-participant analysis sample.")
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

fit_choice_model <- function(stan_file, stan_data, sampling, output_dir, model_name) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  model_file <- file.path(output_dir, paste0(model_name, ".stan"))
  require_condition(file.exists(stan_file), paste("Missing Stan source:", stan_file))
  file.copy(stan_file, model_file, overwrite = TRUE)
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
  refit <- fit_choice_model(PRIMARY_STAN_FILE, sim$data, settings, output_dir, "dm_recovery")
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

run_decision_making <- function(input_dir, output_dir, include_visit = TRUE, include_recovery = TRUE) {
  prepared <- prepare_choice_data(file.path(input_dir, "dm_trials.csv"), include_visit = TRUE)
  behaviour <- prepare_behaviour_data(prepared$trials)
  write_table(prepared$sessions, file.path(output_dir, "processed"), "dm_session_index")
  write_table(prepared$template, file.path(output_dir, "processed"), "dm_offer_template")
  primary_dir <- file.path(output_dir, "dm_primary")
  primary <- fit_choice_model(PRIMARY_STAN_FILE, prepared$stan_data, PRIMARY_SAMPLING, primary_dir, "dm_primary")
  export_choice_results(primary, prepared$sessions, primary_dir)
  parameters <- readr::read_csv(file.path(primary_dir, "session_parameters.csv"), show_col_types = FALSE)
  write_table(behaviour, file.path(output_dir, "processed"), "dm_behaviour")
  visit <- NULL
  if (include_visit) {
    visit_data <- prepared$stan_data
    visit_data$visit2 <- as.integer(prepared$sessions$visit) - 1L
    visit_dir <- file.path(output_dir, "dm_visit")
    visit <- fit_choice_model(VISIT_STAN_FILE, visit_data, VISIT_SAMPLING, visit_dir, "dm_phase_visit")
    export_choice_results(visit, prepared$sessions, visit_dir, include_visit = TRUE)
    fit_behaviour_models(behaviour, file.path(output_dir, "behaviour_visit"))
  }
  recovery <- if (include_recovery) run_recovery(primary, prepared, file.path(output_dir, "dm_recovery")) else NULL
  list(fit = primary, parameters = parameters, prepared = prepared, behaviour = behaviour,
       visit_fit = visit, recovery = recovery)
}

# Posterior variability in the effort phase contrast
export_effort_variability <- function(fit, sessions, output_dir) {
  draws <- fit$draws(variables = c("delta_eff", "sigma_subj", "subj_eff"), format = "matrix")
  write_table(summarise_draws(draws, "sigma_subj[4]"), output_dir, "effort_phase_random_effect_sd")
  ids <- dplyr::distinct(sessions, participant_id, subj_idx)
  individual <- dplyr::bind_rows(lapply(seq_len(nrow(ids)), function(i) {
    x <- as.numeric(draws[, "delta_eff"] + draws[, paste0("subj_eff[4,", ids$subj_idx[i], "]")])
    data.frame(participant_id = ids$participant_id[i], mean = mean(x),
      lower = unname(stats::quantile(x, .025)), upper = unname(stats::quantile(x, .975)),
      prob_positive = mean(x > 0))
  }))
  write_table(individual, output_dir, "participant_effort_phase_contrasts")
  write_table(data.frame(n = nrow(individual), positive_mean = sum(individual$mean > 0),
    positive_probability95 = sum(individual$prob_positive > .95)), output_dir, "participant_effort_phase_counts")
}
