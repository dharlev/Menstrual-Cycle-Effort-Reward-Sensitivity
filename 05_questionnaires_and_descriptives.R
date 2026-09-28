# Paired questionnaire comparisons and symptom moderation
run_questionnaires <- function(inputs, dm_parameters, output_dir) {
  model_dir <- file.path(output_dir, "models")
  table_dir <- file.path(output_dir, "tables")
  models <- list()
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
  list(models = models, table = q_table)
}

export_moderator_results <- function(models, output_dir) {
  table_dir <- file.path(output_dir, "tables")
  write_table(dplyr::bind_rows(lapply(models, `[[`, "coefficients")), table_dir, "all_mixed_model_coefficients")
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
  invisible(models)
}

# Task descriptives use all available session summaries; paired tests use complete pairs.
descriptive_phase_summary <- function(data, variable, label = variable) {
  result <- paired_summary(data, variable, label)
  all_values <- data[[variable]][is.finite(data[[variable]])]
  result$overall_mean <- mean(all_values)
  result$overall_sd <- stats::sd(all_values)
  for (ph in PHASES) {
    values <- data[[variable]][data$phase == ph & is.finite(data[[variable]])]
    prefix <- if (ph == "Follicular") "LF" else "ML"
    result[[paste0(prefix, "_mean")]] <- mean(values)
    result[[paste0(prefix, "_sd")]] <- stats::sd(values)
    result[[paste0(prefix, "_min")]] <- min(values)
    result[[paste0(prefix, "_max")]] <- max(values)
  }
  result
}

run_descriptives <- function(inputs, prepared, dm_parameters, input_dir, output_dir) {
  table_dir <- file.path(output_dir, "tables")
  rt <- read_analysis_table(input_dir, "dm_trials.csv",
    c("participant_id", "phase", "trial_num", "choice_rt"), c("participant_id", "phase", "trial_num"))
  x <- dplyr::left_join(prepared$trials, rt, by = c("participant_id", "phase", "trial_num"),
                        relationship = "one-to-one") |>
    dplyr::mutate(completed = tolower(as.character(surpassed_goal)) %in% c("true", "1", "yes"),
      hard = ifelse(option_chosen == 0, effort1 > effort2, effort2 > effort1),
      option1_dominant = reward1 > reward2 & effort1 < effort2,
      option2_dominant = reward2 > reward1 & effort2 < effort1,
      is_catch = option1_dominant | option2_dominant,
      catch_correct = dplyr::case_when(option1_dominant ~ option_chosen == 0,
                                     option2_dominant ~ option_chosen == 1, TRUE ~ NA))
  # Phase descriptives use all test trials; visit models use recorded execution trials.
  dm_phase <- x |> dplyr::group_by(participant_id, phase) |>
    dplyr::summarise(success_pct = 100 * mean(completed, na.rm = TRUE),
      hard_pct = 100 * mean(hard, na.rm = TRUE), choice_rt = mean(choice_rt, na.rm = TRUE),
      catch_accuracy_pct = 100 * mean(catch_correct[is_catch], na.rm = TRUE), .groups = "drop")
  ep_phase <- inputs$ep_phase |>
    dplyr::mutate(completed = tolower(as.character(surpassed_goal)) %in% c("true", "1", "yes")) |>
    dplyr::group_by(participant_id, phase) |>
    dplyr::summarise(success_pct = 100 * mean(completed, na.rm = TRUE),
      mean_rating = mean(scale_rating, na.rm = TRUE), rating_rt = mean(vas_rt, na.rm = TRUE), .groups = "drop")
  hrv_phase <- read_analysis_table(input_dir, "hrv_descriptive_sessions.csv",
    c("participant_id", "phase", "hrv_baseline", "hrv_mvc"))
  ema_phase <- aggregate_ema(inputs$ema, unique(inputs$ep_phase$participant_id))
  rows <- list()
  add_rows <- function(data, variables, family) {
    dplyr::bind_rows(lapply(variables, function(v) descriptive_phase_summary(data, v))) |>
      dplyr::mutate(analysis = family)
  }
  rows$DM <- add_rows(dm_phase, c("success_pct", "hard_pct", "catch_accuracy_pct", "choice_rt"), "DM")
  rows$EP <- add_rows(ep_phase, c("success_pct", "mean_rating", "rating_rt"), "EP")
  rows$HRV <- add_rows(hrv_phase, c("hrv_baseline", "hrv_mvc"), "HRV")
  rows$EMA <- add_rows(ema_phase, EMA_VARIABLES, "EMA")
  write_table(dplyr::bind_rows(rows), table_dir, "table1_task_hrv_ema")
  write_table(rows$DM, table_dir, "dm_behaviour_phase_comparisons")
  write_table(dm_phase, file.path(output_dir, "processed"), "dm_phase_descriptive_sessions")
  write_table(ep_phase, file.path(output_dir, "processed"), "ep_phase_descriptive_sessions")
  people <- read_analysis_table(input_dir, "participants.csv",
    c("participant_id", "age", "education_years"), "participant_id")
  demographics <- dplyr::bind_rows(lapply(c("age", "education_years"), function(v) {
    values <- people[[v]][is.finite(people[[v]])]
    data.frame(measure = v, n = length(values), mean = mean(values), sd = stats::sd(values),
               min = min(values), max = max(values))
  }))
  write_table(demographics, table_dir, "table1_demographics")
  q <- dplyr::semi_join(inputs$questionnaires, dm_parameters, by = c("participant_id", "phase"))
  screening <- dplyr::bind_rows(lapply(c("PHQ9_Score", "GAD7_Score"), function(v) {
    people <- q |> dplyr::group_by(participant_id) |>
      dplyr::summarise(n_valid = sum(is.finite(.data[[v]])),
                       above = any(.data[[v]] >= 10, na.rm = TRUE), .groups = "drop")
    data.frame(measure = v, observed_n = sum(people$n_valid > 0), above_n = sum(people$above),
               percent = 100 * sum(people$above) / sum(people$n_valid > 0))
  }))
  write_table(screening, table_dir, "questionnaire_screening_thresholds")
  list(dm_phase = dm_phase, ep_phase = ep_phase, tables = dplyr::bind_rows(rows))
}
