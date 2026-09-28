# Three-day pre-session EMA summaries and moderator analyses
# Means and sample SDs are calculated across prompts on days -3 to -1.
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

run_ema_and_moderators <- function(inputs, ep, dm, dm_parameters, output_dir) {
  model_dir <- file.path(output_dir, "models")
  models <- list()
  ema_ep <- aggregate_ema(inputs$ema, inputs$ep_reference$participant_id)
  ema_dm <- aggregate_ema(inputs$ema)
  ep_ema <- dplyr::inner_join(ep, ema_ep, by = c("participant_id", "phase"), relationship = "many-to-one")
  dm_ema <- dplyr::left_join(dm_parameters, ema_dm, by = c("participant_id", "phase"), relationship = "one-to-one")
  for (v in paste0(EMA_VARIABLES, "_z")) {
    name <- paste0("EP_", v)
    d <- dplyr::filter(ep_ema, is.finite(.data[[v]]))
    f <- stats::as.formula(paste("rating_z ~ force_z * phase *", v, "+ (1 | participant_id)"))
    models[[name]] <- fit_mixed_model(f, d, name, model_dir)
  }
  d <- ep |> dplyr::filter(is.finite(AMI_Score)) |> dplyr::mutate(ami_z = z_score(AMI_Score))
  models$EP_AMI <- fit_mixed_model(rating_z ~ force_z * ami_z + trial_z + (1 | participant_id),
                                  d, "EP_AMI", model_dir)
  for (y in c("effSens", "rewSens")) {
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
  write_table(ema_ep, file.path(output_dir, "processed"), "ep_ema_session_summary")
  write_table(ema_dm, file.path(output_dir, "processed"), "dm_ema_session_summary")
  write_table(dm_ema, file.path(output_dir, "processed"), "dm_sessions_with_ema")
  write_table(ep, file.path(output_dir, "processed"), "ep_moderator_trials")
  list(models = models, ep_trials = ep_ema, dm = dm_ema)
}
