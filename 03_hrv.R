# HRV models use previously derived and validated RMSSD measurements.
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

run_hrv <- function(inputs, ep, dm_parameters, output_dir) {
  model_dir <- file.path(output_dir, "models")
  models <- list()
  hrv <- prepare_ep_hrv(inputs$ep_hrv, unique(ep$participant_id))
  ep_hrv <- dplyr::inner_join(ep, hrv, by = c("participant_id", "phase"), relationship = "many-to-one")
  for (v in c("hrv_baseline_z", "hrv_mvc_z")) {
    name <- paste0("EP_", v)
    f <- stats::as.formula(paste("rating_z ~ force_z * phase *", v, "+ trial_z * phase + (1 | participant_id)"))
    models[[name]] <- fit_mixed_model(f, ep_hrv, name, model_dir)
  }
  covariates <- inputs$dm_covariates
  require_condition(nrow(dplyr::anti_join(covariates, dm_parameters,
    by = c("participant_id", "phase"))) == 0L, "DM covariates contain sessions outside the fitted sample.")
  dm <- covariates |> dplyr::mutate(hrv_reactivity = log(hrv_mvc / hrv_baseline))
  for (y in c("effSens", "rewSens")) for (v in c("hrv_baseline", "hrv_mvc", "hrv_reactivity")) {
    d <- dm |> dplyr::filter(is.finite(.data[[v]]), is.finite(.data[[y]]))
    if (v != "hrv_reactivity") d <- dplyr::filter(d, .data[[v]] > 0)
    d <- d |> dplyr::mutate(outcome_z = z_score(.data[[y]]),
      pred_z = z_score(if (v == "hrv_reactivity") .data[[v]] else log(.data[[v]])))
    name <- paste("DM", y, v, sep = "_")
    models[[name]] <- fit_mixed_model(outcome_z ~ phase * pred_z + (1 | participant_id),
                                     d, name, model_dir, reml = FALSE)
  }
  write_table(hrv, file.path(output_dir, "processed"), "ep_hrv_sessions")
  list(models = models, sessions = hrv, ep_trials = ep_hrv, dm = dm)
}
