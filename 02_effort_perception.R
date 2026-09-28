# Session-level effort differentiation and paired phase comparison
ep_session_slopes <- function(trials) {
  trials |>
    dplyr::group_by(participant_id, phase) |>
    dplyr::summarise(effort_slope = unname(stats::coef(stats::lm(scale_rating ~ force))[2]),
                      .groups = "drop")
}

run_effort_perception <- function(inputs, output_dir) {
  ep <- inputs$ep_mod |>
    dplyr::mutate(rating_z = z_score(scale_rating), force_z = z_score(force), trial_z = z_score(trial_num))
  slopes <- ep_session_slopes(inputs$ep_phase)
  write_table(slopes, file.path(output_dir, "processed"), "ep_session_slopes")
  write_table(paired_summary(slopes, "effort_slope"), file.path(output_dir, "tables"), "ep_phase_comparison")
  list(trials = ep, slopes = slopes)
}
