#!/usr/bin/env Rscript
# Run from the repository root.
for (script in c("00_setup.R", "01_decision_making.R", "02_effort_perception.R", "03_hrv.R",
                 "04_ema_and_moderators.R", "05_questionnaires_and_descriptives.R", "06_figures.R",
                 "07_moderation_sensitivity.R")) {
  source(script)
}

run_all <- function(input_dir = INPUT_DIR, output_dir = OUTPUT_DIR,
                    include_visit = TRUE, include_recovery = TRUE, include_mdes = TRUE) {
  packages <- check_packages()
  inputs <- read_secondary_inputs(input_dir)
  dm <- run_decision_making(input_dir, output_dir, include_visit, include_recovery)
  ep <- run_effort_perception(inputs, output_dir)
  hrv <- run_hrv(inputs, ep$trials, dm$parameters, output_dir)
  ema <- run_ema_and_moderators(inputs, ep$trials, hrv$dm, dm$parameters, output_dir)
  questionnaires <- run_questionnaires(inputs, dm$parameters, output_dir)
  models <- c(hrv$models, ema$models, questionnaires$models)
  export_moderator_results(models, output_dir)
  descriptives <- run_descriptives(inputs, dm$prepared, dm$parameters, input_dir, output_dir)
  export_effort_variability(dm$fit, dm$prepared$sessions, file.path(output_dir, "tables"))
  secondary <- list(models = models, slopes = ep$slopes, ep_trials = ep$trials,
                    ep_hrv = hrv$ep_trials, ep_ema = ema$ep_trials, dm = hrv$dm,
                    dm_ema = ema$dm, questionnaires = questionnaires$table)
  write_table(data.frame(analysis = c("DM", "EP_phase", "EP_moderators"),
    participants = c(dm$prepared$stan_data$N_subj, dplyr::n_distinct(inputs$ep_phase$participant_id),
                     dplyr::n_distinct(inputs$ep_mod$participant_id)),
    sessions = c(dm$prepared$stan_data$N_sess, nrow(dplyr::distinct(inputs$ep_phase, participant_id, phase)),
                 nrow(dplyr::distinct(inputs$ep_mod, participant_id, phase))),
    trials = c(nrow(dm$prepared$trials), nrow(inputs$ep_phase), nrow(inputs$ep_mod))),
    file.path(output_dir, "tables"), "analysis_counts")
  make_figures(inputs, secondary, dm$fit, dm$prepared, descriptives$dm_phase, dm$recovery,
               file.path(output_dir, "figures"))
  if (include_mdes) run_moderation_sensitivity(models, output_dir)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  write_table(data.frame(component = c("CmdStan", packages),
    version = c(as.character(cmdstanr::cmdstan_version()),
      vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)))),
    output_dir, "software_versions")
  invisible(list(dm = dm, secondary = secondary))
}

if (sys.nframe() == 0L) run_all()
