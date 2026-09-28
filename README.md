# Menstrual cycle phase selectively alters anticipated effort costs in naturally cycling women

R and Stan code for the decision-making, effort psychophysics, HRV, EMA and
questionnaire analyses, chronological-visit sensitivity analysis, parameter
recovery, moderation sensitivity simulations and analytical figures. The pipeline
starts from analysis-ready tables with task quality control and session assignment
completed. Participant-level data are not yet included and will be made available
according to the paper's Data Availability statement.

## Files and outputs

| File | Purpose and outputs |
|---|---|
| `run_all.R` | Runs the pipeline and records software versions. |
| `00_setup.R` | Configuration, input checks and shared statistical functions. |
| `01_decision_making.R` | Primary and visit-adjusted Stan fits, posterior contrasts, session parameters, diagnostics and parameter recovery. |
| `02_effort_perception.R` | Session effort-differentiation slopes and their paired phase comparison. |
| `03_hrv.R` | Resting, active-MVC and reactivity moderator models. |
| `04_ema_and_moderators.R` | Three-day EMA summaries, EMA interactions and apathy models. |
| `05_questionnaires_and_descriptives.R` | Table 1 summaries, questionnaire comparisons, symptom models and moderator simple slopes. |
| `06_figures.R` | Figures 2, 3, S2 and S3, plus a posterior predictive plot, as PDF/PNG. |
| `07_moderation_sensitivity.R` | Simulation-based minimum detectable effects for the reported HRV/EMA interactions. |
| `dm_hierarchical_phase.stan` | Primary decision-making model. |
| `dm_hierarchical_phase_visit.stan` | Separate chronological-visit sensitivity model. |
| `DATA_SCHEMA.md` | Required input tables, columns and units. |

## Software

The manuscript reports R 4.3.2 and MATLAB R2023b. This R pipeline was checked
with R 4.4.2 and CmdStan 2.37.0; it does not require MATLAB. Required R packages
are `cmdstanr`, `posterior`, `dplyr`, `tidyr`, `readr`, `lme4`, `lmerTest`,
`broom.mixed`, `ggplot2`, `patchwork` and `png`. CmdStan requires a working C++
toolchain. Configure its location with `cmdstanr::set_cmdstan_path()` or the
`CMDSTAN` environment variable. Installed software versions are written to
`outputs/software_versions.csv` and `outputs/sessionInfo.txt`.

## Run

Place the ten CSV inputs listed in `DATA_SCHEMA.md` in `data/`. From the
repository root, run:

```sh
Rscript run_all.R
```

Fits, diagnostics, processed summaries, numerical tables and figures are written
under `outputs/`. The full run includes parameter recovery and 2,200 simulated
datasets for each of 18 moderation models; these simulations can take substantial
time. To omit these simulations while running the observed-data analyses:

```r
source("run_all.R")
run_all(include_recovery = FALSE, include_mdes = FALSE)
```

## Analysis specification

The primary DM input contains 44 participants, 88 sessions and 44 choices per
session, including both catch trials. Reward is divided by the larger reward
within each offered pair; effort is divided by 100 and squared in the value
function. Choice of option 2 follows a Bernoulli-logit likelihood. LF is coded 0
and ML 1. The two Stan files specify all priors and participant effects; the visit
model adds population visit-2 effects separately from the primary analysis.
Sampling settings are in `00_setup.R`.

EP differentiation is the within-session slope of rating on force. Its primary
phase comparison is paired. Trial-level EP moderator models use participant
intercepts and REML; DM second-stage models use participant intercepts and ML.
DM HRV/apathy models use the supplied session parameter table, while DM EMA and
symptom models use parameters extracted from the primary fit. These second-stage
regressions treat parameter estimates as observed outcomes.

EMA valence is mood + (1 − anxiety), and arousal is energy + (1 − fatigue).
Pre-session means and sample SDs use prompts on days −3 to −1. HRV analyses
start from derived RMSSD measurements; raw ECG preprocessing is not performed.
Questionnaire overall statistics use each participant's two-phase average, with
BH correction across PHQ-9, GAD-7, AMI and CIRENS. Other Table 1 overall SDs use
the pooled session summaries. Data requirements and analysis-specific input
samples are specified in `DATA_SCHEMA.md`.

## Preregistrations

- [Decision-making task](https://aspredicted.org/32gz-pw5p.pdf)
- [Effort psychophysics task](https://aspredicted.org/33m9rz.pdf)
