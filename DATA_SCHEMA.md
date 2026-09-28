# Analysis-ready input tables

All files are CSVs in `data/`; no participant rows or identifiers are distributed
with this code. Use a consistent anonymous `participant_id` across tables.
`phase` is `LF`/`ML` or `Follicular`/`Luteal`; missing values are `NA`. Session
tables have one row per participant and phase, and trial tables one row per
participant, phase and trial. Input tables already reflect their respective
analysis samples and physical-visit matching. The code does not infer phase
from calendar dates or reconstruct task exclusions from raw recordings.

## `dm_trials.csv`

The primary sample: 44 participants, two sessions each, 44 test choices per
session. Include both catch trials and exclude practice trials. Anonymous IDs
determine the model's lexicographic participant order; retain the supplied ID
ordering when reproducing a fit with a fixed seed.

| Column | Meaning |
|---|---|
| `participant_id`, `phase` | Participant and fixed analysis phase. |
| `trial_num` | Integer 1–44; the offer pair at each trial number is the same across sessions. |
| `option_chosen` | 0 = option 1; 1 = option 2. |
| `reward1`, `reward2` | Original reward points, before normalization. |
| `effort1`, `effort2` | Required force as percent MVC, before division by 100. |
| `visit` | Chronological physical visit, 1 or 2. |
| `is_force_exerted` | Recorded execution flag: 0/1 or FALSE/TRUE; preserve missing flags as NA. |
| `surpassed_goal` | Recorded completion flag: 0/1 or FALSE/TRUE. |
| `choice_rt` | Choice response time in seconds, used for descriptive summaries only. |

Catch trials are identified by strict dominance: higher reward and lower force.
All their choices enter the primary likelihood. Phase success descriptives use
all test-trial completion flags; only a recorded true flag counts as completion.
Visit-adjusted execution-success models use trials with recorded execution flag
equal to 1 and require known completion on those trials. No choice is removed
from the primary model because its execution flag is missing.

## `ep_phase_trials.csv`

The 41-participant paired EP sample, with 16 test trials in each phase.

| Column | Meaning |
|---|---|
| `participant_id`, `phase`, `trial_num` | Trial key. |
| `force` | Required force: 10, 30, 50 or 70 percent MVC. |
| `scale_rating` | Subjective effort, 0–100. |
| `vas_rt` | Rating response time in seconds. |
| `surpassed_goal` | Recorded completion flag, as above. |

## `ep_moderator_trials.csv`

The EP trial records included in the reported available-case moderator analyses.
This input is distinct from the paired primary EP table. Required columns are
`participant_id`, `phase`, `trial_num`, `force`, `scale_rating` and `AMI_Score`
(summed Apathy Motivation Index score). Force, rating and trial number are
standardized before covariate-specific complete-case selection.

## `ep_hrv_sessions.csv`

Derived HRV records for the EP moderator sample.

| Column | Meaning |
|---|---|
| `participant_id`, `phase` | Session key. |
| `hrv_baseline`, `hrv_mvc` | Resting and active-MVC RMSSD, milliseconds. |
| `rr_n_baseline`, `rr_n_mvc` | Retained RR-interval counts. |
| `baseline_ok`, `mvc_ok` | Logical recording-quality flags. |

EP HRV models require both recordings to pass quality control, positive finite
RMSSD, at least 30 resting and 10 MVC RR intervals, and both phases. Log RMSSD
is standardized over eligible sessions before joining the trial table.

## `dm_moderator_sessions.csv`

Required columns: `participant_id`, `phase`, `effSens`, `rewSens`,
`hrv_baseline`, `hrv_mvc`, `AMI_Score`. Sensitivities are the positive-scale
session estimates used in the reported HRV/apathy regressions; their estimation
precedes this input. They are not replaced by a different fitted parameter table.
RMSSD is in milliseconds. HRV reactivity is calculated as
`log(hrv_mvc / hrv_baseline)`.

## `ema_observations.csv`

The available prompt-level EMA sample before task-specific joins.

| Column | Meaning |
|---|---|
| `participant_id`, `observation_id` | Prompt key; observation ID is unique within participant. |
| `day_before_LF`, `day_before_ML` | Day relative to the corresponding laboratory session; preceding days are negative. |
| `state_energy`, `state_fatigue`, `state_mood`, `state_anxiety` | Item scores on a 0–1 scale. |

Composites range from 0 to 2. A composite observation requires both items.
A session mean requires one valid observation and an SD requires two.
Summaries use individual prompts on days −3 to −1, without daily averaging or
imputation. DM standardization uses all available session windows.

## `ep_ema_reference.csv`

One column, `participant_id`, with one row per participant in the EP EMA
standardization reference sample. Standardization occurs before joining the
EP moderator trials. Table 1 EMA summaries use the paired primary EP sample.

## `questionnaires.csv`

Required columns: `participant_id`, `phase`, `PHQ9_Score`, `GAD7_Score`,
`AMI_Score`, `CIRENS_Score`. These are questionnaire totals linked to the same
physical visit as the corresponding DM phase session. PHQ-9, GAD-7 and AMI
use the original 44-person sample; CIRENS uses its available complete pairs.
Each measure's paired comparison reports its own n. Raw and four-test
BH-adjusted p values are both exported.

## `participants.csv`

One row per participant in the overall behavioral cohort used for the Table 1
demographics. Columns: `participant_id`, `age` (years), `education_years`.
This table does not determine inclusion in the primary DM or EP models.

## `hrv_descriptive_sessions.csv`

The eligible paired HRV session records used in Table 1. Columns:
`participant_id`, `phase`, `hrv_baseline`, `hrv_mvc` (RMSSD in milliseconds).
Quality control and the descriptive sample assignment precede this input.
This table is used only for descriptive statistics and paired phase comparisons;
HRV moderator models use their respective tables above.
