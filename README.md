# Daytime and nighttime temperature analysis

This repository contains five R scripts for risk estimation, burden allocation,
heterogeneity testing, sensitivity analysis, and result collection across six
cardiocerebrovascular mortality outcomes.
The analysis starts from prepared, time-stratified case-crossover data. It does
not download mortality records, construct matched referent dates, or extract
gridded environmental exposures.

Individual records are not distributed with this repository. Running the study
analysis requires authorized access to the input datasets. A synthetic integration
test is provided and can run without the study data.

## Files and workflow

| Script | Input | Main output |
|---|---|---|
| `01_main_risk_models.R` | Matched temperature/humidity RDS files | Continuous temperature models, common exposure definitions, nine-state risks, covariance matrices, counts, support and exclusion diagnostics |
| `02_main_burden.R` | Current primary risk tables, covariance, annual counts and metadata from script 01 | Nine-state burden, four-component Shapley allocation, Monte Carlo intervals, diagnostic figures |
| `03_main_heterogeneity.R` | Primary outputs from scripts 01 and 02 | Global risk tests, within-outcome component tests, between-disease profile tests and share intervals |
| `04_sensitivity_analysis.R` | Matched inputs, primary calibration/models, and pollutant attachments for pollution analyses | Sensitivity risk estimates, burden, inference, sample-flow tables and run status |
| `05_collect_results.R` | Current primary and sensitivity result directories | One workbook, with one result type per worksheet and primary comparator rows |

Run scripts 01, 02 and 03 in order. Then run 04 and 05. Script 04 fits the
sensitivity risk models and calculates their burden allocation, uncertainty
intervals, and heterogeneity tests. It also produces the descriptive differences
between pollutant-adjusted models and their same-sample comparators. Keep
`RUN_HETEROGENEITY = TRUE` to include sensitivity inference. Script 05 collects
the completed primary and sensitivity results into one workbook; it does not
perform additional statistical calculations.
Publication-specific Word formatting and combined figure-layout scripts are
outside this repository. The numerical source tables needed for those displays
are retained.

## Software

Use R with `data.table`, `survival`, `splines`, `ggplot2`, and `openxlsx`.
`splines` is included with R. For example:

```r
install.packages(c("data.table", "survival", "ggplot2", "openxlsx"))
```

Model objects default to RDS. The optional `qs` package is needed only when
`OBJECT_FORMAT = "qs"` is selected or existing QS objects are read. RDS and QS
are different formats; changing a filename extension does not convert an object.
The study outputs record R 3.6.3; this distribution was tested with R 4.4.2.
See [VALIDATION.md](VALIDATION.md) for the tested environment and scope.

## Configure input locations

Run commands from the repository root. The defaults are `data/` for inputs and
`results/` for outputs. Edit `config/input_files.csv` to map source codes to local
files. The `file` column is relative to `data/`, unless an absolute path is supplied.
Keep the `kind` and `code` identifiers unchanged. See
[INPUT_SCHEMA.md](INPUT_SCHEMA.md) for the required columns and units.

```text
repository/
  01_main_risk_models.R
  ...
  config/input_files.csv
  data/
    matched/I00_I52.rds
    matched/I60_I69.rds
    matched/I10_I15.rds
    matched/I20_I25.rds
    matched/I50.rds
    matched/I60_I62.rds
    matched/I63.rds
    pollution/I00_I52_I60_I69.rds
    pollution/I10_I15.rds
    pollution/I20_I25.rds
    pollution/I50.rds
    pollution/I60_I62.rds
    pollution/I63.rds
```

Environment variables can override `DAYNIGHT_PROJECT_ROOT`, `DAYNIGHT_DATA_ROOT`,
and `DAYNIGHT_RESULTS_ROOT`. The project root contains the scripts and `config/`.
Use absolute paths for environment overrides. No machine-specific paths are
embedded in the distributed scripts. Input files and generated results are
excluded by `.gitignore`.

```sh
Rscript 01_main_risk_models.R
Rscript 02_main_burden.R
Rscript 03_main_heterogeneity.R
Rscript 04_sensitivity_analysis.R
Rscript 05_collect_results.R
```

The complete study has millions of matched sets and requires substantial memory
and runtime. Run sequentially unless sufficient memory is available. Each script
has a settings section before its statistical functions. For function-level use,
the header specifies the option that disables automatic execution on `source()`.

## Select outcomes and specifications

Scripts 01 and 02 have six `RUN_*` switches. Set unwanted outcomes to `FALSE`.
Script 03 uses `RUN_OUTCOMES`. Overall mortality calibration is needed even when
fitting only one specific cause. Script 01 fits or reuses that calibration.

Script 04 supports:

- `RUN_MODE = "all"`: every outcome and specification;
- `RUN_MODE = "selected"`: the Cartesian product of `SELECTED_OUTCOMES` and
  `SELECTED_SENSITIVITIES`;
- `RUN_MODE = "custom"`: explicit pairs in `CUSTOM_RUNS`.

For all outcomes and all non-pollution sensitivity analyses:

```r
RUN_MODE <- "selected"
SELECTED_OUTCOMES <- names(OUTCOME_NAMES)
SELECTED_SENSITIVITIES <- c("S2_mmt_ci", "S3_outcome_lags", "S4_fixed_lags")
```

The four sensitivity specifications are:

| Identifier | Definition |
|---|---|
| `S1_pollution` | Separate adjustment for PM2.5, NO2 or O3_8H, and joint adjustment for all three; each has an unadjusted comparator using exactly the same retained sample |
| `S2_mmt_ci` | Reference intervals use the overall outcome's daytime and nighttime MMT 95% intervals; primary lag windows are retained |
| `S3_outcome_lags` | Outcome-specific automatically selected windows; primary reference intervals are retained |
| `S4_fixed_lags` | Cold lags 0-20 and heat lags 0-2 for both periods; primary reference intervals are retained |

`POLLUTION_MODELS_TO_RUN` selects pollutant sets within S1. With all options,
script 04 requests 24 outcome/specification combinations and fits 66 models:
48 pollution models including comparators, plus 18 models for S2-S4.
Pollutant attachments are unnecessary when S1 is not selected. Their mapping
rows can remain in the configuration file without the files being present.

`CONTINUE_ON_ERROR = TRUE` allows independent combinations to finish after a
failure. Always inspect `_summary/tables/run_status.csv` and
`heterogeneity_status.csv`; a finished R process alone does not imply every
requested analysis succeeded. Unavailable results retain an explicit status.

## Statistical definitions

Separate daytime and nighttime distributed lag nonlinear models use natural
cubic temperature splines with four degrees of freedom and lag splines with
three internal knots equally spaced on the log scale, including an intercept.
The lag domain is 0-20. Models adjust for holiday and a natural cubic spline
with three degrees of freedom for mean relative humidity over day/night lags
0-2. Conditional logistic regression uses the Efron method.

For each period, the primary reference interval is the connected range containing
the overall-outcome minimum mortality temperature with cumulative risk no more
than 1% above that minimum. The search is restricted to the first through 99th
temperature percentiles. A truncated interval stops by default. Automatic windows
extend from lag 0 through the last lag with lower 95% RR bound above 1 for the
2.5th-percentile cold or 97.5th-percentile heat contrast against the MMT.
Intervening nonsignificant lags are included. If none qualifies, the configured
fallback is 21 days for cold and 3 days for heat. Manual windows count lag 0.

Each period has a cold indicator based on its cold-window mean below the lower
reference bound and a heat indicator based on its heat-window mean above the
upper bound. With both indicators absent, the period is classified as reference.
Because the windows can differ, the two indicators can overlap. Such records
are excluded, followed by retention of sets with one case and at least one
referent. Day and night classifications form nine mutually exclusive states.
All nine enter burden calculations, including the two discordant cold/heat states.
Statistical significance is not an inclusion criterion for burden allocation.

Let `beta_s` be the joint-state log RR and `n_sy` its death count in year `y`.
The state burden is `AN_sy = n_sy * (1 - exp(-beta_s))`. Reference-state burden
is zero. For a compound state with components `a` and `b`, define
`v_s = 1 - exp(-beta_s)` and allocate
`phi_a = (v_a + v_ab - v_b)/2`, `phi_b = (v_b + v_ab - v_a)/2`.
A single active component receives the whole state value. Summing these
allocations weighted by observed deaths gives four component ANs that sum to
the joint AN. Negative contributions are retained.

AN denotes attributable number and AF attributable fraction. Annual AF divides
annual AN by the final included death count in that year. Mean annual AN and AF
are arithmetic means over years. Component shares divide mean annual component
AN by mean annual joint AN; they are not normalized mean AFs. This decomposition
describes the allocation of the estimated burden under the fitted state model.
Component shares do not identify effects of separate heating or cooling interventions.

Correlated coefficient draws propagate uncertainty conditionally on the supplied
covariance, reference intervals, selected windows and counts. Primary burden
uses 1,000 draws with seed 20260523; MMT intervals use 1,000 draws with seed
20260612. Heterogeneity and share intervals use 10,000 draws with base seed
20260925 and outcome offsets. Sensitivity calculations use base seed 20260926.
Intervals are empirical 2.5th and 97.5th percentiles (R quantile type 8).

Within-outcome tests assess equality of four component ANs using three contrasts.
Between-disease tests compare complete four-component share profiles using three
independent coordinates: 12 degrees of freedom for five causes and 3 for each
disease pair when full rank is available. Overall mortality is excluded from
between-disease tests because it overlaps its constituent causes. Cross-cause
coefficient errors are treated as approximately independent. Separate Holm
families apply to global risk tests, within-outcome tests and ten disease-pair
tests. Pairwise tests follow a significant omnibus test. Profile inference
requires stable nonzero total burden and adequate covariance rank. S3 does not
perform between-disease tests because exposure definitions differ across outcomes.

## Analysis outputs

Each primary outcome uses
`results/main/<outcome>/joint_analysis_<configuration>/` with four subdirectories:
`tables/`, `figures/`, `models/`, and `plot_data/`. An optional pipeline index in
`results/` locates outputs. It is not a file-integrity lock. Scripts 02 and 03
can discover configured result directories without that index. Select a run
with `ANALYSIS_TAG` or `ANALYSIS_DIRS` when several configurations coexist.

Important primary outputs are:

| Output | Use |
|---|---|
| `table_00b_joint_analysis_sample_flow.csv` and `table_00d_overlap_exclusion_distribution.csv` | Sample retention and overlap exclusions |
| `table_03_joint_state_counts.csv` | Death and referent counts by state |
| `table_03c_matched_state_comparison_support.csv` and `table_03d_case_referent_state_comparisons.csv` | Within-set comparison support |
| `table_04_joint_9_state_rr.csv` | Nine-state risks; source for selected-state risk figures |
| `table_04a_joint_logrr_covariance.csv` | Full nonreference coefficient covariance |
| `table_09_exact_9_state_burden_annual.csv` and `table_10_exact_9_state_burden_mean_annual.csv` | Nine-state burden |
| `table_11_shapley_4_component_burden_annual.csv` and `table_12_shapley_4_component_burden_mean_annual.csv` | Annual and mean annual component burden and composition |
| `plotdata_cumulative_curves_single_logknots.csv` in `plot_data/` | Continuous exposure-response curves |
| `table_05_component_burden_and_profile.csv` in the heterogeneity output | Component estimates with stability-screened share intervals |

Primary inference is under `results/main/_core_heterogeneity_inference/`.
Sensitivity outputs are under
`results/sensitivity/<specification>/<outcome>/<variant>/`;
inference and run summaries are under `_heterogeneity/` and `_summary/`.

Primary downstream calculations read current risk tables. They require internally
consistent `beta`, `se`, RR and RR limits; they do not infer a covariance matrix
from RR values alone. Script 02 defaults to the most recent supplied covariance
CSV/RDS and rescales its standard deviations to the risk table's SEs while
preserving correlations. Script 04 passes the covariance from each fitted
sensitivity model to its burden and inference calculations. A refitted model
should supply its complete covariance matrix with its risks. Rescaling is
a conditional covariance assumption, not a reconstruction of a refitted model.

Script 01 caches overall calibration for matching inputs and calibration settings;
an input/settings mismatch triggers recomputation. This cache is separate from
the direct current-table inputs used by subsequent calculations.

## Combined workbook

Script 05 writes
`results/sensitivity/_summary/tables/Main_and_sensitivity_all_results.xlsx`.
Each result type has one worksheet containing all available analyses and outcomes.
Primary rows have `analysis_id = "Main"`. Data headers are on row 4; the model
key is `analysis_id`, `outcome_code`, `model_variant`.

Keep `INCLUDE_ADDITIONAL_MODEL_TABLES = TRUE` to retain coefficient covariances,
annual state counts, exposure definitions, and allocation diagnostics alongside
the reported estimates. `Joint_state_RR` contains the risks, `Component_mean_annual`
contains burden estimates, and `Component_profile_inference` contains
stability-screened component-share intervals. Overall and pairwise disease
comparisons are stored separately in `Between_disease_global` and
`Between_disease_pairwise`.

Script 05 rebuilds the workbook from the current result directories and replaces
the previous workbook. To update sensitivity estimates, rerun the required
models in script 04 and then rerun script 05. Editing workbook cells does not
update the fitted models or propagate changes to other results. Refer to the
directory outputs and their recorded settings for computational provenance.

## Synthetic test

```sh
Rscript tests/smoke_test.R
```

The test creates synthetic matched and pollutant data in a temporary directory,
runs all six primary outcomes, all 66 sensitivity models, and result collection.
It checks allocation identities, paired samples, design settings, and the
coverage of the combined workbook. It uses reduced simulation counts and image
resolution for speed; sensitivity figure export is disabled in the test. The
test directory is printed and retained for inspection. These are software checks,
not estimates from the study population.
