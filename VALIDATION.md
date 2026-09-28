# Validation scope

## Environment

Validation used R 4.4.2 on Windows, `data.table` 1.16.4, `survival` 3.7-0,
`splines` 4.4.2, `ggplot2` 3.5.2, and `openxlsx` 4.2.8. RDS serialization was
tested. Optional QS serialization was not tested in this environment.
The archived study outputs were generated under R 3.6.3 on Linux.

## Checks performed

- All five analysis scripts parsed and loaded with automatic execution disabled.
- The distributed synthetic test completed all six primary outcomes and all
  66 sensitivity model runs, including all four pollutant adjustment sets and
  their same-sample comparators. All 11 sensitivity inference groups completed.
- Result collection completed for the six primary outcomes and all 66
  sensitivity models. Risk estimates, annual counts, and burden results were
  included in the combined workbook.
- State-to-component additivity, final sample counts, connected reference-band
  selection, shared reference definitions, fixed windows, and the S3
  between-disease exclusion were checked.
- Refitting after the final configuration/metadata cleanup reproduced the
  synthetic overall and selected sensitivity risk coefficients. Core burden,
  simulation, Wald-test and profile-test function bodies were also compared with
  their analysis counterparts.
- Main profile inference was checked to report unavailable between-disease
  comparisons when actual reference definitions differ across outcomes.
- No personal/server paths or individual-level study records are included in
  this distribution. The two publication-formatting scripts are excluded.

## Recalculation from study summaries

Primary calculations were rerun in an isolated copy of the six study outcomes'
current risk tables, full covariance matrices, annual state counts and analysis
metadata. Original study files were not overwritten. Annual and mean annual
nine-state and four-component AN, AF and component-share point estimates
agreed with the supplied results to numerical precision (maximum absolute AN
difference approximately 1.1e-12 deaths). Global risk estimates and between-disease
share differences also agreed. The significance decisions at 0.05 for the
global, within-outcome and between-disease tests were unchanged, including Holm
adjustment where applicable.

This is not a claim of bit-for-bit replication of all Monte Carlo outputs.
Some interval endpoints and tests based on simulated transformed covariances
differed. Across the validation runs, stage-three share-interval endpoints
differed by at most approximately 0.59 percentage points from the supplied
output; the between-disease omnibus statistic was approximately 68.13 versus
68.23, with P < 0.001 in both. Some stage-two ratio interval endpoints differed
more (up to approximately 7.67 percentage points), especially in annual/small
burden calculations. Stage-three stability-screened profile intervals are the
designated output for component-share inference.

Finite simulation, floating-point covariance operations, matrix square-root
decompositions and software environments can affect realized draws even with
the same integer seed. Retain the recorded software environment, coefficient
order, covariance, simulation settings and final analysis outputs when exact
numerical reporting is required. The distribution preserves the underlying
statistical estimators and does not substitute stored manuscript numbers into
the calculations.

## Limits of verification

The full individual-level mortality models were not refitted on the study
population during this packaging check. The complete synthetic pipeline and
the study's downstream summary-based calculations were tested. This does not
independently validate upstream matching, exposure extraction, cause coding,
or pollutant data construction. Reproducing the population risk estimates
requires the authorized prepared datasets described in `INPUT_SCHEMA.md`.

The synthetic test reduces Monte Carlo counts and image resolution and disables
sensitivity figure export. Its results are software test outputs, not study
estimates. The main scripts retain their full analysis defaults.
