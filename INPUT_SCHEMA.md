# Input data specification

## Matched mortality and meteorological records

Each file is a genuine RDS object containing a data frame or `data.table`.
One row represents one case or referent date. Upstream preparation must define
referents on the same day of week within the same calendar month and year for
each death, with residential exposure linkage. Each `id` identifies one matched
set with one case and at least one referent. IDs need only be unique within a
source file; the overall-outcome importer offsets IDs when combining sources.

| Column | Type and meaning |
|---|---|
| `id` | Nonmissing matched-set identifier; integer/numeric or consistently coercible to an identifier |
| `date` | `Date`, `IDate`, or an ISO `YYYY-MM-DD` date coercible to `data.table::as.IDate` |
| `case` | Numeric/integer 1 for the death date and 0 for a matched referent |
| `holiday` | Numeric/integer 1 for a holiday and 0 otherwise |
| `temp_day_lag01` through `temp_day_lag21` | Numeric daytime temperatures in degrees Celsius |
| `temp_night_lag01` through `temp_night_lag21` | Numeric nighttime temperatures in degrees Celsius |
| `rh_day_lag01` through `rh_day_lag03` | Numeric daytime relative humidity, percent |
| `rh_night_lag01` through `rh_night_lag03` | Numeric nighttime relative humidity, percent |
| Optional geographic columns | Used only for exclusion diagnostics when recognized or configured in `DIAGNOSTIC_GEOGRAPHY_COLS` |

All column names are case sensitive. Temperature/humidity suffixes are zero
padded and one based: suffix `01` means model lag 0, `02` lag 1, and `21` lag 20.
The study uses consecutive daytime-nighttime cycles ending on the morning of
the case/referent date: daytime follows sunrise and nighttime follows sunset
on the preceding day. The supplied exposures use local sunrise/sunset with
one hour removed from each boundary, and China Standard Time (UTC+8).
The code consumes these prepared exposures and cannot verify the upstream
meteorological extraction or infer its timing from column names.

Records with missing required predictors are excluded before retaining valid
matched sets. The input datasets should cover 2013-2019 for reproduction of the
study. The code summarizes the actual supplied years; it does not replace a
missing study year with zero deaths. No names, addresses or direct identifiers
are needed in the analytic files.

## Outcome mapping

| Analytic outcome | Source file code(s) | Underlying-cause group |
|---|---|---|
| `I00_I52_I60_I69` | `I00_I52` and `I60_I69` | Overall cardiocerebrovascular mortality |
| `I10_I15` | `I10_I15` | Hypertensive diseases |
| `I20_I25` | `I20_I25` | Ischemic heart disease |
| `I50` | `I50` | Heart failure |
| `I60_I62` | `I60_I62` | Hemorrhagic stroke |
| `I63` | `I63` | Ischemic stroke |

The two overall source groups are combined once. Do not obtain overall mortality
by concatenating the five selected specific causes: they do not exhaust the
overall definition. Between-disease inference uses the five mutually exclusive
underlying-cause outcomes, excluding the overlapping overall outcome.

## Pollutant attachments

There is one pollutant RDS data frame for each of the six analytic outcomes.
Match keys refer to the source files before IDs are recoded for combined analyses.
Rows can be in any order. Duplicate or ambiguous linkage keys are invalid.

| Column | Definition |
|---|---|
| `source_code` | Original source code, such as `I00_I52`, `I60_I69`, or `I50` |
| `source_id` | The corresponding original matched-set `id`, represented consistently as text |
| `date` | Case/referent date with the same meaning as in the matched input |
| `case` | Original case/referent indicator |
| `PM25_day_lag0` ... `PM25_day_lag2` | Daytime PM2.5 concentrations, micrograms per cubic meter |
| `PM25_night_lag0` ... `PM25_night_lag2` | Nighttime PM2.5 concentrations, same units |
| `NO2_day_lag0` ... `NO2_day_lag2`, `NO2_night_lag0` ... `NO2_night_lag2` | Corresponding NO2 concentrations |
| `O3_8H_day_lag0` ... `O3_8H_day_lag2`, `O3_8H_night_lag0` ... `O3_8H_night_lag2` | Corresponding supplied eight-hour ozone metric |

Pollutant suffixes are zero based and unpadded. Each pollutant adjustment uses
the arithmetic mean of its six day/night values over lags 0-2, entered linearly.
The mean is standardized by default for numerical conditioning. The script
does not calculate the upstream eight-hour ozone metric. Additional pollutants
or lags may be present but are not used by the four configured adjustment sets.
Missing pollutant values lead to matched-set filtering separately for each
adjustment set. Its adjusted and unadjusted models use the same resulting sample.

## Nine-state risk and downstream input tables

State codes are `D0_N0`, `DC_N0`, `DH_N0`, `D0_NC`, `DC_NC`, `DH_NC`, `D0_NH`,
`DC_NH`, and `DH_NH`. D/N denote day/night, C/H cold/heat, and 0 reference.
`D0_N0` is the reference. Code never assumes that the CSV row order is a
coefficient order; state labels determine alignment.

`table_04_joint_9_state_rr.csv` contains nine rows, including the reference.
Essential fields include `state`, `beta`, `se`, `rr`, `rr_low`, `rr_high`,
`estimable`, `n_rows`, `n_cases`, `n_referents`, `n_sets`, and `support_status`.
Keep additional fields exported by script 01. For estimable nonreference states,
`rr = exp(beta)` and the 95% limits are `exp(beta +/- 1.96 * se)`.
The reference has beta and SE equal to zero and all RR values equal to one.
Absent states have zero counts, missing coefficients, and an explicit
non-estimable support status. An observed compound state requires estimable
single-component states for Shapley allocation.

`table_04a_joint_logrr_covariance.csv` contains a row-label column `state` and
one numeric column per estimable nonreference state, forming a symmetric
positive-semidefinite matrix. Its diagonal is the squared coefficient SE.
The RDS counterpart is `models/joint_logrr_covariance.rds`, containing
`covariance_non_reference` and `mu_non_reference`.

`table_00c_annual_state_death_counts.csv` has `year`, `state`, `n_cases`, and
`n_deaths`. There are nine rows per year, with `n_deaths` repeated and equal
to the sum of state death counts. Counts describe the final model sample.
`table_02_joint_state_thresholds.csv` and `models/source_parameters.rds`
describe thresholds, reference method, windows, overlap policy and computation
settings. Preserve these outputs alongside the corresponding risk/count data.

Workbook inputs use the same numerical schemas plus `analysis_id`,
`outcome_code`, and `model_variant`. `Coefficient_covariance` stacks the
state-labeled matrix rows for each model. The first three workbook rows hold
the title, note and spacing; headers are on row 4. Read full numeric cells,
not rounded display text or `rr_label`, for subsequent calculations.
