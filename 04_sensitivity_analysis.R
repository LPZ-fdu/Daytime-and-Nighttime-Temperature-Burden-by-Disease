#!/usr/bin/env Rscript
# =============================================================================
# Sensitivity risk models, burden allocation, and inference
#
# Purpose: repeat the nine-state risk and Shapley analyses under four specifications:
# S1_pollution, separate and joint PM25/NO2/O3_8H adjustment with an unadjusted fit
# on the identical pollutant-complete sample for each adjustment set;
# S2_mmt_ci, overall-outcome MMT confidence intervals as reference intervals;
# S3_outcome_lags, outcome-specific windows with the common primary references;
# S4_fixed_lags, cold lags 0:20 and heat lags 0:2 in both periods.
#
# Inputs: matched RDS frames with the schema used by 01_main_risk_models.R; primary
# calibration and model objects if available; and RDS pollutant attachments for S1.
# Attachments contain source_code, source_id, date, case, and numeric
# <pollutant>_<day|night>_lag<0|1|2> columns. Source identifiers refer to the original
# matched files. O3_8H is the supplied eight-hour metric. See README.md for linkage.
#
# Outputs: results/sensitivity/<specification>/<outcome>/<variant>/ contains tables,
# figures, models, and plot_data; _heterogeneity contains profile inference;
# _summary contains run status, configuration, and collected results.
#
# RUN_MODE selects all, selected Cartesian combinations, or CUSTOM_RUNS. The four
# pollution adjustment sets each have a same-sample comparator. Between-disease
# tests are omitted for S3 because exposure definitions differ across outcomes.
# Burden and inference use 1,000 and 10,000 coefficient draws, respectively.
# The statistical engines are included below so this script can also run alone.
# Functions only: options(daynight.sensitivity.run = FALSE); source(this_file).
# =============================================================================

# Local paths are configured through environment variables; defaults are relative.
PROJECT_ROOT <- normalizePath(Sys.getenv("DAYNIGHT_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
DATA_ROOT <- Sys.getenv("DAYNIGHT_DATA_ROOT", unset = file.path(PROJECT_ROOT, "data"))
RESULTS_BASE <- Sys.getenv("DAYNIGHT_RESULTS_ROOT", unset = file.path(PROJECT_ROOT, "results"))
WORK_DIR <- PROJECT_ROOT
INPUT_MAP_FILE <- file.path(PROJECT_ROOT, "config", "input_files.csv")
read_input_paths <- function(kind) {
  x <- read.csv(INPUT_MAP_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("kind", "code", "file") %in% names(x)))
    stop("Input mapping requires kind, code, and file columns.")
  x <- x[x$kind == kind, , drop = FALSE]
  if (!nrow(x) || anyNA(x) || anyDuplicated(x$code) || any(!nzchar(x$file)))
    stop("Invalid input mapping for: ", kind)
  absolute <- grepl("^(/|[A-Za-z]:[/\\\\])", x$file)
  setNames(ifelse(absolute, x$file, file.path(DATA_ROOT, x$file)), x$code)
}
MATCHED_FILES <- read_input_paths("matched")

OUTPUT_ROOT <- file.path(RESULTS_BASE, "sensitivity")
PRIMARY_RESULTS_ROOT <- file.path(RESULTS_BASE, "main")

# Set a specific overall-outcome analysis directory if several primary runs exist.
PRIMARY_OVERALL_DIR <- NULL
PRIMARY_CALIBRATION_FILE <- NULL

OUTCOME_NAMES <- c(I00_I52_I60_I69 = "Overall cardiocerebrovascular mortality",
                   I10_I15 = "Hypertensive diseases", I20_I25 = "Ischemic heart disease",
                   I50 = "Heart failure", I60_I62 = "Hemorrhagic stroke", I63 = "Ischemic stroke")
OUTCOME_SOURCES <- list(
  I00_I52_I60_I69 = c(
    I00_I52 = MATCHED_FILES[["I00_I52"]],
    I60_I69 = MATCHED_FILES[["I60_I69"]]),
  I10_I15 = c(I10_I15 = MATCHED_FILES[["I10_I15"]]),
  I20_I25 = c(I20_I25 = MATCHED_FILES[["I20_I25"]]),
  I50 = c(I50 = MATCHED_FILES[["I50"]]),
  I60_I62 = c(I60_I62 = MATCHED_FILES[["I60_I62"]]),
  I63 = c(I63 = MATCHED_FILES[["I63"]]))
POLLUTION_FILES <- read_input_paths("pollution")

SENSITIVITY_NAMES <- c(S1_pollution = "Air-pollutant adjustment",
                       S2_mmt_ci = "Overall MMT confidence-interval reference",
                       S3_outcome_lags = "Outcome-specific lag windows",
                       S4_fixed_lags = "Fixed cold and heat lag windows")

RUN_MODE <- "all"  # "all", "selected", or "custom".
SELECTED_OUTCOMES <- names(OUTCOME_NAMES)
SELECTED_SENSITIVITIES <- names(SENSITIVITY_NAMES)

CUSTOM_RUNS <- data.frame(
  outcome_code = c("I50", "I63"),
  sensitivity = c("S1_pollution", "S3_outcome_lags"), stringsAsFactors = FALSE)

POLLUTION_MODELS <- list(PM25 = "PM25", NO2 = "NO2", O3_8H = "O3_8H",
                         PM25_NO2_O3_8H = c("PM25", "NO2", "O3_8H"))
POLLUTION_MODELS_TO_RUN <- "all"           # e.g. "PM25_NO2_O3_8H"
POLLUTION_LAGS <- 0:2
STANDARDIZE_POLLUTION <- TRUE

PRIMARY_RR_TOLERANCE <- 0.01
ALLOW_TRUNCATED_REFERENCE_BAND <- FALSE
NO_SIGNIFICANT_LAG_ACTION <- "fallback"
FALLBACK_MA_DAYS <- c(cold = 21L, heat = 3L)
FIXED_MA_DAYS <- c(day_cold = 21L, day_heat = 3L, night_cold = 21L, night_heat = 3L)
N_MMT_DRAWS <- 1000L
N_BURDEN_DRAWS <- 1000L
N_HETEROGENEITY_DRAWS <- 10000L
SEED <- 20260926L
FIGURE_DPI <- 600L
WRITE_FIGURES <- TRUE
WRITE_SUMMARY_WORKBOOK <- TRUE
RUN_HETEROGENEITY <- TRUE
REUSE_SAVED_PRIMARY_MODELS <- TRUE
CONTINUE_ON_ERROR <- TRUE

required_packages <- c("data.table", "survival", "splines", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(survival); library(splines); library(ggplot2)})

# 2. Selection, data linkage and nuisance covariates -------------------------------
sens_log <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                  " | ", paste(..., collapse = " "))
clone_object <- function(x) unserialize(serialize(x, NULL))

read_saved_object <- function(path) {
  if (grepl("[.]qs$", path, ignore.case = TRUE)) {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Install qs to read: ", path)
    return(qs::qread(path))
  }
  readRDS(path)
}

find_object <- function(stem) {
  files <- paste0(stem, c(".rds", ".qs"))
  hit <- files[file.exists(files)]
  if (length(hit)) hit[1L] else NULL
}

select_jobs <- function() {
  if (!RUN_MODE %in% c("all", "selected", "custom")) stop("Invalid RUN_MODE.")
  if (RUN_MODE == "custom") jobs <- as.data.table(copy(CUSTOM_RUNS)) else {
    oc <- if (RUN_MODE == "all") names(OUTCOME_NAMES) else SELECTED_OUTCOMES
    ss <- if (RUN_MODE == "all") names(SENSITIVITY_NAMES) else SELECTED_SENSITIVITIES
    jobs <- CJ(outcome_code = oc, sensitivity = ss, sorted = FALSE)
  }
  if (!all(c("outcome_code", "sensitivity") %in% names(jobs)) || !nrow(jobs) ||
      anyNA(jobs) || any(!jobs$outcome_code %in% names(OUTCOME_NAMES)) ||
      any(!jobs$sensitivity %in% names(SENSITIVITY_NAMES))) stop("Invalid outcome/specification selection.")
  unique(jobs[, .(outcome_code, sensitivity)])
}

selected_pollution_models <- function() {
  take <- if (identical(POLLUTION_MODELS_TO_RUN, "all")) names(POLLUTION_MODELS) else POLLUTION_MODELS_TO_RUN
  if (!length(take) || any(!take %in% names(POLLUTION_MODELS)) || anyDuplicated(take) ||
      any(!grepl("^[A-Za-z0-9_]+$", take))) stop("Invalid pollution-model selection.")
  ans <- POLLUTION_MODELS[take]
  allowed <- c("PM25", "PM10", "O3_8H", "O3", "NO2", "SO2", "CO")
  if (any(vapply(ans, function(p) !length(p) || anyNA(p) || anyDuplicated(p) ||
                 any(!p %in% allowed), logical(1)))) stop("Invalid pollutant names/combinations.")
  if (!length(POLLUTION_LAGS) || anyNA(POLLUTION_LAGS) || any(POLLUTION_LAGS < 0) ||
      any(POLLUTION_LAGS != floor(POLLUTION_LAGS)) || anyDuplicated(POLLUTION_LAGS))
    stop("POLLUTION_LAGS must contain distinct nonnegative integers.")
  ans
}

strict_key <- function(x, columns, label) {
  if (!all(columns %in% names(x)) || anyNA(x[, ..columns]) || anyDuplicated(x, by = columns))
    stop("Missing or duplicate join keys in ", label, ": ", paste(columns, collapse = ", "))
}

read_pollution_summaries <- function(code, pollutants) {
  path <- POLLUTION_FILES[[code]]
  if (!file.exists(path)) stop("Pollution attachment not found: ", path)
  z <- as.data.table(read_saved_object(path))
  required <- c("source_code", "source_id", "date", "case")
  columns <- unlist(lapply(pollutants, function(p) paste0(p, "_",
                                                          rep(c("day", "night"), each = length(POLLUTION_LAGS)), "_lag", rep(POLLUTION_LAGS, 2))))
  if (length(setdiff(c(required, columns), names(z)))) stop("Missing pollution columns in ", path,
                                                            ": ", paste(setdiff(c(required, columns), names(z)), collapse = ", "))
  z[, `:=`(source_code = as.character(source_code), source_id = as.character(source_id), date = as.IDate(date))]
  strict_key(z, required, path)
  out <- z[, ..required]
  for (p in pollutants) {
    cc <- paste0(p, "_", rep(c("day", "night"), each = length(POLLUTION_LAGS)),
                 "_lag", rep(POLLUTION_LAGS, 2))
    if (any(!vapply(z[, ..cc], is.numeric, logical(1)))) stop("Pollution values must be numeric: ", p)
    v <- rowMeans(as.matrix(z[, ..cc]), na.rm = FALSE)
    v[!is.finite(v)] <- NA_real_
    out[, (paste0("adj_", p)) := v]
  }
  out
}

read_sensitivity_data <- function(e, code, pollutants = character()) {
  paths <- OUTCOME_SOURCES[[code]]
  if (any(!file.exists(paths))) stop("Missing outcome input(s): ", paste(paths[!file.exists(paths)], collapse = " | "))
  poll <- if (length(pollutants)) read_pollution_summaries(code, pollutants) else NULL
  extra <- paste0("adj_", pollutants)
  required <- c("id", "date", "case", "holiday", sprintf("temp_day_lag%02d", e$LAG_DAYS),
                sprintf("temp_night_lag%02d", e$LAG_DAYS), sprintf("rh_day_lag%02d", 1:3), sprintf("rh_night_lag%02d", 1:3))
  pieces <- audits <- vector("list", length(paths)); offset <- 0
  for (i in seq_along(paths)) {
    sens_log("Reading:", code, names(paths)[i])
    z <- as.data.table(read_saved_object(paths[i]))
    if (length(setdiff(required, names(z)))) stop("Missing model variables: ", paste(setdiff(required, names(z)), collapse = ", "))
    keep <- unique(c(required, e$diagnostic_geography_columns(names(z))))
    drop <- setdiff(names(z), keep); if (length(drop)) z[, (drop) := NULL]
    z[, date := as.IDate(date)]
    z[, case := e$coerce_binary01(case, "case")]
    z[, holiday := e$coerce_binary01(holiday, "holiday")]
    strict_key(z, c("id", "date"), paths[i])
    if (anyNA(z$case)) stop("Invalid case coding in ", paths[i])
    matched <- nrow(z)
    if (length(pollutants)) {
      lookup <- data.table(source_code = names(paths)[i], source_id = as.character(z$id), date = z$date, case = z$case)
      at <- poll[lookup, on = .(source_code, source_id, date, case), which = TRUE]
      if (anyNA(at)) stop("Unmatched source keys in pollution attachment: ", code, "/", names(paths)[i],
                          " (", sum(is.na(at)), " rows). Use the attachment for this complete source dataset.")
      for (p in extra) set(z, j = p, value = poll[[p]][at])
      matched <- sum(!is.na(at))
    }
    old_ids <- unique(z$id)
    z[, id := as.double(match(id, old_ids)) + offset]
    offset <- offset + length(old_ids)
    numeric_cols <- setdiff(required, c("id", "date", "case", "holiday"))
    if (any(!vapply(z[, ..numeric_cols], is.numeric, logical(1)))) stop("Expected numeric temperature/humidity.")
    for (v in numeric_cols) if (any(is.infinite(z[[v]]))) stop("Infinite values in ", v)
    audits[[i]] <- data.table(outcome_code = code, source_code = names(paths)[i], input_file = paths[i],
                              n_rows = nrow(z), n_cases = sum(z$case), n_sets = length(old_ids), n_pollution_keys_matched = matched,
                              pollution_join = if (length(pollutants)) "source_code + source_id + date + case" else "Not applicable")
    pieces[[i]] <- z
    rm(z, old_ids); invisible(gc())
  }
  out <- rbindlist(pieces, fill = TRUE)
  attr(out, "input_audit") <- rbindlist(audits)
  out
}

build_pollution_adjustment <- function(dt, variables) {
  if (!length(variables)) return(list(data = data.table(), summary = data.table()))
  output <- data.table(); rows <- list()
  for (v in variables) {
    x <- dt[[v]]
    if (is.null(x) || any(!is.finite(x))) stop("Incomplete pollution covariate: ", v)
    within_sd <- dt[, .(varies = uniqueN(get(v)) > 1L), by = id]
    scale <- if (STANDARDIZE_POLLUTION) sd(x) else 1
    centre <- if (STANDARDIZE_POLLUTION) mean(x) else 0
    if (!is.finite(scale) || scale <= 0 || !any(within_sd$varies))
      stop("Pollution covariate has no usable within-set variation: ", v)
    output[, (v) := (x - centre) / scale]
    rows[[v]] <- data.table(covariate = v, mean = mean(x), sd = sd(x), minimum = min(x),
                            maximum = max(x), model_center = centre, model_scale = scale,
                            matched_sets_with_variation = sum(within_sd$varies), model_function = "Linear")
  }
  list(data = output, summary = rbindlist(rows))
}

configure_engine <- function(e) {
  for (name in c("DIR_TABLE", "DIR_MODEL", "DIR_PLOT_DATA", "DIR_FIG", "DIR_LOG", "DIR_CONTRIB", "LOG_FILE"))
    assign(name, NULL, envir = e)
  e$OBJECT_FORMAT <- "rds"; e$N_SIM <- N_BURDEN_DRAWS; e$SIM_SEED <- SEED
  e$N_MMT_CI_SIM <- N_MMT_DRAWS; e$FIG_DPI <- FIGURE_DPI
  e$CAP_NEGATIVE_AF <- FALSE; e$SAVE_CASE_CONTRIB <- FALSE
  e$PIPELINE_INDEX_FILE <- file.path(OUTPUT_ROOT, "_summary", "models", "unused_directory_index.rds")
  e$EXTRA_COVARIATES <- character(); e$ADJUSTMENT_COVARIATES <- character()
  if (!WRITE_FIGURES) {e$ggsave_png <- function(...) invisible(NULL); e$ggsave_out <- function(...) invisible(NULL)}
  e
}

set_context <- function(e, code, root, design = NULL) {
  e$OUTCOME_CODE <- code; e$OUTCOME_NAME <- unname(OUTCOME_NAMES[[code]])
  e$OUT_ROOT <- root; e$OUT_BASE <- dirname(root); e$DATA_PATH <- OUTCOME_SOURCES[[code]]
  e$set_output_dirs(root, "sensitivity_analysis_log.txt", FALSE)
  if (!is.null(design)) {
    e$CURRENT_SHARED_DESIGN <- design; e$STATE_THRESHOLD_STRATEGY <- design$thresholds$strategy
    e$apply_auto_lag_windows(design$lags)
  }
}

# 3. Primary calibration and sensitivity designs --------------------------------
resolve_primary_directory <- function() {
  if (!is.null(PRIMARY_OVERALL_DIR)) {
    if (!file.exists(file.path(PRIMARY_OVERALL_DIR, "models", "source_parameters.rds")))
      stop("PRIMARY_OVERALL_DIR has no models/source_parameters.rds: ", PRIMARY_OVERALL_DIR)
    return(PRIMARY_OVERALL_DIR)
  }
  base <- file.path(PRIMARY_RESULTS_ROOT, "I00_I52_I60_I69")
  if (!dir.exists(base)) return(NULL)
  candidates <- list.files(base, pattern = "^source_parameters[.]rds$", recursive = TRUE, full.names = TRUE)
  suitable <- vapply(candidates, function(p) {
    m <- readRDS(p); d <- m$shared_design
    identical(d$settings$reference_method, "rr_tolerance") &&
      isTRUE(all.equal(d$settings$rr_tolerance, PRIMARY_RR_TOLERANCE)) &&
      identical(d$settings$overlap_policy, "exclude")
  }, logical(1))
  candidates <- candidates[suitable]
  if (length(candidates) > 1L) stop("Several primary analyses match. Set PRIMARY_OVERALL_DIR:\n",
                                    paste(dirname(dirname(candidates)), collapse = "\n"))
  if (length(candidates)) dirname(dirname(candidates)) else NULL
}

primary_outcome_directory <- function(overall_dir, code) {
  if (is.null(overall_dir)) return(NULL)
  path <- file.path(dirname(dirname(overall_dir)), code, basename(overall_dir))
  if (dir.exists(path)) path else NULL
}

fit_single_calibration <- function(e, code, root, raw = NULL) {
  set_context(e, code, root)
  if (is.null(raw)) raw <- read_sensitivity_data(e, code)
  e$save_dt(attr(raw, "input_audit"), file.path(e$DIR_TABLE, "table_00_input_source_audit.csv"))
  e$run_single_temperature_models(raw)
}

get_primary_context <- function() {
  e <- configure_engine(make_risk_engine())
  overall_dir <- if (REUSE_SAVED_PRIMARY_MODELS) resolve_primary_directory() else NULL
  meta <- if (!is.null(overall_dir)) readRDS(file.path(overall_dir, "models", "source_parameters.rds")) else NULL
  design <- if (!is.null(meta)) clone_object(meta$shared_design) else NULL
  if (!is.null(design) && (!identical(design$settings$reference_method, "rr_tolerance") ||
                           !isTRUE(all.equal(design$settings$rr_tolerance, PRIMARY_RR_TOLERANCE)) ||
                           !identical(design$settings$overlap_policy, "exclude")))
    stop("The selected primary design must use the configured RR-tolerance interval and exclusion policy.")
  source <- NULL; calibration <- NULL
  if (REUSE_SAVED_PRIMARY_MODELS) {
    candidates <- c(PRIMARY_CALIBRATION_FILE,
                    if (!is.null(design)) design$calibration_file,
                    file.path(PRIMARY_RESULTS_ROOT, "_shared_overall_calibration", c("overall_dlnm_calibration.rds", "overall_dlnm_calibration.qs")))
    candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
    for (path in candidates[file.exists(candidates)]) {
      z <- read_saved_object(path)
      if (!is.null(z$single_results) && (is.null(design) || identical(z$calibration_id, design$calibration_id))) {
        calibration <- z; source <- path; break
      }
    }
  }
  if (is.null(calibration)) {
    path <- if (!is.null(overall_dir)) find_object(file.path(overall_dir, "models", "analysis_results_single_logknots")) else NULL
    if (!is.null(path)) {single <- read_saved_object(path); source <- path} else {
      source <- file.path(OUTPUT_ROOT, "_calibration", "primary_overall")
      single <- fit_single_calibration(e, "I00_I52_I60_I69", source)
    }
    calibration <- list(single_results = single, source_code = "I00_I52_I60_I69",
                        source_name = OUTCOME_NAMES[["I00_I52_I60_I69"]], signature = list(source = source),
                        calibration_id = if (!is.null(design)) design$calibration_id else "sensitivity_overall_primary_calibration",
                        output_root = source, cache_file = source)
  }
  if (is.null(design)) design <- e$make_shared_design(calibration, "rr_tolerance", PRIMARY_RR_TOLERANCE,
                                                      "auto", FIXED_MA_DAYS, NULL, NO_SIGNIFICANT_LAG_ACTION, FALLBACK_MA_DAYS, ALLOW_TRUNCATED_REFERENCE_BAND)
  saved_settings <- calibration$signature$settings
  inherited <- intersect(names(saved_settings), c("LAG_DAYS", "LAG_INDEX", "VAR_DF", "RH_DF",
                                                  "LOG_LAG_N_INTERNAL_KNOTS", "LOG_LAG_INTERCEPT", "MMT_RANGE_PROBS", "CURVE_BY", "CLOGIT_METHOD"))
  for (name in inherited) assign(name, saved_settings[[name]], envir = e)
  e$LAG_N_INTERNAL_KNOTS <- e$LOG_LAG_N_INTERNAL_KNOTS
  e$LAG_INTERCEPT <- e$LOG_LAG_INTERCEPT
  e$LAG_DF <- e$LOG_LAG_N_INTERNAL_KNOTS + 1L + as.integer(e$LOG_LAG_INTERCEPT)
  keys <- c("day_cold_threshold", "day_heat_threshold", "night_cold_threshold", "night_heat_threshold")
  if (length(unlist(design$thresholds[keys])) != 4L || any(!is.finite(unlist(design$thresholds[keys]))))
    stop("The primary design has invalid temperature thresholds.")
  snapshot <- file.path(OUTPUT_ROOT, "_calibration", "models")
  dir.create(snapshot, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(design = design, source = source, primary_overall_directory = overall_dir),
          file.path(snapshot, "primary_design_used.rds"))
  sens_log("Primary calibration source:", source)
  list(engine = e, calibration = calibration, design = design, overall_dir = overall_dir, source = source)
}

get_outcome_single <- function(context, code, raw = NULL) {
  if (code == "I00_I52_I60_I69") return(context$calibration$single_results)
  root <- primary_outcome_directory(context$overall_dir, code)
  path <- if (!is.null(root) && REUSE_SAVED_PRIMARY_MODELS)
    find_object(file.path(root, "models", "analysis_results_single_logknots")) else NULL
  if (!is.null(path)) return(read_saved_object(path))
  fit_single_calibration(context$engine, code, file.path(OUTPUT_ROOT, "_calibration", "outcome_lags", code), raw)
}

build_sensitivity_design <- function(context, specification, code, single = NULL) {
  e <- context$engine; d <- clone_object(context$design)
  if (specification == "S2_mmt_ci") {
    bands <- e$make_reference_bands(context$calibration, "mmt_ci", PRIMARY_RR_TOLERANCE, FALSE)
    d$bands <- bands; d$reference_tag <- "MMT_CI95"
    d$settings$reference_method <- "mmt_ci"; d$settings$rr_tolerance <- NA_real_
    d$thresholds$reference_method <- "mmt_ci"; d$thresholds$rr_tolerance <- NA_real_
    d$thresholds$strategy <- "overall_mmt_ci"; d$thresholds$strategy_tag <- "MMT_CI95"
    for (period in c("day", "night")) {
      b <- bands[exposure == if (period == "day") "daytime" else "nighttime"]
      d$thresholds[[paste0(period, "_cold_threshold")]] <- b$reference_low
      d$thresholds[[paste0(period, "_heat_threshold")]] <- b$reference_high
    }
  }
  if (specification %in% c("S3_outcome_lags", "S4_fixed_lags")) {
    mode <- if (specification == "S3_outcome_lags") "auto" else "manual"
    d$lags <- e$resolve_shared_lags(if (is.null(single)) context$calibration$single_results else single,
                                    mode, FIXED_MA_DAYS, NULL, NO_SIGNIFICANT_LAG_ACTION, FALLBACK_MA_DAYS)
    d$settings$lag_mode <- mode
    d$settings$manual_days <- if (mode == "manual") FIXED_MA_DAYS else NULL
    d$settings$manual_overrides <- NULL
    d$settings$no_sig_action <- NO_SIGNIFICANT_LAG_ACTION
    d$settings$fallback_days <- FALLBACK_MA_DAYS
    d$lags$summary[, lag_source_outcome := if (mode == "auto") code else "Fixed across outcomes"]
    d$lag_tag <- paste0(mode, "_", paste(vapply(d$lags[c("day_cold", "day_heat", "night_cold", "night_heat")],
                                                function(x) length(x$lag_days), integer(1)), collapse = "_"))
  }
  d$settings$sensitivity_specification <- specification
  d$settings$n_sim <- N_BURDEN_DRAWS; d$settings$sim_seed <- SEED; d$settings$cap_negative_af <- FALSE
  d$analysis_tag <- specification
  d
}

# 4. Model execution and matched-sample comparisons -------------------------------
pollution_sample_flow <- function(e, raw, covariates, design) {
  e$EXTRA_COVARIATES <- character(); e$apply_auto_lag_windows(design$lags)
  primary <- e$prepare_joint_model_data(raw)$dt
  keys <- primary[, .(id, date, case)]
  at <- raw[keys, on = .(id, date, case), which = TRUE]
  values <- raw[at, ..covariates]
  keep <- complete.cases(values) & rowSums(!is.finite(as.matrix(values))) == 0L
  complete <- primary[keep]
  good_ids <- complete[, .(n_case = sum(case == 1L), n_control = sum(case == 0L)), by = id][n_case == 1L & n_control >= 1L, id]
  final <- complete[id %in% good_ids]
  describe <- function(x, label) data.table(step = label, n_rows = nrow(x), n_cases = sum(x$case == 1L),
                                            n_referents = sum(x$case == 0L), n_sets = uniqueN(x$id))
  rbindlist(list(describe(raw, "Source records"), describe(primary, "Primary complete valid matched sets before state exclusions"),
                 describe(complete, "Pollutant-complete records"), describe(final, "Pollutant-complete valid matched sets before state exclusions")))
}

write_configuration <- function(e, code, specification, variant, design, pollutants) {
  settings <- list(outcome_code = code, outcome = OUTCOME_NAMES[[code]], specification = specification,
                   specification_label = SENSITIVITY_NAMES[[specification]], model_variant = variant,
                   pollutants = pollutants, pollution_lags = if (length(pollutants)) POLLUTION_LAGS else integer(),
                   pollution_summary = "Arithmetic mean of daytime and nighttime values over selected lags; complete measurements required",
                   reference_source = design$source_code, reference_method = design$settings$reference_method,
                   reference_thresholds = design$thresholds[c("day_cold_threshold", "day_heat_threshold", "night_cold_threshold", "night_heat_threshold")],
                   lag_windows = lapply(design$lags[c("day_cold", "day_heat", "night_cold", "night_heat")], `[[`, "lag_days"),
                   negative_components = "Retained", uncertainty = "Conditional joint-coefficient uncertainty; reference, lag windows and counts fixed",
                   burden_draws = N_BURDEN_DRAWS, seed = SEED)
  tab <- rbindlist(lapply(names(settings), function(k) data.table(setting = k,
                                                                  value = paste(capture.output(dput(settings[[k]])), collapse = " "))))
  fwrite(tab, file.path(e$DIR_TABLE, "table_sensitivity_configuration.csv"))
  saveRDS(settings, file.path(e$DIR_MODEL, "sensitivity_configuration.rds"))
}

fit_model_and_burden <- function(context, raw, code, specification, variant, design,
                                 extra = character(), adjusted = character(), pollutants = character()) {
  e <- context$engine
  root <- file.path(OUTPUT_ROOT, specification, code, variant)
  set_context(e, code, root, design)
  e$EXTRA_COVARIATES <- extra; e$ADJUSTMENT_COVARIATES <- adjusted
  e$AUTO_LAG_NO_SIG_ACTION <- design$settings$no_sig_action
  e$AUTO_LAG_FALLBACK_COLD <- seq_len(design$settings$fallback_days[["cold"]]) - 1L
  e$AUTO_LAG_FALLBACK_HEAT <- seq_len(design$settings$fallback_days[["heat"]]) - 1L
  marker <- file.path(e$DIR_MODEL, "sensitivity_completed.rds")
  if (file.exists(marker)) unlink(marker)
  write_configuration(e, code, specification, variant, design, pollutants)
  e$save_dt(attr(raw, "input_audit"), file.path(e$DIR_TABLE, "table_00_input_source_audit.csv"))
  e$export_shared_design(design)
  e$run_stage1_internal_checks()
  risk <- e$run_joint_temperature_model(raw, context$calibration$single_results)
  b <- configure_engine(make_burden_engine())
  burden <- b$run_stage2_outcome(code, rr_file = file.path(root, "tables", "table_04_joint_9_state_rr.csv"),
                                 n_sim = N_BURDEN_DRAWS, sim_seed = SEED, cap_negative_af = FALSE, covariance_source = "rds")
  saveRDS(list(outcome_code = code, specification = specification, variant = variant,
               completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), marker)
  list(root = root, risk = risk$rr_table, burden = burden$shapley_mean_out,
       population = risk$analysis_meta[c("n_rows_model", "n_cases", "n_matched_sets")])
}

run_combination <- function(context, code, specification, record_result) {
  e <- context$engine
  if (specification == "S1_pollution") {
    models <- selected_pollution_models()
    raw <- read_sensitivity_data(e, code, unique(unlist(models, use.names = FALSE)))
    design <- build_sensitivity_design(context, specification, code)
    for (model in names(models)) {
      pollutants <- models[[model]]; covariates <- paste0("adj_", pollutants)
      root <- file.path(OUTPUT_ROOT, specification, code, paste0(model, "_adjusted"))
      set_context(e, code, root, design)
      flow <- pollution_sample_flow(e, raw, covariates, design)
      pair <- list()
      for (role in c("same_sample_unadjusted", "adjusted")) {
        variant <- paste0(model, "_", role)
        pair[role] <- list(record_result(code, specification, variant, function() {
          z <- fit_model_and_burden(context, raw, code, specification, variant, design,
                                    extra = covariates, adjusted = if (role == "adjusted") covariates else character(), pollutants = pollutants)
          fwrite(flow, file.path(z$root, "tables", "table_pollution_sample_flow.csv"))
          z
        }))
      }
      if (all(vapply(pair, function(x) !is.null(x), logical(1)))) {
        if (!identical(pair$adjusted$population, pair$same_sample_unadjusted$population))
          stop("Adjusted and same-sample comparator populations differ.")
        diff <- merge(pair$adjusted$burden[, .(component, adjusted_AN = AN, adjusted_AF_percent = AF_percent,
                                               adjusted_share_percent = share_percent)], pair$same_sample_unadjusted$burden[, .(component,
                                                                                                                                unadjusted_AN = AN, unadjusted_AF_percent = AF_percent, unadjusted_share_percent = share_percent)], by = "component")
        diff[, `:=`(AN_difference = adjusted_AN - unadjusted_AN,
                    AF_difference_pp = adjusted_AF_percent - unadjusted_AF_percent,
                    share_difference_pp = adjusted_share_percent - unadjusted_share_percent)]
        fwrite(diff, file.path(pair$adjusted$root, "tables", "table_pollution_same_sample_descriptive_comparison.csv"))
      }
    }
  } else {
    raw <- read_sensitivity_data(e, code)
    single <- if (specification == "S3_outcome_lags") get_outcome_single(context, code, raw) else NULL
    design <- build_sensitivity_design(context, specification, code, single)
    record_result(code, specification, "main", function()
      fit_model_and_burden(context, raw, code, specification, "main", design))
  }
  invisible(gc())
}

# 5. Heterogeneity and combined reporting ----------------------------------------
run_sensitivity_heterogeneity <- function(completed) {
  if (!RUN_HETEROGENEITY || !nrow(completed)) return(invisible(NULL))
  h <- make_heterogeneity_engine()
  groups <- unique(completed[, .(sensitivity, variant)])
  inference_status <- list()
  for (i in seq_len(nrow(groups))) {
    specification <- groups$sensitivity[i]; model_variant <- groups$variant[i]
    rows <- completed[sensitivity == specification & variant == model_variant]
    root <- file.path(OUTPUT_ROOT, "_heterogeneity", specification, model_variant)
    h$prepare_analysis_directory(root)
    marker <- file.path(root, "models", "heterogeneity_completed.rds")
    if (file.exists(marker)) unlink(marker)
    error <- tryCatch({
      results <- list(); bundles <- list()
      index <- list(outcomes = setNames(lapply(rows$output_root, function(p) list(output_root = p)), rows$outcome_code),
                    active_design = list(analysis_tag = specification))
      for (code in rows$outcome_code) {
        bundle <- h$read_outcome_bundle(code, index)
        bundles[[code]] <- bundle
        results[[code]] <- h$analyse_outcome(bundle, OUTCOME_NAMES[[code]], N_HETEROGENEITY_DRAWS,
                                             as.integer(SEED + 100L * match(code, names(OUTCOME_NAMES))))
      }
      disease_codes <- setdiff(names(OUTCOME_NAMES), "I00_I52_I60_I69")
      common <- specification != "S3_outcome_lags"
      if (common && all(disease_codes %in% names(bundles))) {
        key <- function(b) b$design[c("reference_method", "rr_tolerance", "thresholds", "lag_column_indices", "overlap_policy")]
        common <- all(vapply(bundles[disease_codes], function(b) isTRUE(all.equal(key(b), key(bundles[[disease_codes[1]]]))), logical(1)))
      }
      between <- h$between_disease_tests(if (common) results else list())
      if (!common) {
        between$omnibus[, status := "Not performed: exposure definitions differ across outcomes"]
        between$pairwise[, status := "Not performed: exposure definitions differ across outcomes"]
      }
      tables <- list(risk = rbindlist(lapply(results, `[[`, "risk")),
                     within = rbindlist(lapply(results, `[[`, "within")),
                     components = rbindlist(lapply(results, `[[`, "components")),
                     diagnostics = rbindlist(lapply(results, `[[`, "diagnostics")),
                     validation = rbindlist(lapply(results, `[[`, "validation")),
                     between = between$omnibus, pairwise = between$pairwise)
      for (name in c("risk", "within")) {
        tables[[name]][, p_holm_across_outcomes := p.adjust(p_value, "holm", n = length(OUTCOME_NAMES))]
        tables[[name]][, multiplicity_family := "Six planned outcomes within this specification/model variant"]
      }
      for (name in names(tables)) fwrite(tables[[name]], file.path(root, "tables", paste0("heterogeneity_", name, ".csv")))
      saveRDS(list(results = results, tables = tables, outcomes_in_this_run = rows$outcome_code,
                   between_disease_common_definition = common, n_sim = N_HETEROGENEITY_DRAWS,
                   multiplicity = "Separate six-outcome Holm families within each specification/variant; ten gated disease pairs",
                   uncertainty = "Conditional coefficient uncertainty; approximate cross-disease independence"),
              file.path(root, "models", "heterogeneity_results.rds"))
      saveRDS(list(run_id = rows$run_id[1], completed = TRUE), marker)
      NULL
    }, error = conditionMessage)
    inference_status[[i]] <- data.table(sensitivity = specification, variant = model_variant,
                                        status = if (is.null(error)) "completed" else "failed", error = if (is.null(error)) "" else error)
    if (!is.null(error)) sens_log("Heterogeneity failed:", specification, model_variant, "|", error)
  }
  rbindlist(inference_status)
}

write_summary <- function(status, jobs, context) {
  root <- file.path(OUTPUT_ROOT, "_summary")
  for (part in c("tables", "figures", "models", "plot_data")) dir.create(file.path(root, part), recursive = TRUE, showWarnings = FALSE)
  fwrite(status, file.path(root, "tables", "run_status.csv"))
  fwrite(jobs, file.path(root, "tables", "requested_combinations.csv"))
  tables <- list(Run_status = status, Requested_combinations = jobs)
  completed <- status[status == "completed"]
  file_map <- c(Joint_state_RR = "table_04_joint_9_state_rr.csv",
                Component_burden = "table_12_shapley_4_component_burden_mean_annual.csv",
                State_counts = "table_03_joint_state_counts.csv",
                Matched_state_support = "table_03c_matched_state_comparison_support.csv",
                Exclusion_distribution = "table_00d_overlap_exclusion_distribution.csv",
                Sample_flow = "table_00b_joint_analysis_sample_flow.csv",
                Pollution_sample_flow = "table_pollution_sample_flow.csv",
                Same_sample_comparison = "table_pollution_same_sample_descriptive_comparison.csv",
                Reference_intervals = "table_02a_shared_reference_intervals.csv",
                Selected_windows = "table_03a_auto_selected_lag_windows.csv")
  for (name in names(file_map)) {
    parts <- lapply(seq_len(nrow(completed)), function(i) {
      file <- file.path(completed$output_root[i], "tables", file_map[[name]])
      if (!file.exists(file)) return(NULL)
      z <- fread(file); z[, `:=`(outcome_code = completed$outcome_code[i],
                                 sensitivity = completed$sensitivity[i], model_variant = completed$variant[i])]
      z
    })
    z <- rbindlist(parts, fill = TRUE)
    if (!ncol(z)) z <- data.table(status = "not_available", note = "No completed applicable models in this invocation.")
    tables[[name]] <- z; fwrite(z, file.path(root, "tables", paste0(name, ".csv")))
  }
  completed_groups <- unique(completed[, .(sensitivity, variant)])
  for (type in c("risk", "within", "components", "diagnostics", "between", "pairwise")) {
    parts <- lapply(seq_len(nrow(completed_groups)), function(i) {
      g <- completed_groups[i]
      marker <- file.path(OUTPUT_ROOT, "_heterogeneity", g$sensitivity, g$variant, "models", "heterogeneity_completed.rds")
      if (!RUN_HETEROGENEITY || !file.exists(marker) || !identical(readRDS(marker)$run_id, status$run_id[1])) return(NULL)
      path <- file.path(OUTPUT_ROOT, "_heterogeneity", g$sensitivity, g$variant, "tables", paste0("heterogeneity_", type, ".csv"))
      if (!file.exists(path)) return(NULL)
      z <- fread(path); z[, `:=`(sensitivity = g$sensitivity, model_variant = g$variant)]; z
    })
    z <- rbindlist(parts, fill = TRUE)
    if (!ncol(z)) z <- data.table(status = "not_available", note = "No completed applicable inference in this invocation.")
    tables[[paste0("Heterogeneity_", type)]] <- z; fwrite(z, file.path(root, "tables", paste0("Heterogeneity_", type, ".csv")))
  }
  notes <- data.table(section = c("Run selection", "Burden", "Uncertainty", "Pollution", "Paired comparison", "Multiplicity"),
                      note = c("Only combinations completed in this invocation enter summaries; inspect run_status.csv for failures.",
                               "AN and AF describe the model-based burden in the final included matched population; component shares use mean annual AN.",
                               "Intervals condition on temperature reference intervals, lag windows and counts; negative contributions are retained. Use heterogeneity component tables for stability-screened share intervals.",
                               "Pollutant terms average complete daytime and nighttime values over the configured lags; O3_8H retains its upstream harmonization assumptions.",
                               "Adjusted and unadjusted pollution models use the same sample. Differences are descriptive and have no independence-based P value.",
                               "Within each specification/variant, risk and within-outcome tests use separate planned six-outcome Holm families; disease pairs use a gated ten-pair family."))
  tables <- c(list(Notes = notes), tables)
  if (WRITE_SUMMARY_WORKBOOK) {
    wb <- openxlsx::createWorkbook()
    title_style <- openxlsx::createStyle(fontSize = 12, textDecoration = "bold", wrapText = TRUE)
    header_style <- openxlsx::createStyle(textDecoration = "bold", border = "bottom", wrapText = TRUE)
    for (name in names(tables)) {
      z <- tables[[name]]; sheet <- substr(name, 1, 31)
      openxlsx::addWorksheet(wb, sheet)
      openxlsx::writeData(wb, sheet, paste("Sensitivity analysis:", gsub("_", " ", name)), startRow = 1, colNames = FALSE)
      openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = 1)
      openxlsx::writeData(wb, sheet, z, startRow = 3, headerStyle = header_style, keepNA = TRUE)
      openxlsx::freezePane(wb, sheet, firstActiveRow = 4, firstActiveCol = 1)
      openxlsx::setColWidths(wb, sheet, cols = seq_len(max(1L, ncol(z))), widths = 18)
    }
    openxlsx::saveWorkbook(wb, file.path(root, "tables", "Sensitivity_analysis_results.xlsx"), overwrite = TRUE)
  }
  saveRDS(list(configuration = list(run_mode = RUN_MODE, jobs = jobs, primary_source = context$source),
               status = status), file.path(root, "models", "sensitivity_run_summary.rds"))
  writeLines(capture.output(sessionInfo()), file.path(root, "models", "sessionInfo.txt"))
  invisible(tables)
}

run_sensitivity_analysis <- function() {
  jobs <- select_jobs()
  if (WRITE_SUMMARY_WORKBOOK && !requireNamespace("openxlsx", quietly = TRUE)) stop("Install openxlsx or set WRITE_SUMMARY_WORKBOOK = FALSE.")
  if (N_BURDEN_DRAWS < 2L || N_MMT_DRAWS < 2L || N_HETEROGENEITY_DRAWS < 1000L)
    stop("Invalid simulation counts; heterogeneity requires at least 1000 draws.")
  if (any(jobs$sensitivity == "S1_pollution")) selected_pollution_models()
  dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
  status_dir <- file.path(OUTPUT_ROOT, "_summary", "tables")
  dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(jobs, file.path(status_dir, "requested_combinations.csv"))
  sens_log("Requested outcome/specification combinations:", nrow(jobs))
  context <- get_primary_context()
  run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%OS6"), "_", Sys.getpid())
  status_rows <- list()
  record_result <- function(code, specification, variant, work) {
    started <- Sys.time(); warnings <- character(); error <- NULL
    sens_log("Starting:", code, specification, variant)
    result <- tryCatch(withCallingHandlers(work(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      if (grepl("coefficient may be infinite|did not converge|ran out of iterations|Loglik converged before", conditionMessage(w), ignore.case = TRUE))
        stop("Model estimation failed: ", conditionMessage(w), call. = FALSE)
    }), error = function(x) {error <<- conditionMessage(x); NULL})
    root <- file.path(OUTPUT_ROOT, specification, code, variant)
    status_rows[[length(status_rows) + 1L]] <<- data.table(run_id = run_id, outcome_code = code, sensitivity = specification,
                                                           variant = variant, status = if (is.null(error)) "completed" else "failed", output_root = root,
                                                           elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
                                                           error = if (is.null(error)) "" else error, warnings = paste(unique(warnings), collapse = " | "))
    fwrite(rbindlist(status_rows, fill = TRUE), file.path(status_dir, "run_status.csv"))
    if (!is.null(error)) {sens_log("Failed:", code, specification, variant, "|", error); if (!CONTINUE_ON_ERROR) stop(error)}
    result
  }
  for (i in seq_len(nrow(jobs))) {
    code <- jobs$outcome_code[i]; specification <- jobs$sensitivity[i]
    tryCatch(run_combination(context, code, specification, record_result), error = function(x) {
      record_result(code, specification, "combination_setup", function() stop(conditionMessage(x), call. = FALSE))
    })
    invisible(gc())
  }
  status <- rbindlist(status_rows, fill = TRUE)
  completed <- status[status == "completed"]
  inference_status <- tryCatch(run_sensitivity_heterogeneity(completed),
                               error = function(x) data.table(status = "failed", error = conditionMessage(x)))
  if (is.null(inference_status)) inference_status <- data.table(status = if (RUN_HETEROGENEITY)
    "no_completed_models" else "not_requested", error = "")
  fwrite(inference_status, file.path(status_dir, "heterogeneity_status.csv"))
  write_summary(status, jobs, context)
  sens_log("Finished. Completed model runs:", nrow(completed), "| failed model/setup runs:", sum(status$status == "failed"))
  sens_log("Summary:", file.path(OUTPUT_ROOT, "_summary", "tables",
                                 if (WRITE_SUMMARY_WORKBOOK) "Sensitivity_analysis_results.xlsx" else "run_status.csv"))
  invisible(status)
}

# 6. Statistical engines ---------------------------------------------------------
make_risk_engine <- function() {
  e <- new.env(parent = environment())
  evalq({
    analysis_file <- function(root, name) {
      section <- if (grepl("^plotdata_", name)) 
        "plot_data"
      else if (grepl("^(input_manifest|output_manifest|source_file_manifest|output_file_manifest)[.]", name)) 
        "models"
      else if (grepl("[.](csv|xlsx)$", name, ignore.case = TRUE)) 
        "tables"
      else if (grepl("[.](png|pdf|svg|tif|tiff|jpg|jpeg)$", name, ignore.case = TRUE)) 
        "figures"
      else "models"
      file.path(root, section, name)
    }
    
    prepare_analysis_directory <- function(root) {
      dir.create(root, recursive = TRUE, showWarnings = FALSE)
      root <- normalizePath(root, winslash = "/", mustWork = TRUE)
      for (section in c("tables", "figures", "models", "plot_data"))
        dir.create(file.path(root, section), showWarnings = FALSE)
      invisible(root)
    }
    
    current_file <- function(root, name, required = TRUE) {
      target <- analysis_file(root, name)
      if (file.exists(target)) return(target)
      if (required) stop("Required input not found: ", target, call. = FALSE)
      NULL
    }
    
    OBJECT_FORMAT <- "rds"
    
    OVERALL_CODE <- "I00_I52_I60_I69"
    
    OVERALL_NAME <- "Overall cardiocerebrovascular mortality"
    
    RUN_REFERENCE_METHOD <- "rr_tolerance"
    
    RUN_RR_TOLERANCE <- 0.01
    
    RUN_ALLOW_TRUNCATED_RR_BAND <- FALSE
    
    RUN_LAG_MODE <- "auto"
    
    RUN_MANUAL_MA_DAYS <- c(day_cold = 21L, day_heat = 3L, night_cold = 21L, night_heat = 3L)
    
    RUN_MANUAL_OVERRIDES <- NULL
    
    RUN_NO_SIG_ACTION <- "fallback"
    
    RUN_FALLBACK_MA_DAYS <- c(cold = 21L, heat = 3L)
    
    CURRENT_SHARED_DESIGN <- NULL
    
    STATE_THRESHOLD_STRATEGY <- NULL
    
    DATA_PATH <- NULL
    
    OUT_BASE <- NULL
    
    OUT_ROOT <- NULL
    
    OUTCOME_CODE <- NULL
    
    OUTCOME_NAME <- NULL
    
    LAG_DAYS <- 1:21
    
    LAG_INDEX <- 0:(length(LAG_DAYS) - 1)
    
    VAR_DF <- 4
    
    RH_DF <- 3
    
    LOG_LAG_N_INTERNAL_KNOTS <- 3
    
    LOG_LAG_INTERCEPT <- TRUE
    
    LAG_DF <- LOG_LAG_N_INTERNAL_KNOTS + 1L + as.integer(LOG_LAG_INTERCEPT)
    
    LAG_KNOT_TYPE <- "log"
    
    LAG_N_INTERNAL_KNOTS <- LOG_LAG_N_INTERNAL_KNOTS
    
    LAG_INTERCEPT <- LOG_LAG_INTERCEPT
    
    MMT_RANGE_PROBS <- c(0.01, 0.99)
    
    CURVE_BY <- 0.1
    
    N_MMT_CI_SIM <- 1000
    
    MMT_CI_SEED <- 20260612
    
    CLOGIT_METHOD <- "efron"
    
    TEST_N <- Inf
    
    AUTO_LAG_ALPHA <- 0.05
    
    AUTO_LAG_SELECTION_RULE <- "zero_to_last"
    
    AUTO_LAG_NO_SIG_ACTION <- RUN_NO_SIG_ACTION
    
    AUTO_LAG_FALLBACK_HEAT <- seq_len(RUN_FALLBACK_MA_DAYS[["heat"]]) - 1L
    
    AUTO_LAG_FALLBACK_COLD <- seq_len(RUN_FALLBACK_MA_DAYS[["cold"]]) - 1L
    
    DAY_COLD_CLASS_LAG_DAYS <- NULL
    
    DAY_HEAT_CLASS_LAG_DAYS <- NULL
    
    NIGHT_COLD_CLASS_LAG_DAYS <- NULL
    
    NIGHT_HEAT_CLASS_LAG_DAYS <- NULL
    
    AUTO_LAG_SELECTION_SUMMARY <- NULL
    
    MODEL_LABEL <- "daynight_3x3_9state"
    
    REF_STATE <- "D0_N0"
    
    BOTH_ANOMALY_POLICY <- "exclude"
    
    N_SIM <- 1000
    
    SIM_SEED <- 20260523
    
    CAP_NEGATIVE_AF <- FALSE
    
    SAVE_CASE_CONTRIB <- FALSE
    
    DIAGNOSTIC_GEOGRAPHY_COLS <- NULL
    
    FONT_FAMILY <- "serif"
    
    FIG_DPI <- 600
    
    SAVE_PDF <- FALSE
    
    DAY_COLOUR <- "#D55E00"
    
    NIGHT_COLOUR <- "#0072B2"
    
    COLD_COLOUR <- "#2166AC"
    
    HEAT_COLOUR <- "#B2182B"
    
    set_output_dirs <- function(out_root, log_file_name, include_contrib = FALSE) {
      prepare_analysis_directory(out_root)
      DIR_TABLE <<- file.path(out_root, "tables")
      DIR_MODEL <<- file.path(out_root, "models")
      DIR_PLOT_DATA <<- file.path(out_root, "plot_data")
      DIR_FIG <<- file.path(out_root, "figures")
      DIR_LOG <<- file.path(out_root, "models")
      DIR_CONTRIB <<- file.path(out_root, "tables")
      dirs <- c(DIR_TABLE, DIR_MODEL, DIR_PLOT_DATA, DIR_FIG, DIR_LOG)
      if (isTRUE(include_contrib)) 
        dirs <- c(dirs, DIR_CONTRIB)
      for (dd in dirs) dir.create(dd, recursive = TRUE, showWarnings = FALSE)
      LOG_FILE <<- file.path(DIR_LOG, log_file_name)
    }
    
    log_msg <- function(...) {
      msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
      cat(msg, "\n")
      cat(msg, "\n", file = LOG_FILE, append = TRUE)
    }
    
    save_dt <- function(dt, file) {
      fwrite(dt, file)
      log_msg("Saved table:", file)
    }
    
    object_path <- function(path) {
      if (identical(OBJECT_FORMAT, "rds")) 
        sub("[.]qs$", ".rds", path)
      else path
    }
    
    write_model_object <- function(obj, file) {
      file <- object_path(file)
      if (OBJECT_FORMAT == "rds") 
        saveRDS(obj, file)
      else {
        if (!requireNamespace("qs", quietly = TRUE)) 
          stop("Install qs or select OBJECT_FORMAT = 'rds'.")
        qs::qsave(obj, file, preset = "fast")
      }
      invisible(file)
    }
    
    read_model_object <- function(file) {
      file <- object_path(file)
      if (tolower(tools::file_ext(file)) == "rds") 
        readRDS(file)
      else {
        if (!requireNamespace("qs", quietly = TRUE)) 
          stop("Reading qs objects requires package qs.")
        qs::qread(file)
      }
    }
    
    save_obj <- function(obj, file) {
      file <- write_model_object(obj, file)
      log_msg("Saved object:", file)
    }
    
    coerce_binary01 <- function(x, variable_name) {
      if (is.factor(x)) 
        x <- as.character(x)
      if (is.logical(x)) 
        x <- as.integer(x)
      if (is.character(x)) {
        x_trim <- trimws(tolower(x))
        mapped <- rep(NA_integer_, length(x_trim))
        mapped[x_trim %in% c("0", "false", "no")] <- 0L
        mapped[x_trim %in% c("1", "true", "yes")] <- 1L
        x <- mapped
      }
      else {
        x <- suppressWarnings(as.numeric(x))
      }
      finite_values <- sort(unique(x[is.finite(x)]))
      if (!all(finite_values %in% c(0, 1))) {
        stop(variable_name, " must contain only binary values coded as 0 and 1. Observed finite values: ", paste(finite_values, 
                                                                                                                 collapse = ", "))
      }
      as.integer(x)
    }
    
    validate_unique_key <- function(dt, key_cols, object_name) {
      if (nrow(dt) == 0L) 
        stop(object_name, " is empty.")
      duplicate_rows <- dt[, .N, by = key_cols][N > 1L]
      if (nrow(duplicate_rows) > 0L) {
        stop(object_name, " contains duplicate rows for key: ", paste(key_cols, collapse = ", "))
      }
      invisible(TRUE)
    }
    
    validate_finite_named_vector <- function(x, expected_names, object_name) {
      if (is.null(names(x))) 
        stop(object_name, " must be a named vector.")
      missing_names <- setdiff(expected_names, names(x))
      if (length(missing_names) > 0L) {
        stop(object_name, " is missing required elements: ", paste(missing_names, collapse = ", "))
      }
      bad_names <- expected_names[!is.finite(x[expected_names])]
      if (length(bad_names) > 0L) {
        stop(object_name, " contains non-finite values for: ", paste(bad_names, collapse = ", "))
      }
      invisible(TRUE)
    }
    
    safe_empirical_quantile <- function(x, probability) {
      x <- x[is.finite(x)]
      if (length(x) == 0L) 
        return(NA_real_)
      as.numeric(quantile(x, probs = probability, na.rm = TRUE, type = 8))
    }
    
    theme_pub <- function(base_size = 10) {
      theme_classic(base_size = base_size, base_family = FONT_FAMILY) + theme(text = element_text(family = FONT_FAMILY, 
                                                                                                  colour = "black"), axis.title = element_text(size = base_size + 1, colour = "black"), axis.text = element_text(size = base_size, 
                                                                                                                                                                                                                 colour = "black"), axis.line = element_line(size = 0.35, colour = "black"), axis.ticks = element_line(size = 0.35, 
                                                                                                                                                                                                                                                                                                                       colour = "black"), strip.background = element_rect(fill = "grey95", colour = "grey70", size = 0.25), strip.text = element_text(size = base_size, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                      colour = "black"), legend.title = element_text(size = base_size, colour = "black"), legend.text = element_text(size = base_size - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       1, colour = "black"), panel.grid = element_blank(), plot.title = element_blank())
    }
    
    theme_nature <- function(base_size = 10) {
      theme_classic(base_size = base_size, base_family = FONT_FAMILY) + theme(text = element_text(family = FONT_FAMILY, 
                                                                                                  colour = "black"), axis.title = element_text(size = base_size + 1, colour = "black"), axis.text = element_text(size = base_size, 
                                                                                                                                                                                                                 colour = "black"), axis.line = element_line(size = 0.35, colour = "black"), axis.ticks = element_line(size = 0.35, 
                                                                                                                                                                                                                                                                                                                       colour = "black"), legend.title = element_text(size = base_size, colour = "black"), legend.text = element_text(size = base_size - 
                                                                                                                                                                                                                                                                                                                                                                                                                                        1, colour = "black"), panel.border = element_blank(), panel.grid = element_blank(), plot.title = element_blank())
    }
    
    ggsave_png <- function(p, filename_base, width, height) {
      png_file <- file.path(DIR_FIG, paste0(filename_base, ".png"))
      ggsave(png_file, p, width = width, height = height, dpi = FIG_DPI, units = "in")
      log_msg("Saved figure:", png_file)
    }
    
    ggsave_out <- function(p, filename_base, width, height) {
      png_file <- file.path(DIR_FIG, paste0(filename_base, ".png"))
      ggsave(png_file, p, width = width, height = height, dpi = FIG_DPI, units = "in")
      log_msg("Saved figure:", png_file)
      if (isTRUE(SAVE_PDF)) {
        pdf_file <- file.path(DIR_FIG, paste0(filename_base, ".pdf"))
        ggsave(pdf_file, p, width = width, height = height, units = "in")
        log_msg("Saved figure:", pdf_file)
      }
    }
    
    make_ns_spec <- function(x, df) {
      x <- x[is.finite(x)]
      if (length(x) == 0) 
        stop("No finite values for spline basis.")
      if (df <= 1) 
        stop("df must be > 1 for ns basis.")
      n_knots <- df - 1
      probs <- seq(0, 1, length.out = n_knots + 2)[-c(1, n_knots + 2)]
      list(type = "ns_quantile_knots", df = df, n_internal_knots = n_knots, knots = as.numeric(quantile(x, probs = probs, 
                                                                                                        na.rm = TRUE, type = 8)), Boundary.knots = range(x, na.rm = TRUE), intercept = FALSE)
    }
    
    make_lag_spec_logknots <- function(lag_index, n_internal_knots = 3, intercept = TRUE) {
      lag_index <- sort(unique(as.numeric(lag_index)))
      lag_index <- lag_index[is.finite(lag_index)]
      if (length(lag_index) < 3) 
        stop("lag_index should contain at least 3 values.")
      if (any(lag_index < 0)) 
        stop("lag_index should be non-negative for log-knots.")
      lag_max <- max(lag_index)
      if (!is.finite(lag_max) || lag_max <= 1) 
        stop("Maximum lag should be > 1 for log-knots.")
      knots <- exp(seq(log(1), log(lag_max), length.out = n_internal_knots + 2))[-c(1, n_internal_knots + 2)]
      list(type = "ns_log_knots", df = NULL, n_internal_knots = n_internal_knots, knots = as.numeric(knots), Boundary.knots = range(lag_index), 
           intercept = intercept)
    }
    
    predict_ns_basis <- function(x, spec) {
      as.matrix(ns(x, knots = spec$knots, Boundary.knots = spec$Boundary.knots, intercept = spec$intercept))
    }
    
    constant_history <- function(x, n_lag = length(LAG_DAYS)) {
      matrix(rep(x, each = n_lag), nrow = length(x), ncol = n_lag, byrow = TRUE)
    }
    
    build_cb_matrix <- function(xmat, var_spec, lag_spec, prefix) {
      xmat <- as.matrix(xmat)
      n <- nrow(xmat)
      L <- ncol(xmat)
      lag_basis <- predict_ns_basis(LAG_INDEX[seq_len(L)], lag_spec)
      p <- ncol(predict_ns_basis(xmat[seq_len(min(10, n)), 1], var_spec))
      q <- ncol(lag_basis)
      out <- matrix(0, nrow = n, ncol = p * q)
      col_names <- character(p * q)
      idx <- 1L
      for (a in seq_len(p)) {
        for (b in seq_len(q)) {
          col_names[idx] <- paste0(prefix, "_v", a, "_l", b)
          idx <- idx + 1L
        }
      }
      colnames(out) <- col_names
      for (ll in seq_len(L)) {
        xb <- predict_ns_basis(xmat[, ll], var_spec)
        idx <- 1L
        for (a in seq_len(p)) {
          for (b in seq_len(q)) {
            out[, idx] <- out[, idx] + xb[, a] * lag_basis[ll, b]
            idx <- idx + 1L
          }
        }
      }
      out
    }
    
    build_cb_constant <- function(x, var_spec, lag_spec, prefix) {
      build_cb_matrix(constant_history(x), var_spec, lag_spec, prefix)
    }
    
    build_cb_one_lag <- function(x, lag_position, var_spec, lag_spec, prefix) {
      x <- as.numeric(x)
      if (lag_position < 1 || lag_position > length(LAG_DAYS)) 
        stop("Invalid lag_position.")
      lag_basis_all <- predict_ns_basis(LAG_INDEX, lag_spec)
      lag_basis <- matrix(lag_basis_all[lag_position, ], nrow = 1)
      xb <- predict_ns_basis(x, var_spec)
      p <- ncol(xb)
      q <- ncol(lag_basis)
      out <- matrix(0, nrow = length(x), ncol = p * q)
      col_names <- character(p * q)
      idx <- 1L
      for (a in seq_len(p)) {
        for (b in seq_len(q)) {
          col_names[idx] <- paste0(prefix, "_v", a, "_l", b)
          out[, idx] <- xb[, a] * lag_basis[1, b]
          idx <- idx + 1L
        }
      }
      colnames(out) <- col_names
      out
    }
    
    fit_clogit_from_design <- function(dt, design_mat, model_name) {
      log_msg("Preparing model data:", model_name)
      design_dt <- as.data.table(design_mat)
      model_dt <- data.table(case = dt$case, id = dt$id, holiday = dt$holiday, rh_lag01_03 = dt$rh_lag01_03)
      rh_basis <- as.data.table(ns(model_dt$rh_lag01_03, df = RH_DF))
      setnames(rh_basis, paste0("rh_ns", seq_len(ncol(rh_basis))))
      model_dt <- cbind(model_dt, rh_basis, design_dt)
      rhs_terms <- c(names(design_dt), names(rh_basis), "holiday")
      form <- as.formula(paste0("case ~ ", paste(rhs_terms, collapse = " + "), " + strata(id)"))
      formula_text <- paste(deparse(form), collapse = " ")
      log_msg("Fitting model:", model_name)
      fit <- clogit(form, data = model_dt, method = CLOGIT_METHOD, model = FALSE, x = FALSE, y = FALSE)
      coefficients <- coef(fit)
      covariance <- vcov(fit)
      fitted_exposure_coefficients <- coefficients[names(design_dt)]
      if (length(fitted_exposure_coefficients) != ncol(design_dt) || any(!is.finite(fitted_exposure_coefficients))) {
        stop("The model ", model_name, " did not return finite coefficients for every exposure-basis term.")
      }
      compact_fit <- list(coefficients = coefficients, covariance = covariance, formula_text = formula_text, exposure_terms = names(design_dt), 
                          rh_terms = names(rh_basis), model_name = model_name, storage_mode = "compact_coefficients_and_covariance_only")
      rm(fit, model_dt, rh_basis, design_dt, form)
      invisible(gc(full = TRUE))
      log_msg("Finished model:", model_name)
      compact_fit
    }
    
    predict_contrast_from_design <- function(fit_obj, x_design, ref_design) {
      terms <- colnames(x_design)
      beta <- fit_obj$coefficients[terms]
      vc <- fit_obj$covariance[terms, terms, drop = FALSE]
      xdiff <- sweep(x_design, 2, ref_design[1, ], "-")
      logrr <- as.numeric(xdiff %*% beta)
      se <- sqrt(rowSums((xdiff %*% vc) * xdiff))
      data.table(logrr = logrr, se = se, rr = exp(logrr), rr_low = exp(logrr - 1.96 * se), rr_high = exp(logrr + 
                                                                                                           1.96 * se))
    }
    
    find_mmt_from_single_model <- function(fit_obj, temp_seq, var_spec, lag_spec, prefix, temp_for_median) {
      terms <- fit_obj$exposure_terms
      b <- fit_obj$coefficients[terms]
      x <- cumulative_design_fast(temp_seq, var_spec, lag_spec, prefix)[, terms, drop = FALSE]
      eta <- function(t) as.numeric(cumulative_design_fast(t, var_spec, lag_spec, prefix)[, terms, drop = FALSE] %*% 
                                      b)
      refine_mmt(temp_seq, as.numeric(x %*% b), eta)
    }
    
    rmvnorm_eigen <- function(n, mu, Sigma) {
      mu <- as.numeric(mu)
      Sigma <- as.matrix(Sigma)
      p <- length(mu)
      if (n < 1L || p < 1L) 
        stop("The simulation size and coefficient dimension must be positive.")
      if (!all(dim(Sigma) == c(p, p))) 
        stop("The covariance matrix has incompatible dimensions.")
      if (any(!is.finite(mu)) || any(!is.finite(Sigma))) {
        stop("Non-finite values were found in the simulation mean or covariance matrix.")
      }
      Sigma <- (Sigma + t(Sigma))/2
      eg <- eigen(Sigma, symmetric = TRUE)
      tolerance <- max(1, max(abs(eg$values))) * 1e-08
      if (min(eg$values) < -tolerance) {
        warning("The covariance matrix had negative eigenvalues; negative values were truncated to zero for simulation.")
      }
      vals <- pmax(eg$values, 0)
      A <- eg$vectors %*% diag(sqrt(vals), nrow = p)
      Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
      sweep(Z %*% t(A), 2, mu, "+")
    }
    
    simulate_mmt_ci_from_single_model <- function(fit_obj, temp_seq, var_spec, lag_spec, prefix, n_sim = N_MMT_CI_SIM, 
                                                  seed = MMT_CI_SEED) {
      terms <- fit_obj$exposure_terms
      x_curve <- build_cb_constant(temp_seq, var_spec, lag_spec, prefix)
      beta <- fit_obj$coefficients[terms]
      vc <- fit_obj$covariance[terms, terms, drop = FALSE]
      if (length(beta) == 0 || any(!is.finite(beta)) || any(!is.finite(vc))) {
        warning("Non-finite beta/vcov found. MMT CI will be returned as NA.")
        return(list(mmt_low = NA_real_, mmt_high = NA_real_, mmt_sim = numeric(0), n_sim = 0L))
      }
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        old_seed <- .Random.seed
      }
      else {
        old_seed <- NULL
      }
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        } else {
          .Random.seed <<- old_seed
        }
      }, add = TRUE)
      set.seed(seed)
      beta_draw <- rmvnorm_eigen(n_sim, beta, vc)
      colnames(beta_draw) <- terms
      eta_draw <- x_curve[, terms, drop = FALSE] %*% t(beta_draw)
      mmt_sim <- vapply(seq_len(n_sim), function(i) {
        eta <- function(t) as.numeric(cumulative_design_fast(t, var_spec, lag_spec, prefix)[, terms, drop = FALSE] %*% 
                                        beta_draw[i, ])
        refine_mmt(temp_seq, eta_draw[, i], eta)
      }, numeric(1))
      list(mmt_low = as.numeric(quantile(mmt_sim, 0.025, na.rm = TRUE, type = 8)), mmt_high = as.numeric(quantile(mmt_sim, 
                                                                                                                  0.975, na.rm = TRUE, type = 8)), mmt_sim = mmt_sim, n_sim = length(mmt_sim))
    }
    
    make_cumulative_curve <- function(fit_obj, temp_seq, mmt, var_spec, lag_spec, prefix) {
      x_curve <- build_cb_constant(temp_seq, var_spec, lag_spec, prefix)
      x_ref <- build_cb_constant(mmt, var_spec, lag_spec, prefix)
      out <- predict_contrast_from_design(fit_obj, x_curve, x_ref)
      out[, `:=`(temperature, temp_seq)]
      out[]
    }
    
    make_lag_curve <- function(fit_obj, temp_value, ref_value, var_spec, lag_spec, prefix) {
      ans <- rbindlist(lapply(seq_along(LAG_DAYS), function(ii) {
        x_lag <- build_cb_one_lag(temp_value, ii, var_spec, lag_spec, prefix)
        ref_lag <- build_cb_one_lag(ref_value, ii, var_spec, lag_spec, prefix)
        tmp <- predict_contrast_from_design(fit_obj, x_lag, ref_lag)
        tmp[, `:=`(lag_day, LAG_INDEX[ii])]
        tmp
      }))
      ans[]
    }
    
    make_mmt_annot <- function(curve_dt, mmt_value) {
      x_min <- min(curve_dt$temperature, na.rm = TRUE)
      x_max <- max(curve_dt$temperature, na.rm = TRUE)
      x_rng <- x_max - x_min
      y_max <- max(curve_dt$rr_high, na.rm = TRUE)
      y_min <- min(curve_dt$rr_low, na.rm = TRUE)
      y_rng <- y_max - y_min
      if (mmt_value > x_min + 0.75 * x_rng) {
        x_lab <- mmt_value - 0.03 * x_rng
        hjust_lab <- 1
      }
      else {
        x_lab <- mmt_value + 0.03 * x_rng
        hjust_lab <- 0
      }
      data.table(x = x_lab, y = y_max - 0.08 * y_rng, label = sprintf("MMT = %.1f °C", mmt_value), hjust = hjust_lab)
    }
    
    plot_cumulative_curve <- function(curve_dt, aa, method_suffix) {
      mmt_value <- unique(curve_dt$mmt)
      mmt_low <- unique(curve_dt$mmt_low)
      mmt_high <- unique(curve_dt$mmt_high)
      p_curve <- ggplot(curve_dt, aes(x = temperature, y = rr)) + geom_ribbon(aes(ymin = rr_low, ymax = rr_high), 
                                                                              fill = aa$colour, alpha = 0.22) + geom_line(colour = aa$colour, size = 0.65) + geom_hline(yintercept = 1, 
                                                                                                                                                                        linetype = "dashed", size = 0.35) + geom_vline(xintercept = mmt_value, linetype = "dotted", size = 0.45, 
                                                                                                                                                                                                                       colour = aa$colour) + geom_vline(xintercept = c(mmt_low, mmt_high), linetype = "dashed", size = 0.35, colour = aa$colour, 
                                                                                                                                                                                                                                                        alpha = 0.85) + geom_vline(xintercept = unique(c(curve_dt$p025, curve_dt$p975)), linetype = "longdash", 
                                                                                                                                                                                                                                                                                   size = 0.3, colour = "grey35") + labs(x = aa$xlab, y = "Cumulative relative risk") + theme_pub()
      mmt_lab <- make_mmt_annot(curve_dt, mmt_value)
      if (is.finite(mmt_low) && is.finite(mmt_high)) {
        mmt_text <- sprintf("MMT = %.1f °C\n95%% CI: %.1f–%.1f °C", mmt_value, mmt_low, mmt_high)
      }
      else {
        mmt_text <- sprintf("MMT = %.1f °C", mmt_value)
      }
      p_curve <- p_curve + annotate("text", x = mmt_lab$x, y = mmt_lab$y, label = mmt_text, hjust = mmt_lab$hjust, 
                                    vjust = 1, size = 3, family = FONT_FAMILY)
      ggsave_png(p_curve, paste0("fig_01_cumulative_curve_", aa$exposure, "_single_", method_suffix), width = 3.8, 
                 height = 3.2)
    }
    
    plot_lag_curve <- function(lag_dt, aa, method_suffix) {
      p_lag <- ggplot(lag_dt, aes(x = lag_day, y = rr, colour = contrast_label, fill = contrast_label)) + geom_ribbon(aes(ymin = rr_low, 
                                                                                                                          ymax = rr_high), alpha = 0.18, colour = NA) + geom_line(size = 0.65) + geom_hline(yintercept = 1, linetype = "dashed", 
                                                                                                                                                                                                            size = 0.35) + scale_colour_manual(values = c(`P2.5 vs MMT` = COLD_COLOUR, `P97.5 vs MMT` = HEAT_COLOUR), 
                                                                                                                                                                                                                                               breaks = c("P2.5 vs MMT", "P97.5 vs MMT"), labels = c(expression(P[2.5] ~ "vs MMT"), expression(P[97.5] ~ 
                                                                                                                                                                                                                                                                                                                                                 "vs MMT")), name = NULL) + scale_fill_manual(values = c(`P2.5 vs MMT` = COLD_COLOUR, `P97.5 vs MMT` = HEAT_COLOUR), 
                                                                                                                                                                                                                                                                                                                                                                                              breaks = c("P2.5 vs MMT", "P97.5 vs MMT"), labels = c(expression(P[2.5] ~ "vs MMT"), expression(P[97.5] ~ 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                "vs MMT")), name = NULL) + labs(x = "Lag days", y = "Relative risk") + theme_pub() + theme(legend.position = "top")
      ggsave_png(p_lag, paste0("fig_02_lag_effects_", aa$exposure, "_single_", method_suffix), width = 4, height = 3.2)
    }
    
    prepare_single_model_data <- function(dt_raw) {
      required_base <- c("id", "date", "case", "holiday")
      day_cols <- sprintf("temp_day_lag%02d", LAG_DAYS)
      night_cols <- sprintf("temp_night_lag%02d", LAG_DAYS)
      rh_day_cols <- sprintf("rh_day_lag%02d", 1:3)
      rh_night_cols <- sprintf("rh_night_lag%02d", 1:3)
      required_cols <- c(required_base, day_cols, night_cols, rh_day_cols, rh_night_cols)
      missing_cols <- setdiff(required_cols, names(dt_raw))
      if (length(missing_cols) > 0) 
        stop("Missing required columns for the single-temperature analysis: ", paste(missing_cols, collapse = ", "))
      dt <- copy(dt_raw[, ..required_cols])
      if (is.finite(TEST_N)) {
        log_msg("TEST_N is finite. Keeping complete sets whose ids occur in the first", TEST_N, "rows.")
        test_ids <- unique(dt$id[seq_len(min(TEST_N, nrow(dt)))])
        dt <- dt[id %in% test_ids]
      }
      dt[, `:=`(date, as.IDate(date))]
      dt[, `:=`(case, coerce_binary01(case, "case"))]
      dt[, `:=`(holiday, coerce_binary01(holiday, "holiday"))]
      rh_short_cols <- c(sprintf("rh_day_lag%02d", 1:3), sprintf("rh_night_lag%02d", 1:3))
      dt[, `:=`(rh_lag01_03, rowMeans(.SD, na.rm = FALSE)), .SDcols = rh_short_cols]
      dt[, `:=`((rh_short_cols), NULL)]
      model_vars_for_complete <- c("id", "date", "case", "holiday", "rh_lag01_03", day_cols, night_cols)
      log_msg("Filtering complete cases for the single-temperature models...")
      dt <- dt[dt[, complete.cases(.SD), .SDcols = model_vars_for_complete]]
      log_msg("Filtering valid matched sets for the single-temperature models...")
      valid_id <- dt[, .(n_case = sum(case == 1), n_ref = sum(case == 0)), by = id][n_case == 1 & n_ref >= 1, id]
      dt <- dt[id %in% valid_id]
      rm(valid_id)
      gc()
      if (!nrow(dt)) 
        stop("No complete valid matched sets remain for the single-temperature models.")
      log_msg("Rows after filtering:", nrow(dt))
      log_msg("Cases:", dt[, sum(case == 1)])
      log_msg("Matched sets:", dt[, uniqueN(id)])
      list(dt = dt, day_cols = day_cols, night_cols = night_cols, rh_day_cols = rh_day_cols, rh_night_cols = rh_night_cols)
    }
    
    diagnostic_geography_columns <- function(column_names) {
      if (!is.null(DIAGNOSTIC_GEOGRAPHY_COLS)) {
        if (!is.character(DIAGNOSTIC_GEOGRAPHY_COLS) || anyNA(DIAGNOSTIC_GEOGRAPHY_COLS)) 
          stop("DIAGNOSTIC_GEOGRAPHY_COLS must be NULL or a character vector of source column names.")
        return(intersect(unique(DIAGNOSTIC_GEOGRAPHY_COLS), column_names))
      }
      candidates <- c("province", "province_code", "province_name", "prov", "prov_code", "city", "city_code", "city_name", 
                      "prefecture", "prefecture_code", "county", "county_code", "county_name", "district", "district_code", "region", 
                      "region_code", "region_name", "site", "site_id", "location", "location_id", "\u7701", "\u7701\u4efd", "\u7701\u4ee3\u7801", 
                      "\u7701\u7f16\u7801", "\u5e02", "\u57ce\u5e02", "\u5730\u5e02", "\u5e02\u4ee3\u7801", "\u5e02\u7f16\u7801", "\u53bf", "\u533a\u53bf", "\u53bf\u4ee3\u7801", "\u53bf\u7f16\u7801", 
                      "\u5730\u533a")
      column_names[tolower(column_names) %in% candidates]
    }
    
    summarize_overlap_exclusions <- function(dt) {
      geo_cols <- diagnostic_geography_columns(names(dt))
      keep <- unique(c("id", "date", "case", "joint_state", "D_C_raw", "D_H_raw", "N_C_raw", "N_H_raw", geo_cols))
      z <- copy(dt[, ..keep])
      z[, `:=`(day_overlap = D_C_raw == 1L & D_H_raw == 1L, night_overlap = N_C_raw == 1L & N_H_raw == 1L, directly_excluded = is.na(joint_state))]
      if (anyNA(z$day_overlap) || anyNA(z$night_overlap)) 
        stop("Overlap diagnostics require complete temperature indicators.")
      sets <- z[, .(n_case_before = sum(case == 1L), n_referents_before = sum(case == 0L), n_case_after = sum(case == 
                                                                                                                1L & !directly_excluded), n_referents_after = sum(case == 0L & !directly_excluded)), by = id]
      if (any(sets$n_case_before != 1L | sets$n_referents_before < 1L)) 
        stop("Exclusion diagnostics must start from valid matched sets.")
      z[sets, on = "id", `:=`(case_absent = i.n_case_after != 1L, referents_absent = i.n_referents_after < 1L)]
      z[, `:=`(indirect_case_absent = !directly_excluded & case_absent, indirect_referents_absent = !directly_excluded & 
                 !case_absent & referents_absent, retained = !directly_excluded & !case_absent & !referents_absent, record_role = fifelse(case == 
                                                                                                                                            1L, "case", "referent"))]
      if (any(z[, directly_excluded + indirect_case_absent + indirect_referents_absent + retained] != 1L)) 
        stop("Exclusion categories do not partition the pre-exclusion records.")
      count_records <- function(x) {
        list(n_records_before = nrow(x), n_sets_before = uniqueN(x$id), n_no_overlap = sum(!x$day_overlap & !x$night_overlap), 
             n_day_only_overlap = sum(x$day_overlap & !x$night_overlap), n_night_only_overlap = sum(!x$day_overlap & 
                                                                                                      x$night_overlap), n_both_periods_overlap = sum(x$day_overlap & x$night_overlap), n_any_overlap = sum(x$day_overlap | 
                                                                                                                                                                                                             x$night_overlap), n_directly_excluded = sum(x$directly_excluded), n_indirectly_excluded_case_absent = sum(x$indirect_case_absent), 
             n_indirectly_excluded_all_referents_absent = sum(x$indirect_referents_absent), n_excluded_total = sum(!x$retained), 
             n_retained = sum(x$retained), n_sets_with_retained_records = uniqueN(x$id[x$retained]))
      }
      summarize_partition <- function(values, partition, source_column) {
        values <- as.character(values)
        values[is.na(values) | !nzchar(trimws(values))] <- "(Missing)"
        z[, `:=`(diagnostic_group, values)]
        roles <- z[, count_records(.SD), by = .(group_value = diagnostic_group, record_role)]
        totals <- z[, count_records(.SD), by = .(group_value = diagnostic_group)]
        totals[, `:=`(record_role, "all_records")]
        out <- rbindlist(list(totals, roles), use.names = TRUE)
        out[, `:=`(group_type = partition, grouping_variable = source_column)]
        out
      }
      month <- as.integer(format(z$date, "%m"))
      seasons <- c("Winter_DJF", "Spring_MAM", "Summer_JJA", "Autumn_SON")
      ans <- list(summarize_partition(rep("All", nrow(z)), "overall", "none"), summarize_partition(format(z$date, 
                                                                                                          "%Y"), "year", "date"), summarize_partition(sprintf("%02d", month), "month", "date"), summarize_partition(format(z$date, 
                                                                                                                                                                                                                           "%Y-%m"), "year_month", "date"), summarize_partition(seasons[(month%%12L)%/%3L + 1L], "season", "date"))
      for (field in geo_cols) ans[[length(ans) + 1L]] <- summarize_partition(z[[field]], "geography", field)
      ans <- rbindlist(ans, use.names = TRUE)
      rate_counts <- c("n_no_overlap", "n_day_only_overlap", "n_night_only_overlap", "n_both_periods_overlap", "n_any_overlap", 
                       "n_directly_excluded", "n_indirectly_excluded_case_absent", "n_indirectly_excluded_all_referents_absent", 
                       "n_excluded_total", "n_retained")
      for (field in rate_counts) ans[, `:=`((paste0(sub("^n_", "", field), "_percent")), 100 * get(field)/n_records_before)]
      ans[, `:=`(outcome = OUTCOME_NAME, outcome_code = OUTCOME_CODE, both_anomaly_policy = BOTH_ANOMALY_POLICY, 
                 denominator_population = "Complete-case valid matched records before joint-state exclusion", percentage_denominator = "n_records_before within the same group and record_role", 
                 geography_columns_used = if (length(geo_cols)) 
                   paste(geo_cols, collapse = "|")
                 else "none")]
      setcolorder(ans, c("outcome", "outcome_code", "group_type", "grouping_variable", "group_value", "record_role"))
      setorderv(ans, c("group_type", "grouping_variable", "group_value", "record_role"))
      list(table = ans, valid_ids = sets[n_case_after == 1L & n_referents_after >= 1L, id])
    }
    
    summarize_matched_state_support <- function(dt, state_levels = make_state_levels(), ref_state = REF_STATE) {
      z <- dt[, .(id, case, state = as.character(joint_state))]
      if (!nrow(z) || anyNA(z$state) || any(!z$state %in% state_levels)) 
        stop("Matched-state support requires a nonempty final sample with valid states.")
      sets <- z[, .(n_case = sum(case == 1L), n_referents = sum(case == 0L), n_states = uniqueN(state)), by = id]
      if (any(sets$n_case != 1L | sets$n_referents < 1L)) 
        stop("Matched-state support requires one case and at least one referent per set.")
      cases <- z[case == 1L, .(id, case_state = state)]
      pairs <- merge(z[case == 0L, .(id, referent_state = state)], cases, by = "id", all.x = TRUE, sort = FALSE)
      observed_pairs <- pairs[, .(n_case_referent_pairs = .N, n_matched_sets = uniqueN(id)), by = .(case_state, referent_state)]
      matrix_table <- merge(CJ(case_state = state_levels, referent_state = state_levels, sorted = FALSE), observed_pairs, 
                            by = c("case_state", "referent_state"), all.x = TRUE, sort = FALSE)
      matrix_table[is.na(n_matched_sets), `:=`(n_case_referent_pairs = 0L, n_matched_sets = 0L)]
      matrix_table[, `:=`(discordant_state_pair, case_state != referent_state)]
      matrix_table[, `:=`(case_order = match(case_state, state_levels), referent_order = match(referent_state, state_levels))]
      setorder(matrix_table, case_order, referent_order)
      matrix_table[, `:=`(c("case_order", "referent_order"), NULL)]
      presence <- unique(z[, .(id, state)])
      ref_ids <- presence[state == ref_state, id]
      observed_states <- unique(presence$state)
      edges <- observed_pairs[case_state != referent_state & n_matched_sets > 0L]
      connected <- intersect(ref_state, observed_states)
      repeat {
        expanded <- union(connected, unique(c(edges[case_state %in% connected, referent_state], edges[referent_state %in% 
                                                                                                        connected, case_state])))
        if (length(expanded) == length(connected)) 
          break
        connected <- expanded
      }
      support <- rbindlist(lapply(state_levels, function(st) {
        ids <- presence[state == st, id]
        n_forward <- uniqueN(pairs[case_state == st & referent_state != st, id])
        n_reverse <- uniqueN(pairs[case_state != st & referent_state == st, id])
        n_concordant <- sum(sets$id %in% ids & sets$n_states == 1L)
        n_informative <- n_forward + n_reverse
        if (length(ids) != n_concordant + n_informative) 
          stop("Within-state matched support counts do not reconcile.")
        data.table(state = st, n_cases = sum(cases$case_state == st), n_referent_records = sum(pairs$referent_state == 
                                                                                                 st), n_sets_containing_state = length(ids), n_sets_case_state_with_other_referent = n_forward, n_sets_other_case_with_state_referent = n_reverse, 
                   n_sets_with_state_contrast = n_informative, n_sets_entirely_in_state = n_concordant, n_sets_cooccurring_with_reference = if (st == 
                                                                                                                                                ref_state) 
                     NA_integer_
                   else sum(ids %in% ref_ids), connected_to_reference = if (length(ids)) 
                     st %in% connected
                   else NA, comparison_support = if (!length(ids)) 
                     "not_observed"
                   else if (!n_informative) 
                     "no_within_set_state_contrast"
                   else if (!st %in% connected) 
                     "disconnected_from_reference"
                   else if (!n_forward || !n_reverse) 
                     "one_case_referent_direction_only"
                   else "both_directions_observed")
      }))
      support[, `:=`(n_matched_sets_total = nrow(sets), n_sets_with_any_state_contrast = sum(sets$n_states > 1L), 
                     n_sets_without_state_contrast = sum(sets$n_states == 1L), percent_sets_without_state_contrast = 100 * mean(sets$n_states == 
                                                                                                                                  1L), interpretation = "Observed comparison support; connectivity alone does not rule out separation or imprecision")]
      for (out in list(support, matrix_table)) {
        out[, `:=`(outcome = OUTCOME_NAME, outcome_code = OUTCOME_CODE, population = "Final included valid matched sets", 
                   reference_state = ref_state)]
        setcolorder(out, c("outcome", "outcome_code"))
      }
      if (sum(matrix_table$n_case_referent_pairs) != sum(sets$n_referents)) 
        stop("Case-referent comparison counts do not reconcile.")
      list(summary = support, comparisons = matrix_table)
    }
    
    prepare_joint_model_data <- function(dt_raw) {
      required_base <- c("id", "date", "case", "holiday")
      day_cols <- sprintf("temp_day_lag%02d", LAG_DAYS)
      night_cols <- sprintf("temp_night_lag%02d", LAG_DAYS)
      rh_day_cols <- sprintf("rh_day_lag%02d", 1:3)
      rh_night_cols <- sprintf("rh_night_lag%02d", 1:3)
      heat_day_cols <- sprintf("temp_day_lag%02d", DAY_HEAT_CLASS_LAG_DAYS)
      cold_day_cols <- sprintf("temp_day_lag%02d", DAY_COLD_CLASS_LAG_DAYS)
      heat_night_cols <- sprintf("temp_night_lag%02d", NIGHT_HEAT_CLASS_LAG_DAYS)
      cold_night_cols <- sprintf("temp_night_lag%02d", NIGHT_COLD_CLASS_LAG_DAYS)
      required_cols <- c(required_base, day_cols, night_cols, rh_day_cols, rh_night_cols)
      missing_cols <- setdiff(required_cols, names(dt_raw))
      if (length(missing_cols) > 0) 
        stop("Missing required columns for the joint analysis: ", paste(missing_cols, collapse = ", "))
      retained_cols <- unique(c(required_cols, EXTRA_COVARIATES, diagnostic_geography_columns(names(dt_raw))))
      dt <- copy(dt_raw[, ..retained_cols])
      if (is.finite(TEST_N)) {
        log_msg("TEST_N is finite. Keeping complete sets whose ids occur in the first", TEST_N, "rows.")
        test_ids <- unique(dt$id[seq_len(min(TEST_N, nrow(dt)))])
        dt <- dt[id %in% test_ids]
      }
      dt[, `:=`(date, as.IDate(date))]
      dt[, `:=`(year, as.integer(format(date, "%Y")))]
      dt[, `:=`(case, coerce_binary01(case, "case"))]
      dt[, `:=`(holiday, coerce_binary01(holiday, "holiday"))]
      rh_short_cols <- c(sprintf("rh_day_lag%02d", 1:3), sprintf("rh_night_lag%02d", 1:3))
      dt[, `:=`(rh_lag01_03, rowMeans(.SD, na.rm = FALSE)), .SDcols = rh_short_cols]
      dt[, `:=`((rh_short_cols), NULL)]
      model_vars_for_complete <- c("id", "case", "holiday", "rh_lag01_03", "year", day_cols, night_cols, EXTRA_COVARIATES)
      log_msg("Filtering complete cases for the joint model...")
      cc <- dt[, complete.cases(.SD), .SDcols = model_vars_for_complete]
      dt <- dt[cc]
      rm(cc)
      gc()
      log_msg("Filtering valid matched sets for the joint model...")
      valid_id <- dt[, .(n_case = sum(case == 1), n_ref = sum(case == 0), n = .N), by = id][n_case == 1 & n_ref >= 
                                                                                              1, id]
      dt <- dt[id %in% valid_id]
      rm(valid_id)
      gc()
      if (!nrow(dt)) 
        stop("No complete valid matched sets remain for the joint model.")
      log_msg("Joint-analysis rows after filtering:", nrow(dt))
      log_msg("Joint-analysis cases:", dt[, sum(case == 1)])
      log_msg("Joint-analysis matched sets:", dt[, uniqueN(id)])
      list(dt = dt, day_cols = day_cols, night_cols = night_cols, rh_day_cols = rh_day_cols, rh_night_cols = rh_night_cols, 
           heat_day_cols = heat_day_cols, cold_day_cols = cold_day_cols, heat_night_cols = heat_night_cols, cold_night_cols = cold_night_cols)
    }
    
    run_single_temperature_models <- function(dt_raw) {
      log_msg("Single-temperature analysis started: single daytime and nighttime DLNM analyses.")
      log_msg("Output root:", OUT_ROOT)
      prep <- prepare_single_model_data(dt_raw)
      dt <- prep$dt
      day_cols <- prep$day_cols
      night_cols <- prep$night_cols
      sample_desc <- data.table(outcome = OUTCOME_NAME, n_rows = nrow(dt), n_cases = dt[, sum(case == 1)], n_referents = dt[, 
                                                                                                                            sum(case == 0)], n_matched_sets = dt[, uniqueN(id)], date_min = as.character(min(dt$date)), date_max = as.character(max(dt$date)), 
                                lag_column_suffix_min = min(LAG_DAYS), lag_column_suffix_max = max(LAG_DAYS), lag_index_min = min(LAG_INDEX), 
                                lag_index_max = max(LAG_INDEX))
      save_dt(sample_desc, file.path(DIR_TABLE, "table_00_analysis_sample.csv"))
      log_msg("Building exposure matrices and spline specifications...")
      day_mat <- as.matrix(dt[, ..day_cols])
      night_mat <- as.matrix(dt[, ..night_cols])
      temperature_history_cols <- c(day_cols, night_cols)
      dt[, `:=`((temperature_history_cols), NULL)]
      rm(temperature_history_cols)
      gc()
      day_spec <- make_ns_spec(as.vector(day_mat), VAR_DF)
      night_spec <- make_ns_spec(as.vector(night_mat), VAR_DF)
      lag_spec <- make_lag_spec_logknots(lag_index = LAG_INDEX, n_internal_knots = LOG_LAG_N_INTERNAL_KNOTS, intercept = LOG_LAG_INTERCEPT)
      method_suffix <- "logknots"
      method_label <- "Log-knots lag basis: natural spline with log-spaced lag knots and intercept"
      basis_specs <- list(method_suffix = method_suffix, method_label = method_label, day_spec = day_spec, night_spec = night_spec, 
                          lag_spec = lag_spec, VAR_DF = VAR_DF, LAG_DF = LAG_DF, RH_DF = RH_DF, LAG_DAYS = LAG_DAYS, LAG_INDEX = LAG_INDEX, 
                          lag_note = "LAG_DAYS is used for column suffixes; LAG_INDEX is used for modeling and plotting lag days.")
      save_obj(basis_specs, file.path(DIR_MODEL, "basis_specs_single_logknots.qs"))
      log_msg("Fitting daytime single-index DLNM: logknots")
      x_day <- build_cb_matrix(day_mat, day_spec, lag_spec, "day_single")
      fit_day <- fit_clogit_from_design(dt, x_day, "day_single_dlnm_logknots")
      save_obj(fit_day, file.path(DIR_MODEL, "model_day_single_dlnm_logknots.qs"))
      rm(x_day)
      gc()
      log_msg("Fitting nighttime single-index DLNM: logknots")
      x_night <- build_cb_matrix(night_mat, night_spec, lag_spec, "night_single")
      fit_night <- fit_clogit_from_design(dt, x_night, "night_single_dlnm_logknots")
      save_obj(fit_night, file.path(DIR_MODEL, "model_night_single_dlnm_logknots.qs"))
      rm(x_night)
      gc()
      analysis_list <- list(list(exposure = "daytime", exposure_label = "Daytime", prefix = "day_single", fit = fit_day, 
                                 mat = day_mat, spec = day_spec, colour = DAY_COLOUR, xlab = "Daytime temperature (°C)"), list(exposure = "nighttime", 
                                                                                                                               exposure_label = "Nighttime", prefix = "night_single", fit = fit_night, mat = night_mat, spec = night_spec, 
                                                                                                                               colour = NIGHT_COLOUR, xlab = "Nighttime temperature (°C)"))
      curve_all <- list()
      lag_all <- list()
      rr_summary_all <- list()
      mmt_summary_all <- list()
      for (aa in analysis_list) {
        log_msg("Producing single-temperature results for:", aa$exposure)
        temp_vec <- as.vector(aa$mat)
        temp_range <- quantile(temp_vec, MMT_RANGE_PROBS, na.rm = TRUE)
        temp_seq <- sort(unique(c(seq(temp_range[1], temp_range[2], by = CURVE_BY), unname(temp_range[2]))))
        mmt <- find_mmt_from_single_model(fit_obj = aa$fit, temp_seq = temp_seq, var_spec = aa$spec, lag_spec = lag_spec, 
                                          prefix = aa$prefix, temp_for_median = temp_vec)
        temp_seq <- sort(unique(c(temp_seq, mmt)))
        mmt_ci <- simulate_mmt_ci_from_single_model(fit_obj = aa$fit, temp_seq = temp_seq, var_spec = aa$spec, 
                                                    lag_spec = lag_spec, prefix = aa$prefix, n_sim = N_MMT_CI_SIM, seed = MMT_CI_SEED + ifelse(aa$exposure == 
                                                                                                                                                 "daytime", 1L, 2L))
        p025 <- as.numeric(quantile(temp_vec, 0.025, na.rm = TRUE, type = 8))
        p975 <- as.numeric(quantile(temp_vec, 0.975, na.rm = TRUE, type = 8))
        curve_dt <- make_cumulative_curve(aa$fit, temp_seq, mmt, aa$spec, lag_spec, aa$prefix)
        curve_dt[, `:=`(outcome = OUTCOME_NAME, model = "single_index", lag_method = method_suffix, lag_method_label = method_label, 
                        exposure = aa$exposure, exposure_label = aa$exposure_label, mmt = mmt, mmt_low = mmt_ci$mmt_low, mmt_high = mmt_ci$mmt_high, 
                        mmt_n_sim = mmt_ci$n_sim, p025 = p025, p975 = p975)]
        setcolorder(curve_dt, c("outcome", "model", "lag_method", "lag_method_label", "exposure", "exposure_label", 
                                "temperature", "mmt", "mmt_low", "mmt_high", "mmt_n_sim", "p025", "p975", "rr", "rr_low", "rr_high", 
                                "logrr", "se"))
        curve_all[[aa$exposure]] <- curve_dt
        plot_cumulative_curve(curve_dt, aa, method_suffix)
        lag_cold <- make_lag_curve(aa$fit, p025, mmt, aa$spec, lag_spec, aa$prefix)
        lag_cold[, `:=`(outcome = OUTCOME_NAME, model = "single_index", lag_method = method_suffix, lag_method_label = method_label, 
                        exposure = aa$exposure, exposure_label = aa$exposure_label, contrast = "P2.5_vs_MMT", contrast_label = "P2.5 vs MMT", 
                        temperature = p025, reference_temperature = mmt)]
        lag_heat <- make_lag_curve(aa$fit, p975, mmt, aa$spec, lag_spec, aa$prefix)
        lag_heat[, `:=`(outcome = OUTCOME_NAME, model = "single_index", lag_method = method_suffix, lag_method_label = method_label, 
                        exposure = aa$exposure, exposure_label = aa$exposure_label, contrast = "P97.5_vs_MMT", contrast_label = "P97.5 vs MMT", 
                        temperature = p975, reference_temperature = mmt)]
        lag_dt <- rbindlist(list(lag_cold, lag_heat), use.names = TRUE)
        setcolorder(lag_dt, c("outcome", "model", "lag_method", "lag_method_label", "exposure", "exposure_label", 
                              "contrast", "contrast_label", "lag_day", "temperature", "reference_temperature", "rr", "rr_low", "rr_high", 
                              "logrr", "se"))
        lag_all[[aa$exposure]] <- lag_dt
        plot_lag_curve(lag_dt, aa, method_suffix)
        x_ext <- build_cb_constant(c(p025, p975), aa$spec, lag_spec, aa$prefix)
        x_ref <- build_cb_constant(mmt, aa$spec, lag_spec, aa$prefix)
        rr_sum <- predict_contrast_from_design(aa$fit, x_ext, x_ref)
        rr_sum[, `:=`(outcome = OUTCOME_NAME, model = "single_index", lag_method = method_suffix, lag_method_label = method_label, 
                      exposure = aa$exposure, exposure_label = aa$exposure_label, contrast = c("P2.5_vs_MMT", "P97.5_vs_MMT"), 
                      contrast_label = c("P2.5 vs MMT", "P97.5 vs MMT"), temperature = c(p025, p975), reference_temperature = mmt, 
                      mmt = mmt, mmt_low = mmt_ci$mmt_low, mmt_high = mmt_ci$mmt_high, mmt_n_sim = mmt_ci$n_sim)]
        setcolorder(rr_sum, c("outcome", "model", "lag_method", "lag_method_label", "exposure", "exposure_label", 
                              "contrast", "contrast_label", "temperature", "reference_temperature", "mmt", "rr", "rr_low", "rr_high", 
                              "logrr", "se"))
        rr_summary_all[[aa$exposure]] <- rr_sum
        mmt_summary_all[[aa$exposure]] <- data.table(outcome = OUTCOME_NAME, model = "single_index", lag_method = method_suffix, 
                                                     lag_method_label = method_label, exposure = aa$exposure, exposure_label = aa$exposure_label, mmt = mmt, 
                                                     mmt_low = mmt_ci$mmt_low, mmt_high = mmt_ci$mmt_high, mmt_n_sim = mmt_ci$n_sim, mmt_percentile = mean(temp_vec <= 
                                                                                                                                                             mmt, na.rm = TRUE) * 100, p025 = p025, p975 = p975, var_df = VAR_DF, lag_df_original_setting = LAG_DF, 
                                                     lag_spec_type = lag_spec$type, lag_internal_knots = paste(round(lag_spec$knots, 6), collapse = ";"), 
                                                     lag_intercept = lag_spec$intercept, rh_df = RH_DF, lag_column_suffix_min = min(LAG_DAYS), lag_column_suffix_max = max(LAG_DAYS), 
                                                     lag_index_min = min(LAG_INDEX), lag_index_max = max(LAG_INDEX))
      }
      curve_all_dt <- rbindlist(curve_all, use.names = TRUE)
      lag_all_dt <- rbindlist(lag_all, use.names = TRUE)
      rr_summary_dt <- rbindlist(rr_summary_all, use.names = TRUE)
      mmt_summary_dt <- rbindlist(mmt_summary_all, use.names = TRUE)
      save_obj(curve_all_dt, file.path(DIR_PLOT_DATA, "plotdata_cumulative_curves_single_logknots.qs"))
      save_dt(curve_all_dt, file.path(DIR_PLOT_DATA, "plotdata_cumulative_curves_single_logknots.csv"))
      save_obj(lag_all_dt, file.path(DIR_PLOT_DATA, "plotdata_lag_effects_single_logknots.qs"))
      save_dt(lag_all_dt, file.path(DIR_PLOT_DATA, "plotdata_lag_effects_single_logknots.csv"))
      save_dt(rr_summary_dt, file.path(DIR_TABLE, "table_01_cumulative_rr_p025_p975_vs_mmt_single_logknots.csv"))
      save_dt(lag_all_dt, file.path(DIR_TABLE, "table_02_lag_specific_rr_p025_p975_vs_mmt_single_logknots.csv"))
      save_dt(mmt_summary_dt, file.path(DIR_TABLE, "table_03_mmt_summary_single_logknots.csv"))
      analysis_meta <- list(outcome = OUTCOME_NAME, data_path = DATA_PATH, output_root = OUT_ROOT, method_suffix = method_suffix, 
                            method_label = method_label, n_rows_model = nrow(dt), n_cases = dt[, sum(case == 1)], settings = list(VAR_DF = VAR_DF, 
                                                                                                                                  LAG_DF = LAG_DF, RH_DF = RH_DF, LAG_DAYS = LAG_DAYS, LAG_INDEX = LAG_INDEX, LOG_LAG_N_INTERNAL_KNOTS = LOG_LAG_N_INTERNAL_KNOTS, 
                                                                                                                                  LOG_LAG_INTERCEPT = LOG_LAG_INTERCEPT, N_MMT_CI_SIM = N_MMT_CI_SIM, lag_spec = lag_spec, lag_note = "LAG_DAYS is used for column suffixes; LAG_INDEX is used for modeling and plotting lag days."), 
                            mmt_summary = mmt_summary_dt, cumulative_rr_summary = rr_summary_dt)
      save_obj(analysis_meta, file.path(DIR_MODEL, "analysis_metadata_single_logknots.qs"))
      result <- list(curve_all_dt = curve_all_dt, lag_all_dt = lag_all_dt, rr_summary_dt = rr_summary_dt, mmt_summary_dt = mmt_summary_dt, 
                     fit_day = fit_day, fit_night = fit_night, basis_specs = basis_specs, day_spec = day_spec, night_spec = night_spec, 
                     lag_spec = lag_spec)
      save_obj(result, file.path(DIR_MODEL, "analysis_results_single_logknots.qs"))
      log_msg("Single daytime and nighttime analyses finished successfully.")
      log_msg("Single-exposure cumulative and lag-response outputs saved.")
      invisible(result)
    }
    
    format_lag_window <- function(lag_days) {
      lag_days <- sort(unique(as.integer(lag_days)))
      if (length(lag_days) == 0L) 
        return(NA_character_)
      if (length(lag_days) == 1L) 
        return(paste0("lag", lag_days))
      paste0("lag", min(lag_days), "-lag", max(lag_days))
    }
    
    format_lag_window_from_suffix <- function(lag_suffix) {
      format_lag_window(as.integer(lag_suffix) - 1L)
    }
    
    validate_fallback_window <- function(x, name) {
      ok <- is.numeric(x) && length(x) > 0L && !anyNA(x) && all(is.finite(x)) && all(x == as.integer(x)) && all(x %in% 
                                                                                                                  LAG_INDEX) && identical(as.integer(x), seq.int(min(x), max(x)))
      if (!ok) {
        stop(name, " must be ordered consecutive integer model lags within 0:20.")
      }
      invisible(TRUE)
    }
    
    select_one_lag_window <- function(lag_dt, exposure_value, contrast_value, component, thermal, alpha = AUTO_LAG_ALPHA, 
                                      selection_rule = AUTO_LAG_SELECTION_RULE, no_sig_action = AUTO_LAG_NO_SIG_ACTION, fallback_heat = AUTO_LAG_FALLBACK_HEAT, 
                                      fallback_cold = AUTO_LAG_FALLBACK_COLD) {
      if (!identical(selection_rule, "zero_to_last")) {
        stop("selection_rule must be 'zero_to_last'.")
      }
      if (!no_sig_action %in% c("fallback", "stop")) {
        stop("no_sig_action must be either 'fallback' or 'stop'.")
      }
      if (!thermal %in% c("cold", "heat")) {
        stop("thermal must be either 'cold' or 'heat'.")
      }
      if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
        stop("alpha must be one finite value between 0 and 1.")
      }
      validate_fallback_window(fallback_heat, "fallback_heat")
      validate_fallback_window(fallback_cold, "fallback_cold")
      required_cols <- c("exposure", "contrast", "lag_day", "rr", "rr_low", "rr_high", "logrr", "se", "temperature", 
                         "reference_temperature")
      missing_cols <- setdiff(required_cols, names(lag_dt))
      if (length(missing_cols) > 0L) {
        stop("Lag-response results are missing columns: ", paste(missing_cols, collapse = ", "))
      }
      x <- copy(lag_dt[exposure == exposure_value & contrast == contrast_value])
      setorder(x, lag_day)
      if (nrow(x) != length(LAG_INDEX) || anyDuplicated(x$lag_day) || !identical(as.integer(x$lag_day), as.integer(LAG_INDEX))) {
        stop("Expected exactly ", length(LAG_INDEX), " lag estimates (", min(LAG_INDEX), ":", max(LAG_INDEX), ") for ", 
             component, ".")
      }
      if (any(!is.finite(x$rr)) || any(!is.finite(x$rr_low)) || any(!is.finite(x$rr_high)) || any(!is.finite(x$logrr)) || 
          any(!is.finite(x$se)) || any(x$rr_low <= 0) || any(x$rr_low > x$rr) || any(x$rr > x$rr_high)) {
        stop("Non-finite or invalid lag estimates/confidence limits for ", component, ".")
      }
      if (uniqueN(x$temperature) != 1L || uniqueN(x$reference_temperature) != 1L || any(!is.finite(x$temperature)) || 
          any(!is.finite(x$reference_temperature))) {
        stop("Invalid temperature contrast for ", component, ".")
      }
      if (thermal == "cold" && x$temperature[1] >= x$reference_temperature[1]) {
        stop("Cold contrast temperature is not below MMT for ", component, ".")
      }
      if (thermal == "heat" && x$temperature[1] <= x$reference_temperature[1]) {
        stop("Heat contrast temperature is not above MMT for ", component, ".")
      }
      zcrit <- qnorm(1 - alpha/2)
      x[, `:=`(selection_ci_low, exp(logrr - zcrit * se))]
      x[, `:=`(selection_ci_high, exp(logrr + zcrit * se))]
      x[, `:=`(significant_harmful, selection_ci_low > 1)]
      sig_lags <- as.integer(x[significant_harmful == TRUE, lag_day])
      fallback_used <- length(sig_lags) == 0L
      if (!fallback_used) {
        selected_lags <- seq.int(0L, max(sig_lags))
        selection_reason <- "zero_to_last_significant_harmful_lag"
      }
      else {
        if (identical(no_sig_action, "stop")) {
          stop("No harmful significant lag (two-sided CI lower bound > 1) was found for ", component, ". Single-temperature results have been saved; the joint model was not run.")
        }
        selected_lags <- if (thermal == "heat") 
          fallback_heat
        else fallback_cold
        selected_lags <- as.integer(selected_lags)
        selection_reason <- paste0("prespecified_fallback_", thermal, "_", format_lag_window(selected_lags), "_no_significant_harmful_lag")
      }
      selected_lags <- as.integer(selected_lags)
      selected_suffix <- selected_lags + 1L
      selected_flag <- x$lag_day %in% selected_lags
      selected_non_significant <- selected_flag & !x$significant_harmful
      summary <- data.table(outcome = OUTCOME_NAME, exposure = exposure_value, thermal = thermal, contrast = contrast_value, 
                            component = component, significance_rule = paste0(round((1 - alpha) * 100, 1), "% CI lower bound > 1"), 
                            lag_selection_rule = "lag0-to-last-significant continuous window", status = if (fallback_used) 
                              "prespecified_fallback"
                            else "auto_selected", n_significant_harmful_lags = length(sig_lags), significant_harmful_lags = if (length(sig_lags)) 
                              paste(sig_lags, collapse = ",")
                            else "", selected_lag_start = min(selected_lags), selected_lag_end = max(selected_lags), selected_lag_days = paste(selected_lags, 
                                                                                                                                               collapse = ","), selected_data_suffix = paste(selected_suffix, collapse = ","), selected_window = format_lag_window(selected_lags), 
                            n_selected_lags = length(selected_lags), n_selected_non_significant = sum(selected_non_significant), significant_at_lag20 = any(sig_lags == 
                                                                                                                                                              max(LAG_INDEX)), fallback_used = fallback_used, selection_reason = selection_reason, temperature = x$temperature[1], 
                            reference_temperature = x$reference_temperature[1])
      diagnostics <- copy(x)
      diagnostics[, `:=`(outcome = OUTCOME_NAME, thermal = thermal, component = component, selected = lag_day %in% 
                           selected_lags, selected_non_significant = (lag_day %in% selected_lags) & !significant_harmful, fallback_used = fallback_used, 
                         selection_reason = selection_reason)]
      list(lag_days = selected_lags, data_suffix = selected_suffix, summary = summary, diagnostics = diagnostics)
    }
    
    apply_auto_lag_windows <- function(auto_lag) {
      AUTO_LAG_SELECTION_SUMMARY <<- copy(auto_lag$summary)
      DAY_COLD_CLASS_LAG_DAYS <<- auto_lag$day_cold$data_suffix
      DAY_HEAT_CLASS_LAG_DAYS <<- auto_lag$day_heat$data_suffix
      NIGHT_COLD_CLASS_LAG_DAYS <<- auto_lag$night_cold$data_suffix
      NIGHT_HEAT_CLASS_LAG_DAYS <<- auto_lag$night_heat$data_suffix
      all_suffix <- c(DAY_COLD_CLASS_LAG_DAYS, DAY_HEAT_CLASS_LAG_DAYS, NIGHT_COLD_CLASS_LAG_DAYS, NIGHT_HEAT_CLASS_LAG_DAYS)
      if (any(!all_suffix %in% LAG_DAYS)) {
        stop("At least one selected lag falls outside the available lag columns.")
      }
      invisible(TRUE)
    }
    
    make_state_levels <- function() {
      c("D0_N0", "DC_N0", "DH_N0", "D0_NC", "DC_NC", "DH_NC", "D0_NH", "DC_NH", "DH_NH")
    }
    
    state_dt_from_levels <- function(state_levels) {
      dt_state <- data.table(state = state_levels)
      dt_state[, `:=`(day_code, sub("_.*$", "", state))]
      dt_state[, `:=`(night_code, sub("^.*_", "", state))]
      dt_state[, `:=`(D_C, as.integer(day_code == "DC"))]
      dt_state[, `:=`(D_H, as.integer(day_code == "DH"))]
      dt_state[, `:=`(N_C, as.integer(night_code == "NC"))]
      dt_state[, `:=`(N_H, as.integer(night_code == "NH"))]
      dt_state[, `:=`(day_status, fifelse(day_code == "D0", "No daytime anomaly", fifelse(day_code == "DC", "Daytime cold", 
                                                                                          "Daytime heat")))]
      dt_state[, `:=`(night_status, fifelse(night_code == "N0", "No nighttime anomaly", fifelse(night_code == "NC", 
                                                                                                "Nighttime cold", "Nighttime heat")))]
      dt_state[, `:=`(day_axis, factor(day_status, levels = c("No daytime anomaly", "Daytime cold", "Daytime heat")))]
      dt_state[, `:=`(night_axis, factor(night_status, levels = c("No nighttime anomaly", "Nighttime cold", "Nighttime heat")))]
      dt_state[, `:=`(day_order, match(day_code, c("D0", "DC", "DH")) - 1L)]
      dt_state[, `:=`(night_order, match(night_code, c("N0", "NC", "NH")) - 1L)]
      dt_state[, `:=`(mask, day_order + 3L * night_order)]
      dt_state[, `:=`(component_mask, D_C * 1L + D_H * 2L + N_C * 4L + N_H * 8L)]
      dt_state[, `:=`(c("day_order", "night_order"), NULL)]
      dt_state[]
    }
    
    COMPONENTS <- c("D_C", "D_H", "N_C", "N_H")
    
    COMPONENT_BITS <- c(D_C = 1L, D_H = 2L, N_C = 4L, N_H = 8L)
    
    COMPONENT_INFO <- data.table(component = c(COMPONENTS, "joint_total"), component_order = 1:5, component_label = c("Daytime cold", 
                                                                                                                      "Daytime heat", "Nighttime cold", "Nighttime heat", "Joint total"), component_group = c("daytime_cold", "daytime_heat", 
                                                                                                                                                                                                              "nighttime_cold", "nighttime_heat", "joint_total"))
    
    state_from_mask <- function(mask) {
      day_levels <- c("D0", "DC", "DH")
      night_levels <- c("N0", "NC", "NH")
      day_i <- mask%%3L
      night_i <- mask%/%3L
      paste0(day_levels[day_i + 1L], "_", night_levels[night_i + 1L])
    }
    
    mask_from_daynight_codes <- function(day_code, night_code) {
      day_i <- match(day_code, c("D0", "DC", "DH")) - 1L
      night_i <- match(night_code, c("N0", "NC", "NH")) - 1L
      as.integer(day_i + 3L * night_i)
    }
    
    resolve_period_status <- function(cold, heat, period_name, policy = BOTH_ANOMALY_POLICY) {
      code <- rep(NA_character_, length(cold))
      code[!cold & !heat] <- if (period_name == "day") 
        "D0"
      else "N0"
      code[cold & !heat] <- if (period_name == "day") 
        "DC"
      else "NC"
      code[!cold & heat] <- if (period_name == "day") 
        "DH"
      else "NH"
      both <- cold & heat
      if (any(both, na.rm = TRUE)) {
        if (policy == "exclude") {
          code[both] <- NA_character_
        }
        else if (policy == "heat_priority") {
          code[both] <- if (period_name == "day") 
            "DH"
          else "NH"
        }
        else if (policy == "cold_priority") {
          code[both] <- if (period_name == "day") 
            "DC"
          else "NC"
        }
        else {
          stop("Unknown BOTH_ANOMALY_POLICY: ", policy)
        }
      }
      code
    }
    
    build_joint_indicators <- function(dt, day_cold_threshold, day_heat_threshold, night_cold_threshold, night_heat_threshold) {
      threshold_values <- c(day_cold_threshold, day_heat_threshold, night_cold_threshold, night_heat_threshold)
      if (any(!is.finite(threshold_values))) {
        stop("Joint-state reference thresholds must all be finite.")
      }
      out <- copy(dt)
      out[, `:=`(D_C_raw, as.integer(day_cold_temp < day_cold_threshold))]
      out[, `:=`(D_H_raw, as.integer(day_heat_temp > day_heat_threshold))]
      out[, `:=`(N_C_raw, as.integer(night_cold_temp < night_cold_threshold))]
      out[, `:=`(N_H_raw, as.integer(night_heat_temp > night_heat_threshold))]
      out[, `:=`(day_code, resolve_period_status(D_C_raw == 1, D_H_raw == 1, "day", BOTH_ANOMALY_POLICY))]
      out[, `:=`(night_code, resolve_period_status(N_C_raw == 1, N_H_raw == 1, "night", BOTH_ANOMALY_POLICY))]
      out[, `:=`(D_C, as.integer(day_code == "DC"))]
      out[, `:=`(D_H, as.integer(day_code == "DH"))]
      out[, `:=`(N_C, as.integer(night_code == "NC"))]
      out[, `:=`(N_H, as.integer(night_code == "NH"))]
      out[, `:=`(joint_state, fifelse(is.na(day_code) | is.na(night_code), NA_character_, paste0(day_code, "_", night_code)))]
      out[, `:=`(joint_mask, mask_from_daynight_codes(day_code, night_code))]
      attr(out, "state_thresholds") <- list(strategy = STATE_THRESHOLD_STRATEGY, strategy_tag = reference_strategy_tag(), 
                                            reference_method = CURRENT_SHARED_DESIGN$settings$reference_method, reference_source_outcome = CURRENT_SHARED_DESIGN$source_name, 
                                            reference_source_code = CURRENT_SHARED_DESIGN$source_code, calibration_id = CURRENT_SHARED_DESIGN$calibration_id, 
                                            rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance, day_cold_threshold = day_cold_threshold, day_heat_threshold = day_heat_threshold, 
                                            night_cold_threshold = night_cold_threshold, night_heat_threshold = night_heat_threshold)
      out[]
    }
    
    make_beta_state_vector <- function(coefficients, covariance, state_levels, ref_state) {
      cf <- coefficients
      vc <- covariance
      beta_state <- setNames(rep(NA_real_, length(state_levels)), state_levels)
      beta_state[ref_state] <- 0
      state_se <- setNames(rep(NA_real_, length(state_levels)), state_levels)
      state_se[ref_state] <- 0
      for (st in setdiff(state_levels, ref_state)) {
        term <- paste0("joint_state", st)
        if (term %in% names(cf)) {
          beta_state[st] <- unname(cf[term])
          state_se[st] <- sqrt(vc[term, term])
        }
      }
      list(beta_state = beta_state, state_se = state_se)
    }
    
    make_rr_table <- function(beta_state, state_se, state_info) {
      out <- copy(state_info)
      out[, `:=`(beta, beta_state[state])]
      out[, `:=`(se, state_se[state])]
      out[, `:=`(rr, fifelse(is.finite(beta), exp(beta), NA_real_))]
      out[, `:=`(rr_low, fifelse(is.finite(beta) & is.finite(se), exp(beta - 1.96 * se), NA_real_))]
      out[, `:=`(rr_high, fifelse(is.finite(beta) & is.finite(se), exp(beta + 1.96 * se), NA_real_))]
      out[state == REF_STATE, `:=`(beta = 0, se = 0, rr = 1, rr_low = 1, rr_high = 1)]
      out[, `:=`(estimable, is.finite(beta) & is.finite(se))]
      out[, `:=`(rr_label, fifelse(state == REF_STATE, "1.00\nReference", fifelse(estimable, sprintf("%.2f\n(%.2f, %.2f)", 
                                                                                                     rr, rr_low, rr_high), "NE\nNot observed")))]
      out[]
    }
    
    run_joint_temperature_model <- function(dt_raw, single_logknots_results) {
      log_msg("Joint analysis started: 3 x 3 theoretical joint-state framework with dynamic support and Shapley burden decomposition.")
      log_msg("Output root:", OUT_ROOT)
      log_msg("Joint-state threshold strategy:", STATE_THRESHOLD_STRATEGY)
      log_msg("Reference source:", CURRENT_SHARED_DESIGN$source_name)
      log_msg("Reference definition:", reference_strategy_tag())
      log_msg("BOTH_ANOMALY_POLICY:", BOTH_ANOMALY_POLICY)
      prep <- prepare_joint_model_data(dt_raw)
      dt <- prep$dt
      day_cols <- prep$day_cols
      night_cols <- prep$night_cols
      heat_day_cols <- prep$heat_day_cols
      cold_day_cols <- prep$cold_day_cols
      heat_night_cols <- prep$heat_night_cols
      cold_night_cols <- prep$cold_night_cols
      desc_basic <- data.table(outcome = OUTCOME_NAME, model_label = MODEL_LABEL, n_rows = nrow(dt), n_cases = dt[, 
                                                                                                                  sum(case == 1)], n_referents = dt[, sum(case == 0)], n_matched_sets = dt[, uniqueN(id)], date_min = as.character(min(dt$date)), 
                               date_max = as.character(max(dt$date)), day_heat_lag_window = format_lag_window_from_suffix(DAY_HEAT_CLASS_LAG_DAYS), 
                               day_cold_lag_window = format_lag_window_from_suffix(DAY_COLD_CLASS_LAG_DAYS), night_heat_lag_window = format_lag_window_from_suffix(NIGHT_HEAT_CLASS_LAG_DAYS), 
                               night_cold_lag_window = format_lag_window_from_suffix(NIGHT_COLD_CLASS_LAG_DAYS), state_threshold_strategy = STATE_THRESHOLD_STRATEGY, 
                               threshold_strategy_tag = reference_strategy_tag(), reference_method = CURRENT_SHARED_DESIGN$settings$reference_method, 
                               reference_source_outcome = CURRENT_SHARED_DESIGN$source_name, reference_source_code = CURRENT_SHARED_DESIGN$source_code, 
                               calibration_id = CURRENT_SHARED_DESIGN$calibration_id, rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance, 
                               both_anomaly_policy = BOTH_ANOMALY_POLICY, 
                               mmt_lag_knot_type = LAG_KNOT_TYPE, mmt_lag_n_internal_knots = LAG_N_INTERNAL_KNOTS, mmt_lag_intercept = LAG_INTERCEPT, 
                               n_sim = N_SIM, cap_negative_af = CAP_NEGATIVE_AF)
      save_dt(desc_basic, file.path(DIR_TABLE, "table_00_basic_analysis_sample_before_joint_state_filter.csv"))
      single_mmt_summary_logknots <- single_logknots_results$mmt_summary_dt
      if (is.null(single_mmt_summary_logknots)) {
        stop("MMT summary not found in the single-temperature results.")
      }
      basis_specs <- single_logknots_results$basis_specs
      day_spec <- single_logknots_results$day_spec
      night_spec <- single_logknots_results$night_spec
      lag_spec <- single_logknots_results$lag_spec
      if (is.null(basis_specs) || is.null(day_spec) || is.null(night_spec) || is.null(lag_spec)) {
        stop("The single-temperature analysis did not return complete spline specifications.")
      }
      save_obj(basis_specs, file.path(DIR_MODEL, "basis_specs_for_mmt.qs"))
      mmt_day <- single_mmt_summary_logknots[exposure == "daytime", mmt][1]
      mmt_night <- single_mmt_summary_logknots[exposure == "nighttime", mmt][1]
      mmt_day_low <- single_mmt_summary_logknots[exposure == "daytime", mmt_low][1]
      mmt_day_high <- single_mmt_summary_logknots[exposure == "daytime", mmt_high][1]
      mmt_night_low <- single_mmt_summary_logknots[exposure == "nighttime", mmt_low][1]
      mmt_night_high <- single_mmt_summary_logknots[exposure == "nighttime", mmt_high][1]
      if (any(!is.finite(c(mmt_day, mmt_night, mmt_day_low, mmt_day_high, mmt_night_low, mmt_night_high)))) {
        stop("At least one MMT or MMT confidence limit is not finite.")
      }
      mmt_summary <- data.table(outcome = OUTCOME_NAME, model_label = MODEL_LABEL, day_mmt = mmt_day, day_mmt_low = mmt_day_low, 
                                day_mmt_high = mmt_day_high, night_mmt = mmt_night, night_mmt_low = mmt_night_low, night_mmt_high = mmt_night_high, 
                                state_threshold_strategy = STATE_THRESHOLD_STRATEGY, threshold_strategy_tag = reference_strategy_tag(), 
                                reference_method = CURRENT_SHARED_DESIGN$settings$reference_method, reference_source_outcome = CURRENT_SHARED_DESIGN$source_name, 
                                reference_source_code = CURRENT_SHARED_DESIGN$source_code, calibration_id = CURRENT_SHARED_DESIGN$calibration_id, 
                                rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance, var_df = VAR_DF, rh_df = RH_DF, lag_knot_type = LAG_KNOT_TYPE, 
                                lag_n_internal_knots = LAG_N_INTERNAL_KNOTS, lag_intercept = LAG_INTERCEPT, lag_knots = paste(round(lag_spec$knots, 
                                                                                                                                    4), collapse = ","), lag_boundary_knots = paste(round(lag_spec$Boundary.knots, 4), collapse = ","), 
                                mmt_source = paste("Overall-outcome single-temperature log-knots DLNM:", CURRENT_SHARED_DESIGN$source_name))
      save_dt(mmt_summary, file.path(DIR_TABLE, "table_01_mmt_summary.csv"))
      save_obj(mmt_summary, file.path(DIR_MODEL, "mmt_summary.qs"))
      daytime_curve_for_mmt <- copy(single_logknots_results$curve_all_dt[exposure == "daytime"])
      nighttime_curve_for_mmt <- copy(single_logknots_results$curve_all_dt[exposure == "nighttime"])
      daytime_curve_for_mmt[, `:=`(exposure, "daytime_temperature")]
      nighttime_curve_for_mmt[, `:=`(exposure, "nighttime_temperature")]
      save_dt(daytime_curve_for_mmt, file.path(DIR_PLOT_DATA, "plotdata_daytime_temperature_single_curve_for_mmt.csv"))
      save_obj(daytime_curve_for_mmt, file.path(DIR_PLOT_DATA, "plotdata_daytime_temperature_single_curve_for_mmt.qs"))
      save_dt(nighttime_curve_for_mmt, file.path(DIR_PLOT_DATA, "plotdata_nighttime_temperature_single_curve_for_mmt.csv"))
      save_obj(nighttime_curve_for_mmt, file.path(DIR_PLOT_DATA, "plotdata_nighttime_temperature_single_curve_for_mmt.qs"))
      log_msg("Daytime MMT:", round(mmt_day, 3), "95% CI:", round(mmt_day_low, 3), "to", round(mmt_day_high, 3))
      log_msg("Nighttime MMT:", round(mmt_night, 3), "95% CI:", round(mmt_night_low, 3), "to", round(mmt_night_high, 
                                                                                                     3))
      log_msg("Joint-state threshold strategy:", STATE_THRESHOLD_STRATEGY)
      log_msg("Reference source:", CURRENT_SHARED_DESIGN$source_name)
      log_msg("Reference definition:", reference_strategy_tag())
      log_msg("Building heat and cold lag-window summaries...")
      dt[, `:=`(day_heat_temp, rowMeans(.SD, na.rm = FALSE)), .SDcols = heat_day_cols]
      dt[, `:=`(day_cold_temp, rowMeans(.SD, na.rm = FALSE)), .SDcols = cold_day_cols]
      dt[, `:=`(night_heat_temp, rowMeans(.SD, na.rm = FALSE)), .SDcols = heat_night_cols]
      dt[, `:=`(night_cold_temp, rowMeans(.SD, na.rm = FALSE)), .SDcols = cold_night_cols]
      temperature_history_cols <- c(day_cols, night_cols)
      dt[, `:=`((temperature_history_cols), NULL)]
      rm(temperature_history_cols)
      gc()
      STATE_LEVELS <- make_state_levels()
      STATE_INFO <- state_dt_from_levels(STATE_LEVELS)
      state_thresholds <- compute_joint_state_thresholds()
      log_msg("Daytime cold threshold:", round(state_thresholds$day_cold_threshold, 3), "| lower shared reference bound; applied to", 
              format_lag_window_from_suffix(DAY_COLD_CLASS_LAG_DAYS), "mean temperature")
      log_msg("Daytime heat threshold:", round(state_thresholds$day_heat_threshold, 3), "| upper shared reference bound; applied to", 
              format_lag_window_from_suffix(DAY_HEAT_CLASS_LAG_DAYS), "mean temperature")
      log_msg("Nighttime cold threshold:", round(state_thresholds$night_cold_threshold, 3), "| lower shared reference bound; applied to", 
              format_lag_window_from_suffix(NIGHT_COLD_CLASS_LAG_DAYS), "mean temperature")
      log_msg("Nighttime heat threshold:", round(state_thresholds$night_heat_threshold, 3), "| upper shared reference bound; applied to", 
              format_lag_window_from_suffix(NIGHT_HEAT_CLASS_LAG_DAYS), "mean temperature")
      dt <- build_joint_indicators(dt, day_cold_threshold = state_thresholds$day_cold_threshold, day_heat_threshold = state_thresholds$day_heat_threshold, 
                                   night_cold_threshold = state_thresholds$night_cold_threshold, night_heat_threshold = state_thresholds$night_heat_threshold)
      threshold_attr <- attr(dt, "state_thresholds")
      n_before_joint_filter <- nrow(dt)
      n_case_before_joint_filter <- dt[, sum(case == 1)]
      n_excluded_both_anomaly <- dt[, sum(is.na(joint_state))]
      overlap_diagnostics <- data.table(outcome = OUTCOME_NAME, threshold_strategy_tag = state_thresholds$strategy_tag, 
                                        reference_method = CURRENT_SHARED_DESIGN$settings$reference_method, reference_source_outcome = CURRENT_SHARED_DESIGN$source_name, 
                                        reference_source_code = CURRENT_SHARED_DESIGN$source_code, calibration_id = CURRENT_SHARED_DESIGN$calibration_id, 
                                        rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance, n_rows = nrow(dt), n_cases = dt[, sum(case == 1)], day_cold_heat_overlap_rows = dt[, 
                                                                                                                                                                    sum(D_C_raw == 1L & D_H_raw == 1L)], day_cold_heat_overlap_percent = dt[, mean(D_C_raw == 1L & D_H_raw == 
                                                                                                                                                                                                                                                     1L) * 100], day_cold_heat_overlap_cases = dt[case == 1L, sum(D_C_raw == 1L & D_H_raw == 1L)], day_cold_heat_overlap_case_percent = dt[case == 
                                                                                                                                                                                                                                                                                                                                                                                             1L, mean(D_C_raw == 1L & D_H_raw == 1L) * 100], night_cold_heat_overlap_rows = dt[, sum(N_C_raw == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       1L & N_H_raw == 1L)], night_cold_heat_overlap_percent = dt[, mean(N_C_raw == 1L & N_H_raw == 1L) * 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    100], night_cold_heat_overlap_cases = dt[case == 1L, sum(N_C_raw == 1L & N_H_raw == 1L)], night_cold_heat_overlap_case_percent = dt[case == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          1L, mean(N_C_raw == 1L & N_H_raw == 1L) * 100], excluded_joint_state_rows = n_excluded_both_anomaly)
      save_dt(overlap_diagnostics, file.path(DIR_TABLE, "table_00a_cold_heat_overlap_diagnostics.csv"))
      exclusion_details <- summarize_overlap_exclusions(dt)
      exclusion_details$table[, `:=`(calibration_id = CURRENT_SHARED_DESIGN$calibration_id, reference_method = CURRENT_SHARED_DESIGN$settings$reference_method)]
      save_dt(exclusion_details$table, file.path(DIR_TABLE, "table_00d_overlap_exclusion_distribution.csv"))
      log_msg("Exclusion diagnostic geography fields:", exclusion_details$table$geography_columns_used[1L])
      ref_diagnostics <- dt[joint_state == REF_STATE, .(n_ref_state_rows = .N, n_ref_state_cases = sum(case == 1), 
                                                        n_ref_state_referents = sum(case == 0), n_ref_state_sets = uniqueN(id))]
      reference_diagnostics <- data.table(state_threshold_strategy = STATE_THRESHOLD_STRATEGY, threshold_strategy_tag = state_thresholds$strategy_tag, 
                                          reference_method = CURRENT_SHARED_DESIGN$settings$reference_method, reference_source_outcome = CURRENT_SHARED_DESIGN$source_name, 
                                          reference_source_code = CURRENT_SHARED_DESIGN$source_code, calibration_id = CURRENT_SHARED_DESIGN$calibration_id, 
                                          rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance, day_cold_threshold = state_thresholds$day_cold_threshold, day_heat_threshold = state_thresholds$day_heat_threshold, 
                                          night_cold_threshold = state_thresholds$night_cold_threshold, night_heat_threshold = state_thresholds$night_heat_threshold, 
                                          n_rows_after_joint_state_filter = n_before_joint_filter - n_excluded_both_anomaly, n_excluded_both_anomaly = n_excluded_both_anomaly, 
                                          n_ref_state_rows = ref_diagnostics$n_ref_state_rows, n_ref_state_cases = ref_diagnostics$n_ref_state_cases, 
                                          n_ref_state_referents = ref_diagnostics$n_ref_state_referents, n_ref_state_sets = ref_diagnostics$n_ref_state_sets)
      save_dt(reference_diagnostics, file.path(DIR_TABLE, "table_00_reference_state_diagnostics.csv"))
      dt <- dt[!is.na(joint_state)]
      valid_id2 <- exclusion_details$valid_ids
      dt <- dt[id %in% valid_id2]
      rm(valid_id2, exclusion_details)
      gc()
      log_msg("Rows before removing within-period cold+heat records:", n_before_joint_filter)
      log_msg("Cases before removing within-period cold+heat records:", n_case_before_joint_filter)
      if (!nrow(dt) || !any(dt$case == 1L)) 
        stop("No valid matched sets remain after joint-state exclusion.")
      save_dt(data.table(outcome = OUTCOME_NAME, source_reference = CURRENT_SHARED_DESIGN$source_name, n_rows_before_joint_filter = n_before_joint_filter, 
                         n_cases_before_joint_filter = n_case_before_joint_filter, n_rows_removed_for_cold_heat_overlap = n_excluded_both_anomaly, 
                         n_rows_final = nrow(dt), n_cases_final = sum(dt$case == 1L), n_cases_lost_total = n_case_before_joint_filter - 
                           sum(dt$case == 1L), case_retention_percent = 100 * sum(dt$case == 1L)/n_case_before_joint_filter, burden_denominator = "Final included cases in valid matched sets"), 
              file.path(DIR_TABLE, "table_00b_joint_analysis_sample_flow.csv"))
      log_msg("Rows after joint-state filter and matched-set refilter:", nrow(dt))
      log_msg("Cases after joint-state filter and matched-set refilter:", dt[, sum(case == 1)])
      dt[, `:=`(joint_state, factor(joint_state, levels = STATE_LEVELS))]
      dt[, `:=`(joint_state, relevel(joint_state, ref = REF_STATE))]
      threshold_summary <- data.table(outcome = OUTCOME_NAME, model_label = MODEL_LABEL, both_anomaly_policy = BOTH_ANOMALY_POLICY, 
                                      state_threshold_strategy = STATE_THRESHOLD_STRATEGY, threshold_strategy_tag = state_thresholds$strategy_tag, 
                                      reference_method = CURRENT_SHARED_DESIGN$settings$reference_method, reference_source_outcome = CURRENT_SHARED_DESIGN$source_name, 
                                      reference_source_code = CURRENT_SHARED_DESIGN$source_code, calibration_id = CURRENT_SHARED_DESIGN$calibration_id, 
                                      rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance, day_cold_threshold = state_thresholds$day_cold_threshold, day_heat_threshold = state_thresholds$day_heat_threshold, 
                                      night_cold_threshold = state_thresholds$night_cold_threshold, night_heat_threshold = state_thresholds$night_heat_threshold, 
                                      threshold_distribution = state_thresholds$threshold_distribution, day_mmt = mmt_day, day_mmt_low = mmt_day_low, 
                                      day_mmt_high = mmt_day_high, night_mmt = mmt_night, night_mmt_low = mmt_night_low, night_mmt_high = mmt_night_high, 
                                      day_heat_lag_days_data_suffix = paste(DAY_HEAT_CLASS_LAG_DAYS, collapse = ","), day_cold_lag_days_data_suffix = paste(DAY_COLD_CLASS_LAG_DAYS, 
                                                                                                                                                            collapse = ","), night_heat_lag_days_data_suffix = paste(NIGHT_HEAT_CLASS_LAG_DAYS, collapse = ","), 
                                      night_cold_lag_days_data_suffix = paste(NIGHT_COLD_CLASS_LAG_DAYS, collapse = ","), day_heat_lag_window = format_lag_window_from_suffix(DAY_HEAT_CLASS_LAG_DAYS), 
                                      day_cold_lag_window = format_lag_window_from_suffix(DAY_COLD_CLASS_LAG_DAYS), night_heat_lag_window = format_lag_window_from_suffix(NIGHT_HEAT_CLASS_LAG_DAYS), 
                                      night_cold_lag_window = format_lag_window_from_suffix(NIGHT_COLD_CLASS_LAG_DAYS), mmt_lag_knot_type = LAG_KNOT_TYPE, 
                                      mmt_lag_n_internal_knots = LAG_N_INTERNAL_KNOTS, mmt_lag_intercept = LAG_INTERCEPT, mmt_lag_knots = paste(round(lag_spec$knots, 
                                                                                                                                                      4), collapse = ","), mmt_lag_boundary_knots = paste(round(lag_spec$Boundary.knots, 4), collapse = ","))
      save_dt(threshold_summary, file.path(DIR_TABLE, "table_02_joint_state_thresholds.csv"))
      save_obj(threshold_summary, file.path(DIR_MODEL, "joint_state_thresholds.qs"))
      state_counts <- dt[, .(n_rows = .N, n_cases = sum(case == 1), n_referents = sum(case == 0), n_sets = uniqueN(id)), 
                         by = .(joint_state)]
      setnames(state_counts, "joint_state", "state")
      state_counts[, `:=`(state, as.character(state))]
      state_counts <- merge(STATE_INFO, state_counts, by = "state", all.x = TRUE, sort = FALSE)
      state_counts[is.na(n_rows), `:=`(n_rows = 0L, n_cases = 0L, n_referents = 0L, n_sets = 0L)]
      state_counts[, `:=`(support_status, fifelse(n_rows == 0L, "not_observed", fifelse(n_cases == 0L | n_referents == 
                                                                                          0L, "observed_but_not_estimable", "observed_estimable")))]
      setorder(state_counts, mask)
      save_dt(state_counts, file.path(DIR_TABLE, "table_03_joint_state_counts.csv"))
      matched_support <- summarize_matched_state_support(dt, STATE_LEVELS, REF_STATE)
      save_dt(matched_support$summary, file.path(DIR_TABLE, "table_03c_matched_state_comparison_support.csv"))
      save_dt(matched_support$comparisons, file.path(DIR_TABLE, "table_03d_case_referent_state_comparisons.csv"))
      log_msg("Matched sets with a within-set state contrast:", matched_support$summary$n_sets_with_any_state_contrast[1L], 
              "of", matched_support$summary$n_matched_sets_total[1L])
      limited_support <- matched_support$summary[!comparison_support %in% c("not_observed", "both_directions_observed"), 
                                                 state]
      if (length(limited_support)) 
        log_msg("Inspect within-set support diagnostics for:", paste(limited_support, collapse = ", "))
      rm(matched_support, limited_support)
      log_msg("Joint-state counts before model fitting:")
      for (ii in seq_len(nrow(state_counts))) {
        log_msg("  ", state_counts$state[ii], "| rows:", state_counts$n_rows[ii], "| cases:", state_counts$n_cases[ii], 
                "| referents:", state_counts$n_referents[ii], "| sets:", state_counts$n_sets[ii], "| support:", state_counts$support_status[ii])
      }
      if (state_counts[state == REF_STATE, n_rows] == 0L) {
        stop("The reference state D0_N0 is not observed and the joint model cannot be fitted.")
      }
      unobserved_states <- state_counts[support_status == "not_observed", state]
      unsupported_observed_states <- state_counts[support_status == "observed_but_not_estimable", state]
      if (length(unsupported_observed_states) > 0L) {
        stop("The following states contain observations but have no cases or no referents and cannot be safely omitted: ", 
             paste(unsupported_observed_states, collapse = ", "), ". Inspect table_03_joint_state_counts.csv.")
      }
      observed_state_levels <- STATE_LEVELS[!STATE_LEVELS %in% unobserved_states]
      estimable_non_ref_states <- setdiff(observed_state_levels, REF_STATE)
      if (length(unobserved_states) > 0L) {
        log_msg("Structurally unobserved joint states will remain in reporting tables as NE but will not be included in the conditional logistic regression:", 
                paste(unobserved_states, collapse = ", "))
      }
      log_msg("Joint model will estimate", length(estimable_non_ref_states), "non-reference state coefficients; global Wald df =", 
              length(estimable_non_ref_states))
      log_msg("Fitting joint categorical model:", MODEL_LABEL)
      joint_model_dt <- data.table(case = dt$case, id = dt$id, joint_state = factor(as.character(dt$joint_state), 
                                                                                    levels = observed_state_levels), holiday = dt$holiday, rh_lag01_03 = dt$rh_lag01_03)
      joint_model_dt[, `:=`(joint_state, relevel(joint_state, ref = REF_STATE))]
      rh_basis <- as.data.table(ns(joint_model_dt$rh_lag01_03, df = RH_DF))
      setnames(rh_basis, paste0("rh_ns", seq_len(ncol(rh_basis))))
      joint_model_dt <- cbind(joint_model_dt, rh_basis)
      pollution <- build_pollution_adjustment(dt, ADJUSTMENT_COVARIATES)
      if (length(ADJUSTMENT_COVARIATES)) {
        joint_model_dt <- cbind(joint_model_dt, pollution$data)
        save_dt(pollution$summary, file.path(DIR_TABLE, "table_pollution_covariate_summary.csv"))
      }
      rhs_terms <- c("joint_state", names(rh_basis), "holiday", ADJUSTMENT_COVARIATES)
      form_joint <- as.formula(paste0("case ~ ", paste(rhs_terms, collapse = " + "), " + strata(id)"))
      formula_joint_text <- paste(deparse(form_joint), collapse = " ")
      fit_joint <- clogit(form_joint, data = joint_model_dt, method = CLOGIT_METHOD, model = FALSE, x = FALSE, y = FALSE)
      coef_joint <- coef(fit_joint)
      if (length(ADJUSTMENT_COVARIATES) && (any(!ADJUSTMENT_COVARIATES %in% names(coef_joint)) || any(!is.finite(coef_joint[ADJUSTMENT_COVARIATES])))) 
        stop("At least one selected pollution coefficient is not estimable.")
      vcov_joint <- vcov(fit_joint)
      if (length(ADJUSTMENT_COVARIATES)) {
        nuisance_table <- data.table(covariate = ADJUSTMENT_COVARIATES, beta = unname(coef_joint[ADJUSTMENT_COVARIATES]), 
                                     se = sqrt(diag(vcov_joint)[ADJUSTMENT_COVARIATES]))
        save_dt(nuisance_table, file.path(DIR_TABLE, "table_pollution_adjustment_coefficients.csv"))
      }
      required_joint_terms <- paste0("joint_state", estimable_non_ref_states)
      names(required_joint_terms) <- estimable_non_ref_states
      missing_joint_terms <- setdiff(required_joint_terms, names(coef_joint))
      available_joint_terms <- intersect(required_joint_terms, names(coef_joint))
      nonfinite_joint_terms <- available_joint_terms[!is.finite(coef_joint[available_joint_terms])]
      if (length(missing_joint_terms) > 0L || length(nonfinite_joint_terms) > 0L) {
        failed_joint_compact <- list(formula_text = formula_joint_text, model_label = MODEL_LABEL, state_levels = STATE_LEVELS, 
                                     observed_state_levels = observed_state_levels, estimable_non_reference_states = estimable_non_ref_states, 
                                     unobserved_states = unobserved_states, ref_state = REF_STATE, coefficients = coef_joint, covariance = vcov_joint, 
                                     missing_joint_terms = missing_joint_terms, nonfinite_joint_terms = nonfinite_joint_terms, state_counts = state_counts, 
                                     threshold_summary = threshold_summary, storage_mode = "compact_failed_joint_model_diagnostics")
        save_obj(failed_joint_compact, file.path(DIR_MODEL, paste0("model_joint_", MODEL_LABEL, "_FAILED_COMPACT.qs")))
        rm(fit_joint, joint_model_dt, rh_basis, form_joint)
        invisible(gc(full = TRUE))
        stop("At least one OBSERVED joint state was not estimable after removing structurally unobserved states. Missing terms: ", 
             paste(missing_joint_terms, collapse = ", "), "; non-finite terms: ", paste(nonfinite_joint_terms, collapse = ", "), 
             ". This indicates separation or insufficient within-stratum support rather than a structurally empty state.")
      }
      beta_info <- make_beta_state_vector(coefficients = coef_joint, covariance = vcov_joint, state_levels = STATE_LEVELS, 
                                          ref_state = REF_STATE)
      beta_state <- beta_info$beta_state
      state_se <- beta_info$state_se
      covariance_state <- matrix(NA_real_, nrow = length(STATE_LEVELS), ncol = length(STATE_LEVELS), dimnames = list(STATE_LEVELS, 
                                                                                                                     STATE_LEVELS))
      covariance_state[REF_STATE, REF_STATE] <- 0
      if (length(estimable_non_ref_states) > 0L) {
        estimable_terms <- required_joint_terms[estimable_non_ref_states]
        covariance_state[estimable_non_ref_states, estimable_non_ref_states] <- vcov_joint[estimable_terms, estimable_terms, 
                                                                                           drop = FALSE]
      }
      mu_non_reference <- setNames(as.numeric(coef_joint[required_joint_terms]), estimable_non_ref_states)
      covariance_non_reference <- vcov_joint[required_joint_terms, required_joint_terms, drop = FALSE]
      dimnames(covariance_non_reference) <- list(estimable_non_ref_states, estimable_non_ref_states)
      joint_model_compact <- list(coefficients = coef_joint, covariance = vcov_joint, beta_state = beta_state, covariance_state = covariance_state, 
                                  mu_non_reference = mu_non_reference, covariance_non_reference = covariance_non_reference, theoretical_state_levels = STATE_LEVELS, 
                                  state_levels = STATE_LEVELS, observed_state_levels = observed_state_levels, estimable_non_reference_states = estimable_non_ref_states, 
                                  unobserved_states = unobserved_states, n_estimable_non_reference = length(estimable_non_ref_states), state_terms = required_joint_terms, 
                                  formula_text = formula_joint_text, model_label = MODEL_LABEL, both_anomaly_policy = BOTH_ANOMALY_POLICY, 
                                  ref_state = REF_STATE, state_info = STATE_INFO, state_counts = state_counts, threshold_summary = threshold_summary, 
                                  state_threshold_strategy = STATE_THRESHOLD_STRATEGY, threshold_strategy_tag = state_thresholds$strategy_tag, 
                                  reference_method = CURRENT_SHARED_DESIGN$settings$reference_method, reference_source_outcome = CURRENT_SHARED_DESIGN$source_name, 
                                  reference_source_code = CURRENT_SHARED_DESIGN$source_code, calibration_id = CURRENT_SHARED_DESIGN$calibration_id, 
                                  rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance, state_thresholds = state_thresholds, mmt_lag_spec = lag_spec, 
                                  storage_mode = "compact_coefficients_and_covariance_only_dynamic_state_support")
      save_obj(joint_model_compact, file.path(DIR_MODEL, paste0("model_joint_", MODEL_LABEL, ".qs")))
      rm(fit_joint, joint_model_dt, rh_basis, form_joint)
      invisible(gc(full = TRUE))
      log_msg("Finished joint categorical model; retained compact coefficients/covariance only.")
      required_states <- unique(as.character(dt[case == 1, joint_state]))
      missing_required <- required_states[!is.finite(beta_state[required_states])]
      if (length(missing_required) > 0) {
        stop("Some states required for burden calculation are not estimable in the joint model: ", paste(unique(missing_required), 
                                                                                                         collapse = ", "), ". Check sparse states or the selected shared reference thresholds.")
      }
      rr_table <- make_rr_table(beta_state, state_se, STATE_INFO)
      rr_table <- merge(rr_table, state_counts[, .(state, n_rows, n_cases, n_referents, n_sets, support_status)], 
                        by = "state", all.x = TRUE, sort = FALSE)
      setorder(rr_table, mask)
      save_dt(rr_table, file.path(DIR_TABLE, "table_04_joint_9_state_rr.csv"))
      save_obj(rr_table, file.path(DIR_MODEL, "joint_9_state_rr.qs"))
      case_dt <- dt[case == 1, .(id, date, year, joint_state = as.character(joint_state), joint_mask)]
      years <- sort(unique(case_dt$year))
      counts_by_year <- setNames(vector("list", length(years)), as.character(years))
      n_deaths_by_year <- setNames(vector("list", length(years)), as.character(years))
      for (yy in years) {
        tmp <- case_dt[year == yy, .N, by = joint_state]
        counts <- setNames(rep(0, length(STATE_LEVELS)), STATE_LEVELS)
        counts[tmp$joint_state] <- tmp$N
        counts_by_year[[as.character(yy)]] <- counts
        n_deaths_by_year[[as.character(yy)]] <- nrow(case_dt[year == yy])
      }
      annual_counts <- rbindlist(lapply(names(counts_by_year), function(y) {
        data.table(year = as.integer(y), state = STATE_LEVELS, n_cases = as.numeric(counts_by_year[[y]][STATE_LEVELS]), 
                   n_deaths = as.numeric(n_deaths_by_year[[y]]))
      }))
      save_dt(annual_counts, file.path(DIR_TABLE, "table_00c_annual_state_death_counts.csv"))
      export_covariance_csv(covariance_non_reference, file.path(DIR_TABLE, "table_04a_joint_logrr_covariance.csv"))
      analysis_meta <- list(outcome = OUTCOME_NAME, data_path = DATA_PATH, output_root = OUT_ROOT, model_label = MODEL_LABEL, 
                            both_anomaly_policy = BOTH_ANOMALY_POLICY, n_rows_model = nrow(dt), n_cases = dt[, sum(case == 1)], n_matched_sets = dt[, 
                                                                                                                                                    uniqueN(id)], mmt_summary = mmt_summary, threshold_summary = threshold_summary, state_levels = STATE_LEVELS, 
                            observed_state_levels = observed_state_levels, estimable_non_reference_states = estimable_non_ref_states, 
                            unobserved_states = unobserved_states, n_estimable_non_reference = length(estimable_non_ref_states), reference_state = REF_STATE, 
                            state_threshold_strategy = STATE_THRESHOLD_STRATEGY, threshold_strategy_tag = state_thresholds$strategy_tag, 
                            reference_method = CURRENT_SHARED_DESIGN$settings$reference_method, reference_source_outcome = CURRENT_SHARED_DESIGN$source_name, 
                            reference_source_code = CURRENT_SHARED_DESIGN$source_code, calibration_id = CURRENT_SHARED_DESIGN$calibration_id, 
                            rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance, state_thresholds = state_thresholds, overlap_diagnostics = overlap_diagnostics, 
                            lag_response_setting_for_mmt = list(knot_type = LAG_KNOT_TYPE, n_internal_knots = LAG_N_INTERNAL_KNOTS, 
                                                                intercept = LAG_INTERCEPT, knots = lag_spec$knots, boundary_knots = lag_spec$Boundary.knots, note = "Natural cubic spline with intercept and three internal knots placed at equally spaced values on the log scale."), 
                            day_heat_class_lag_days = DAY_HEAT_CLASS_LAG_DAYS, day_cold_class_lag_days = DAY_COLD_CLASS_LAG_DAYS, night_heat_class_lag_days = NIGHT_HEAT_CLASS_LAG_DAYS, 
                            night_cold_class_lag_days = NIGHT_COLD_CLASS_LAG_DAYS, auto_lag_selection = list(alpha = AUTO_LAG_ALPHA, 
                                                                                                             selection_rule = AUTO_LAG_SELECTION_RULE, no_significant_lag_action = AUTO_LAG_NO_SIG_ACTION, fallback_heat_model_lags = AUTO_LAG_FALLBACK_HEAT, 
                                                                                                             fallback_cold_model_lags = AUTO_LAG_FALLBACK_COLD, selection_summary = AUTO_LAG_SELECTION_SUMMARY), 
                            n_sim = N_SIM, sim_seed = SIM_SEED, cap_negative_af = CAP_NEGATIVE_AF, shared_design = CURRENT_SHARED_DESIGN, 
                            uncertainty_scope = "Joint-model coefficient Monte Carlo; conditional on shared thresholds and lag selection; case counts treated as fixed", 
                            mean_annual_AF_definition = "Arithmetic mean of annual AF", burden_population = "Final included cases after complete-case, state-overlap and matched-set filtering", 
                            analysis_stage = "risk_estimation", n_years = length(years))
      ready <- write_stage1_input_bundle(rr_table, annual_counts, case_dt, joint_model_compact, analysis_meta)
      log_msg("Stage-one risk estimation and input preparation completed:", OUT_ROOT)
      invisible(list(rr_table = rr_table, joint_model_compact = joint_model_compact, annual_counts = annual_counts, 
                     threshold_summary = threshold_summary, analysis_meta = analysis_meta, ready = ready))
    }
    
    validate_scalar_choice <- function(x, allowed, label) {
      if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% allowed) 
        stop(label, " must be one of: ", paste(allowed, collapse = ", "))
    }
    
    validate_ma_days <- function(x, required_names, allow_partial = FALSE) {
      if (allow_partial && is.null(x)) 
        return(invisible(TRUE))
      if (!is.numeric(x) || !length(x) || is.null(names(x)) || anyNA(x) || any(!is.finite(x)) || any(x != floor(x)) || 
          any(x < 1L | x > length(LAG_INDEX)) || anyDuplicated(names(x)) || any(!names(x) %in% required_names) || 
          (!allow_partial && !setequal(names(x), required_names))) {
        stop("Moving-average days must be named integers between 1 and ", length(LAG_INDEX), "; required names: ", 
             paste(required_names, collapse = ", "))
      }
      invisible(TRUE)
    }
    
    validate_run_settings <- function(reference_method, rr_tolerance, lag_mode, manual_days, overrides, no_sig_action, 
                                      fallback_days, allow_truncated) {
      validate_scalar_choice(reference_method, c("mmt_ci", "rr_tolerance"), "reference_method")
      validate_scalar_choice(lag_mode, c("auto", "manual"), "lag_mode")
      validate_scalar_choice(no_sig_action, c("fallback", "stop"), "no_sig_action")
      validate_scalar_choice(BOTH_ANOMALY_POLICY, c("exclude", "heat_priority", "cold_priority"), "BOTH_ANOMALY_POLICY")
      if (!is.numeric(rr_tolerance) || length(rr_tolerance) != 1L || !is.finite(rr_tolerance) || rr_tolerance <= 
          0) 
        stop("rr_tolerance must be positive, e.g. 0.01.")
      if (!is.logical(allow_truncated) || length(allow_truncated) != 1L || is.na(allow_truncated)) 
        stop("allow_truncated_rr_band must be TRUE or FALSE.")
      keys <- c("day_cold", "day_heat", "night_cold", "night_heat")
      validate_ma_days(overrides, keys, TRUE)
      if (lag_mode == "manual") 
        validate_ma_days(manual_days, keys)
      validate_ma_days(fallback_days, c("cold", "heat"))
      for (v in list(N_SIM, N_MMT_CI_SIM)) if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 2 || v != 
                                               floor(v)) 
        stop("N_SIM and N_MMT_CI_SIM must be integers >= 2.")
      if (!is.numeric(TEST_N) || length(TEST_N) != 1L || is.na(TEST_N) || TEST_N <= 0 || (is.finite(TEST_N) && TEST_N != 
                                                                                          floor(TEST_N))) 
        stop("TEST_N must be a positive integer or Inf.")
      invisible(TRUE)
    }
    
    cumulative_design_fast <- function(x, var_spec, lag_spec, prefix) {
      xb <- predict_ns_basis(x, var_spec)
      lb <- colSums(predict_ns_basis(LAG_INDEX, lag_spec))
      out <- do.call(cbind, lapply(seq_len(ncol(xb)), function(a) outer(xb[, a], lb)))
      colnames(out) <- unlist(lapply(seq_len(ncol(xb)), function(a) paste0(prefix, "_v", a, "_l", seq_along(lb))))
      out
    }
    
    refine_mmt <- function(temp_seq, eta_grid, eta_function) {
      n <- length(temp_seq)
      if (n < 3L || any(!is.finite(eta_grid))) 
        stop("Invalid MMT search grid.")
      mids <- 2L:(n - 1L)
      candidates <- mids[eta_grid[mids] <= eta_grid[mids - 1L] & eta_grid[mids] <= eta_grid[mids + 1L]]
      roots <- vapply(candidates, function(i) optimize(eta_function, c(temp_seq[i - 1L], temp_seq[i + 1L]), tol = 1e-07)$minimum, 
                      numeric(1))
      candidates_t <- unique(c(temp_seq[c(1L, n)], roots))
      values <- vapply(candidates_t, eta_function, numeric(1))
      candidates_t[which.min(values)]
    }
    
    reference_strategy_tag <- function() {
      if (is.null(CURRENT_SHARED_DESIGN)) 
        stop("Shared design has not been initialized.")
      CURRENT_SHARED_DESIGN$reference_tag
    }
    
    compute_joint_state_thresholds <- function(dt = NULL) {
      if (is.null(CURRENT_SHARED_DESIGN)) 
        stop("Shared design has not been initialized.")
      CURRENT_SHARED_DESIGN$thresholds
    }
    
    rr_tolerance_band <- function(curve, mmt, epsilon, allow_truncated = FALSE) {
      x <- copy(curve[, .(temperature, logrr)])
      setorder(x, temperature)
      if (anyDuplicated(x$temperature) || any(!is.finite(unlist(x)))) 
        stop("Invalid cumulative curve.")
      k <- which.min(abs(x$temperature - mmt))
      if (abs(x$temperature[k] - mmt) > 1e-06) 
        stop("MMT must be present in the prediction grid.")
      target <- log1p(epsilon)
      inside <- x$logrr <= target + 1e-12
      if (!inside[k]) 
        stop("MMT is outside RR-tolerance band.")
      left <- right <- k
      while (left > 1L && inside[left - 1L]) left <- left - 1L
      while (right < nrow(x) && inside[right + 1L]) right <- right + 1L
      cross <- function(a, b) {
        x$temperature[a] + (target - x$logrr[a]) * (x$temperature[b] - x$temperature[a])/(x$logrr[b] - x$logrr[a])
      }
      low <- if (left == 1L) 
        x$temperature[1L]
      else cross(left - 1L, left)
      high <- if (right == nrow(x)) 
        x$temperature[nrow(x)]
      else cross(right, right + 1L)
      truncated <- c(left = left == 1L, right = right == nrow(x))
      if (any(truncated) && !allow_truncated) {
        stop("RR-tolerance interval reaches the MMT search-domain boundary. Inspect the overall curve; ", "RUN_ALLOW_TRUNCATED_RR_BAND = TRUE explicitly permits a truncated band.")
      }
      list(low = low, high = high, truncated_left = truncated[1L], truncated_right = truncated[2L])
    }
    
    make_reference_bands <- function(calibration, method, epsilon, allow_truncated) {
      validate_scalar_choice(method, c("mmt_ci", "rr_tolerance"), "reference_method")
      if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0) 
        stop("RR tolerance must be one positive finite number, e.g. 0.01.")
      sr <- calibration$single_results
      bands <- lapply(c("daytime", "nighttime"), function(period) {
        m <- sr$mmt_summary_dt[exposure == period]
        curve <- sr$curve_all_dt[exposure == period]
        if (nrow(m) != 1L) 
          stop("Expected one overall MMT per period.")
        band <- if (method == "mmt_ci") {
          list(low = m$mmt_low, high = m$mmt_high, truncated_left = m$mmt_low <= min(curve$temperature) + 1e-07, 
               truncated_right = m$mmt_high >= max(curve$temperature) - 1e-07)
        }
        else rr_tolerance_band(curve, m$mmt, epsilon, allow_truncated)
        if (any(!is.finite(c(band$low, band$high))) || band$low >= band$high) 
          stop("Empty or invalid overall reference band for ", period, "; inspect the overall MMT/curve.")
        data.table(source_outcome = calibration$source_name, source_code = calibration$source_code, calibration_id = calibration$calibration_id, 
                   exposure = period, method = method, rr_tolerance = if (method == "rr_tolerance") 
                     epsilon
                   else NA_real_, mmt = m$mmt, mmt_low = m$mmt_low, mmt_high = m$mmt_high, reference_low = band$low, reference_high = band$high, 
                   contains_point_mmt = band$low <= m$mmt && band$high >= m$mmt, lower_search_limit = min(curve$temperature), 
                   upper_search_limit = max(curve$temperature), truncated_left = as.logical(band$truncated_left), truncated_right = as.logical(band$truncated_right), 
                   definition = if (method == "mmt_ci") 
                     "Operational band from overall MMT 95% CI"
                   else "Connected set containing overall MMT with fitted cumulative RR <= 1 + epsilon")
      })
      rbindlist(bands)
    }
    
    resolve_shared_lags <- function(single_results, lag_mode, manual_days, overrides, no_sig_action, fallback_days) {
      keys <- c("day_cold", "day_heat", "night_cold", "night_heat")
      validate_scalar_choice(lag_mode, c("auto", "manual"), "lag_mode")
      validate_scalar_choice(no_sig_action, c("fallback", "stop"), "no_sig_action")
      validate_ma_days(fallback_days, c("cold", "heat"))
      validate_ma_days(overrides, keys, allow_partial = TRUE)
      if (lag_mode == "manual") 
        validate_ma_days(manual_days, keys)
      specs <- list(day_cold = c("daytime", "P2.5_vs_MMT", "Daytime cold", "cold"), day_heat = c("daytime", "P97.5_vs_MMT", 
                                                                                                 "Daytime heat", "heat"), night_cold = c("nighttime", "P2.5_vs_MMT", "Nighttime cold", "cold"), night_heat = c("nighttime", 
                                                                                                                                                                                                               "P97.5_vs_MMT", "Nighttime heat", "heat"))
      ans <- setNames(vector("list", length(keys)), keys)
      for (key in keys) {
        z <- specs[[key]]
        forced <- if (lag_mode == "manual") 
          manual_days[[key]]
        else NULL
        if (key %in% names(overrides)) 
          forced <- overrides[[key]]
        if (is.null(forced)) {
          a <- select_one_lag_window(single_results$lag_all_dt, z[1L], z[2L], z[3L], z[4L], alpha = AUTO_LAG_ALPHA, 
                                     selection_rule = "zero_to_last", no_sig_action = no_sig_action, fallback_heat = seq_len(fallback_days[["heat"]]) - 
                                       1L, fallback_cold = seq_len(fallback_days[["cold"]]) - 1L)
        }
        else {
          d <- copy(single_results$lag_all_dt[exposure == z[1L] & contrast == z[2L]])
          setorder(d, lag_day)
          if (nrow(d) != length(LAG_INDEX)) 
            stop("Missing overall lag diagnostics for ", key)
          d[, `:=`(selection_ci_low, exp(logrr - qnorm(1 - AUTO_LAG_ALPHA/2) * se))]
          d[, `:=`(selection_ci_high, exp(logrr + qnorm(1 - AUTO_LAG_ALPHA/2) * se))]
          d[, `:=`(significant_harmful, selection_ci_low > 1)]
          lags <- seq_len(forced) - 1L
          reason <- if (lag_mode == "manual") 
            "manual_moving_average_days"
          else "manual_component_override"
          d[, `:=`(thermal = z[4L], component = z[3L], selected = lag_day %in% lags, selected_non_significant = (lag_day %in% 
                                                                                                                   lags) & !significant_harmful, fallback_used = FALSE, selection_reason = reason)]
          sig <- d[significant_harmful == TRUE, lag_day]
          s <- data.table(outcome = OVERALL_NAME, exposure = z[1L], thermal = z[4L], contrast = z[2L], component = z[3L], 
                          significance_rule = "Not used for manual selection", lag_selection_rule = "manual lag0 to lag(days - 1)", 
                          status = "manual", n_significant_harmful_lags = length(sig), significant_harmful_lags = paste(sig, 
                                                                                                                        collapse = ","), selected_lag_start = 0L, selected_lag_end = max(lags), selected_lag_days = paste(lags, 
                                                                                                                                                                                                                          collapse = ","), selected_data_suffix = paste(lags + 1L, collapse = ","), selected_window = format_lag_window(lags), 
                          n_selected_lags = length(lags), n_selected_non_significant = sum(d$selected_non_significant), significant_at_lag20 = max(LAG_INDEX) %in% 
                            sig, fallback_used = FALSE, selection_reason = reason, temperature = d$temperature[1L], reference_temperature = d$reference_temperature[1L])
          a <- list(lag_days = lags, data_suffix = lags + 1L, summary = s, diagnostics = d)
        }
        a$summary[, `:=`(outcome = OVERALL_NAME, source_outcome = OVERALL_NAME, component_key = key)]
        a$diagnostics[, `:=`(outcome = OVERALL_NAME, source_outcome = OVERALL_NAME, component_key = key)]
        ans[[key]] <- a
      }
      ans$summary <- rbindlist(lapply(ans[keys], `[[`, "summary"), fill = TRUE)
      ans$diagnostics <- rbindlist(lapply(ans[keys], `[[`, "diagnostics"), fill = TRUE)
      ans
    }
    
    make_shared_design <- function(calibration, reference_method, rr_tolerance, lag_mode, manual_days, manual_overrides, 
                                   no_sig_action, fallback_days, allow_truncated) {
      bands <- make_reference_bands(calibration, reference_method, rr_tolerance, allow_truncated)
      lags <- resolve_shared_lags(calibration$single_results, lag_mode, manual_days, manual_overrides, no_sig_action, 
                                  fallback_days)
      ref_tag <- if (reference_method == "mmt_ci") 
        "MMT_CI95"
      else paste0("RRtol_", gsub("[.]", "p", format(rr_tolerance, scientific = FALSE, trim = TRUE)))
      keys <- c("day_cold", "day_heat", "night_cold", "night_heat")
      lag_sizes <- vapply(lags[keys], function(x) length(x$lag_days), integer(1))
      selection_modes <- vapply(lags[keys], function(x) x$summary$status[1L], character(1))
      lag_tag <- paste0(lag_mode, "_", paste(paste0(c("DC", "DH", "NC", "NH"), lag_sizes), collapse = "_"))
      if (!is.null(manual_overrides)) 
        lag_tag <- paste0(lag_tag, "_override_", paste(names(manual_overrides), collapse = "-"))
      if (any(selection_modes == "prespecified_fallback")) 
        lag_tag <- paste0(lag_tag, "_fallback")
      day <- bands[exposure == "daytime"]
      night <- bands[exposure == "nighttime"]
      thresholds <- list(strategy = paste0("overall_", reference_method), strategy_tag = ref_tag, reference_method = reference_method, 
                         rr_tolerance = if (reference_method == "rr_tolerance") rr_tolerance else NA_real_, source_outcome = calibration$source_name, 
                         source_code = calibration$source_code, calibration_id = calibration$calibration_id, cold_quantile_prob = NA_real_, 
                         heat_quantile_prob = NA_real_, day_cold_threshold = day$reference_low, day_heat_threshold = day$reference_high, 
                         night_cold_threshold = night$reference_low, night_heat_threshold = night$reference_high, threshold_distribution = "Overall-outcome full-lag cumulative DLNM; same thresholds applied to the four selected moving-average exposures")
      settings <- list(reference_method = reference_method, rr_tolerance = rr_tolerance, lag_mode = lag_mode, manual_days = if (lag_mode == 
                                                                                                                                "manual") manual_days else NULL, manual_overrides = manual_overrides, no_sig_action = no_sig_action, fallback_days = fallback_days, 
                       alpha = AUTO_LAG_ALPHA, allow_truncated = allow_truncated, overlap_policy = BOTH_ANOMALY_POLICY, n_sim = N_SIM, 
                       sim_seed = SIM_SEED, cap_negative_af = CAP_NEGATIVE_AF)
      list(reference_tag = ref_tag, lag_tag = lag_tag, analysis_tag = paste(ref_tag, lag_tag, BOTH_ANOMALY_POLICY, 
                                                                            sep = "__"), bands = bands, lags = lags, thresholds = thresholds, settings = settings, calibration_id = calibration$calibration_id, 
           calibration_signature = calibration$signature, source_code = calibration$source_code, source_name = calibration$source_name, 
           calibration_file = calibration$cache_file)
    }
    
    export_shared_design <- function(design) {
      save_dt(design$bands, file.path(DIR_TABLE, "table_02a_shared_reference_intervals.csv"))
      s <- copy(design$lags$summary)
      d <- copy(design$lags$diagnostics)
      s[, `:=`(applied_outcome, OUTCOME_NAME)]
      d[, `:=`(applied_outcome, OUTCOME_NAME)]
      save_dt(s, file.path(DIR_TABLE, "table_03a_auto_selected_lag_windows.csv"))
      save_dt(d, file.path(DIR_TABLE, "table_03b_auto_lag_selection_diagnostics.csv"))
      save_obj(s, file.path(DIR_MODEL, "auto_selected_lag_windows.qs"))
      save_obj(d, file.path(DIR_MODEL, "auto_lag_selection_diagnostics.qs"))
      save_obj(design, file.path(DIR_MODEL, "shared_analysis_design.qs"))
    }
    
    export_covariance_csv <- function(cv, file) {
      cols <- setNames(lapply(seq_len(ncol(cv)), function(j) formatC(cv[, j], digits = 17L, format = "g")), colnames(cv))
      save_dt(as.data.table(c(list(state = rownames(cv)), cols)), file)
    }
    
    write_stage1_input_bundle <- function(rr_table, annual_counts, case_dt, compact, metadata) {
      folder <- DIR_TABLE
      required_tables <- c("table_04_joint_9_state_rr.csv", "table_03_joint_state_counts.csv", "table_02_joint_state_thresholds.csv", 
                           "table_02a_shared_reference_intervals.csv", "table_03a_auto_selected_lag_windows.csv", "table_00a_cold_heat_overlap_diagnostics.csv", 
                           "table_00b_joint_analysis_sample_flow.csv", "table_00d_overlap_exclusion_distribution.csv", "table_03c_matched_state_comparison_support.csv", 
                           "table_03d_case_referent_state_comparisons.csv")
      for (nm in required_tables) {
        if (!file.exists(file.path(folder, nm))) 
          stop("Required stage-one table missing: ", nm)
      }
      save_dt(annual_counts, file.path(folder, "table_00c_annual_state_death_counts.csv"))
      cv <- compact$covariance_non_reference
      export_covariance_csv(cv, file.path(folder, "table_04a_joint_logrr_covariance.csv"))
      saveRDS(list(covariance_non_reference = cv, mu_non_reference = compact$mu_non_reference), file.path(DIR_MODEL, 
                                                                                                          "joint_logrr_covariance.rds"))
      metadata$outcome_code <- OUTCOME_CODE
      metadata$stage1_id <- paste0(OUTCOME_CODE, "_", format(Sys.time(), "%Y%m%dT%H%M%OS6"), "_", Sys.getpid())
      metadata$pipeline_protocol <- "daynight_two_stage_v1"
      metadata$object_format <- OBJECT_FORMAT
      metadata$save_case_contrib <- SAVE_CASE_CONTRIB
      metadata$plot_settings <- list(font_family = FONT_FAMILY, fig_dpi = FIG_DPI, save_pdf = SAVE_PDF)
      metadata$pipeline_index_file <- PIPELINE_INDEX_FILE
      metadata$stage1_output_root <- normalizePath(OUT_ROOT, winslash = "/", mustWork = TRUE)
      metadata$stage1_input_directory <- normalizePath(folder, winslash = "/", mustWork = TRUE)
      saveRDS(metadata, file.path(DIR_MODEL, "source_parameters.rds"))
      case_path <- file.path(DIR_MODEL, "case_states.rds")
      if (SAVE_CASE_CONTRIB) {
        saveRDS(case_dt, case_path)
      }
      else if (file.exists(case_path) && !file.remove(case_path)) {
        stop("Cannot clear an obsolete stage-one case-state file.")
      }
      writeLines(c("Current analysis outputs and downstream calculation inputs.", "table_04 contains state-specific risk estimates; table_00c contains annual final-sample death counts.", 
                   "joint_logrr_covariance.rds stores the full-precision covariance; the CSV provides an inspection copy.", 
                   "source_parameters.rds stores reference definitions, lag windows, sample information and computational settings.", 
                   "table_00d describes overlap patterns, direct exclusions, subsequent matched-set losses and retention.", 
                   "Its percentages use pre-exclusion records within each group and record role; missing geography is retained.", 
                   "Calendar distributions use record dates; winter is December-February. Different partitions must not be added together.", 
                   "table_03c reports final-sample within-set state contrast support and reference connectivity.", "table_03d gives case-state by referent-state pair counts and distinct matched sets in each cell.", 
                   "Referents are records, not independent people; matched sets can contribute to several comparison cells.", 
                   "Comparison support and connectivity do not establish absence of separation, confounding or imprecision.", 
                   "Burden decomposition is completed by this sensitivity script after risk estimation.", "Both stages write tables, figures and model objects into this analysis directory.", 
                   "A specific analysis can be selected with run_stage2_outcome(rr_file='tables/table_04_joint_9_state_rr.csv').", 
                   "Downstream scripts read current files from tables, models and plot_data."), file.path(DIR_MODEL, "README.txt"), 
                 useBytes = TRUE)
      ready <- list(protocol = "daynight_two_stage_v1", stage1_id = metadata$stage1_id, outcome_code = OUTCOME_CODE, 
                    outcome_name = OUTCOME_NAME, output_root = metadata$stage1_output_root, input_directory = metadata$stage1_input_directory, 
                    calibration_id = CURRENT_SHARED_DESIGN$calibration_id, analysis_tag = CURRENT_SHARED_DESIGN$analysis_tag, 
                    n_cases = sum(annual_counts$n_cases), completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
      save_obj(metadata, file.path(DIR_MODEL, "analysis_metadata_stage1.qs"))
      saveRDS(ready, file.path(DIR_MODEL, "stage1_completed.rds"))
      ready
    }
    
    run_stage1_internal_checks <- function() {
      states <- make_state_levels()
      info <- state_dt_from_levels(states)
      if (length(states) != 9L || uniqueN(states) != 9L || REF_STATE != states[1L] || uniqueN(info$mask) != 9L || 
          uniqueN(info$component_mask) != 9L) 
        stop("The nine-state definition failed the stage-one consistency check.")
      invisible(TRUE)
    }
    
  }, envir = e)
  e
}

make_burden_engine <- function() {
  e <- new.env(parent = environment())
  evalq({
    analysis_file <- function(root, name) {
      section <- if (grepl("^plotdata_", name)) 
        "plot_data"
      else if (grepl("^(input_manifest|output_manifest|source_file_manifest|output_file_manifest)[.]", name)) 
        "models"
      else if (grepl("[.](csv|xlsx)$", name, ignore.case = TRUE)) 
        "tables"
      else if (grepl("[.](png|pdf|svg|tif|tiff|jpg|jpeg)$", name, ignore.case = TRUE)) 
        "figures"
      else "models"
      file.path(root, section, name)
    }
    
    prepare_analysis_directory <- function(root) {
      dir.create(root, recursive = TRUE, showWarnings = FALSE)
      root <- normalizePath(root, winslash = "/", mustWork = TRUE)
      for (section in c("tables", "figures", "models", "plot_data"))
        dir.create(file.path(root, section), showWarnings = FALSE)
      invisible(root)
    }
    
    current_file <- function(root, name, required = TRUE) {
      target <- analysis_file(root, name)
      if (file.exists(target)) return(target)
      if (required) stop("Required input not found: ", target, call. = FALSE)
      NULL
    }
    
    RUN_COVARIANCE_SOURCE <- "auto"
    
    RUN_N_SIM <- NULL
    
    RUN_SIM_SEED <- NULL
    
    RUN_CAP_NEGATIVE_AF <- NULL
    
    REF_STATE <- "D0_N0"
    
    MODEL_LABEL <- "daynight_3x3_9state"
    
    N_SIM <- 1000L
    
    SIM_SEED <- 20260523L
    
    CAP_NEGATIVE_AF <- FALSE
    
    OBJECT_FORMAT <- "rds"
    
    SAVE_CASE_CONTRIB <- FALSE
    
    FONT_FAMILY <- "serif"
    
    FIG_DPI <- 600
    
    SAVE_PDF <- FALSE
    
    OUTCOME_CODE <- OUTCOME_NAME <- OUT_ROOT <- NULL
    
    set_output_dirs <- function(out_root, log_file_name, include_contrib = FALSE) {
      prepare_analysis_directory(out_root)
      DIR_TABLE <<- file.path(out_root, "tables")
      DIR_MODEL <<- file.path(out_root, "models")
      DIR_PLOT_DATA <<- file.path(out_root, "plot_data")
      DIR_FIG <<- file.path(out_root, "figures")
      DIR_LOG <<- file.path(out_root, "models")
      DIR_CONTRIB <<- file.path(out_root, "tables")
      dirs <- c(DIR_TABLE, DIR_MODEL, DIR_PLOT_DATA, DIR_FIG, DIR_LOG)
      if (isTRUE(include_contrib)) 
        dirs <- c(dirs, DIR_CONTRIB)
      for (dd in dirs) dir.create(dd, recursive = TRUE, showWarnings = FALSE)
      LOG_FILE <<- file.path(DIR_LOG, log_file_name)
    }
    
    log_msg <- function(...) {
      msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
      cat(msg, "\n")
      cat(msg, "\n", file = LOG_FILE, append = TRUE)
    }
    
    save_dt <- function(dt, file) {
      fwrite(dt, file)
      log_msg("Saved table:", file)
    }
    
    save_obj <- function(obj, file) {
      if (OBJECT_FORMAT == "rds") {
        file <- sub("[.]qs$", ".rds", file)
        saveRDS(obj, file)
      }
      else if (OBJECT_FORMAT == "qs") {
        if (!requireNamespace("qs", quietly = TRUE)) 
          stop("Saving .qs requires package qs.")
        qs::qsave(obj, file, preset = "fast")
      }
      else stop("OBJECT_FORMAT must be qs or rds.")
      log_msg("Saved object:", file)
    }
    
    coerce_binary01 <- function(x, variable_name) {
      if (is.factor(x)) 
        x <- as.character(x)
      if (is.logical(x)) 
        x <- as.integer(x)
      if (is.character(x)) {
        x_trim <- trimws(tolower(x))
        mapped <- rep(NA_integer_, length(x_trim))
        mapped[x_trim %in% c("0", "false", "no")] <- 0L
        mapped[x_trim %in% c("1", "true", "yes")] <- 1L
        x <- mapped
      }
      else {
        x <- suppressWarnings(as.numeric(x))
      }
      finite_values <- sort(unique(x[is.finite(x)]))
      if (!all(finite_values %in% c(0, 1))) {
        stop(variable_name, " must contain only binary values coded as 0 and 1. Observed finite values: ", paste(finite_values, 
                                                                                                                 collapse = ", "))
      }
      as.integer(x)
    }
    
    validate_unique_key <- function(dt, key_cols, object_name) {
      if (nrow(dt) == 0L) 
        stop(object_name, " is empty.")
      duplicate_rows <- dt[, .N, by = key_cols][N > 1L]
      if (nrow(duplicate_rows) > 0L) {
        stop(object_name, " contains duplicate rows for key: ", paste(key_cols, collapse = ", "))
      }
      invisible(TRUE)
    }
    
    validate_finite_named_vector <- function(x, expected_names, object_name) {
      if (is.null(names(x))) 
        stop(object_name, " must be a named vector.")
      missing_names <- setdiff(expected_names, names(x))
      if (length(missing_names) > 0L) {
        stop(object_name, " is missing required elements: ", paste(missing_names, collapse = ", "))
      }
      bad_names <- expected_names[!is.finite(x[expected_names])]
      if (length(bad_names) > 0L) {
        stop(object_name, " contains non-finite values for: ", paste(bad_names, collapse = ", "))
      }
      invisible(TRUE)
    }
    
    safe_empirical_quantile <- function(x, probability) {
      x <- x[is.finite(x)]
      if (length(x) == 0L) 
        return(NA_real_)
      as.numeric(quantile(x, probs = probability, na.rm = TRUE, type = 8))
    }
    
    theme_pub <- function(base_size = 10) {
      theme_classic(base_size = base_size, base_family = FONT_FAMILY) + theme(text = element_text(family = FONT_FAMILY, 
                                                                                                  colour = "black"), axis.title = element_text(size = base_size + 1, colour = "black"), axis.text = element_text(size = base_size, 
                                                                                                                                                                                                                 colour = "black"), axis.line = element_line(size = 0.35, colour = "black"), axis.ticks = element_line(size = 0.35, 
                                                                                                                                                                                                                                                                                                                       colour = "black"), strip.background = element_rect(fill = "grey95", colour = "grey70", size = 0.25), strip.text = element_text(size = base_size, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                      colour = "black"), legend.title = element_text(size = base_size, colour = "black"), legend.text = element_text(size = base_size - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       1, colour = "black"), panel.grid = element_blank(), plot.title = element_blank())
    }
    
    theme_nature <- function(base_size = 10) {
      theme_classic(base_size = base_size, base_family = FONT_FAMILY) + theme(text = element_text(family = FONT_FAMILY, 
                                                                                                  colour = "black"), axis.title = element_text(size = base_size + 1, colour = "black"), axis.text = element_text(size = base_size, 
                                                                                                                                                                                                                 colour = "black"), axis.line = element_line(size = 0.35, colour = "black"), axis.ticks = element_line(size = 0.35, 
                                                                                                                                                                                                                                                                                                                       colour = "black"), legend.title = element_text(size = base_size, colour = "black"), legend.text = element_text(size = base_size - 
                                                                                                                                                                                                                                                                                                                                                                                                                                        1, colour = "black"), panel.border = element_blank(), panel.grid = element_blank(), plot.title = element_blank())
    }
    
    save_plot_device <- function(p, file, width, height, pdf = FALSE) {
      destination <- normalizePath(dirname(file), winslash = "/", mustWork = TRUE)
      previous <- getwd()
      on.exit(setwd(previous), add = TRUE)
      setwd(destination)
      if (pdf) 
        grDevices::pdf(basename(file), width = width, height = height)
      else grDevices::png(basename(file), width = width, height = height, units = "in", res = FIG_DPI, bg = "white")
      on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
      print(p)
    }
    
    ggsave_png <- function(p, filename_base, width, height) {
      png_file <- file.path(DIR_FIG, paste0(filename_base, ".png"))
      save_plot_device(p, png_file, width, height)
      log_msg("Saved figure:", png_file)
    }
    
    ggsave_out <- function(p, filename_base, width, height) {
      png_file <- file.path(DIR_FIG, paste0(filename_base, ".png"))
      save_plot_device(p, png_file, width, height)
      log_msg("Saved figure:", png_file)
      if (isTRUE(SAVE_PDF)) {
        pdf_file <- file.path(DIR_FIG, paste0(filename_base, ".pdf"))
        save_plot_device(p, pdf_file, width, height, pdf = TRUE)
        log_msg("Saved figure:", pdf_file)
      }
    }
    
    rmvnorm_eigen <- function(n, mu, Sigma) {
      mu <- as.numeric(mu)
      Sigma <- as.matrix(Sigma)
      p <- length(mu)
      if (n < 1L || p < 1L) 
        stop("The simulation size and coefficient dimension must be positive.")
      if (!all(dim(Sigma) == c(p, p))) 
        stop("The covariance matrix has incompatible dimensions.")
      if (any(!is.finite(mu)) || any(!is.finite(Sigma))) {
        stop("Non-finite values were found in the simulation mean or covariance matrix.")
      }
      Sigma <- (Sigma + t(Sigma))/2
      eg <- eigen(Sigma, symmetric = TRUE)
      tolerance <- max(1, max(abs(eg$values))) * 1e-08
      if (min(eg$values) < -tolerance) {
        warning("The covariance matrix had negative eigenvalues; negative values were truncated to zero for simulation.")
      }
      vals <- pmax(eg$values, 0)
      A <- eg$vectors %*% diag(sqrt(vals), nrow = p)
      Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
      sweep(Z %*% t(A), 2, mu, "+")
    }
    
    make_state_levels <- function() {
      c("D0_N0", "DC_N0", "DH_N0", "D0_NC", "DC_NC", "DH_NC", "D0_NH", "DC_NH", "DH_NH")
    }
    
    state_dt_from_levels <- function(state_levels) {
      dt_state <- data.table(state = state_levels)
      dt_state[, `:=`(day_code, sub("_.*$", "", state))]
      dt_state[, `:=`(night_code, sub("^.*_", "", state))]
      dt_state[, `:=`(D_C, as.integer(day_code == "DC"))]
      dt_state[, `:=`(D_H, as.integer(day_code == "DH"))]
      dt_state[, `:=`(N_C, as.integer(night_code == "NC"))]
      dt_state[, `:=`(N_H, as.integer(night_code == "NH"))]
      dt_state[, `:=`(day_status, fifelse(day_code == "D0", "No daytime anomaly", fifelse(day_code == "DC", "Daytime cold", 
                                                                                          "Daytime heat")))]
      dt_state[, `:=`(night_status, fifelse(night_code == "N0", "No nighttime anomaly", fifelse(night_code == "NC", 
                                                                                                "Nighttime cold", "Nighttime heat")))]
      dt_state[, `:=`(day_axis, factor(day_status, levels = c("No daytime anomaly", "Daytime cold", "Daytime heat")))]
      dt_state[, `:=`(night_axis, factor(night_status, levels = c("No nighttime anomaly", "Nighttime cold", "Nighttime heat")))]
      dt_state[, `:=`(day_order, match(day_code, c("D0", "DC", "DH")) - 1L)]
      dt_state[, `:=`(night_order, match(night_code, c("N0", "NC", "NH")) - 1L)]
      dt_state[, `:=`(mask, day_order + 3L * night_order)]
      dt_state[, `:=`(component_mask, D_C * 1L + D_H * 2L + N_C * 4L + N_H * 8L)]
      dt_state[, `:=`(c("day_order", "night_order"), NULL)]
      dt_state[]
    }
    
    COMPONENTS <- c("D_C", "D_H", "N_C", "N_H")
    
    COMPONENT_BITS <- c(D_C = 1L, D_H = 2L, N_C = 4L, N_H = 8L)
    
    COMPONENT_INFO <- data.table(component = c(COMPONENTS, "joint_total"), component_order = 1:5, component_label = c("Daytime cold", 
                                                                                                                      "Daytime heat", "Nighttime cold", "Nighttime heat", "Joint total"), component_group = c("daytime_cold", "daytime_heat", 
                                                                                                                                                                                                              "nighttime_cold", "nighttime_heat", "joint_total"))
    
    state_from_mask <- function(mask) {
      day_levels <- c("D0", "DC", "DH")
      night_levels <- c("N0", "NC", "NH")
      day_i <- mask%%3L
      night_i <- mask%/%3L
      paste0(day_levels[day_i + 1L], "_", night_levels[night_i + 1L])
    }
    
    mask_from_daynight_codes <- function(day_code, night_code) {
      day_i <- match(day_code, c("D0", "DC", "DH")) - 1L
      night_i <- match(night_code, c("N0", "NC", "NH")) - 1L
      as.integer(day_i + 3L * night_i)
    }
    
    make_beta_state_vector <- function(coefficients, covariance, state_levels, ref_state) {
      cf <- coefficients
      vc <- covariance
      beta_state <- setNames(rep(NA_real_, length(state_levels)), state_levels)
      beta_state[ref_state] <- 0
      state_se <- setNames(rep(NA_real_, length(state_levels)), state_levels)
      state_se[ref_state] <- 0
      for (st in setdiff(state_levels, ref_state)) {
        term <- paste0("joint_state", st)
        if (term %in% names(cf)) {
          beta_state[st] <- unname(cf[term])
          state_se[st] <- sqrt(vc[term, term])
        }
      }
      list(beta_state = beta_state, state_se = state_se)
    }
    
    make_rr_table <- function(beta_state, state_se, state_info) {
      out <- copy(state_info)
      out[, `:=`(beta, beta_state[state])]
      out[, `:=`(se, state_se[state])]
      out[, `:=`(rr, fifelse(is.finite(beta), exp(beta), NA_real_))]
      out[, `:=`(rr_low, fifelse(is.finite(beta) & is.finite(se), exp(beta - 1.96 * se), NA_real_))]
      out[, `:=`(rr_high, fifelse(is.finite(beta) & is.finite(se), exp(beta + 1.96 * se), NA_real_))]
      out[state == REF_STATE, `:=`(beta = 0, se = 0, rr = 1, rr_low = 1, rr_high = 1)]
      out[, `:=`(estimable, is.finite(beta) & is.finite(se))]
      out[, `:=`(rr_label, fifelse(state == REF_STATE, "1.00\nReference", fifelse(estimable, sprintf("%.2f\n(%.2f, %.2f)", 
                                                                                                     rr, rr_low, rr_high), "NE\nNot observed")))]
      out[]
    }
    
    safe_af_from_beta <- function(beta) {
      af <- 1 - exp(-beta)
      if (CAP_NEGATIVE_AF) 
        af <- pmax(af, 0)
      af
    }
    
    compute_state_burden_tables <- function(beta_state, counts_by_year, n_deaths_by_year, state_info) {
      state_order <- state_info$state
      if (is.null(names(beta_state)) || !all(state_order %in% names(beta_state))) {
        stop("The joint-state coefficient vector does not contain all theoretical state names.")
      }
      if (!identical(names(counts_by_year), names(n_deaths_by_year))) {
        stop("Annual state counts and annual death totals must have identical year names and order.")
      }
      total_counts <- setNames(rep(0, length(state_order)), state_order)
      for (y in names(counts_by_year)) {
        counts_y <- counts_by_year[[y]]
        if (is.null(names(counts_y)) || !all(state_order %in% names(counts_y))) {
          stop("The state-count vector for year ", y, " does not contain all theoretical states.")
        }
        total_counts <- total_counts + as.numeric(counts_y[state_order])
      }
      non_estimable_states <- state_order[!is.finite(beta_state[state_order])]
      non_estimable_with_cases <- non_estimable_states[total_counts[non_estimable_states] > 0]
      if (length(non_estimable_with_cases) > 0L) {
        stop("At least one non-estimable joint state has observed cases and therefore cannot be assigned burden: ", 
             paste(non_estimable_with_cases, collapse = ", "))
      }
      af_state <- safe_af_from_beta(beta_state)
      af_state[REF_STATE] <- 0
      annual <- rbindlist(lapply(names(counts_by_year), function(y) {
        counts <- counts_by_year[[y]][state_order]
        n_deaths <- as.numeric(n_deaths_by_year[[y]])
        if (length(n_deaths) != 1L || !is.finite(n_deaths) || n_deaths <= 0) {
          stop("The total number of deaths for year ", y, " must be one positive finite value.")
        }
        if (any(!is.finite(counts)) || any(counts < 0)) {
          stop("The state counts for year ", y, " must be finite and non-negative.")
        }
        if (abs(sum(counts) - n_deaths) > 1e-08) {
          stop("State counts do not sum to the total number of deaths for year ", y, ".")
        }
        dt_y <- data.table(year = as.integer(y), state = state_order, n_cases = as.numeric(counts), n_deaths = n_deaths, 
                           state_af_value = af_state[state_order])
        dt_y[, `:=`(AN, fifelse(n_cases == 0, 0, n_cases * state_af_value))]
        if (any(!is.finite(dt_y$AN))) {
          bad <- dt_y[!is.finite(AN), state]
          stop("Non-finite attributable deaths were obtained for states: ", paste(bad, collapse = ", "))
        }
        joint_AN <- sum(dt_y$AN)
        dt_y[, `:=`(AF, AN/n_deaths)]
        dt_y[, `:=`(AF_percent, AF * 100)]
        if (!is.finite(joint_AN) || joint_AN == 0) {
          dt_y[, `:=`(share_percent, NA_real_)]
        }
        else {
          dt_y[, `:=`(share_percent, AN/joint_AN * 100)]
        }
        dt_y[]
      }))
      mean_annual <- annual[, .(n_cases = mean(n_cases, na.rm = TRUE), n_deaths = mean(n_deaths, na.rm = TRUE), AN = mean(AN, 
                                                                                                                          na.rm = TRUE), AF = mean(AF, na.rm = TRUE), AF_percent = mean(AF_percent, na.rm = TRUE)), by = state]
      mean_joint_AN <- mean_annual[, sum(AN)]
      if (!is.finite(mean_joint_AN) || mean_joint_AN == 0) {
        mean_annual[, `:=`(share_percent, NA_real_)]
      }
      else {
        mean_annual[, `:=`(share_percent, AN/mean_joint_AN * 100)]
      }
      mean_annual[, `:=`(period, "mean_annual")]
      annual <- merge(annual, state_info, by = "state", all.x = TRUE, sort = FALSE)
      mean_annual <- merge(mean_annual, state_info, by = "state", all.x = TRUE, sort = FALSE)
      list(annual = annual, mean_annual = mean_annual)
    }
    
    compute_state_shapley_values <- function(beta_state, state_info) {
      if (!all(state_info$state %in% names(beta_state))) {
        stop("The state coefficient vector does not contain all theoretical joint-state names.")
      }
      af_state <- safe_af_from_beta(beta_state)
      af_state[REF_STATE] <- 0
      mask_to_state <- setNames(state_info$state, as.character(state_info$component_mask))
      phi <- matrix(NA_real_, nrow = nrow(state_info), ncol = length(COMPONENTS), dimnames = list(state_info$state, 
                                                                                                  COMPONENTS))
      phi[REF_STATE, ] <- 0
      for (ii in seq_len(nrow(state_info))) {
        st <- state_info$state[ii]
        if (st == REF_STATE) 
          next
        if (!is.finite(af_state[st])) 
          next
        active <- COMPONENTS[as.integer(unlist(state_info[ii, ..COMPONENTS], use.names = FALSE)) == 1L]
        phi[st, ] <- 0
        if (length(active) == 0L) 
          next
        if (length(active) == 1L) {
          phi[st, active] <- af_state[st]
          next
        }
        if (length(active) != 2L) {
          stop("Each admissible joint exposure pattern must contain at most two active components.")
        }
        component_a <- active[1]
        component_b <- active[2]
        state_a <- unname(mask_to_state[as.character(COMPONENT_BITS[component_a])])
        state_b <- unname(mask_to_state[as.character(COMPONENT_BITS[component_b])])
        if (length(state_a) != 1L || length(state_b) != 1L || is.na(state_a) || is.na(state_b)) {
          stop("A required single-component state is unavailable for Shapley decomposition.")
        }
        value_a <- af_state[state_a]
        value_b <- af_state[state_b]
        value_ab <- af_state[st]
        if (any(!is.finite(c(value_a, value_b, value_ab)))) {
          stop("An observed compound state requires a non-estimable single-component counterfactual for Shapley decomposition: ", 
               st, ". Required single-component states: ", state_a, ", ", state_b, ".")
        }
        phi[st, component_a] <- 0.5 * (value_a + value_ab - value_b)
        phi[st, component_b] <- 0.5 * (value_b + value_ab - value_a)
      }
      finite_rows <- apply(phi, 1, function(z) all(is.finite(z))) & is.finite(af_state[rownames(phi)])
      if (any(finite_rows)) {
        additivity_error <- rowSums(phi[finite_rows, , drop = FALSE]) - af_state[rownames(phi)[finite_rows]]
        if (max(abs(additivity_error), na.rm = TRUE) > 1e-10) {
          stop("State-specific Shapley values failed the additivity check.")
        }
      }
      phi
    }
    
    compute_shapley_burden_tables <- function(beta_state, counts_by_year, n_deaths_by_year, state_info) {
      state_order <- state_info$state
      if (is.null(names(beta_state)) || !all(state_order %in% names(beta_state))) {
        stop("The joint-state coefficient vector does not contain all theoretical state names.")
      }
      if (!identical(names(counts_by_year), names(n_deaths_by_year))) {
        stop("Annual state counts and annual death totals must have identical year names and order.")
      }
      if (length(counts_by_year) == 0L) 
        stop("No annual death counts were available.")
      phi_state <- compute_state_shapley_values(beta_state, state_info)
      annual <- rbindlist(lapply(names(counts_by_year), function(y) {
        counts_named <- counts_by_year[[y]]
        if (is.null(names(counts_named)) || !all(state_order %in% names(counts_named))) {
          stop("The state-count vector for year ", y, " does not contain all theoretical states.")
        }
        counts <- as.numeric(counts_named[state_order])
        names(counts) <- state_order
        n_deaths <- as.numeric(n_deaths_by_year[[y]])
        if (length(n_deaths) != 1L || !is.finite(n_deaths) || n_deaths <= 0) {
          stop("The total number of deaths for year ", y, " must be one positive finite value.")
        }
        if (any(!is.finite(counts)) || any(counts < 0)) {
          stop("The state counts for year ", y, " must be finite and non-negative.")
        }
        if (abs(sum(counts) - n_deaths) > 1e-08) {
          stop("State counts do not sum to the total number of deaths for year ", y, ".")
        }
        positive_states <- state_order[counts > 0]
        if (length(positive_states) > 0L) {
          finite_phi <- apply(phi_state[positive_states, COMPONENTS, drop = FALSE], 1, function(z) all(is.finite(z)))
          if (any(!finite_phi)) {
            bad_states <- positive_states[!finite_phi]
            stop("Observed case states lack estimable Shapley values in year ", y, ": ", paste(bad_states, 
                                                                                               collapse = ", "))
          }
          component_an <- as.numeric(crossprod(counts[positive_states], phi_state[positive_states, COMPONENTS, 
                                                                                  drop = FALSE]))
        }
        else {
          component_an <- rep(0, length(COMPONENTS))
        }
        names(component_an) <- COMPONENTS
        joint_an <- sum(component_an)
        out_y <- data.table(year = as.integer(y), component = COMPONENTS, n_deaths = n_deaths, AN = component_an)
        out_y[, `:=`(AF, AN/n_deaths)]
        out_y[, `:=`(AF_percent, AF * 100)]
        if (!is.finite(joint_an) || joint_an == 0) {
          out_y[, `:=`(share_percent, NA_real_)]
        }
        else {
          out_y[, `:=`(share_percent, AN/joint_an * 100)]
        }
        joint_y <- data.table(year = as.integer(y), component = "joint_total", n_deaths = n_deaths, AN = joint_an, 
                              AF = joint_an/n_deaths, AF_percent = joint_an/n_deaths * 100, share_percent = 100)
        rbindlist(list(out_y, joint_y), use.names = TRUE, fill = TRUE)
      }))
      mean_annual <- annual[, .(n_deaths = mean(n_deaths, na.rm = TRUE), AN = mean(AN, na.rm = TRUE), AF = mean(AF, 
                                                                                                                na.rm = TRUE), AF_percent = mean(AF_percent, na.rm = TRUE)), by = component]
      mean_joint_an <- mean_annual[component == "joint_total", AN]
      if (length(mean_joint_an) != 1L || !is.finite(mean_joint_an) || mean_joint_an == 0) {
        mean_annual[, `:=`(share_percent, NA_real_)]
      }
      else {
        mean_annual[component != "joint_total", `:=`(share_percent, AN/mean_joint_an * 100)]
        mean_annual[component == "joint_total", `:=`(share_percent, 100)]
      }
      mean_annual[, `:=`(period, "mean_annual")]
      if (is.finite(mean_joint_an) && mean_joint_an != 0) {
        component_share_sum <- mean_annual[component %in% COMPONENTS, if (any(is.finite(share_percent))) 
          sum(share_percent, na.rm = TRUE)
          else NA_real_]
        if (!is.finite(component_share_sum) || abs(component_share_sum - 100) > 1e-08) {
          stop("Mean-annual Shapley component shares failed the 100% sum check.")
        }
      }
      annual <- merge(annual, COMPONENT_INFO, by = "component", all.x = TRUE, sort = FALSE)
      mean_annual <- merge(mean_annual, COMPONENT_INFO, by = "component", all.x = TRUE, sort = FALSE)
      setorder(annual, year, component_order)
      setorder(mean_annual, component_order)
      list(annual = annual, mean_annual = mean_annual, phi_by_state = phi_state)
    }
    
    add_ci <- function(obs, sims, by_cols, value_cols) {
      validate_unique_key(obs, by_cols, "Observed burden table")
      if (nrow(sims) == 0L) 
        stop("The simulation table is empty.")
      missing_obs <- setdiff(c(by_cols, value_cols), names(obs))
      missing_sims <- setdiff(c(by_cols, value_cols), names(sims))
      if (length(missing_obs) > 0L) {
        stop("Observed burden table is missing columns: ", paste(missing_obs, collapse = ", "))
      }
      if (length(missing_sims) > 0L) {
        stop("Simulation burden table is missing columns: ", paste(missing_sims, collapse = ", "))
      }
      ci <- sims[, c(setNames(lapply(.SD, safe_empirical_quantile, probability = 0.025), paste0(value_cols, "_low")), 
                     setNames(lapply(.SD, safe_empirical_quantile, probability = 0.975), paste0(value_cols, "_high"))), by = by_cols, 
                 .SDcols = value_cols]
      validate_unique_key(ci, by_cols, "Simulation confidence-interval table")
      merge(obs, ci, by = by_cols, all.x = TRUE, sort = FALSE)
    }
    
    run_internal_checks <- function() {
      states <- make_state_levels()
      if (length(states) != 9L || uniqueN(states) != 9L || REF_STATE != states[1]) {
        stop("The joint-state definition must contain 9 unique states with D0_N0 as the reference.")
      }
      state_info <- state_dt_from_levels(states)
      if (nrow(state_info) != 9L || uniqueN(state_info$mask) != 9L || uniqueN(state_info$component_mask) != 9L) {
        stop("The joint-state metadata failed the uniqueness check.")
      }
      beta_test <- setNames(c(0, 0.04, 0.06, 0.03, 0.08, 0.1, 0.05, 0.09, 0.12), states)
      counts_test <- setNames(c(100, 20, 15, 25, 10, 8, 18, 7, 5), states)
      deaths_test <- sum(counts_test)
      counts_by_year_test <- list(`2013` = counts_test, `2014` = counts_test + 1)
      n_deaths_by_year_test <- list(`2013` = deaths_test, `2014` = sum(counts_test + 1))
      shapley_test <- compute_shapley_burden_tables(beta_state = beta_test, counts_by_year = counts_by_year_test, 
                                                    n_deaths_by_year = n_deaths_by_year_test, state_info = state_info)
      if (nrow(shapley_test$annual) != 10L || nrow(shapley_test$mean_annual) != 5L) {
        stop("The Shapley burden table dimensions failed the internal check.")
      }
      share_sum <- shapley_test$mean_annual[component %in% COMPONENTS, sum(share_percent)]
      if (!is.finite(share_sum) || abs(share_sum - 100) > 1e-08) {
        stop("The Shapley burden shares failed the internal 100% sum check.")
      }
      phi_sum <- rowSums(shapley_test$phi_by_state)
      af_test <- safe_af_from_beta(beta_test)
      af_test[REF_STATE] <- 0
      if (max(abs(phi_sum - af_test[rownames(shapley_test$phi_by_state)])) > 1e-10) {
        stop("The state-specific Shapley values failed the internal additivity check.")
      }
      invisible(TRUE)
    }
    
    read_result_object <- function(path) {
      if (!file.exists(path)) 
        stop("Result object not found: ", path)
      ext <- tolower(tools::file_ext(path))
      if (ext == "rds") 
        return(readRDS(path))
      if (ext == "qs") {
        if (!requireNamespace("qs", quietly = TRUE)) 
          stop("Reading .qs requires package qs: ", path)
        return(qs::qread(path))
      }
      stop("Expected .qs or .rds object: ", path)
    }
    
    load_rr_input <- function(path) {
      rr <- if (is.data.frame(path)) 
        copy(as.data.table(path))
      else fread(path)
      need <- c("state", "rr", "rr_low", "rr_high", "n_cases", "n_rows", "n_referents", "n_sets", "support_status")
      missing <- setdiff(need, names(rr))
      if (length(missing)) 
        stop("RR table missing columns from the original table_04: ", paste(missing, collapse = ", "))
      states <- make_state_levels()
      if (nrow(rr) != length(states) || anyDuplicated(rr$state) || !setequal(rr$state, states)) 
        stop("table_04 must contain exactly the nine unique original state names.")
      rr <- rr[match(states, state)]
      if (!"beta" %in% names(rr)) 
        rr[, `:=`(beta, log(rr))]
      if (!"se" %in% names(rr)) 
        rr[, `:=`(se, (log(rr_high) - log(rr_low))/(2 * 1.96))]
      for (v in c("beta", "se", "rr", "rr_low", "rr_high", "n_cases", "n_rows", "n_referents", "n_sets")) if (!is.numeric(rr[[v]])) 
        stop("Non-numeric RR table column: ", v)
      for (v in c("n_cases", "n_rows", "n_referents", "n_sets")) if (any(!is.finite(rr[[v]])) || any(rr[[v]] < 0) || 
                                                                     any(rr[[v]] != floor(rr[[v]]))) 
        stop("Invalid counts in RR table column: ", v)
      if (any(rr$n_rows != rr$n_cases + rr$n_referents)) 
        stop("RR n_rows differs from cases + referents.")
      valid_support <- c("observed_estimable", "not_observed")
      if (any(!rr$support_status %in% valid_support)) 
        stop("Observed but non-estimable states cannot be resumed.")
      observed <- rr$n_rows > 0
      if (any((rr$support_status == "observed_estimable") != observed)) 
        stop("Inconsistent RR state support flags/counts.")
      r <- rr[observed]
      if (any(!is.finite(r$beta)) || any(!is.finite(r$se)) || any(r$se < 0) || any(!is.finite(r$rr)) || any(!is.finite(r$rr_low)) || 
          any(!is.finite(r$rr_high)) || any(r$rr <= 0 | r$rr_low <= 0 | r$rr_high <= 0)) 
        stop("Invalid observed-state RR/beta/SE/CI.")
      close <- function(x, y) all(abs(x - y) <= 1e-06 * (1 + abs(y)))
      if (!close(log(r$rr), r$beta) || !close(log(r$rr_low), r$beta - 1.96 * r$se) || !close(log(r$rr_high), r$beta + 
                                                                                             1.96 * r$se)) {
        stop("RR, beta, SE and CI columns are inconsistent in the supplied table. Supply one internally consistent set of estimates.")
      }
      ref <- rr[state == REF_STATE]
      if (ref$n_rows == 0 || abs(ref$beta) > 1e-12 || abs(ref$se) > 1e-12 || abs(ref$rr - 1) > 1e-12) 
        stop("D0_N0 must be an observed reference with beta=0, SE=0, RR=1.")
      if (any(is.finite(rr[!observed, beta])) || any(rr[!observed, n_cases] != 0)) 
        stop("Unobserved states must have zero cases and non-estimable beta.")
      info <- state_dt_from_levels(states)
      canonical <- make_rr_table(setNames(rr$beta, rr$state), setNames(rr$se, rr$state), info)
      csv_risk_values <- rr[, .(rr, rr_low, rr_high)]
      canonical[, `:=`(rr = csv_risk_values$rr, rr_low = csv_risk_values$rr_low, rr_high = csv_risk_values$rr_high)]
      canonical <- merge(canonical, rr[, .(state, n_rows, n_cases, n_referents, n_sets, support_status)], by = "state", 
                         all.x = TRUE, sort = FALSE)
      setorder(canonical, mask)
      list(table = canonical, original = rr, beta = setNames(rr$beta, rr$state), se = setNames(rr$se, rr$state), 
           estimable = setdiff(rr[observed, state], REF_STATE), unobserved = rr[!observed, state])
    }
    
    normalise_annual_counts <- function(x, risks) {
      x <- copy(as.data.table(x))
      cols <- c("year", "state", "n_cases", "n_deaths")
      if (!all(cols %in% names(x))) 
        stop("Annual counts require year, state, n_cases, n_deaths columns.")
      x <- x[, ..cols]
      validate_unique_key(x, c("year", "state"), "Annual state death counts")
      for (v in c("year", "n_cases", "n_deaths")) if (!is.numeric(x[[v]]) || any(!is.finite(x[[v]])) || any(x[[v]] != 
                                                                                                            floor(x[[v]]))) 
        stop("Invalid integer annual count/year column: ", v)
      if (any(x$n_cases < 0) || any(x$n_deaths <= 0)) 
        stop("Invalid annual deaths.")
      states <- make_state_levels()
      for (yr in unique(x$year)) {
        y <- x[year == yr]
        if (nrow(y) != 9L || !setequal(y$state, states) || uniqueN(y$n_deaths) != 1L || sum(y$n_cases) != y$n_deaths[1L]) 
          stop("Annual counts must contain all nine states and sum to n_deaths for year ", yr)
      }
      totals <- x[, .(total = sum(n_cases)), by = state]
      expected <- risks$original$n_cases[match(totals$state, risks$original$state)]
      if (any(totals$total != expected)) 
        stop("Annual state counts do not match n_cases in table_04. Sample/threshold/window/outcome mismatch.")
      x[, `:=`(state_order, match(state, states))]
      setorder(x, year, state_order)
      x[, `:=`(state_order, NULL)]
      x
    }
    
    read_stage1_package <- function(outcome_code = NULL, rr_file = NULL, pipeline_index_file = PIPELINE_INDEX_FILE, 
                                    skip_if_missing = FALSE) {
      index <- NULL
      if (is.null(rr_file)) {
        index <- load_directory_index(pipeline_index_file)
        entry <- index$outcomes[[outcome_code]]
        if (is.null(entry)) {
          if (skip_if_missing) {
            message("Skipping unconfigured outcome: ", outcome_code)
            return(NULL)
          }
          stop("No input directory is configured for: ", outcome_code)
        }
        root <- normalizePath(entry$output_root, winslash = "/", mustWork = TRUE)
      }
      else {
        root <- normalizePath(dirname(rr_file), winslash = "/", mustWork = TRUE)
        if (basename(root) == "tables") 
          root <- dirname(root)
      }
      required <- c("table_04_joint_9_state_rr.csv", "table_00c_annual_state_death_counts.csv", "source_parameters.rds", 
                    "table_02_joint_state_thresholds.csv")
      invisible(lapply(required, function(name) current_file(root, name)))
      for (name in c("joint_logrr_covariance.rds", "table_04a_joint_logrr_covariance.csv", "case_states.rds")) current_file(root, 
                                                                                                                            name, required = FALSE)
      params <- readRDS(current_file(root, "source_parameters.rds"))
      if (is.null(outcome_code)) 
        outcome_code <- params$outcome_code
      ready <- list(protocol = "daynight_two_stage_v1", outcome_code = outcome_code, outcome_name = if (is.null(params$outcome)) outcome_code else params$outcome, 
                    stage1_id = if (is.null(params$stage1_id)) outcome_code else params$stage1_id, calibration_id = params$calibration_id, 
                    analysis_tag = if (is.null(params$shared_design$analysis_tag)) basename(root) else params$shared_design$analysis_tag, 
                    input_directory = file.path(root, "tables"), output_root = root)
      list(root = root, folder = file.path(root, "tables"), rr_file = current_file(root, "table_04_joint_9_state_rr.csv"), 
           ready = ready, metadata = params, index = index, pipeline_index_file = pipeline_index_file)
    }
    
    read_working_covariance <- function(package, risks, source = RUN_COVARIANCE_SOURCE) {
      source <- match.arg(source, c("auto", "rescale", "csv", "rds"))
      states <- risks$estimable
      read_matrix <- function(file) {
        if (tolower(tools::file_ext(file)) == "csv") {
          x <- fread(file)
          labels <- intersect(c("state", "term"), names(x))
          if (length(labels) != 1L) 
            stop("Covariance CSV requires one state or term column.")
          row_names <- x[[labels]]
          x[, `:=`((labels), NULL)]
          if (!all(vapply(x, is.numeric, logical(1)))) 
            stop("Covariance entries must be numeric.")
          V <- as.matrix(x)
          rownames(V) <- row_names
        }
        else {
          x <- readRDS(file)
          if (is.matrix(x)) 
            V <- x
          else {
            V <- x$covariance_non_reference
            if (is.null(V)) 
              V <- x$covariance
            if (is.null(V)) 
              V <- x$covariance_state
          }
        }
        if (!is.matrix(V) || !is.numeric(V) || is.null(rownames(V)) || is.null(colnames(V))) 
          stop("The covariance must be a named numeric matrix.")
        rownames(V) <- sub("^joint_state", "", rownames(V))
        colnames(V) <- sub("^joint_state", "", colnames(V))
        if (anyDuplicated(rownames(V)) || anyDuplicated(colnames(V)) || !all(states %in% rownames(V)) || !all(states %in% 
                                                                                                              colnames(V))) 
          stop("The covariance does not cover the estimable states.")
        V <- V[states, states, drop = FALSE]
        if (any(!is.finite(V))) 
          stop("The covariance contains non-finite entries.")
        if (length(states)) {
          scale <- max(abs(V), .Machine$double.eps)
          if (max(abs(V - t(V))) > 1e-08 * scale) 
            stop("The covariance is not symmetric.")
          V <- (V + t(V))/2
          if (min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < -1e-08 * scale) 
            stop("The covariance is not positive semidefinite.")
        }
        V
      }
      rds_name <- "joint_logrr_covariance.rds"
      csv_name <- "table_04a_joint_logrr_covariance.csv"
      can_rescale <- source %in% c("auto", "rescale")
      if (source == "auto") {
        candidates <- c(csv = analysis_file(package$root, csv_name), rds = analysis_file(package$root, rds_name))
        candidates <- candidates[file.exists(candidates)]
        if (!length(candidates)) 
          stop("A full coefficient covariance matrix is required.")
        source <- names(candidates)[which.max(as.numeric(file.info(candidates)$mtime))]
      }
      V <- read_matrix(analysis_file(package$root, if (source == "csv") 
        csv_name
        else rds_name))
      se <- risks$se[states]
      if (can_rescale && length(states)) {
        old_se <- sqrt(diag(V))
        if (any(old_se == 0 & se != 0)) 
          stop("Cannot recover correlations for a zero-variance coefficient; provide a full covariance.")
        ratio <- ifelse(old_se == 0, 1, se/old_se)
        V <- V * outer(ratio, ratio)
      }
      if (length(states) && any(abs(diag(V) - se^2) > 1e-07 * pmax(1e-12, se^2))) 
        stop("The supplied full covariance diagonal disagrees with SE squared in the risk table.")
      list(matrix = V, model = NULL, source = source)
    }
    
    save_working_inputs <- function(package, risks, covariance) {
      root <- package$root
      states <- risks$estimable
      saveRDS(list(mu_non_reference = risks$beta[states], covariance_non_reference = covariance$matrix), analysis_file(root, 
                                                                                                                       "joint_logrr_covariance.rds"))
      fwrite(cbind(data.table(state = states), as.data.table(covariance$matrix)), analysis_file(root, "table_04a_joint_logrr_covariance.csv"))
      package$metadata$stage1_input_directory <- file.path(root, "tables")
      package$metadata$output_root <- root
      saveRDS(package$metadata, analysis_file(root, "source_parameters.rds"))
      package
    }
    
    run_stage2_outcome <- function(outcome_code = NULL, rr_file = NULL, pipeline_index_file = PIPELINE_INDEX_FILE, 
                                   skip_if_missing = FALSE, n_sim = RUN_N_SIM, sim_seed = RUN_SIM_SEED, cap_negative_af = RUN_CAP_NEGATIVE_AF, 
                                   covariance_source = RUN_COVARIANCE_SOURCE) {
      package <- read_stage1_package(outcome_code, rr_file, pipeline_index_file, skip_if_missing)
      if (is.null(package)) 
        return(invisible(NULL))
      metadata <- package$metadata
      ready <- package$ready
      bundle <- list(directory = package$folder, rr_file = package$rr_file, covariance_file = analysis_file(package$root, 
                                                                                                            "joint_logrr_covariance.rds"), annual_counts_file = file.path(package$folder, "table_00c_annual_state_death_counts.csv"), 
                     metadata_file = analysis_file(package$root, "source_parameters.rds"), threshold_file = file.path(package$folder, 
                                                                                                                      "table_02_joint_state_thresholds.csv"))
      risks <- load_rr_input(bundle$rr_file)
      covariance <- read_working_covariance(package, risks, covariance_source)
      annual <- normalise_annual_counts(fread(bundle$annual_counts_file), risks)
      n_sim <- if (is.null(n_sim)) 
        metadata$n_sim
      else n_sim
      sim_seed <- if (is.null(sim_seed)) 
        metadata$sim_seed
      else sim_seed
      cap_negative_af <- if (is.null(cap_negative_af)) 
        metadata$cap_negative_af
      else cap_negative_af
      if (!is.numeric(n_sim) || length(n_sim) != 1L || !is.finite(n_sim) || n_sim < 2 || n_sim != floor(n_sim)) 
        stop("n_sim must be an integer >= 2.")
      if (!is.numeric(sim_seed) || length(sim_seed) != 1L || !is.finite(sim_seed) || sim_seed < 0 || sim_seed > .Machine$integer.max || 
          sim_seed != floor(sim_seed)) 
        stop("Invalid Monte Carlo seed.")
      if (!is.logical(cap_negative_af) || length(cap_negative_af) != 1L || is.na(cap_negative_af)) 
        stop("Invalid negative-AF setting.")
      if (!metadata$object_format %in% c("qs", "rds")) 
        stop("Invalid stage-one object format.")
      if (metadata$object_format == "qs" && !requireNamespace("qs", quietly = TRUE)) 
        stop("Install qs to retain stage one's object output format.")
      package <- save_working_inputs(package, risks, covariance)
      metadata <- package$metadata
      ready <- package$ready
      OUTCOME_CODE <<- ready$outcome_code
      OUTCOME_NAME <<- ready$outcome_name
      OUT_ROOT <<- package$root
      N_SIM <<- as.integer(n_sim)
      SIM_SEED <<- as.integer(sim_seed)
      CAP_NEGATIVE_AF <<- cap_negative_af
      OBJECT_FORMAT <<- metadata$object_format
      SAVE_CASE_CONTRIB <<- isTRUE(metadata$save_case_contrib)
      FONT_FAMILY <<- metadata$plot_settings$font_family
      FIG_DPI <<- metadata$plot_settings$fig_dpi
      SAVE_PDF <<- isTRUE(metadata$plot_settings$save_pdf)
      set_output_dirs(OUT_ROOT, "stage2_burden_analysis_log.txt", SAVE_CASE_CONTRIB)
      log_msg("Stage two started; input bundle:", bundle$directory)
      log_msg("Shared output directory:", OUT_ROOT)
      owned_markers <- file.path(DIR_MODEL, c("stage2_completed.rds", "run_completed.qs", "run_completed.rds"))
      for (p in owned_markers[file.exists(owned_markers)]) if (!file.remove(p)) 
        stop("Cannot invalidate stage-two completion marker: ", p)
      run_internal_checks()
      threshold_summary <- fread(bundle$threshold_file)
      STATE_LEVELS <- make_state_levels()
      STATE_INFO <- state_dt_from_levels(STATE_LEVELS)
      rr_table <- risks$table
      beta_state <- risks$beta
      state_se <- risks$se
      estimable_non_ref_states <- risks$estimable
      unobserved_states <- risks$unobserved
      observed_state_levels <- setdiff(STATE_LEVELS, unobserved_states)
      coef_joint <- setNames(beta_state[estimable_non_ref_states], paste0("joint_state", estimable_non_ref_states))
      vcov_joint <- covariance$matrix
      dimnames(vcov_joint) <- list(names(coef_joint), names(coef_joint))
      compact_file <- file.path(DIR_MODEL, paste0("model_joint_", MODEL_LABEL, ".", OBJECT_FORMAT))
      compact <- if (file.exists(compact_file)) 
        read_result_object(compact_file)
      else list()
      compact$coefficients <- coef_joint
      compact$covariance <- vcov_joint
      compact$beta_state <- beta_state
      compact$mu_non_reference <- risks$beta[estimable_non_ref_states]
      compact$covariance_non_reference <- covariance$matrix
      compact$covariance_state <- matrix(NA_real_, length(STATE_LEVELS), length(STATE_LEVELS), dimnames = list(STATE_LEVELS, 
                                                                                                               STATE_LEVELS))
      compact$covariance_state[REF_STATE, REF_STATE] <- 0
      compact$covariance_state[estimable_non_ref_states, estimable_non_ref_states] <- covariance$matrix
      compact$storage_mode <- "compact_state_coefficients_and_covariance"
      save_obj(compact, file.path(DIR_MODEL, paste0("model_joint_", MODEL_LABEL, ".qs")))
      save_obj(risks$original, file.path(DIR_MODEL, "joint_9_state_rr.qs"))
      save_obj(metadata, file.path(DIR_MODEL, "analysis_metadata_stage1.qs"))
      years <- sort(unique(annual$year))
      counts_by_year <- setNames(lapply(years, function(yr) {
        y <- annual[year == yr]
        setNames(y$n_cases, y$state)[STATE_LEVELS]
      }), as.character(years))
      n_deaths_by_year <- setNames(lapply(years, function(yr) unique(annual[year == yr, n_deaths])), as.character(years))
      case_dt <- NULL
      if (SAVE_CASE_CONTRIB) {
        if (!file.exists(analysis_file(package$root, "case_states.rds"))) 
          stop("Stage one did not save required case-level input.")
        case_dt <- as.data.table(readRDS(analysis_file(package$root, "case_states.rds")))
        if (!all(c("id", "date", "year", "joint_state", "joint_mask") %in% names(case_dt)) || nrow(case_dt) != 
            sum(annual$n_cases) || anyNA(case_dt$year) || any(!case_dt$joint_state %in% STATE_LEVELS) || anyDuplicated(case_dt, 
                                                                                                                       by = c("id", "date"))) 
          stop("Invalid saved case-state input.")
        actual <- case_dt[, .(n_cases = .N), by = .(year, state = as.character(joint_state))]
        compare <- merge(annual[, .(year, state, expected = n_cases)], actual, by = c("year", "state"), all.x = TRUE)
        compare[is.na(n_cases), `:=`(n_cases, 0L)]
        if (any(compare$expected != compare$n_cases)) 
          stop("Case-level and annual counts disagree.")
        case_dt <- case_dt[, .(id, date, year, joint_state, joint_mask)]
      }
      plot_rr <- copy(rr_table)
      finite_rr <- plot_rr[is.finite(rr), rr]
      fill_limits <- range(finite_rr, na.rm = TRUE)
      if (!is.finite(fill_limits[1]) || !is.finite(fill_limits[2]) || fill_limits[1] == fill_limits[2]) {
        fill_limits <- c(0.9, 1.1)
      }
      p_heatmap <- ggplot(plot_rr, aes(x = day_axis, y = night_axis, fill = rr)) + geom_tile(colour = "white", size = 0.35) + 
        geom_text(aes(label = rr_label), size = 2.6, family = FONT_FAMILY, lineheight = 0.9) + scale_fill_gradientn(colours = c("#2166AC", 
                                                                                                                                "#F7F7F7", "#B2182B"), limits = fill_limits, na.value = "white", name = "RR") + labs(x = "Daytime state", 
                                                                                                                                                                                                                     y = "Nighttime state") + coord_equal() + theme_nature(base_size = 9) + theme(axis.text.x = element_text(angle = 35, 
                                                                                                                                                                                                                                                                                                                             hjust = 1, vjust = 1), legend.position = "right")
      save_dt(plot_rr, file.path(DIR_PLOT_DATA, "plotdata_joint_3x3_rr_heatmap.csv"))
      save_obj(plot_rr, file.path(DIR_PLOT_DATA, "plotdata_joint_3x3_rr_heatmap.qs"))
      ggsave_out(p_heatmap, "fig_01_joint_3x3_rr_heatmap", width = 5.6, height = 4.8)
      log_msg("Computing observed exact nine-state attributable burden...")
      state_obs <- compute_state_burden_tables(beta_state, counts_by_year, n_deaths_by_year, STATE_INFO)
      rr_cols <- rr_table[, .(state, rr, rr_low, rr_high, support_status)]
      state_obs$annual <- merge(state_obs$annual, rr_cols, by = "state", all.x = TRUE, sort = FALSE)
      state_obs$mean_annual <- merge(state_obs$mean_annual, rr_cols, by = "state", all.x = TRUE, sort = FALSE)
      log_msg("Computing observed four-component Shapley burden from the nine-state model...")
      shapley_obs <- compute_shapley_burden_tables(beta_state, counts_by_year, n_deaths_by_year, STATE_INFO)
      phi_observed <- as.data.table(shapley_obs$phi_by_state, keep.rownames = "state")
      missing_phi_columns <- setdiff(COMPONENTS, names(phi_observed))
      if (length(missing_phi_columns) > 0L) {
        stop("The state-specific Shapley matrix is missing component columns: ", paste(missing_phi_columns, collapse = ", "))
      }
      phi_columns <- paste0("phi_", COMPONENTS)
      setnames(phi_observed, COMPONENTS, phi_columns)
      phi_observed <- merge(phi_observed, STATE_INFO, by = "state", all.x = TRUE, sort = FALSE)
      phi_observed[, `:=`(state_af_value, safe_af_from_beta(beta_state[state]))]
      phi_observed[state == REF_STATE, `:=`(state_af_value, 0)]
      phi_observed[, `:=`(shapley_sum, rowSums(.SD)), .SDcols = phi_columns]
      phi_observed[, `:=`(additivity_error, shapley_sum - state_af_value)]
      setorder(phi_observed, mask)
      save_dt(phi_observed, file.path(DIR_TABLE, "table_07_state_specific_shapley_values.csv"))
      log_msg("Starting Monte Carlo uncertainty estimation for state and Shapley burden. N_SIM:", N_SIM)
      set.seed(SIM_SEED)
      non_ref_states <- setdiff(STATE_LEVELS, REF_STATE)
      term_names_all <- paste0("joint_state", non_ref_states)
      names(term_names_all) <- non_ref_states
      estimable_states <- estimable_non_ref_states
      estimable_terms <- term_names_all[estimable_states]
      if (length(unobserved_states) > 0L) {
        log_msg("Monte Carlo simulation excludes structurally unobserved coefficient(s):", paste(unobserved_states, 
                                                                                                 collapse = ", "))
      }
      mu <- coef_joint[estimable_terms]
      Sigma <- vcov_joint[estimable_terms, estimable_terms, drop = FALSE]
      if (length(mu) > 0 && (any(!is.finite(mu)) || any(!is.finite(Sigma)))) {
        stop("Non-finite coefficient or covariance values were found in the joint model.")
      }
      if (length(mu) > 0) {
        draw_terms <- rmvnorm_eigen(N_SIM, mu, Sigma)
        colnames(draw_terms) <- estimable_states
      }
      else {
        draw_terms <- matrix(nrow = N_SIM, ncol = 0)
      }
      sim_state_annual_list <- vector("list", N_SIM)
      sim_state_mean_list <- vector("list", N_SIM)
      sim_shapley_annual_list <- vector("list", N_SIM)
      sim_shapley_mean_list <- vector("list", N_SIM)
      for (ss in seq_len(N_SIM)) {
        if (ss%%100 == 0) 
          log_msg("Monte Carlo simulation:", ss, "of", N_SIM)
        beta_draw <- setNames(rep(NA_real_, length(STATE_LEVELS)), STATE_LEVELS)
        beta_draw[REF_STATE] <- 0
        beta_draw[estimable_states] <- draw_terms[ss, estimable_states]
        state_sim <- compute_state_burden_tables(beta_draw, counts_by_year, n_deaths_by_year, STATE_INFO)
        shapley_sim <- compute_shapley_burden_tables(beta_draw, counts_by_year, n_deaths_by_year, STATE_INFO)
        sim_state_annual_list[[ss]] <- state_sim$annual[, .(sim = ss, year, state, AN, AF_percent, share_percent)]
        sim_state_mean_list[[ss]] <- state_sim$mean_annual[, .(sim = ss, period, state, AN, AF_percent, share_percent)]
        sim_shapley_annual_list[[ss]] <- shapley_sim$annual[, .(sim = ss, year, component, AN, AF_percent, share_percent)]
        sim_shapley_mean_list[[ss]] <- shapley_sim$mean_annual[, .(sim = ss, period, component, AN, AF_percent, 
                                                                   share_percent)]
      }
      sim_state_annual <- rbindlist(sim_state_annual_list, use.names = TRUE, fill = TRUE)
      sim_state_mean <- rbindlist(sim_state_mean_list, use.names = TRUE, fill = TRUE)
      sim_shapley_annual <- rbindlist(sim_shapley_annual_list, use.names = TRUE, fill = TRUE)
      sim_shapley_mean <- rbindlist(sim_shapley_mean_list, use.names = TRUE, fill = TRUE)
      rm(sim_state_annual_list, sim_state_mean_list, sim_shapley_annual_list, sim_shapley_mean_list, draw_terms)
      gc()
      state_annual_ci <- add_ci(obs = state_obs$annual, sims = sim_state_annual, by_cols = c("year", "state"), value_cols = c("AN", 
                                                                                                                              "AF_percent", "share_percent"))
      state_mean_ci <- add_ci(obs = state_obs$mean_annual, sims = sim_state_mean, by_cols = c("period", "state"), 
                              value_cols = c("AN", "AF_percent", "share_percent"))
      shapley_annual_ci <- add_ci(obs = shapley_obs$annual, sims = sim_shapley_annual, by_cols = c("year", "component"), 
                                  value_cols = c("AN", "AF_percent", "share_percent"))
      shapley_mean_ci <- add_ci(obs = shapley_obs$mean_annual, sims = sim_shapley_mean, by_cols = c("period", "component"), 
                                value_cols = c("AN", "AF_percent", "share_percent"))
      setorder(state_annual_ci, year, mask)
      setorder(state_mean_ci, mask)
      setorder(shapley_annual_ci, year, component_order)
      setorder(shapley_mean_ci, component_order)
      for (object_name in c("state_annual_ci", "state_mean_ci", "shapley_annual_ci", "shapley_mean_ci")) {
        object_value <- get(object_name)
        object_value[, `:=`(AN_CI = sprintf("%.3f (%.3f, %.3f)", AN, AN_low, AN_high), AF_percent_CI = sprintf("%.3f (%.3f, %.3f)", 
                                                                                                               AF_percent, AF_percent_low, AF_percent_high), share_percent_CI = sprintf("%.3f (%.3f, %.3f)", share_percent, 
                                                                                                                                                                                        share_percent_low, share_percent_high))]
        assign(object_name, object_value)
      }
      state_share_check_annual <- state_annual_ci[, .(state_share_sum = if (any(is.finite(share_percent))) 
        sum(share_percent, na.rm = TRUE)
        else NA_real_), by = year]
      state_share_check_mean <- state_mean_ci[, .(state_share_sum = if (any(is.finite(share_percent))) 
        sum(share_percent, na.rm = TRUE)
        else NA_real_), by = period]
      shapley_share_check_annual <- shapley_annual_ci[component %in% COMPONENTS, .(component_share_sum = if (any(is.finite(share_percent))) 
        sum(share_percent, na.rm = TRUE)
        else NA_real_), by = year]
      shapley_share_check_mean <- shapley_mean_ci[component %in% COMPONENTS, .(component_share_sum = if (any(is.finite(share_percent))) 
        sum(share_percent, na.rm = TRUE)
        else NA_real_), by = period]
      shapley_additivity_annual <- merge(state_annual_ci[, .(state_joint_AN = sum(AN, na.rm = TRUE)), by = year], 
                                         shapley_annual_ci[component == "joint_total", .(year, shapley_joint_AN = AN)], by = "year", all = TRUE)
      shapley_additivity_annual[, `:=`(difference, shapley_joint_AN - state_joint_AN)]
      shapley_additivity_mean <- data.table(period = "mean_annual", state_joint_AN = state_mean_ci[, sum(AN, na.rm = TRUE)], 
                                            shapley_joint_AN = shapley_mean_ci[component == "joint_total", AN])
      shapley_additivity_mean[, `:=`(difference, shapley_joint_AN - state_joint_AN)]
      if (any(abs(state_share_check_annual$state_share_sum - 100) > 1e-06, na.rm = TRUE) || any(abs(state_share_check_mean$state_share_sum - 
                                                                                                    100) > 1e-06, na.rm = TRUE)) {
        stop("The nine-state burden shares failed the 100% sum check.")
      }
      if (any(abs(shapley_share_check_annual$component_share_sum - 100) > 1e-06, na.rm = TRUE) || any(abs(shapley_share_check_mean$component_share_sum - 
                                                                                                          100) > 1e-06, na.rm = TRUE)) {
        stop("The four-component Shapley burden shares failed the 100% sum check.")
      }
      if (any(abs(shapley_additivity_annual$difference) > 1e-08, na.rm = TRUE) || any(abs(shapley_additivity_mean$difference) > 
                                                                                      1e-08, na.rm = TRUE)) {
        stop("The Shapley burden totals failed the additivity check against the nine-state burden.")
      }
      save_dt(state_share_check_annual, file.path(DIR_TABLE, "table_05_state_share_sum_check_annual.csv"))
      save_dt(state_share_check_mean, file.path(DIR_TABLE, "table_06_state_share_sum_check_mean_annual.csv"))
      save_dt(shapley_share_check_annual, file.path(DIR_TABLE, "table_08_shapley_share_sum_check_annual.csv"))
      save_dt(shapley_share_check_mean, file.path(DIR_TABLE, "table_08_shapley_share_sum_check_mean_annual.csv"))
      save_dt(shapley_additivity_annual, file.path(DIR_TABLE, "table_08_shapley_additivity_check_annual.csv"))
      save_dt(shapley_additivity_mean, file.path(DIR_TABLE, "table_08_shapley_additivity_check_mean_annual.csv"))
      state_annual_out <- state_annual_ci[, .(year, state, mask, component_mask, D_C, D_H, N_C, N_H, day_code, night_code, 
                                              day_status, night_status, support_status, n_cases, n_deaths, rr, rr_low, rr_high, AN, AN_low, AN_high, 
                                              AN_CI, AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI, share_percent, share_percent_low, share_percent_high, 
                                              share_percent_CI)]
      state_mean_out <- state_mean_ci[, .(period, state, mask, component_mask, D_C, D_H, N_C, N_H, day_code, night_code, 
                                          day_status, night_status, support_status, n_cases, n_deaths, rr, rr_low, rr_high, AN, AN_low, AN_high, 
                                          AN_CI, AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI, share_percent, share_percent_low, share_percent_high, 
                                          share_percent_CI)]
      shapley_annual_out <- shapley_annual_ci[, .(year, component, component_order, component_label, component_group, 
                                                  n_deaths, AN, AN_low, AN_high, AN_CI, AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI, share_percent, 
                                                  share_percent_low, share_percent_high, share_percent_CI)]
      shapley_mean_out <- shapley_mean_ci[, .(period, component, component_order, component_label, component_group, 
                                              n_deaths, AN, AN_low, AN_high, AN_CI, AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI, share_percent, 
                                              share_percent_low, share_percent_high, share_percent_CI)]
      save_dt(state_annual_out, file.path(DIR_TABLE, "table_09_exact_9_state_burden_annual.csv"))
      save_dt(state_mean_out, file.path(DIR_TABLE, "table_10_exact_9_state_burden_mean_annual.csv"))
      save_dt(shapley_annual_out, file.path(DIR_TABLE, "table_11_shapley_4_component_burden_annual.csv"))
      save_dt(shapley_mean_out, file.path(DIR_TABLE, "table_12_shapley_4_component_burden_mean_annual.csv"))
      save_obj(state_annual_out, file.path(DIR_MODEL, "exact_9_state_burden_annual.qs"))
      save_obj(state_mean_out, file.path(DIR_MODEL, "exact_9_state_burden_mean_annual.qs"))
      save_obj(shapley_annual_out, file.path(DIR_MODEL, "shapley_4_component_burden_annual.qs"))
      save_obj(shapley_mean_out, file.path(DIR_MODEL, "shapley_4_component_burden_mean_annual.qs"))
      state_summary <- state_mean_out[, .(result_type = "exact_9_state", period, item = state, item_label = paste(day_status, 
                                                                                                                  "|", night_status), AN, AN_low, AN_high, AN_CI, AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI, 
                                          share_percent, share_percent_low, share_percent_high, share_percent_CI)]
      shapley_summary <- shapley_mean_out[, .(result_type = "shapley_4_component", period, item = component, item_label = component_label, 
                                              AN, AN_low, AN_high, AN_CI, AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI, share_percent, 
                                              share_percent_low, share_percent_high, share_percent_CI)]
      mean_summary_requested <- rbindlist(list(state_summary, shapley_summary), use.names = TRUE, fill = TRUE)
      save_dt(mean_summary_requested, file.path(DIR_TABLE, "table_13_requested_mean_annual_summary.csv"))
      save_obj(mean_summary_requested, file.path(DIR_MODEL, "requested_mean_annual_summary.qs"))
      shapley_plot_data <- copy(shapley_mean_out[component %in% COMPONENTS])
      shapley_plot_data[, `:=`(component_label, factor(component_label, levels = COMPONENT_INFO[component %in% COMPONENTS, 
                                                                                                component_label]))]
      save_dt(shapley_plot_data, file.path(DIR_PLOT_DATA, "plotdata_shapley_4_component_burden_share.csv"))
      save_obj(shapley_plot_data, file.path(DIR_PLOT_DATA, "plotdata_shapley_4_component_burden_share.qs"))
      component_colours <- c(`Daytime cold` = "#2C7FB8", `Daytime heat` = "#F4A582", `Nighttime cold` = "#7FCDBB", 
                             `Nighttime heat` = "#C51B1D")
      if (all(is.finite(shapley_plot_data$share_percent)) && all(shapley_plot_data$share_percent >= 0) && sum(shapley_plot_data$share_percent) > 
          0) {
        p_shapley <- ggplot(shapley_plot_data, aes(x = "", y = share_percent, fill = component_label)) + geom_col(width = 1, 
                                                                                                                  colour = "white", size = 0.35) + coord_polar(theta = "y") + geom_text(aes(label = sprintf("%.1f%%", 
                                                                                                                                                                                                            share_percent)), position = position_stack(vjust = 0.5), size = 3, family = FONT_FAMILY) + scale_fill_manual(values = component_colours, 
                                                                                                                                                                                                                                                                                                                         name = NULL) + labs(x = NULL, y = NULL) + theme_void(base_family = FONT_FAMILY) + theme(legend.position = "top", 
                                                                                                                                                                                                                                                                                                                                                                                                                 legend.text = element_text(size = 9, colour = "black"), plot.margin = margin(5, 5, 5, 5))
        ggsave_out(p_shapley, "fig_02_shapley_4_component_burden_share", width = 5.2, height = 4.5)
      }
      else {
        for (extension in c("png", "pdf")) {
          old_figure <- file.path(DIR_FIG, paste0("fig_02_shapley_4_component_burden_share.", extension))
          if (file.exists(old_figure) && !file.remove(old_figure)) 
            stop("Cannot remove an obsolete share figure: ", old_figure)
        }
        log_msg("The Shapley share figure was not generated because one or more component shares were negative or non-finite.")
      }
      rm(sim_state_annual, sim_state_mean, sim_shapley_annual, sim_shapley_mean)
      gc()
      if (SAVE_CASE_CONTRIB) {
        log_msg("Saving individual case-level state AF contributions...")
        af_state_obs <- safe_af_from_beta(beta_state)
        af_state_obs[REF_STATE] <- 0
        contrib <- copy(case_dt)
        contrib[, `:=`(state_af, af_state_obs[joint_state])]
        save_obj(contrib, file.path(DIR_MODEL, "case_level_state_af_contributions.qs"))
        save_dt(contrib, file.path(DIR_CONTRIB, "case_level_state_af_contributions.csv"))
      }
      analysis_meta <- metadata
      analysis_meta$output_root <- OUT_ROOT
      analysis_meta$n_sim <- N_SIM
      analysis_meta$sim_seed <- SIM_SEED
      analysis_meta$cap_negative_af <- CAP_NEGATIVE_AF
      analysis_meta$shapley_decomposition <- list(method = "State-specific Shapley allocation derived from the nine-state model", 
                                                  components = COMPONENTS, component_labels = COMPONENT_INFO, phi_by_state = phi_observed)
      analysis_meta$share_checks <- list(state_annual = state_share_check_annual, state_mean = state_share_check_mean, 
                                         shapley_annual = shapley_share_check_annual, shapley_mean = shapley_share_check_mean, shapley_additivity_annual = shapley_additivity_annual, 
                                         shapley_additivity_mean = shapley_additivity_mean)
      analysis_meta$stage2 <- list(stage1_id = ready$stage1_id, input_directory = bundle$directory, no_model_refit = TRUE, 
                                   completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
      save_obj(analysis_meta, file.path(DIR_MODEL, "analysis_metadata.qs"))
      writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, "sessionInfo_stage2.txt"))
      complete <- list(outcome = OUTCOME_CODE, stage1_id = ready$stage1_id, completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
      saveRDS(complete, file.path(DIR_MODEL, "stage2_completed.rds"))
      save_obj(complete, file.path(DIR_MODEL, "run_completed.qs"))
      log_msg("Both stages now form the complete analysis output:", OUT_ROOT)
      invisible(list(output_root = OUT_ROOT, rr_table = rr_table, annual_counts = annual, state_annual_out = state_annual_out, 
                     state_mean_out = state_mean_out, shapley_annual_out = shapley_annual_out, shapley_mean_out = shapley_mean_out, 
                     phi_observed = phi_observed, analysis_meta = analysis_meta))
    }
    
  }, envir = e)
  e
}

make_heterogeneity_engine <- function() {
  e <- new.env(parent = environment())
  evalq({
    analysis_file <- function(root, name) {
      section <- if (grepl("^plotdata_", name)) 
        "plot_data"
      else if (grepl("^(input_manifest|output_manifest|source_file_manifest|output_file_manifest)[.]", name)) 
        "models"
      else if (grepl("[.](csv|xlsx)$", name, ignore.case = TRUE)) 
        "tables"
      else if (grepl("[.](png|pdf|svg|tif|tiff|jpg|jpeg)$", name, ignore.case = TRUE)) 
        "figures"
      else "models"
      file.path(root, section, name)
    }
    
    prepare_analysis_directory <- function(root) {
      dir.create(root, recursive = TRUE, showWarnings = FALSE)
      root <- normalizePath(root, winslash = "/", mustWork = TRUE)
      for (section in c("tables", "figures", "models", "plot_data"))
        dir.create(file.path(root, section), showWarnings = FALSE)
      invisible(root)
    }
    
    current_file <- function(root, name, required = TRUE) {
      target <- analysis_file(root, name)
      if (file.exists(target)) return(target)
      if (required) stop("Required input not found: ", target, call. = FALSE)
      NULL
    }
    
    N_SIM <- 10000L
    
    SIM_SEED <- 20260925L
    
    ALPHA <- 0.05
    
    WRITE_XLSX <- TRUE
    
    SAVE_SIMULATION_DRAWS <- FALSE
    
    EIGEN_REL_TOL <- 1e-10
    
    RECONSTRUCTION_REL_TOL <- 1e-07
    
    TOTAL_AN_REL_EPS <- 1e-08
    
    MAX_TOTAL_SIGN_SWITCH_FRACTION <- 0.01
    
    MAX_PROFILE_SPLIT_COV_REL_DIFF <- 0.5
    
    OUTCOME_CONFIG <- data.frame(outcome_code = c("I00_I52_I60_I69", "I10_I15", "I20_I25", "I50", "I60_I62", "I63"), 
                                 outcome = c("Overall cardiocerebrovascular mortality", "Hypertensive diseases", "Ischemic heart disease", 
                                             "Heart failure", "Hemorrhagic stroke", "Ischemic stroke"), between_disease = c(FALSE, TRUE, TRUE, TRUE, 
                                                                                                                              TRUE, TRUE), stringsAsFactors = FALSE)
    
    RUN_OUTCOMES <- OUTCOME_CONFIG$outcome_code
    
    REF_STATE <- "D0_N0"
    
    STATE_LEVELS <- c("D0_N0", "DC_N0", "DH_N0", "D0_NC", "DC_NC", "DH_NC", "D0_NH", "DC_NH", "DH_NH")
    
    COMPONENTS <- c("D_C", "N_C", "D_H", "N_H")
    
    COMPONENT_LABELS <- c(D_C = "Daytime cold", N_C = "Nighttime cold", D_H = "Daytime heat", N_H = "Nighttime heat")
    
    PROFILE_COORDINATES <- COMPONENTS[1:3]
    
    PROTOCOL <- "daynight_two_stage_v1"
    
    assert <- function(ok, message) {
      if (length(ok) != 1L || is.na(ok) || !ok) 
        stop(message, call. = FALSE)
      invisible(TRUE)
    }
    
    same_numeric <- function(a, b, tolerance = RECONSTRUCTION_REL_TOL) {
      a <- as.numeric(a)
      b <- as.numeric(b)
      length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) && all(abs(a - b) <= tolerance * (1 + pmax(abs(a), 
                                                                                                                  abs(b))))
    }
    
    read_object <- function(path) {
      assert(file.exists(path), paste("Missing input:", path))
      extension <- tolower(tools::file_ext(path))
      if (extension == "rds") 
        return(readRDS(path))
      if (extension == "qs") {
        assert(requireNamespace("qs", quietly = TRUE), "Install qs to read stage-two .qs metadata.")
        return(qs::qread(path))
      }
      stop("Unsupported object format: ", path)
    }
    
    describe_files <- function(paths) {
      paths <- unique(paths[file.exists(paths)])
      paths <- paths[!file.info(paths)$isdir]
      data.table(path = normalizePath(paths, winslash = "/", mustWork = TRUE), bytes = as.numeric(file.info(paths)$size))
    }
    
    check_psd <- function(V, label) {
      assert(is.matrix(V) && is.numeric(V) && nrow(V) == ncol(V) && nrow(V) > 0L && all(is.finite(V)), paste("Invalid covariance:", 
                                                                                                             label))
      scale <- max(abs(V), .Machine$double.eps)
      assert(max(abs(V - t(V))) <= 1e-08 * scale, paste("Asymmetric covariance:", label))
      V <- (V + t(V))/2
      eg <- eigen(V, symmetric = TRUE)
      threshold <- max(abs(eg$values), .Machine$double.eps) * EIGEN_REL_TOL
      assert(min(eg$values) >= -threshold, paste("Non-positive-semidefinite covariance:", label))
      list(matrix = V, eigen = eg, keep = eg$values > threshold)
    }
    
    wald_test <- function(estimate, covariance, required_rank = NULL) {
      estimate <- as.numeric(estimate)
      unavailable <- function(status, rank = NA_integer_) list(statistic = NA_real_, df = rank, p_value = NA_real_, 
                                                               status = status)
      if (any(!is.finite(estimate)) || any(!is.finite(covariance))) 
        return(unavailable("not_estimable_nonfinite"))
      p <- check_psd(as.matrix(covariance), "Wald contrast")
      assert(length(estimate) == nrow(p$matrix), "Wald estimate/covariance size mismatch.")
      rank <- sum(p$keep)
      if (!is.null(required_rank) && rank != required_rank) 
        return(unavailable("not_estimable_rank_deficient", rank))
      if (rank == 0L) 
        return(unavailable("not_estimable_zero_variance", rank))
      if (any(!p$keep)) {
        unsupported <- crossprod(p$eigen$vectors[, !p$keep, drop = FALSE], estimate)
        if (max(abs(unsupported)) > 1e-07 * max(1, sqrt(sum(estimate^2)))) 
          return(unavailable("not_estimable_contrast_outside_covariance_range", rank))
      }
      coordinates <- crossprod(p$eigen$vectors[, p$keep, drop = FALSE], estimate)
      statistic <- sum(as.numeric(coordinates)^2/p$eigen$values[p$keep])
      list(statistic = statistic, df = rank, p_value = pchisq(statistic, rank, lower.tail = FALSE), status = if (rank == 
                                                                                                                 length(estimate)) "estimated" else "estimated_reduced_rank")
    }
    
    inverse_full <- function(V) {
      p <- check_psd(V, "profile")
      assert(all(p$keep), "Profile covariance is rank deficient.")
      sweep(p$eigen$vectors, 2, 1/p$eigen$values, "*") %*% t(p$eigen$vectors)
    }
    
    simulate_coefficients <- function(mu, V, n, seed) {
      p <- check_psd(V, "joint-state coefficients")
      set.seed(seed)
      root <- sweep(p$eigen$vectors, 2, sqrt(pmax(p$eigen$values, 0)), "*")
      Z <- matrix(rnorm(n * length(mu)), n, length(mu))
      out <- sweep(Z %*% t(root), 2, mu, "+")
      colnames(out) <- names(mu)
      out
    }
    
    ci_columns <- function(draws) {
      draws <- as.matrix(draws)
      if (!all(is.finite(draws))) 
        return(matrix(NA_real_, 2, ncol(draws), dimnames = list(c("low", "high"), colnames(draws))))
      out <- vapply(seq_len(ncol(draws)), function(j) quantile(draws[, j], c(0.025, 0.975), type = 8, names = FALSE), 
                    numeric(2))
      dimnames(out) <- list(c("low", "high"), colnames(draws))
      out
    }
    
    design_key <- function(meta) {
      lag_names <- c(day_cold = "day_cold_class_lag_days", day_heat = "day_heat_class_lag_days", night_cold = "night_cold_class_lag_days", 
                     night_heat = "night_heat_class_lag_days")
      thresholds <- c("day_cold_threshold", "day_heat_threshold", "night_cold_threshold", "night_heat_threshold")
      list(calibration_id = meta$calibration_id, analysis_tag = meta$shared_design$analysis_tag, reference_method = meta$reference_method, 
           rr_tolerance = meta$rr_tolerance, thresholds = unlist(meta$state_thresholds[thresholds]), lag_column_indices = lapply(lag_names, 
                                                                                                                                 function(name) as.integer(meta[[name]])), overlap_policy = meta$both_anomaly_policy, reference_source = meta$reference_source_code, 
           settings = meta$shared_design$settings)
    }
    
    read_outcome_bundle <- function(code, index) {
      entry <- index$outcomes[[code]]
      assert(!is.null(entry), paste("No input directory configured for", code))
      root <- normalizePath(entry$output_root, winslash = "/", mustWork = TRUE)
      names <- c("source_parameters.rds", "table_04_joint_9_state_rr.csv", "table_00c_annual_state_death_counts.csv", 
                 "table_12_shapley_4_component_burden_mean_annual.csv")
      paths <- vapply(names, function(name) current_file(root, name), character(1))
      for (name in c("joint_logrr_covariance.rds", "table_04a_joint_logrr_covariance.csv")) current_file(root, name, 
                                                                                                         required = FALSE)
      meta <- readRDS(paths[["source_parameters.rds"]])
      primary <- fread(paths[["table_12_shapley_4_component_burden_mean_annual.csv"]])
      assert(all(c("component", "AN", "AF_percent", "share_percent", "n_deaths") %in% names(primary)) && nrow(primary) == 
               5L && !anyDuplicated(primary$component) && setequal(primary$component, c(COMPONENTS, "joint_total")), paste("Invalid component burden table:", 
                                                                                                                           code))
      ready <- list(stage1_id = meta$stage1_id, calibration_id = meta$calibration_id, analysis_tag = index$active_design$analysis_tag)
      list(code = code, root = root, folder = file.path(root, "tables"), ready = ready, metadata = meta, design = design_key(meta), 
           primary = primary, input_files = describe_files(unname(paths)))
    }
    
    read_risk_counts <- function(bundle) {
      rr <- fread(file.path(bundle$folder, "table_04_joint_9_state_rr.csv"))
      need <- c("state", "beta", "se", "rr", "rr_low", "rr_high", "n_cases", "n_rows", "n_referents", "support_status")
      assert(all(need %in% names(rr)) && nrow(rr) == 9L && !anyDuplicated(rr$state) && setequal(rr$state, STATE_LEVELS), 
             "Invalid nine-state RR table.")
      rr <- rr[match(STATE_LEVELS, state)]
      for (v in c("n_cases", "n_rows", "n_referents")) assert(is.numeric(rr[[v]]) && all(is.finite(rr[[v]])) && all(rr[[v]] >= 
                                                                                                                      0 & rr[[v]] == floor(rr[[v]])), paste("Invalid counts:", v))
      assert(all(rr$n_rows == rr$n_cases + rr$n_referents), "Cases/referents do not sum to rows.")
      observed <- rr$n_rows > 0
      assert(all(rr$support_status == ifelse(observed, "observed_estimable", "not_observed")), "Invalid observed-state support flags.")
      for (v in c("beta", "se", "rr", "rr_low", "rr_high")) assert(is.numeric(rr[[v]]) && all(is.finite(rr[[v]][observed])), 
                                                                   paste("Invalid risk column:", v))
      r <- rr[observed]
      assert(all(r$se >= 0 & r$rr > 0 & r$rr_low > 0 & r$rr_high > 0) && same_numeric(log(r$rr), r$beta) && same_numeric(log(r$rr_low), 
                                                                                                                         r$beta - 1.96 * r$se) && same_numeric(log(r$rr_high), r$beta + 1.96 * r$se), "Inconsistent RR/beta/SE/CI values.")
      reference <- rr[state == REF_STATE]
      assert(reference$n_rows > 0 && reference$beta == 0 && reference$se == 0 && reference$rr == 1, "Invalid common reference state.")
      assert(all(rr[!observed, n_cases] == 0) && all(is.na(rr[!observed, beta])), "Unobserved states must retain non-estimable risks and zero case counts.")
      states <- rr[observed & state != REF_STATE, state]
      assert(length(states) > 0, "No estimable non-reference state.")
      mu <- setNames(rr$beta[match(states, rr$state)], states)
      csv <- file.path(bundle$folder, "table_04a_joint_logrr_covariance.csv")
      rds <- analysis_file(bundle$root, "joint_logrr_covariance.rds")
      candidates <- c(csv, rds)
      candidates <- candidates[file.exists(candidates)]
      assert(length(candidates) > 0, "A full coefficient covariance is required.")
      chosen <- candidates[which.max(as.numeric(file.info(candidates)$mtime))]
      if (tolower(tools::file_ext(chosen)) == "csv") {
        v <- fread(chosen)
        label <- intersect(c("state", "term"), names(v))
        assert(length(label) == 1L, "Covariance CSV requires one state/term column.")
        rn <- v[[label]]
        v[, `:=`((label), NULL)]
        V <- as.matrix(v)
        rownames(V) <- rn
      }
      else {
        saved <- readRDS(chosen)
        V <- if (is.matrix(saved)) 
          saved
        else saved$covariance_non_reference
      }
      assert(is.matrix(V) && is.numeric(V) && !is.null(rownames(V)) && !is.null(colnames(V)), "Invalid covariance matrix.")
      rownames(V) <- sub("^joint_state", "", rownames(V))
      colnames(V) <- sub("^joint_state", "", colnames(V))
      assert(!anyDuplicated(rownames(V)) && !anyDuplicated(colnames(V)) && all(states %in% rownames(V)) && all(states %in% 
                                                                                                                 colnames(V)), "Covariance does not cover the current states.")
      V <- check_psd(V[states, states, drop = FALSE], "current covariance")$matrix
      supplied_se <- rr$se[match(states, rr$state)]
      previous_se <- sqrt(diag(V))
      assert(!any(previous_se == 0 & supplied_se != 0), "A full covariance is needed for nonzero SEs with zero prior variance.")
      ratio <- ifelse(previous_se == 0, 1, supplied_se/previous_se)
      V <- V * outer(ratio, ratio)
      saveRDS(list(mu_non_reference = mu, covariance_non_reference = V), rds)
      fwrite(cbind(data.table(state = states), as.data.table(V)), csv)
      annual <- fread(file.path(bundle$folder, "table_00c_annual_state_death_counts.csv"))
      assert(all(c("year", "state", "n_cases", "n_deaths") %in% names(annual)), "Annual counts lack required columns.")
      for (v in c("year", "n_cases", "n_deaths")) assert(is.numeric(annual[[v]]) && all(is.finite(annual[[v]])) && 
                                                           all(annual[[v]] == floor(annual[[v]])), paste("Invalid annual column:", v))
      assert(all(annual$n_cases >= 0 & annual$n_deaths > 0) && !anyDuplicated(annual[, .(year, state)]), "Invalid or duplicate annual state counts.")
      years <- sort(unique(annual$year))
      assert(length(years) > 0, "No analysis year.")
      for (y in years) {
        z <- annual[year == y]
        assert(nrow(z) == 9L && setequal(z$state, STATE_LEVELS) && uniqueN(z$n_deaths) == 1L && sum(z$n_cases) == 
                 z$n_deaths[1L], paste("Invalid state counts in year", y))
      }
      totals <- annual[, .(total = sum(n_cases)), by = state]
      assert(same_numeric(totals$total, rr$n_cases[match(totals$state, rr$state)], 0), "Annual counts do not match the current risk table.")
      beta <- setNames(rep(NA_real_, 9), STATE_LEVELS)
      beta[REF_STATE] <- 0
      beta[states] <- mu
      list(rr = rr, states = states, mu = mu, covariance = V, beta = beta, annual = annual, years = years)
    }
    
    state_components <- function(state) {
      tokens <- strsplit(state, "_", fixed = TRUE)[[1]]
      mapping <- c(DC = "D_C", DH = "D_H", NC = "N_C", NH = "N_H")
      unname(mapping[tokens[tokens %in% names(mapping)]])
    }
    
    SINGLETON <- c(D_C = "DC_N0", N_C = "D0_NC", D_H = "DH_N0", N_H = "D0_NH")
    
    shapley_operator <- function(weights, estimable_states) {
      assert(setequal(names(weights), STATE_LEVELS) && all(is.finite(weights)) && all(weights >= 0), "Invalid state weights.")
      weights <- weights[STATE_LEVELS]
      B <- matrix(0, 4, 9, dimnames = list(COMPONENTS, STATE_LEVELS))
      for (s in STATE_LEVELS) {
        active <- state_components(s)
        if (weights[[s]] == 0 || length(active) == 0L) 
          next
        needed <- unique(c(s, unname(SINGLETON[active])))
        assert(all(needed %in% estimable_states), paste("Positive-weight state requires unavailable risks:", s, 
                                                        paste(setdiff(needed, estimable_states), collapse = ", ")))
        if (length(active) == 1L) {
          B[active, s] <- B[active, s] + weights[[s]]
        }
        else {
          for (k in active) {
            other <- setdiff(active, k)
            B[k, SINGLETON[[k]]] <- B[k, SINGLETON[[k]]] + weights[[s]]/2
            B[k, s] <- B[k, s] + weights[[s]]/2
            B[k, SINGLETON[[other]]] <- B[k, SINGLETON[[other]]] - weights[[s]]/2
          }
        }
      }
      expected <- weights
      expected[REF_STATE] <- 0
      assert(same_numeric(colSums(B), expected), "Shapley additivity failed.")
      B[, estimable_states, drop = FALSE]
    }
    
    make_operators <- function(risk) {
      annual <- copy(risk$annual)
      summary <- annual[, .(mean_cases = mean(n_cases), mean_fraction = mean(n_cases/n_deaths)), by = state]
      counts <- setNames(summary$mean_cases, summary$state)
      fractions <- setNames(summary$mean_fraction, summary$state)
      assert(abs(sum(fractions) - 1) < 1e-10, "Mean annual state fractions do not sum to one.")
      list(AN = shapley_operator(counts, risk$states), AF_percent = 100 * shapley_operator(fractions, risk$states), 
           mean_cases = counts, mean_fractions = fractions, mean_deaths = mean(unique(annual[, .(year, n_deaths)])$n_deaths))
    }
    
    validate_primary <- function(primary, an, af, mean_deaths) {
      p <- primary[match(c(COMPONENTS, "joint_total"), component)]
      target_an <- c(an, sum(an))
      target_af <- c(af, sum(af))
      target_share <- if (sum(an) == 0) 
        rep(NA_real_, 5)
      else c(an/sum(an) * 100, 100)
      data.table(metric = c("AN", "AF_percent", "share_percent"), max_absolute_difference = c(max(abs(p$AN - target_an)), 
                                                                                              max(abs(p$AF_percent - target_af)), if (all(is.na(target_share))) NA_real_ else max(abs(p$share_percent - 
                                                                                                                                                                                        target_share))))
    }
    
    profile_diagnostics <- function(an, draws) {
      draws <- as.matrix(draws)
      assert(ncol(draws) == 4L && nrow(draws) >= 4L && all(COMPONENTS %in% names(an)), "Profile diagnostics require four named components and at least four draws.")
      if (is.null(colnames(draws))) 
        colnames(draws) <- COMPONENTS
      assert(setequal(colnames(draws), COMPONENTS), "Unexpected component names in profile draws.")
      draws <- draws[, COMPONENTS, drop = FALSE]
      an <- an[COMPONENTS]
      total <- sum(an)
      total_draws <- rowSums(draws)
      finite <- all(is.finite(draws)) && all(is.finite(total_draws))
      total_ci <- if (finite) 
        quantile(total_draws, c(0.025, 0.975), type = 8, names = FALSE)
      else c(NA_real_, NA_real_)
      epsilon <- TOTAL_AN_REL_EPS * max(1, sum(abs(an)))
      point_ok <- is.finite(total) && abs(total) > epsilon
      excludes_zero <- finite && (total_ci[1] > 0 || total_ci[2] < 0)
      sign_switch <- if (finite && point_ok) 
        mean(sign(total_draws) != sign(total))
      else NA_real_
      near_zero <- if (finite) 
        mean(abs(total_draws) <= epsilon)
      else NA_real_
      ratio_draws <- matrix(NA_real_, nrow(draws), 4, dimnames = list(NULL, COMPONENTS))
      V <- matrix(NA_real_, 3, 3, dimnames = list(PROFILE_COORDINATES, PROFILE_COORDINATES))
      reasons <- character()
      if (!finite) 
        reasons <- c(reasons, "nonfinite_AN_draws")
      if (!point_ok) 
        reasons <- c(reasons, "total_AN_near_zero")
      if (!excludes_zero) 
        reasons <- c(reasons, "total_AN_interval_includes_zero_or_unavailable")
      if (is.finite(sign_switch) && sign_switch > MAX_TOTAL_SIGN_SWITCH_FRACTION) 
        reasons <- c(reasons, "excess_total_sign_switching")
      if (is.finite(near_zero) && near_zero > 0) 
        reasons <- c(reasons, "near_zero_simulated_total")
      split_difference <- NA_real_
      rank <- NA_integer_
      max_abs_share <- NA_real_
      if (finite && all(total_draws != 0)) {
        ratio_draws <- sweep(draws, 1, total_draws, "/") * 100
        if (all(is.finite(ratio_draws))) {
          max_abs_share <- max(abs(ratio_draws))
          z <- ratio_draws[, PROFILE_COORDINATES, drop = FALSE]
          V <- cov(z)
          split <- floor(nrow(z)/2)
          V1 <- cov(z[seq_len(split), , drop = FALSE])
          V2 <- cov(z[seq.int(split + 1L, nrow(z)), , drop = FALSE])
          split_difference <- norm(V1 - V2, "F")/max(norm(V, "F"), .Machine$double.eps)
          rank <- sum(check_psd(V, "simulated profile")$keep)
          if (rank != 3L) 
            reasons <- c(reasons, "profile_covariance_rank_not_three")
          if (!is.finite(split_difference) || split_difference > MAX_PROFILE_SPLIT_COV_REL_DIFF) 
            reasons <- c(reasons, "unstable_split_half_profile_covariance")
        }
        else reasons <- c(reasons, "nonfinite_share_draws")
      }
      else reasons <- c(reasons, "share_draws_unavailable")
      stable <- length(reasons) == 0L
      status <- if (stable) 
        "eligible"
      else paste(unique(reasons), collapse = "; ")
      list(eligible = stable, status = status, covariance = V, share_draws = ratio_draws, row = data.table(total_mean_annual_AN = total, 
                                                                                                           total_AN_CI_low = total_ci[1], total_AN_CI_high = total_ci[2], finite_AN_draw_fraction = mean(apply(is.finite(draws), 
                                                                                                                                                                                                               1, all)), total_sign_switch_fraction = sign_switch, near_zero_total_fraction = near_zero, max_absolute_simulated_share_percent = max_abs_share, 
                                                                                                           profile_covariance_rank = rank, split_half_covariance_relative_difference = split_difference, profile_eligible = stable, 
                                                                                                           profile_status = status, n_sim = nrow(draws), n_discarded_draws = 0L))
    }
    
    analyse_outcome <- function(bundle, outcome_name, n_sim, seed) {
      risk <- read_risk_counts(bundle)
      operators <- make_operators(risk)
      f <- 1 - exp(-risk$mu)
      an <- setNames(as.numeric(operators$AN %*% f), COMPONENTS)
      af <- setNames(as.numeric(operators$AF_percent %*% f), COMPONENTS)
      validation <- validate_primary(bundle$primary, an, af, operators$mean_deaths)
      beta_draws <- simulate_coefficients(risk$mu, risk$covariance, n_sim, seed)
      f_draws <- 1 - exp(-beta_draws)
      an_draws <- f_draws %*% t(operators$AN)
      af_draws <- f_draws %*% t(operators$AF_percent)
      colnames(an_draws) <- colnames(af_draws) <- COMPONENTS
      covariance_an <- if (all(is.finite(an_draws))) 
        cov(an_draws)
      else matrix(NA_real_, 4, 4)
      profile <- profile_diagnostics(an, an_draws)
      observed_share <- if (sum(an) == 0) 
        setNames(rep(NA_real_, 4), COMPONENTS)
      else an/sum(an) * 100
      C <- matrix(c(1, -1, 0, 0, 1, 0, -1, 0, 1, 0, 0, -1), nrow = 3, byrow = TRUE, dimnames = list(c("D_C_minus_N_C", 
                                                                                                      "D_C_minus_D_H", "D_C_minus_N_H"), COMPONENTS))
      within <- wald_test(C %*% an, C %*% covariance_an %*% t(C))
      global_risk <- wald_test(risk$mu, risk$covariance)
      peak <- risk$states[which.max(risk$mu)]
      peak_row <- risk$rr[state == peak]
      risk_row <- data.table(outcome_code = bundle$code, outcome = outcome_name, chi_square = global_risk$statistic, 
                             df = global_risk$df, p_value = global_risk$p_value, status = global_risk$status, n_estimable_nonreference_states = length(risk$states), 
                             unobserved_states = paste(setdiff(STATE_LEVELS, c(REF_STATE, risk$states)), collapse = "; "), highest_RR_state_descriptive = peak, 
                             highest_RR = peak_row$rr, highest_RR_CI_low = peak_row$rr_low, highest_RR_CI_high = peak_row$rr_high, null_hypothesis = "All estimable non-reference joint-state log relative risks equal zero", 
                             interpretation = "Global association with the reference included; the highest RR is descriptive")
      within_row <- data.table(outcome_code = bundle$code, outcome = outcome_name, chi_square = within$statistic, 
                               df = within$df, p_value = within$p_value, status = within$status, null_hypothesis = "Four mean annual component attributable numbers are equal", 
                               total_mean_annual_AN = sum(an), profile_eligible = profile$eligible, profile_status = profile$status)
      an_ci <- ci_columns(an_draws)
      af_ci <- ci_columns(af_draws)
      share_ci <- if (profile$eligible) 
        ci_columns(profile$share_draws)
      else matrix(NA_real_, 2, 4, dimnames = list(c("low", "high"), COMPONENTS))
      components <- data.table(outcome_code = bundle$code, outcome = outcome_name, component = COMPONENTS, component_label = unname(COMPONENT_LABELS[COMPONENTS]), 
                               AN = as.numeric(an), AN_CI_low = an_ci[1, ], AN_CI_high = an_ci[2, ], AF_percent = as.numeric(af), AF_percent_CI_low = af_ci[1, 
                               ], AF_percent_CI_high = af_ci[2, ], share_percent = as.numeric(observed_share), share_CI_low = share_ci[1, 
                               ], share_CI_high = share_ci[2, ], profile_status = profile$status)
      diagnostics <- copy(profile$row)
      diagnostics[, `:=`(outcome_code = bundle$code, outcome = outcome_name, simulation_seed = seed, negative_component_count = sum(an < 
                                                                                                                                      0), n_cases = sum(risk$annual$n_cases), years = paste(risk$years, collapse = "; "), max_AN_reconstruction_error = validation[metric == 
                                                                                                                                                                                                                                                                     "AN", max_absolute_difference], max_AF_reconstruction_error = validation[metric == "AF_percent", max_absolute_difference], 
                         max_share_reconstruction_error = validation[metric == "share_percent", max_absolute_difference])]
      validation[, `:=`(outcome_code, bundle$code)]
      list(risk = risk_row, within = within_row, components = components, diagnostics = diagnostics, validation = validation, 
           years = risk$years, AN = an, AF_percent = af, share = observed_share, profile = list(eligible = profile$eligible, 
                                                                                                status = profile$status, z = observed_share[PROFILE_COORDINATES], covariance = profile$covariance), 
           covariance_AN = covariance_an, coefficient_mean = risk$mu, coefficient_covariance = risk$covariance, draws = if (SAVE_SIMULATION_DRAWS) list(AN = an_draws, 
                                                                                                                                                        AF_percent = af_draws, share_percent = profile$share_draws) else NULL)
    }
    
    between_disease_tests <- function(results, alpha = ALPHA) {
      codes <- OUTCOME_CONFIG$outcome_code[OUTCOME_CONFIG$between_disease]
      available <- codes %in% names(results)
      reason <- character()
      if (!all(available)) 
        reason <- c(reason, paste("missing_disease_outcomes", paste(codes[!available], collapse = ",")))
      if (all(available)) {
        unsuitable <- codes[!vapply(results[codes], function(x) x$profile$eligible, logical(1))]
        if (length(unsuitable)) 
          reason <- c(reason, paste("ineligible_profiles", paste(unsuitable, collapse = ",")))
        if (!all(vapply(results[codes], function(x) identical(x$years, results[[codes[1]]]$years), logical(1)))) 
          reason <- c(reason, "different_analysis_years")
      }
      omnibus <- data.table(test = "Five-disease multivariate burden-profile Wald/Q test", n_diseases_planned = 5L, 
                            profile_dimensions = 3L, chi_square = NA_real_, df = 12L, p_value = NA_real_, status = if (length(reason)) 
                              paste(reason, collapse = "; ")
                            else "estimated", covariance_assumption = "Approximately independent disease-specific coefficient errors; shared calibration and counts fixed")
      if (!length(reason)) {
        precision <- lapply(results[codes], function(x) inverse_full(x$profile$covariance))
        pooled <- as.numeric(solve(Reduce(`+`, precision), Reduce(`+`, lapply(seq_along(codes), function(i) precision[[i]] %*% 
                                                                                results[[codes[i]]]$profile$z))))
        Q <- sum(vapply(seq_along(codes), function(i) {
          difference <- results[[codes[i]]]$profile$z - pooled
          as.numeric(crossprod(difference, precision[[i]] %*% difference))
        }, numeric(1)))
        omnibus[, `:=`(chi_square = Q, p_value = pchisq(Q, 12L, lower.tail = FALSE))]
      }
      perform <- is.finite(omnibus$p_value) && omnibus$p_value < alpha
      pairs <- combn(codes, 2, simplify = FALSE)
      rows <- lapply(pairs, function(pair) {
        row <- data.table(disease_1 = pair[1], disease_2 = pair[2], chi_square = NA_real_, df = 3L, p_value = NA_real_, 
                          status = if (perform) 
                            "estimated"
                          else if (length(reason)) 
                            "not_performed_omnibus_unavailable"
                          else "not_performed_omnibus_not_significant")
        for (k in COMPONENTS) row[, `:=`((paste0(k, "_share_difference_pp")), NA_real_)]
        if (perform) {
          x <- results[[pair[1]]]
          y <- results[[pair[2]]]
          test <- wald_test(x$profile$z - y$profile$z, x$profile$covariance + y$profile$covariance, required_rank = 3L)
          row[, `:=`(chi_square = test$statistic, df = test$df, p_value = test$p_value, status = test$status)]
          for (k in COMPONENTS) row[, `:=`((paste0(k, "_share_difference_pp")), x$share[[k]] - y$share[[k]])]
        }
        row
      })
      pairwise <- rbindlist(rows)
      pairwise[, `:=`(p_holm, p.adjust(p_value, method = "holm", n = 10L))]
      pairwise[, `:=`(significant_after_holm, fifelse(is.finite(p_holm), p_holm < alpha, NA))]
      pairwise[, `:=`(disease_1_name = OUTCOME_CONFIG$outcome[match(disease_1, OUTCOME_CONFIG$outcome_code)], disease_2_name = OUTCOME_CONFIG$outcome[match(disease_2, 
                                                                                                                                                            OUTCOME_CONFIG$outcome_code)])]
      list(omnibus = omnibus, pairwise = pairwise)
    }
    
    add_holm_family <- function(x, alpha) {
      x <- copy(x)
      x[, `:=`(p_holm_across_outcomes, p.adjust(p_value, method = "holm", n = nrow(x)))]
      x[, `:=`(significant_after_holm, fifelse(is.finite(p_holm_across_outcomes), p_holm_across_outcomes < alpha, 
                                               NA))]
      x
    }
    
    format_p <- function(p) {
      out <- rep("NE", length(p))
      valid <- is.finite(p)
      out[valid & p < 0.001] <- "<0.001"
      out[valid & p >= 0.001] <- sprintf("%.3f", p[valid & p >= 0.001])
      out
    }
    
  }, envir = e)
  e
}

# 7. Execute the selected analyses ----------------------------------------------
if (isTRUE(getOption("daynight.sensitivity.run", TRUE))) run_sensitivity_analysis()
