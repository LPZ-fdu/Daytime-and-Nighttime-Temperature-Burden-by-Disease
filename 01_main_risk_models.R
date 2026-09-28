#!/usr/bin/env Rscript
# =============================================================================
# Primary risk estimation
#
# Purpose: fit separate daytime and nighttime distributed lag nonlinear models,
# derive shared reference intervals and lag windows from overall mortality, and
# fit a nine-state conditional logistic model for each of six outcomes.
#
# Inputs: RDS data frames listed under kind=matched in config/input_files.csv.
# Each row is a case or referent date. Required columns: id, date, case, holiday,
# temp_day_lag01:21, temp_night_lag01:21, rh_day_lag01:03, rh_night_lag01:03.
# Temperature and humidity suffix 01 represents cycle lag 0. Case and holiday
# are binary; id identifies one death and its matched referents within a file.
# The overall outcome combines the I00_I52 and I60_I69 source files.
# Input construction, units, cycle timing, and source identifiers are in README.md.
#
# Outputs: results/main/<outcome>/joint_analysis_<configuration>/ contains tables,
# figures, models, and plot_data. Tables include continuous curves, shared reference
# definitions, selected windows, nine-state risks and full coefficient covariance,
# annual state death counts, matched-comparison support, and exclusion diagnostics.
# Shared calibration is stored in results/main/_shared_overall_calibration.
#
# Primary settings: connected cumulative RR <= 1.01 reference region; common
# windows from lag 0 through the last significantly harmful lag; overlap exclusion.
# Outcome switches below select individual analyses. Run from the repository root.
# Functions only: options(daynight.stage1.skip_calls = TRUE); source(this_file).
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

analysis_file <- function(root, name) {
  section <- if (grepl("^plotdata_", name)) "plot_data" else if (
    grepl("^(input_manifest|output_manifest|source_file_manifest|output_file_manifest)[.]", name)) "models" else if (
      grepl("[.](csv|xlsx)$", name, ignore.case = TRUE)) "tables" else if (
        grepl("[.](png|pdf|svg|tif|tiff|jpg|jpeg)$", name, ignore.case = TRUE)) "figures" else "models"
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

RESULTS_ROOT <- file.path(RESULTS_BASE, "main")
PIPELINE_INDEX_FILE <- file.path(RESULTS_BASE, "daynight_pipeline_index.rds")
SHARED_DIR <- file.path(RESULTS_ROOT, "_shared_overall_calibration")
OBJECT_FORMAT <- "rds"  # "rds" or "qs" for model/result objects

OVERALL_CODE <- "I00_I52_I60_I69"
OVERALL_NAME <- "Overall cardiocerebrovascular mortality"
DATA_PATHS <- c(
  I00_I52 = MATCHED_FILES[["I00_I52"]],
  I60_I69 = MATCHED_FILES[["I60_I69"]]
)

# Reference interval: "mmt_ci" or "rr_tolerance".
RUN_REFERENCE_METHOD <- "rr_tolerance"
RUN_RR_TOLERANCE <- 0.01
RUN_ALLOW_TRUNCATED_RR_BAND <- FALSE

# Automatic selection: lag 0 through the last harmful significant lag.
# Manual day counts include lag 0: a three-day mean uses lags 0, 1 and 2.
RUN_LAG_MODE <- "auto"  # "auto" or "manual"
RUN_MANUAL_MA_DAYS <- c(day_cold = 21L, day_heat = 3L,
                        night_cold = 21L, night_heat = 3L)
RUN_MANUAL_OVERRIDES <- NULL  # e.g. c(day_heat = 3L) within automatic mode
RUN_NO_SIG_ACTION <- "fallback"  # "fallback" or "stop"
RUN_FALLBACK_MA_DAYS <- c(cold = 21L, heat = 3L)

# Outcomes are called individually at the end of the script.
RUN_OVERALL <- TRUE
RUN_I10_I15 <- TRUE
RUN_I20_I25 <- TRUE
RUN_I50 <- TRUE
RUN_I60_I62 <- TRUE
RUN_I63 <- TRUE

CURRENT_SHARED_DESIGN <- NULL
STATE_THRESHOLD_STRATEGY <- NULL

# ============================================================
# Model settings, data preparation and distributed lag models
# ============================================================

required_packages <- c("data.table", "survival", "splines", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
  library(ggplot2)
})

DATA_PATH <- NULL
OUT_BASE <- NULL
OUT_ROOT <- NULL
OUTCOME_CODE <- NULL
OUTCOME_NAME <- NULL

LAG_DAYS  <- 1:21
LAG_INDEX <- 0:(length(LAG_DAYS) - 1)

VAR_DF <- 4
RH_DF  <- 3

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
BOTH_ANOMALY_POLICY <- "exclude" # supported: "exclude", "heat_priority", "cold_priority"

N_SIM <- 1000
SIM_SEED <- 20260523
CAP_NEGATIVE_AF <- FALSE
SAVE_CASE_CONTRIB <- FALSE

# NULL detects common geography fields; character(0) disables geographic summaries.
# Specify source column names when needed, e.g. c("province_code", "city_code").
DIAGNOSTIC_GEOGRAPHY_COLS <- NULL

FONT_FAMILY <- "serif"
FIG_DPI <- 600
SAVE_PDF <- FALSE
DAY_COLOUR   <- "#D55E00"
NIGHT_COLOUR <- "#0072B2"
COLD_COLOUR  <- "#2166AC"
HEAT_COLOUR  <- "#B2182B"

set_output_dirs <- function(out_root, log_file_name, include_contrib = FALSE) {
  prepare_analysis_directory(out_root)
  DIR_TABLE <<- file.path(out_root, "tables")
  DIR_MODEL <<- file.path(out_root, "models")
  DIR_PLOT_DATA <<- file.path(out_root, "plot_data")
  DIR_FIG <<- file.path(out_root, "figures")
  DIR_LOG <<- file.path(out_root, "models")
  DIR_CONTRIB <<- file.path(out_root, "tables")
  
  dirs <- c(DIR_TABLE, DIR_MODEL, DIR_PLOT_DATA, DIR_FIG, DIR_LOG)
  if (isTRUE(include_contrib)) dirs <- c(dirs, DIR_CONTRIB)
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
  if (identical(OBJECT_FORMAT, "rds")) sub("[.]qs$", ".rds", path) else path
}

write_model_object <- function(obj, file) {
  file <- object_path(file)
  if (OBJECT_FORMAT == "rds") saveRDS(obj, file) else {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Install qs or select OBJECT_FORMAT = 'rds'.")
    qs::qsave(obj, file, preset = "fast")
  }
  invisible(file)
}

read_model_object <- function(file) {
  file <- object_path(file)
  if (tolower(tools::file_ext(file)) == "rds") readRDS(file) else {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Reading qs objects requires package qs.")
    qs::qread(file)
  }
}

save_obj <- function(obj, file) {
  file <- write_model_object(obj, file)
  log_msg("Saved object:", file)
}

coerce_binary01 <- function(x, variable_name) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  
  if (is.character(x)) {
    x_trim <- trimws(tolower(x))
    mapped <- rep(NA_integer_, length(x_trim))
    mapped[x_trim %in% c("0", "false", "no")] <- 0L
    mapped[x_trim %in% c("1", "true", "yes")] <- 1L
    x <- mapped
  } else {
    x <- suppressWarnings(as.numeric(x))
  }
  
  finite_values <- sort(unique(x[is.finite(x)]))
  if (!all(finite_values %in% c(0, 1))) {
    stop(
      variable_name,
      " must contain only binary values coded as 0 and 1. Observed finite values: ",
      paste(finite_values, collapse = ", ")
    )
  }
  
  as.integer(x)
}

validate_unique_key <- function(dt, key_cols, object_name) {
  if (nrow(dt) == 0L) stop(object_name, " is empty.")
  duplicate_rows <- dt[, .N, by = key_cols][N > 1L]
  if (nrow(duplicate_rows) > 0L) {
    stop(
      object_name,
      " contains duplicate rows for key: ",
      paste(key_cols, collapse = ", ")
    )
  }
  invisible(TRUE)
}

validate_finite_named_vector <- function(x, expected_names, object_name) {
  if (is.null(names(x))) stop(object_name, " must be a named vector.")
  missing_names <- setdiff(expected_names, names(x))
  if (length(missing_names) > 0L) {
    stop(
      object_name,
      " is missing required elements: ",
      paste(missing_names, collapse = ", ")
    )
  }
  bad_names <- expected_names[!is.finite(x[expected_names])]
  if (length(bad_names) > 0L) {
    stop(
      object_name,
      " contains non-finite values for: ",
      paste(bad_names, collapse = ", ")
    )
  }
  invisible(TRUE)
}

safe_empirical_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(quantile(x, probs = probability, na.rm = TRUE, type = 8))
}

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY, colour = "black"),
      axis.title = element_text(size = base_size + 1, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.line = element_line(size = 0.35, colour = "black"),
      axis.ticks = element_line(size = 0.35, colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = "grey70", size = 0.25),
      strip.text = element_text(size = base_size, colour = "black"),
      legend.title = element_text(size = base_size, colour = "black"),
      legend.text = element_text(size = base_size - 1, colour = "black"),
      panel.grid = element_blank(),
      plot.title = element_blank()
    )
}

theme_nature <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY, colour = "black"),
      axis.title = element_text(size = base_size + 1, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.line = element_line(size = 0.35, colour = "black"),
      axis.ticks = element_line(size = 0.35, colour = "black"),
      legend.title = element_text(size = base_size, colour = "black"),
      legend.text = element_text(size = base_size - 1, colour = "black"),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank()
    )
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
  if (length(x) == 0) stop("No finite values for spline basis.")
  if (df <= 1) stop("df must be > 1 for ns basis.")
  
  n_knots <- df - 1
  probs <- seq(0, 1, length.out = n_knots + 2)[-c(1, n_knots + 2)]
  
  list(
    type = "ns_quantile_knots",
    df = df,
    n_internal_knots = n_knots,
    knots = as.numeric(quantile(x, probs = probs, na.rm = TRUE, type = 8)),
    Boundary.knots = range(x, na.rm = TRUE),
    intercept = FALSE
  )
}

make_lag_spec_logknots <- function(lag_index, n_internal_knots = 3, intercept = TRUE) {
  lag_index <- sort(unique(as.numeric(lag_index)))
  lag_index <- lag_index[is.finite(lag_index)]
  
  if (length(lag_index) < 3) stop("lag_index should contain at least 3 values.")
  if (any(lag_index < 0)) stop("lag_index should be non-negative for log-knots.")
  
  lag_max <- max(lag_index)
  if (!is.finite(lag_max) || lag_max <= 1) stop("Maximum lag should be > 1 for log-knots.")
  
  knots <- exp(seq(log(1), log(lag_max), length.out = n_internal_knots + 2))[-c(1, n_internal_knots + 2)]
  
  list(
    type = "ns_log_knots",
    df = NULL,
    n_internal_knots = n_internal_knots,
    knots = as.numeric(knots),
    Boundary.knots = range(lag_index),
    intercept = intercept
  )
}

predict_ns_basis <- function(x, spec) {
  as.matrix(ns(
    x,
    knots = spec$knots,
    Boundary.knots = spec$Boundary.knots,
    intercept = spec$intercept
  ))
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
  if (lag_position < 1 || lag_position > length(LAG_DAYS)) stop("Invalid lag_position.")
  
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
  
  model_dt <- data.table(
    case = dt$case,
    id = dt$id,
    holiday = dt$holiday,
    rh_lag01_03 = dt$rh_lag01_03
  )
  
  rh_basis <- as.data.table(ns(model_dt$rh_lag01_03, df = RH_DF))
  setnames(rh_basis, paste0("rh_ns", seq_len(ncol(rh_basis))))
  model_dt <- cbind(model_dt, rh_basis, design_dt)
  
  rhs_terms <- c(names(design_dt), names(rh_basis), "holiday")
  form <- as.formula(paste0("case ~ ", paste(rhs_terms, collapse = " + "), " + strata(id)"))
  formula_text <- paste(deparse(form), collapse = " ")
  
  log_msg("Fitting model:", model_name)
  
  fit <- clogit(
    form,
    data = model_dt,
    method = CLOGIT_METHOD,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  
  coefficients <- coef(fit)
  covariance <- vcov(fit)
  
  fitted_exposure_coefficients <- coefficients[names(design_dt)]
  if (length(fitted_exposure_coefficients) != ncol(design_dt) ||
      any(!is.finite(fitted_exposure_coefficients))) {
    stop(
      "The model ",
      model_name,
      " did not return finite coefficients for every exposure-basis term."
    )
  }
  
  compact_fit <- list(
    coefficients = coefficients,
    covariance = covariance,
    formula_text = formula_text,
    exposure_terms = names(design_dt),
    rh_terms = names(rh_basis),
    model_name = model_name,
    storage_mode = "compact_coefficients_and_covariance_only"
  )
  
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
  
  data.table(
    logrr = logrr,
    se = se,
    rr = exp(logrr),
    rr_low = exp(logrr - 1.96 * se),
    rr_high = exp(logrr + 1.96 * se)
  )
}

find_mmt_from_single_model <- function(fit_obj, temp_seq, var_spec, lag_spec, prefix, temp_for_median) {
  terms <- fit_obj$exposure_terms
  b <- fit_obj$coefficients[terms]
  x <- cumulative_design_fast(temp_seq, var_spec, lag_spec, prefix)[, terms, drop = FALSE]
  eta <- function(t) as.numeric(cumulative_design_fast(t, var_spec, lag_spec, prefix)[, terms, drop = FALSE] %*% b)
  refine_mmt(temp_seq, as.numeric(x %*% b), eta)
}

rmvnorm_eigen <- function(n, mu, Sigma) {
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  p <- length(mu)
  
  if (n < 1L || p < 1L) stop("The simulation size and coefficient dimension must be positive.")
  if (!all(dim(Sigma) == c(p, p))) stop("The covariance matrix has incompatible dimensions.")
  if (any(!is.finite(mu)) || any(!is.finite(Sigma))) {
    stop("Non-finite values were found in the simulation mean or covariance matrix.")
  }
  
  Sigma <- (Sigma + t(Sigma)) / 2
  eg <- eigen(Sigma, symmetric = TRUE)
  tolerance <- max(1, max(abs(eg$values))) * 1e-8
  if (min(eg$values) < -tolerance) {
    warning(
      "The covariance matrix had negative eigenvalues; negative values were truncated to zero for simulation."
    )
  }
  
  vals <- pmax(eg$values, 0)
  A <- eg$vectors %*% diag(sqrt(vals), nrow = p)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  sweep(Z %*% t(A), 2, mu, "+")
}

simulate_mmt_ci_from_single_model <- function(fit_obj, temp_seq, var_spec, lag_spec, prefix,
                                              n_sim = N_MMT_CI_SIM,
                                              seed = MMT_CI_SEED) {
  terms <- fit_obj$exposure_terms
  x_curve <- build_cb_constant(temp_seq, var_spec, lag_spec, prefix)
  
  beta <- fit_obj$coefficients[terms]
  vc <- fit_obj$covariance[terms, terms, drop = FALSE]
  
  if (length(beta) == 0 || any(!is.finite(beta)) || any(!is.finite(vc))) {
    warning("Non-finite beta/vcov found. MMT CI will be returned as NA.")
    return(list(
      mmt_low = NA_real_,
      mmt_high = NA_real_,
      mmt_sim = numeric(0),
      n_sim = 0L
    ))
  }
  
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- .Random.seed
  } else {
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
    eta <- function(t) as.numeric(cumulative_design_fast(t, var_spec, lag_spec, prefix)[, terms, drop = FALSE] %*% beta_draw[i, ])
    refine_mmt(temp_seq, eta_draw[, i], eta)
  }, numeric(1))
  
  list(
    mmt_low = as.numeric(quantile(mmt_sim, 0.025, na.rm = TRUE, type = 8)),
    mmt_high = as.numeric(quantile(mmt_sim, 0.975, na.rm = TRUE, type = 8)),
    mmt_sim = mmt_sim,
    n_sim = length(mmt_sim)
  )
}

make_cumulative_curve <- function(fit_obj, temp_seq, mmt, var_spec, lag_spec, prefix) {
  x_curve <- build_cb_constant(temp_seq, var_spec, lag_spec, prefix)
  x_ref <- build_cb_constant(mmt, var_spec, lag_spec, prefix)
  out <- predict_contrast_from_design(fit_obj, x_curve, x_ref)
  out[, temperature := temp_seq]
  out[]
}

make_lag_curve <- function(fit_obj, temp_value, ref_value, var_spec, lag_spec, prefix) {
  ans <- rbindlist(lapply(seq_along(LAG_DAYS), function(ii) {
    x_lag <- build_cb_one_lag(temp_value, ii, var_spec, lag_spec, prefix)
    ref_lag <- build_cb_one_lag(ref_value, ii, var_spec, lag_spec, prefix)
    tmp <- predict_contrast_from_design(fit_obj, x_lag, ref_lag)
    tmp[, lag_day := LAG_INDEX[ii]]
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
  } else {
    x_lab <- mmt_value + 0.03 * x_rng
    hjust_lab <- 0
  }
  
  data.table(
    x = x_lab,
    y = y_max - 0.08 * y_rng,
    label = sprintf("MMT = %.1f °C", mmt_value),
    hjust = hjust_lab
  )
}

plot_cumulative_curve <- function(curve_dt, aa, method_suffix) {
  mmt_value <- unique(curve_dt$mmt)
  mmt_low <- unique(curve_dt$mmt_low)
  mmt_high <- unique(curve_dt$mmt_high)
  
  p_curve <- ggplot(curve_dt, aes(x = temperature, y = rr)) +
    geom_ribbon(aes(ymin = rr_low, ymax = rr_high), fill = aa$colour, alpha = 0.22) +
    geom_line(colour = aa$colour, size = 0.65) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.35) +
    geom_vline(xintercept = mmt_value, linetype = "dotted", size = 0.45, colour = aa$colour) +
    geom_vline(xintercept = c(mmt_low, mmt_high), linetype = "dashed", size = 0.35, colour = aa$colour, alpha = 0.85) +
    geom_vline(xintercept = unique(c(curve_dt$p025, curve_dt$p975)), linetype = "longdash", size = 0.30, colour = "grey35") +
    labs(x = aa$xlab, y = "Cumulative relative risk") +
    theme_pub()
  
  mmt_lab <- make_mmt_annot(curve_dt, mmt_value)
  if (is.finite(mmt_low) && is.finite(mmt_high)) {
    mmt_text <- sprintf("MMT = %.1f °C\n95%% CI: %.1f–%.1f °C", mmt_value, mmt_low, mmt_high)
  } else {
    mmt_text <- sprintf("MMT = %.1f °C", mmt_value)
  }
  
  p_curve <- p_curve +
    annotate(
      "text",
      x = mmt_lab$x,
      y = mmt_lab$y,
      label = mmt_text,
      hjust = mmt_lab$hjust,
      vjust = 1,
      size = 3,
      family = FONT_FAMILY
    )
  
  ggsave_png(
    p_curve,
    paste0("fig_01_cumulative_curve_", aa$exposure, "_single_", method_suffix),
    width = 3.8,
    height = 3.2
  )
}

plot_lag_curve <- function(lag_dt, aa, method_suffix) {
  p_lag <- ggplot(lag_dt, aes(x = lag_day, y = rr, colour = contrast_label, fill = contrast_label)) +
    geom_ribbon(aes(ymin = rr_low, ymax = rr_high), alpha = 0.18, colour = NA) +
    geom_line(size = 0.65) +
    geom_hline(yintercept = 1, linetype = "dashed", size = 0.35) +
    scale_colour_manual(
      values = c("P2.5 vs MMT" = COLD_COLOUR, "P97.5 vs MMT" = HEAT_COLOUR),
      breaks = c("P2.5 vs MMT", "P97.5 vs MMT"),
      labels = c(expression(P[2.5]~"vs MMT"), expression(P[97.5]~"vs MMT")),
      name = NULL
    ) +
    scale_fill_manual(
      values = c("P2.5 vs MMT" = COLD_COLOUR, "P97.5 vs MMT" = HEAT_COLOUR),
      breaks = c("P2.5 vs MMT", "P97.5 vs MMT"),
      labels = c(expression(P[2.5]~"vs MMT"), expression(P[97.5]~"vs MMT")),
      name = NULL
    ) +
    labs(x = "Lag days", y = "Relative risk") +
    theme_pub() +
    theme(legend.position = "top")
  
  ggsave_png(
    p_lag,
    paste0("fig_02_lag_effects_", aa$exposure, "_single_", method_suffix),
    width = 4.0,
    height = 3.2
  )
}

prepare_single_model_data <- function(dt_raw) {
  required_base <- c("id", "date", "case", "holiday")
  day_cols <- sprintf("temp_day_lag%02d", LAG_DAYS)
  night_cols <- sprintf("temp_night_lag%02d", LAG_DAYS)
  rh_day_cols <- sprintf("rh_day_lag%02d", 1:3)
  rh_night_cols <- sprintf("rh_night_lag%02d", 1:3)
  
  required_cols <- c(required_base, day_cols, night_cols, rh_day_cols, rh_night_cols)
  missing_cols <- setdiff(required_cols, names(dt_raw))
  if (length(missing_cols) > 0) stop("Missing required columns for the single-temperature analysis: ", paste(missing_cols, collapse = ", "))
  
  dt <- copy(dt_raw[, ..required_cols])
  
  if (is.finite(TEST_N)) {
    log_msg("TEST_N is finite. Keeping complete sets whose ids occur in the first", TEST_N, "rows.")
    test_ids <- unique(dt$id[seq_len(min(TEST_N, nrow(dt)))])
    dt <- dt[id %in% test_ids]
  }
  
  dt[, date := as.IDate(date)]
  dt[, case := coerce_binary01(case, "case")]
  dt[, holiday := coerce_binary01(holiday, "holiday")]
  
  rh_short_cols <- c(sprintf("rh_day_lag%02d", 1:3), sprintf("rh_night_lag%02d", 1:3))
  dt[, rh_lag01_03 := rowMeans(.SD, na.rm = FALSE), .SDcols = rh_short_cols]
  dt[, (rh_short_cols) := NULL]
  
  model_vars_for_complete <- c("id", "date", "case", "holiday", "rh_lag01_03", day_cols, night_cols)
  log_msg("Filtering complete cases for the single-temperature models...")
  dt <- dt[dt[, complete.cases(.SD), .SDcols = model_vars_for_complete]]
  
  log_msg("Filtering valid matched sets for the single-temperature models...")
  valid_id <- dt[, .(n_case = sum(case == 1), n_ref = sum(case == 0)), by = id][n_case == 1 & n_ref >= 1, id]
  dt <- dt[id %in% valid_id]
  rm(valid_id)
  gc()
  
  if (!nrow(dt)) stop("No complete valid matched sets remain for the single-temperature models.")
  log_msg("Rows after filtering:", nrow(dt))
  log_msg("Cases:", dt[, sum(case == 1)])
  log_msg("Matched sets:", dt[, uniqueN(id)])
  
  list(
    dt = dt,
    day_cols = day_cols,
    night_cols = night_cols,
    rh_day_cols = rh_day_cols,
    rh_night_cols = rh_night_cols
  )
}

diagnostic_geography_columns <- function(column_names) {
  if (!is.null(DIAGNOSTIC_GEOGRAPHY_COLS)) {
    if (!is.character(DIAGNOSTIC_GEOGRAPHY_COLS) || anyNA(DIAGNOSTIC_GEOGRAPHY_COLS))
      stop("DIAGNOSTIC_GEOGRAPHY_COLS must be NULL or a character vector of source column names.")
    return(intersect(unique(DIAGNOSTIC_GEOGRAPHY_COLS), column_names))
  }
  candidates <- c("province", "province_code", "province_name", "prov", "prov_code",
                  "city", "city_code", "city_name", "prefecture", "prefecture_code",
                  "county", "county_code", "county_name", "district", "district_code",
                  "region", "region_code", "region_name", "site", "site_id", "location", "location_id",
                  "\u7701", "\u7701\u4efd", "\u7701\u4ee3\u7801", "\u7701\u7f16\u7801",
                  "\u5e02", "\u57ce\u5e02", "\u5730\u5e02", "\u5e02\u4ee3\u7801", "\u5e02\u7f16\u7801",
                  "\u53bf", "\u533a\u53bf", "\u53bf\u4ee3\u7801", "\u53bf\u7f16\u7801", "\u5730\u533a")
  column_names[tolower(column_names) %in% candidates]
}

# Denominators are complete-case valid matched records before state exclusion.
# Calendar summaries use each record's date; DJF is the winter season.
summarize_overlap_exclusions <- function(dt) {
  geo_cols <- diagnostic_geography_columns(names(dt))
  keep <- unique(c("id", "date", "case", "joint_state", "D_C_raw", "D_H_raw",
                   "N_C_raw", "N_H_raw", geo_cols))
  z <- copy(dt[, ..keep])
  z[, `:=`(day_overlap = D_C_raw == 1L & D_H_raw == 1L,
           night_overlap = N_C_raw == 1L & N_H_raw == 1L,
           directly_excluded = is.na(joint_state))]
  if (anyNA(z$day_overlap) || anyNA(z$night_overlap))
    stop("Overlap diagnostics require complete temperature indicators.")
  sets <- z[, .(n_case_before = sum(case == 1L), n_referents_before = sum(case == 0L),
                n_case_after = sum(case == 1L & !directly_excluded),
                n_referents_after = sum(case == 0L & !directly_excluded)), by = id]
  if (any(sets$n_case_before != 1L | sets$n_referents_before < 1L))
    stop("Exclusion diagnostics must start from valid matched sets.")
  z[sets, on = "id", `:=`(case_absent = i.n_case_after != 1L,
                          referents_absent = i.n_referents_after < 1L)]
  z[, `:=`(indirect_case_absent = !directly_excluded & case_absent,
           indirect_referents_absent = !directly_excluded & !case_absent & referents_absent,
           retained = !directly_excluded & !case_absent & !referents_absent,
           record_role = fifelse(case == 1L, "case", "referent"))]
  if (any(z[, directly_excluded + indirect_case_absent + indirect_referents_absent + retained] != 1L))
    stop("Exclusion categories do not partition the pre-exclusion records.")
  
  count_records <- function(x) {
    list(n_records_before = nrow(x), n_sets_before = uniqueN(x$id),
         n_no_overlap = sum(!x$day_overlap & !x$night_overlap),
         n_day_only_overlap = sum(x$day_overlap & !x$night_overlap),
         n_night_only_overlap = sum(!x$day_overlap & x$night_overlap),
         n_both_periods_overlap = sum(x$day_overlap & x$night_overlap),
         n_any_overlap = sum(x$day_overlap | x$night_overlap),
         n_directly_excluded = sum(x$directly_excluded),
         n_indirectly_excluded_case_absent = sum(x$indirect_case_absent),
         n_indirectly_excluded_all_referents_absent = sum(x$indirect_referents_absent),
         n_excluded_total = sum(!x$retained), n_retained = sum(x$retained),
         n_sets_with_retained_records = uniqueN(x$id[x$retained]))
  }
  # Summarize cases, referents and all records within each partition.
  summarize_partition <- function(values, partition, source_column) {
    values <- as.character(values)
    values[is.na(values) | !nzchar(trimws(values))] <- "(Missing)"
    z[, diagnostic_group := values]
    roles <- z[, count_records(.SD), by = .(group_value = diagnostic_group, record_role)]
    totals <- z[, count_records(.SD), by = .(group_value = diagnostic_group)]
    totals[, record_role := "all_records"]
    out <- rbindlist(list(totals, roles), use.names = TRUE)
    out[, `:=`(group_type = partition, grouping_variable = source_column)]
    out
  }
  month <- as.integer(format(z$date, "%m"))
  seasons <- c("Winter_DJF", "Spring_MAM", "Summer_JJA", "Autumn_SON")
  ans <- list(
    summarize_partition(rep("All", nrow(z)), "overall", "none"),
    summarize_partition(format(z$date, "%Y"), "year", "date"),
    summarize_partition(sprintf("%02d", month), "month", "date"),
    summarize_partition(format(z$date, "%Y-%m"), "year_month", "date"),
    summarize_partition(seasons[(month %% 12L) %/% 3L + 1L], "season", "date"))
  for (field in geo_cols)
    ans[[length(ans) + 1L]] <- summarize_partition(z[[field]], "geography", field)
  ans <- rbindlist(ans, use.names = TRUE)
  rate_counts <- c("n_no_overlap", "n_day_only_overlap", "n_night_only_overlap",
                   "n_both_periods_overlap", "n_any_overlap", "n_directly_excluded",
                   "n_indirectly_excluded_case_absent", "n_indirectly_excluded_all_referents_absent",
                   "n_excluded_total", "n_retained")
  for (field in rate_counts)
    ans[, (paste0(sub("^n_", "", field), "_percent")) := 100 * get(field) / n_records_before]
  ans[, `:=`(outcome = OUTCOME_NAME, outcome_code = OUTCOME_CODE,
             both_anomaly_policy = BOTH_ANOMALY_POLICY,
             denominator_population = "Complete-case valid matched records before joint-state exclusion",
             percentage_denominator = "n_records_before within the same group and record_role",
             geography_columns_used = if (length(geo_cols)) paste(geo_cols, collapse = "|") else "none")]
  setcolorder(ans, c("outcome", "outcome_code", "group_type", "grouping_variable", "group_value", "record_role"))
  setorderv(ans, c("group_type", "grouping_variable", "group_value", "record_role"))
  list(table = ans, valid_ids = sets[n_case_after == 1L & n_referents_after >= 1L, id])
}

# Pair counts represent matched referent records. Set counts are deduplicated
# within each cell, and cannot be summed across cells to obtain unique sets.
summarize_matched_state_support <- function(dt, state_levels = make_state_levels(), ref_state = REF_STATE) {
  z <- dt[, .(id, case, state = as.character(joint_state))]
  if (!nrow(z) || anyNA(z$state) || any(!z$state %in% state_levels))
    stop("Matched-state support requires a nonempty final sample with valid states.")
  sets <- z[, .(n_case = sum(case == 1L), n_referents = sum(case == 0L),
                n_states = uniqueN(state)), by = id]
  if (any(sets$n_case != 1L | sets$n_referents < 1L))
    stop("Matched-state support requires one case and at least one referent per set.")
  cases <- z[case == 1L, .(id, case_state = state)]
  pairs <- merge(z[case == 0L, .(id, referent_state = state)], cases,
                 by = "id", all.x = TRUE, sort = FALSE)
  observed_pairs <- pairs[, .(n_case_referent_pairs = .N, n_matched_sets = uniqueN(id)),
                          by = .(case_state, referent_state)]
  matrix_table <- merge(CJ(case_state = state_levels, referent_state = state_levels, sorted = FALSE),
                        observed_pairs, by = c("case_state", "referent_state"), all.x = TRUE, sort = FALSE)
  matrix_table[is.na(n_matched_sets), `:=`(n_case_referent_pairs = 0L, n_matched_sets = 0L)]
  matrix_table[, discordant_state_pair := case_state != referent_state]
  matrix_table[, `:=`(case_order = match(case_state, state_levels), referent_order = match(referent_state, state_levels))]
  setorder(matrix_table, case_order, referent_order)
  matrix_table[, c("case_order", "referent_order") := NULL]
  
  presence <- unique(z[, .(id, state)])
  ref_ids <- presence[state == ref_state, id]
  observed_states <- unique(presence$state)
  edges <- observed_pairs[case_state != referent_state & n_matched_sets > 0L]
  connected <- intersect(ref_state, observed_states)
  repeat {
    expanded <- union(connected, unique(c(edges[case_state %in% connected, referent_state],
                                          edges[referent_state %in% connected, case_state])))
    if (length(expanded) == length(connected)) break
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
    data.table(state = st, n_cases = sum(cases$case_state == st),
               n_referent_records = sum(pairs$referent_state == st), n_sets_containing_state = length(ids),
               n_sets_case_state_with_other_referent = n_forward,
               n_sets_other_case_with_state_referent = n_reverse,
               n_sets_with_state_contrast = n_informative,
               n_sets_entirely_in_state = n_concordant,
               n_sets_cooccurring_with_reference = if (st == ref_state) NA_integer_ else sum(ids %in% ref_ids),
               connected_to_reference = if (length(ids)) st %in% connected else NA,
               comparison_support = if (!length(ids)) "not_observed" else if (!n_informative) "no_within_set_state_contrast"
               else if (!st %in% connected) "disconnected_from_reference"
               else if (!n_forward || !n_reverse) "one_case_referent_direction_only" else "both_directions_observed")
  }))
  support[, `:=`(n_matched_sets_total = nrow(sets),
                 n_sets_with_any_state_contrast = sum(sets$n_states > 1L),
                 n_sets_without_state_contrast = sum(sets$n_states == 1L),
                 percent_sets_without_state_contrast = 100 * mean(sets$n_states == 1L),
                 interpretation = "Observed comparison support; connectivity alone does not rule out separation or imprecision")]
  for (out in list(support, matrix_table)) {
    out[, `:=`(outcome = OUTCOME_NAME, outcome_code = OUTCOME_CODE,
               population = "Final included valid matched sets", reference_state = ref_state)]
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
  if (length(missing_cols) > 0) stop("Missing required columns for the joint analysis: ", paste(missing_cols, collapse = ", "))
  
  retained_cols <- unique(c(required_cols, diagnostic_geography_columns(names(dt_raw))))
  dt <- copy(dt_raw[, ..retained_cols])
  
  if (is.finite(TEST_N)) {
    log_msg("TEST_N is finite. Keeping complete sets whose ids occur in the first", TEST_N, "rows.")
    test_ids <- unique(dt$id[seq_len(min(TEST_N, nrow(dt)))])
    dt <- dt[id %in% test_ids]
  }
  
  dt[, date := as.IDate(date)]
  dt[, year := as.integer(format(date, "%Y"))]
  dt[, case := coerce_binary01(case, "case")]
  dt[, holiday := coerce_binary01(holiday, "holiday")]
  
  rh_short_cols <- c(sprintf("rh_day_lag%02d", 1:3), sprintf("rh_night_lag%02d", 1:3))
  dt[, rh_lag01_03 := rowMeans(.SD, na.rm = FALSE), .SDcols = rh_short_cols]
  dt[, (rh_short_cols) := NULL]
  
  model_vars_for_complete <- c("id", "case", "holiday", "rh_lag01_03", "year", day_cols, night_cols)
  
  log_msg("Filtering complete cases for the joint model...")
  cc <- dt[, complete.cases(.SD), .SDcols = model_vars_for_complete]
  dt <- dt[cc]
  rm(cc)
  gc()
  
  log_msg("Filtering valid matched sets for the joint model...")
  valid_id <- dt[, .(n_case = sum(case == 1), n_ref = sum(case == 0), n = .N), by = id][n_case == 1 & n_ref >= 1, id]
  dt <- dt[id %in% valid_id]
  rm(valid_id)
  gc()
  
  if (!nrow(dt)) stop("No complete valid matched sets remain for the joint model.")
  log_msg("Joint-analysis rows after filtering:", nrow(dt))
  log_msg("Joint-analysis cases:", dt[, sum(case == 1)])
  log_msg("Joint-analysis matched sets:", dt[, uniqueN(id)])
  
  list(
    dt = dt,
    day_cols = day_cols,
    night_cols = night_cols,
    rh_day_cols = rh_day_cols,
    rh_night_cols = rh_night_cols,
    heat_day_cols = heat_day_cols,
    cold_day_cols = cold_day_cols,
    heat_night_cols = heat_night_cols,
    cold_night_cols = cold_night_cols
  )
}

run_single_temperature_models <- function(dt_raw) {
  
  log_msg("Single-temperature analysis started: single daytime and nighttime DLNM analyses.")
  log_msg("Output root:", OUT_ROOT)
  
  prep <- prepare_single_model_data(dt_raw)
  dt <- prep$dt
  day_cols <- prep$day_cols
  night_cols <- prep$night_cols
  
  sample_desc <- data.table(
    outcome = OUTCOME_NAME,
    n_rows = nrow(dt),
    n_cases = dt[, sum(case == 1)],
    n_referents = dt[, sum(case == 0)],
    n_matched_sets = dt[, uniqueN(id)],
    date_min = as.character(min(dt$date)),
    date_max = as.character(max(dt$date)),
    lag_column_suffix_min = min(LAG_DAYS),
    lag_column_suffix_max = max(LAG_DAYS),
    lag_index_min = min(LAG_INDEX),
    lag_index_max = max(LAG_INDEX)
  )
  save_dt(sample_desc, file.path(DIR_TABLE, "table_00_analysis_sample.csv"))
  
  log_msg("Building exposure matrices and spline specifications...")
  day_mat <- as.matrix(dt[, ..day_cols])
  night_mat <- as.matrix(dt[, ..night_cols])
  temperature_history_cols <- c(day_cols, night_cols)
  dt[, (temperature_history_cols) := NULL]
  rm(temperature_history_cols)
  gc()
  
  day_spec <- make_ns_spec(as.vector(day_mat), VAR_DF)
  night_spec <- make_ns_spec(as.vector(night_mat), VAR_DF)
  lag_spec <- make_lag_spec_logknots(
    lag_index = LAG_INDEX,
    n_internal_knots = LOG_LAG_N_INTERNAL_KNOTS,
    intercept = LOG_LAG_INTERCEPT
  )
  
  method_suffix <- "logknots"
  method_label <- "Log-knots lag basis: natural spline with log-spaced lag knots and intercept"
  
  basis_specs <- list(
    method_suffix = method_suffix,
    method_label = method_label,
    day_spec = day_spec,
    night_spec = night_spec,
    lag_spec = lag_spec,
    VAR_DF = VAR_DF,
    LAG_DF = LAG_DF,
    RH_DF = RH_DF,
    LAG_DAYS = LAG_DAYS,
    LAG_INDEX = LAG_INDEX,
    lag_note = "LAG_DAYS is used for column suffixes; LAG_INDEX is used for modeling and plotting lag days."
  )
  save_obj(basis_specs, file.path(DIR_MODEL, "basis_specs_single_logknots.qs"))
  
  log_msg("Fitting daytime single-index DLNM: logknots")
  x_day <- build_cb_matrix(day_mat, day_spec, lag_spec, "day_single")
  fit_day <- fit_clogit_from_design(dt, x_day, "day_single_dlnm_logknots")
  save_obj(fit_day, file.path(DIR_MODEL, "model_day_single_dlnm_logknots.qs"))
  rm(x_day); gc()
  
  log_msg("Fitting nighttime single-index DLNM: logknots")
  x_night <- build_cb_matrix(night_mat, night_spec, lag_spec, "night_single")
  fit_night <- fit_clogit_from_design(dt, x_night, "night_single_dlnm_logknots")
  save_obj(fit_night, file.path(DIR_MODEL, "model_night_single_dlnm_logknots.qs"))
  rm(x_night); gc()
  
  analysis_list <- list(
    list(
      exposure = "daytime",
      exposure_label = "Daytime",
      prefix = "day_single",
      fit = fit_day,
      mat = day_mat,
      spec = day_spec,
      colour = DAY_COLOUR,
      xlab = "Daytime temperature (°C)"
    ),
    list(
      exposure = "nighttime",
      exposure_label = "Nighttime",
      prefix = "night_single",
      fit = fit_night,
      mat = night_mat,
      spec = night_spec,
      colour = NIGHT_COLOUR,
      xlab = "Nighttime temperature (°C)"
    )
  )
  
  curve_all <- list()
  lag_all <- list()
  rr_summary_all <- list()
  mmt_summary_all <- list()
  
  for (aa in analysis_list) {
    log_msg("Producing single-temperature results for:", aa$exposure)
    
    temp_vec <- as.vector(aa$mat)
    temp_range <- quantile(temp_vec, MMT_RANGE_PROBS, na.rm = TRUE)
    temp_seq <- sort(unique(c(seq(temp_range[1], temp_range[2], by = CURVE_BY), unname(temp_range[2]))))
    
    mmt <- find_mmt_from_single_model(
      fit_obj = aa$fit,
      temp_seq = temp_seq,
      var_spec = aa$spec,
      lag_spec = lag_spec,
      prefix = aa$prefix,
      temp_for_median = temp_vec
    )
    
    temp_seq <- sort(unique(c(temp_seq, mmt)))
    mmt_ci <- simulate_mmt_ci_from_single_model(
      fit_obj = aa$fit,
      temp_seq = temp_seq,
      var_spec = aa$spec,
      lag_spec = lag_spec,
      prefix = aa$prefix,
      n_sim = N_MMT_CI_SIM,
      seed = MMT_CI_SEED + ifelse(aa$exposure == "daytime", 1L, 2L)
    )
    
    p025 <- as.numeric(quantile(temp_vec, 0.025, na.rm = TRUE, type = 8))
    p975 <- as.numeric(quantile(temp_vec, 0.975, na.rm = TRUE, type = 8))
    
    curve_dt <- make_cumulative_curve(aa$fit, temp_seq, mmt, aa$spec, lag_spec, aa$prefix)
    curve_dt[, `:=`(
      outcome = OUTCOME_NAME,
      model = "single_index",
      lag_method = method_suffix,
      lag_method_label = method_label,
      exposure = aa$exposure,
      exposure_label = aa$exposure_label,
      mmt = mmt,
      mmt_low = mmt_ci$mmt_low,
      mmt_high = mmt_ci$mmt_high,
      mmt_n_sim = mmt_ci$n_sim,
      p025 = p025,
      p975 = p975
    )]
    setcolorder(curve_dt, c(
      "outcome", "model", "lag_method", "lag_method_label",
      "exposure", "exposure_label", "temperature", "mmt", "mmt_low", "mmt_high", "mmt_n_sim", "p025", "p975",
      "rr", "rr_low", "rr_high", "logrr", "se"
    ))
    curve_all[[aa$exposure]] <- curve_dt
    
    plot_cumulative_curve(curve_dt, aa, method_suffix)
    
    lag_cold <- make_lag_curve(aa$fit, p025, mmt, aa$spec, lag_spec, aa$prefix)
    lag_cold[, `:=`(
      outcome = OUTCOME_NAME,
      model = "single_index",
      lag_method = method_suffix,
      lag_method_label = method_label,
      exposure = aa$exposure,
      exposure_label = aa$exposure_label,
      contrast = "P2.5_vs_MMT",
      contrast_label = "P2.5 vs MMT",
      temperature = p025,
      reference_temperature = mmt
    )]
    
    lag_heat <- make_lag_curve(aa$fit, p975, mmt, aa$spec, lag_spec, aa$prefix)
    lag_heat[, `:=`(
      outcome = OUTCOME_NAME,
      model = "single_index",
      lag_method = method_suffix,
      lag_method_label = method_label,
      exposure = aa$exposure,
      exposure_label = aa$exposure_label,
      contrast = "P97.5_vs_MMT",
      contrast_label = "P97.5 vs MMT",
      temperature = p975,
      reference_temperature = mmt
    )]
    
    lag_dt <- rbindlist(list(lag_cold, lag_heat), use.names = TRUE)
    setcolorder(lag_dt, c(
      "outcome", "model", "lag_method", "lag_method_label",
      "exposure", "exposure_label", "contrast", "contrast_label", "lag_day",
      "temperature", "reference_temperature", "rr", "rr_low", "rr_high", "logrr", "se"
    ))
    lag_all[[aa$exposure]] <- lag_dt
    
    plot_lag_curve(lag_dt, aa, method_suffix)
    
    x_ext <- build_cb_constant(c(p025, p975), aa$spec, lag_spec, aa$prefix)
    x_ref <- build_cb_constant(mmt, aa$spec, lag_spec, aa$prefix)
    rr_sum <- predict_contrast_from_design(aa$fit, x_ext, x_ref)
    rr_sum[, `:=`(
      outcome = OUTCOME_NAME,
      model = "single_index",
      lag_method = method_suffix,
      lag_method_label = method_label,
      exposure = aa$exposure,
      exposure_label = aa$exposure_label,
      contrast = c("P2.5_vs_MMT", "P97.5_vs_MMT"),
      contrast_label = c("P2.5 vs MMT", "P97.5 vs MMT"),
      temperature = c(p025, p975),
      reference_temperature = mmt,
      mmt = mmt,
      mmt_low = mmt_ci$mmt_low,
      mmt_high = mmt_ci$mmt_high,
      mmt_n_sim = mmt_ci$n_sim
    )]
    setcolorder(rr_sum, c(
      "outcome", "model", "lag_method", "lag_method_label",
      "exposure", "exposure_label", "contrast", "contrast_label",
      "temperature", "reference_temperature", "mmt",
      "rr", "rr_low", "rr_high", "logrr", "se"
    ))
    rr_summary_all[[aa$exposure]] <- rr_sum
    
    mmt_summary_all[[aa$exposure]] <- data.table(
      outcome = OUTCOME_NAME,
      model = "single_index",
      lag_method = method_suffix,
      lag_method_label = method_label,
      exposure = aa$exposure,
      exposure_label = aa$exposure_label,
      mmt = mmt,
      mmt_low = mmt_ci$mmt_low,
      mmt_high = mmt_ci$mmt_high,
      mmt_n_sim = mmt_ci$n_sim,
      mmt_percentile = mean(temp_vec <= mmt, na.rm = TRUE) * 100,
      p025 = p025,
      p975 = p975,
      var_df = VAR_DF,
      lag_df_original_setting = LAG_DF,
      lag_spec_type = lag_spec$type,
      lag_internal_knots = paste(round(lag_spec$knots, 6), collapse = ";"),
      lag_intercept = lag_spec$intercept,
      rh_df = RH_DF,
      lag_column_suffix_min = min(LAG_DAYS),
      lag_column_suffix_max = max(LAG_DAYS),
      lag_index_min = min(LAG_INDEX),
      lag_index_max = max(LAG_INDEX)
    )
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
  
  analysis_meta <- list(
    outcome = OUTCOME_NAME,
    data_path = DATA_PATH,
    output_root = OUT_ROOT,
    method_suffix = method_suffix,
    method_label = method_label,
    n_rows_model = nrow(dt),
    n_cases = dt[, sum(case == 1)],
    settings = list(
      VAR_DF = VAR_DF,
      LAG_DF = LAG_DF,
      RH_DF = RH_DF,
      LAG_DAYS = LAG_DAYS,
      LAG_INDEX = LAG_INDEX,
      LOG_LAG_N_INTERNAL_KNOTS = LOG_LAG_N_INTERNAL_KNOTS,
      LOG_LAG_INTERCEPT = LOG_LAG_INTERCEPT,
      N_MMT_CI_SIM = N_MMT_CI_SIM,
      lag_spec = lag_spec,
      lag_note = "LAG_DAYS is used for column suffixes; LAG_INDEX is used for modeling and plotting lag days."
    ),
    mmt_summary = mmt_summary_dt,
    cumulative_rr_summary = rr_summary_dt
  )
  save_obj(analysis_meta, file.path(DIR_MODEL, "analysis_metadata_single_logknots.qs"))
  
  result <- list(
    curve_all_dt = curve_all_dt,
    lag_all_dt = lag_all_dt,
    rr_summary_dt = rr_summary_dt,
    mmt_summary_dt = mmt_summary_dt,
    fit_day = fit_day,
    fit_night = fit_night,
    basis_specs = basis_specs,
    day_spec = day_spec,
    night_spec = night_spec,
    lag_spec = lag_spec
  )
  save_obj(result, file.path(DIR_MODEL, "analysis_results_single_logknots.qs"))
  
  log_msg("Single daytime and nighttime analyses finished successfully.")
  log_msg("Single-exposure cumulative and lag-response outputs saved.")
  invisible(result)
}

format_lag_window <- function(lag_days) {
  lag_days <- sort(unique(as.integer(lag_days)))
  if (length(lag_days) == 0L) return(NA_character_)
  if (length(lag_days) == 1L) return(paste0("lag", lag_days))
  paste0("lag", min(lag_days), "-lag", max(lag_days))
}

format_lag_window_from_suffix <- function(lag_suffix) {
  format_lag_window(as.integer(lag_suffix) - 1L)
}

validate_fallback_window <- function(x, name) {
  ok <- is.numeric(x) && length(x) > 0L && !anyNA(x) &&
    all(is.finite(x)) && all(x == as.integer(x)) &&
    all(x %in% LAG_INDEX) &&
    identical(as.integer(x), seq.int(min(x), max(x)))
  if (!ok) {
    stop(name, " must be ordered consecutive integer model lags within 0:20.")
  }
  invisible(TRUE)
}

select_one_lag_window <- function(lag_dt,
                                  exposure_value,
                                  contrast_value,
                                  component,
                                  thermal,
                                  alpha = AUTO_LAG_ALPHA,
                                  selection_rule = AUTO_LAG_SELECTION_RULE,
                                  no_sig_action = AUTO_LAG_NO_SIG_ACTION,
                                  fallback_heat = AUTO_LAG_FALLBACK_HEAT,
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
  
  required_cols <- c(
    "exposure", "contrast", "lag_day", "rr", "rr_low", "rr_high",
    "logrr", "se", "temperature", "reference_temperature"
  )
  missing_cols <- setdiff(required_cols, names(lag_dt))
  if (length(missing_cols) > 0L) {
    stop("Lag-response results are missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  x <- copy(lag_dt[exposure == exposure_value & contrast == contrast_value])
  setorder(x, lag_day)
  
  if (nrow(x) != length(LAG_INDEX) || anyDuplicated(x$lag_day) ||
      !identical(as.integer(x$lag_day), as.integer(LAG_INDEX))) {
    stop(
      "Expected exactly ", length(LAG_INDEX), " lag estimates (",
      min(LAG_INDEX), ":", max(LAG_INDEX), ") for ", component, "."
    )
  }
  
  if (any(!is.finite(x$rr)) || any(!is.finite(x$rr_low)) ||
      any(!is.finite(x$rr_high)) || any(!is.finite(x$logrr)) ||
      any(!is.finite(x$se)) || any(x$rr_low <= 0) ||
      any(x$rr_low > x$rr) || any(x$rr > x$rr_high)) {
    stop("Non-finite or invalid lag estimates/confidence limits for ", component, ".")
  }
  
  if (uniqueN(x$temperature) != 1L || uniqueN(x$reference_temperature) != 1L ||
      any(!is.finite(x$temperature)) || any(!is.finite(x$reference_temperature))) {
    stop("Invalid temperature contrast for ", component, ".")
  }
  if (thermal == "cold" && x$temperature[1] >= x$reference_temperature[1]) {
    stop("Cold contrast temperature is not below MMT for ", component, ".")
  }
  if (thermal == "heat" && x$temperature[1] <= x$reference_temperature[1]) {
    stop("Heat contrast temperature is not above MMT for ", component, ".")
  }
  
  zcrit <- qnorm(1 - alpha / 2)
  x[, selection_ci_low := exp(logrr - zcrit * se)]
  x[, selection_ci_high := exp(logrr + zcrit * se)]
  x[, significant_harmful := selection_ci_low > 1]
  sig_lags <- as.integer(x[significant_harmful == TRUE, lag_day])
  
  fallback_used <- length(sig_lags) == 0L
  if (!fallback_used) {
    selected_lags <- seq.int(0L, max(sig_lags))
    selection_reason <- "zero_to_last_significant_harmful_lag"
  } else {
    if (identical(no_sig_action, "stop")) {
      stop(
        "No harmful significant lag (two-sided CI lower bound > 1) was found for ",
        component,
        ". Single-temperature results have been saved; the joint model was not run."
      )
    }
    selected_lags <- if (thermal == "heat") fallback_heat else fallback_cold
    selected_lags <- as.integer(selected_lags)
    selection_reason <- paste0("prespecified_fallback_", thermal, "_",
                               format_lag_window(selected_lags), "_no_significant_harmful_lag")
  }
  
  selected_lags <- as.integer(selected_lags)
  selected_suffix <- selected_lags + 1L
  selected_flag <- x$lag_day %in% selected_lags
  selected_non_significant <- selected_flag & !x$significant_harmful
  
  summary <- data.table(
    outcome = OUTCOME_NAME,
    exposure = exposure_value,
    thermal = thermal,
    contrast = contrast_value,
    component = component,
    significance_rule = paste0(round((1 - alpha) * 100, 1), "% CI lower bound > 1"),
    lag_selection_rule = "lag0-to-last-significant continuous window",
    status = if (fallback_used) "prespecified_fallback" else "auto_selected",
    n_significant_harmful_lags = length(sig_lags),
    significant_harmful_lags = if (length(sig_lags)) paste(sig_lags, collapse = ",") else "",
    selected_lag_start = min(selected_lags),
    selected_lag_end = max(selected_lags),
    selected_lag_days = paste(selected_lags, collapse = ","),
    selected_data_suffix = paste(selected_suffix, collapse = ","),
    selected_window = format_lag_window(selected_lags),
    n_selected_lags = length(selected_lags),
    n_selected_non_significant = sum(selected_non_significant),
    significant_at_lag20 = any(sig_lags == max(LAG_INDEX)),
    fallback_used = fallback_used,
    selection_reason = selection_reason,
    temperature = x$temperature[1],
    reference_temperature = x$reference_temperature[1]
  )
  
  diagnostics <- copy(x)
  diagnostics[, `:=`(
    outcome = OUTCOME_NAME,
    thermal = thermal,
    component = component,
    selected = lag_day %in% selected_lags,
    selected_non_significant = (lag_day %in% selected_lags) & !significant_harmful,
    fallback_used = fallback_used,
    selection_reason = selection_reason
  )]
  
  list(
    lag_days = selected_lags,
    data_suffix = selected_suffix,
    summary = summary,
    diagnostics = diagnostics
  )
}

apply_auto_lag_windows <- function(auto_lag) {
  AUTO_LAG_SELECTION_SUMMARY <<- copy(auto_lag$summary)
  DAY_COLD_CLASS_LAG_DAYS <<- auto_lag$day_cold$data_suffix
  DAY_HEAT_CLASS_LAG_DAYS <<- auto_lag$day_heat$data_suffix
  NIGHT_COLD_CLASS_LAG_DAYS <<- auto_lag$night_cold$data_suffix
  NIGHT_HEAT_CLASS_LAG_DAYS <<- auto_lag$night_heat$data_suffix
  
  all_suffix <- c(
    DAY_COLD_CLASS_LAG_DAYS,
    DAY_HEAT_CLASS_LAG_DAYS,
    NIGHT_COLD_CLASS_LAG_DAYS,
    NIGHT_HEAT_CLASS_LAG_DAYS
  )
  if (any(!all_suffix %in% LAG_DAYS)) {
    stop("At least one selected lag falls outside the available lag columns.")
  }
  
  invisible(TRUE)
}

make_state_levels <- function() {
  c(
    "D0_N0", "DC_N0", "DH_N0",
    "D0_NC", "DC_NC", "DH_NC",
    "D0_NH", "DC_NH", "DH_NH"
  )
}

state_dt_from_levels <- function(state_levels) {
  dt_state <- data.table(state = state_levels)
  dt_state[, day_code := sub("_.*$", "", state)]
  dt_state[, night_code := sub("^.*_", "", state)]
  
  dt_state[, D_C := as.integer(day_code == "DC")]
  dt_state[, D_H := as.integer(day_code == "DH")]
  dt_state[, N_C := as.integer(night_code == "NC")]
  dt_state[, N_H := as.integer(night_code == "NH")]
  
  dt_state[, day_status := fifelse(
    day_code == "D0", "No daytime anomaly",
    fifelse(day_code == "DC", "Daytime cold", "Daytime heat")
  )]
  
  dt_state[, night_status := fifelse(
    night_code == "N0", "No nighttime anomaly",
    fifelse(night_code == "NC", "Nighttime cold", "Nighttime heat")
  )]
  
  dt_state[, day_axis := factor(day_status, levels = c(
    "No daytime anomaly", "Daytime cold", "Daytime heat"
  ))]
  
  dt_state[, night_axis := factor(night_status, levels = c(
    "No nighttime anomaly", "Nighttime cold", "Nighttime heat"
  ))]
  
  dt_state[, day_order := match(day_code, c("D0", "DC", "DH")) - 1L]
  dt_state[, night_order := match(night_code, c("N0", "NC", "NH")) - 1L]
  dt_state[, mask := day_order + 3L * night_order]
  dt_state[, component_mask := D_C * 1L + D_H * 2L + N_C * 4L + N_H * 8L]
  dt_state[, c("day_order", "night_order") := NULL]
  
  dt_state[]
}

COMPONENTS <- c("D_C", "D_H", "N_C", "N_H")
COMPONENT_BITS <- c(D_C = 1L, D_H = 2L, N_C = 4L, N_H = 8L)

COMPONENT_INFO <- data.table(
  component = c(COMPONENTS, "joint_total"),
  component_order = 1:5,
  component_label = c(
    "Daytime cold",
    "Daytime heat",
    "Nighttime cold",
    "Nighttime heat",
    "Joint total"
  ),
  component_group = c(
    "daytime_cold",
    "daytime_heat",
    "nighttime_cold",
    "nighttime_heat",
    "joint_total"
  )
)

state_from_mask <- function(mask) {
  day_levels <- c("D0", "DC", "DH")
  night_levels <- c("N0", "NC", "NH")
  day_i <- mask %% 3L
  night_i <- mask %/% 3L
  paste0(day_levels[day_i + 1L], "_", night_levels[night_i + 1L])
}

mask_from_daynight_codes <- function(day_code, night_code) {
  day_i <- match(day_code, c("D0", "DC", "DH")) - 1L
  night_i <- match(night_code, c("N0", "NC", "NH")) - 1L
  as.integer(day_i + 3L * night_i)
}

resolve_period_status <- function(cold, heat, period_name, policy = BOTH_ANOMALY_POLICY) {
  code <- rep(NA_character_, length(cold))
  code[!cold & !heat] <- if (period_name == "day") "D0" else "N0"
  code[cold & !heat] <- if (period_name == "day") "DC" else "NC"
  code[!cold & heat] <- if (period_name == "day") "DH" else "NH"
  
  both <- cold & heat
  if (any(both, na.rm = TRUE)) {
    if (policy == "exclude") {
      code[both] <- NA_character_
    } else if (policy == "heat_priority") {
      code[both] <- if (period_name == "day") "DH" else "NH"
    } else if (policy == "cold_priority") {
      code[both] <- if (period_name == "day") "DC" else "NC"
    } else {
      stop("Unknown BOTH_ANOMALY_POLICY: ", policy)
    }
  }
  
  code
}

build_joint_indicators <- function(dt,
                                   day_cold_threshold,
                                   day_heat_threshold,
                                   night_cold_threshold,
                                   night_heat_threshold) {
  threshold_values <- c(
    day_cold_threshold, day_heat_threshold,
    night_cold_threshold, night_heat_threshold
  )
  if (any(!is.finite(threshold_values))) {
    stop("Joint-state reference thresholds must all be finite.")
  }
  
  out <- copy(dt)
  
  out[, D_C_raw := as.integer(day_cold_temp < day_cold_threshold)]
  out[, D_H_raw := as.integer(day_heat_temp > day_heat_threshold)]
  out[, N_C_raw := as.integer(night_cold_temp < night_cold_threshold)]
  out[, N_H_raw := as.integer(night_heat_temp > night_heat_threshold)]
  
  out[, day_code := resolve_period_status(D_C_raw == 1, D_H_raw == 1, "day", BOTH_ANOMALY_POLICY)]
  out[, night_code := resolve_period_status(N_C_raw == 1, N_H_raw == 1, "night", BOTH_ANOMALY_POLICY)]
  
  out[, D_C := as.integer(day_code == "DC")]
  out[, D_H := as.integer(day_code == "DH")]
  out[, N_C := as.integer(night_code == "NC")]
  out[, N_H := as.integer(night_code == "NH")]
  
  out[, joint_state := fifelse(
    is.na(day_code) | is.na(night_code),
    NA_character_,
    paste0(day_code, "_", night_code)
  )]
  out[, joint_mask := mask_from_daynight_codes(day_code, night_code)]
  
  attr(out, "state_thresholds") <- list(
    strategy = STATE_THRESHOLD_STRATEGY,
    strategy_tag = reference_strategy_tag(),
    reference_method = CURRENT_SHARED_DESIGN$settings$reference_method,
    reference_source_outcome = CURRENT_SHARED_DESIGN$source_name,
    reference_source_code = CURRENT_SHARED_DESIGN$source_code,
    calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
    rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance,
    day_cold_threshold = day_cold_threshold,
    day_heat_threshold = day_heat_threshold,
    night_cold_threshold = night_cold_threshold,
    night_heat_threshold = night_heat_threshold
  )
  
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
  out[, beta := beta_state[state]]
  out[, se := state_se[state]]
  
  out[, rr := fifelse(is.finite(beta), exp(beta), NA_real_)]
  out[, rr_low := fifelse(is.finite(beta) & is.finite(se), exp(beta - 1.96 * se), NA_real_)]
  out[, rr_high := fifelse(is.finite(beta) & is.finite(se), exp(beta + 1.96 * se), NA_real_)]
  
  out[state == REF_STATE, `:=`(
    beta = 0,
    se = 0,
    rr = 1,
    rr_low = 1,
    rr_high = 1
  )]
  
  out[, estimable := is.finite(beta) & is.finite(se)]
  out[, rr_label := fifelse(
    state == REF_STATE,
    "1.00\nReference",
    fifelse(
      estimable,
      sprintf("%.2f\n(%.2f, %.2f)", rr, rr_low, rr_high),
      "NE\nNot observed"
    )
  )]
  
  out[]
}

# ============================================================
# Nine-state risk estimation and annual death counts
# ============================================================

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
  
  desc_basic <- data.table(
    outcome = OUTCOME_NAME,
    model_label = MODEL_LABEL,
    n_rows = nrow(dt),
    n_cases = dt[, sum(case == 1)],
    n_referents = dt[, sum(case == 0)],
    n_matched_sets = dt[, uniqueN(id)],
    date_min = as.character(min(dt$date)),
    date_max = as.character(max(dt$date)),
    day_heat_lag_window = format_lag_window_from_suffix(DAY_HEAT_CLASS_LAG_DAYS),
    day_cold_lag_window = format_lag_window_from_suffix(DAY_COLD_CLASS_LAG_DAYS),
    night_heat_lag_window = format_lag_window_from_suffix(NIGHT_HEAT_CLASS_LAG_DAYS),
    night_cold_lag_window = format_lag_window_from_suffix(NIGHT_COLD_CLASS_LAG_DAYS),
    state_threshold_strategy = STATE_THRESHOLD_STRATEGY,
    threshold_strategy_tag = reference_strategy_tag(),
    reference_method = CURRENT_SHARED_DESIGN$settings$reference_method,
    reference_source_outcome = CURRENT_SHARED_DESIGN$source_name,
    reference_source_code = CURRENT_SHARED_DESIGN$source_code,
    calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
    rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance,
    both_anomaly_policy = BOTH_ANOMALY_POLICY,
    mmt_lag_knot_type = LAG_KNOT_TYPE,
    mmt_lag_n_internal_knots = LAG_N_INTERNAL_KNOTS,
    mmt_lag_intercept = LAG_INTERCEPT,
    n_sim = N_SIM,
    cap_negative_af = CAP_NEGATIVE_AF
  )
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
  
  mmt_summary <- data.table(
    outcome = OUTCOME_NAME,
    model_label = MODEL_LABEL,
    day_mmt = mmt_day,
    day_mmt_low = mmt_day_low,
    day_mmt_high = mmt_day_high,
    night_mmt = mmt_night,
    night_mmt_low = mmt_night_low,
    night_mmt_high = mmt_night_high,
    state_threshold_strategy = STATE_THRESHOLD_STRATEGY,
    threshold_strategy_tag = reference_strategy_tag(),
    reference_method = CURRENT_SHARED_DESIGN$settings$reference_method,
    reference_source_outcome = CURRENT_SHARED_DESIGN$source_name,
    reference_source_code = CURRENT_SHARED_DESIGN$source_code,
    calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
    rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance,
    var_df = VAR_DF,
    rh_df = RH_DF,
    lag_knot_type = LAG_KNOT_TYPE,
    lag_n_internal_knots = LAG_N_INTERNAL_KNOTS,
    lag_intercept = LAG_INTERCEPT,
    lag_knots = paste(round(lag_spec$knots, 4), collapse = ","),
    lag_boundary_knots = paste(round(lag_spec$Boundary.knots, 4), collapse = ","),
    mmt_source = paste("Overall-outcome single-temperature log-knots DLNM:", CURRENT_SHARED_DESIGN$source_name)
  )
  save_dt(mmt_summary, file.path(DIR_TABLE, "table_01_mmt_summary.csv"))
  save_obj(mmt_summary, file.path(DIR_MODEL, "mmt_summary.qs"))
  
  daytime_curve_for_mmt <- copy(single_logknots_results$curve_all_dt[exposure == "daytime"])
  nighttime_curve_for_mmt <- copy(single_logknots_results$curve_all_dt[exposure == "nighttime"])
  daytime_curve_for_mmt[, exposure := "daytime_temperature"]
  nighttime_curve_for_mmt[, exposure := "nighttime_temperature"]
  save_dt(daytime_curve_for_mmt, file.path(DIR_PLOT_DATA, "plotdata_daytime_temperature_single_curve_for_mmt.csv"))
  save_obj(daytime_curve_for_mmt, file.path(DIR_PLOT_DATA, "plotdata_daytime_temperature_single_curve_for_mmt.qs"))
  save_dt(nighttime_curve_for_mmt, file.path(DIR_PLOT_DATA, "plotdata_nighttime_temperature_single_curve_for_mmt.csv"))
  save_obj(nighttime_curve_for_mmt, file.path(DIR_PLOT_DATA, "plotdata_nighttime_temperature_single_curve_for_mmt.qs"))
  
  log_msg("Daytime MMT:", round(mmt_day, 3), "95% CI:", round(mmt_day_low, 3), "to", round(mmt_day_high, 3))
  log_msg("Nighttime MMT:", round(mmt_night, 3), "95% CI:", round(mmt_night_low, 3), "to", round(mmt_night_high, 3))
  log_msg("Joint-state threshold strategy:", STATE_THRESHOLD_STRATEGY)
  log_msg("Reference source:", CURRENT_SHARED_DESIGN$source_name)
  log_msg("Reference definition:", reference_strategy_tag())
  
  log_msg("Building heat and cold lag-window summaries...")
  dt[, day_heat_temp := rowMeans(.SD, na.rm = FALSE), .SDcols = heat_day_cols]
  dt[, day_cold_temp := rowMeans(.SD, na.rm = FALSE), .SDcols = cold_day_cols]
  dt[, night_heat_temp := rowMeans(.SD, na.rm = FALSE), .SDcols = heat_night_cols]
  dt[, night_cold_temp := rowMeans(.SD, na.rm = FALSE), .SDcols = cold_night_cols]
  temperature_history_cols <- c(day_cols, night_cols)
  dt[, (temperature_history_cols) := NULL]
  rm(temperature_history_cols)
  gc()
  
  STATE_LEVELS <- make_state_levels()
  STATE_INFO <- state_dt_from_levels(STATE_LEVELS)
  
  state_thresholds <- compute_joint_state_thresholds()
  
  log_msg(
    "Daytime cold threshold:", round(state_thresholds$day_cold_threshold, 3),
    "| lower shared reference bound; applied to",
    format_lag_window_from_suffix(DAY_COLD_CLASS_LAG_DAYS), "mean temperature"
  )
  log_msg(
    "Daytime heat threshold:", round(state_thresholds$day_heat_threshold, 3),
    "| upper shared reference bound; applied to",
    format_lag_window_from_suffix(DAY_HEAT_CLASS_LAG_DAYS), "mean temperature"
  )
  log_msg(
    "Nighttime cold threshold:", round(state_thresholds$night_cold_threshold, 3),
    "| lower shared reference bound; applied to",
    format_lag_window_from_suffix(NIGHT_COLD_CLASS_LAG_DAYS), "mean temperature"
  )
  log_msg(
    "Nighttime heat threshold:", round(state_thresholds$night_heat_threshold, 3),
    "| upper shared reference bound; applied to",
    format_lag_window_from_suffix(NIGHT_HEAT_CLASS_LAG_DAYS), "mean temperature"
  )
  
  dt <- build_joint_indicators(
    dt,
    day_cold_threshold = state_thresholds$day_cold_threshold,
    day_heat_threshold = state_thresholds$day_heat_threshold,
    night_cold_threshold = state_thresholds$night_cold_threshold,
    night_heat_threshold = state_thresholds$night_heat_threshold
  )
  threshold_attr <- attr(dt, "state_thresholds")
  
  n_before_joint_filter <- nrow(dt)
  n_case_before_joint_filter <- dt[, sum(case == 1)]
  n_excluded_both_anomaly <- dt[, sum(is.na(joint_state))]
  
  overlap_diagnostics <- data.table(
    outcome = OUTCOME_NAME,
    threshold_strategy_tag = state_thresholds$strategy_tag,
    reference_method = CURRENT_SHARED_DESIGN$settings$reference_method,
    reference_source_outcome = CURRENT_SHARED_DESIGN$source_name,
    reference_source_code = CURRENT_SHARED_DESIGN$source_code,
    calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
    rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance,
    n_rows = nrow(dt),
    n_cases = dt[, sum(case == 1)],
    day_cold_heat_overlap_rows = dt[, sum(D_C_raw == 1L & D_H_raw == 1L)],
    day_cold_heat_overlap_percent = dt[, mean(D_C_raw == 1L & D_H_raw == 1L) * 100],
    day_cold_heat_overlap_cases = dt[case == 1L, sum(D_C_raw == 1L & D_H_raw == 1L)],
    day_cold_heat_overlap_case_percent = dt[case == 1L, mean(D_C_raw == 1L & D_H_raw == 1L) * 100],
    night_cold_heat_overlap_rows = dt[, sum(N_C_raw == 1L & N_H_raw == 1L)],
    night_cold_heat_overlap_percent = dt[, mean(N_C_raw == 1L & N_H_raw == 1L) * 100],
    night_cold_heat_overlap_cases = dt[case == 1L, sum(N_C_raw == 1L & N_H_raw == 1L)],
    night_cold_heat_overlap_case_percent = dt[case == 1L, mean(N_C_raw == 1L & N_H_raw == 1L) * 100],
    excluded_joint_state_rows = n_excluded_both_anomaly
  )
  save_dt(overlap_diagnostics, file.path(DIR_TABLE, "table_00a_cold_heat_overlap_diagnostics.csv"))
  
  exclusion_details <- summarize_overlap_exclusions(dt)
  exclusion_details$table[, `:=`(calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
                                 reference_method = CURRENT_SHARED_DESIGN$settings$reference_method)]
  save_dt(exclusion_details$table, file.path(DIR_TABLE, "table_00d_overlap_exclusion_distribution.csv"))
  log_msg("Exclusion diagnostic geography fields:", exclusion_details$table$geography_columns_used[1L])
  
  ref_diagnostics <- dt[joint_state == REF_STATE, .(
    n_ref_state_rows = .N,
    n_ref_state_cases = sum(case == 1),
    n_ref_state_referents = sum(case == 0),
    n_ref_state_sets = uniqueN(id)
  )]
  
  reference_diagnostics <- data.table(
    state_threshold_strategy = STATE_THRESHOLD_STRATEGY,
    threshold_strategy_tag = state_thresholds$strategy_tag,
    reference_method = CURRENT_SHARED_DESIGN$settings$reference_method,
    reference_source_outcome = CURRENT_SHARED_DESIGN$source_name,
    reference_source_code = CURRENT_SHARED_DESIGN$source_code,
    calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
    rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance,
    day_cold_threshold = state_thresholds$day_cold_threshold,
    day_heat_threshold = state_thresholds$day_heat_threshold,
    night_cold_threshold = state_thresholds$night_cold_threshold,
    night_heat_threshold = state_thresholds$night_heat_threshold,
    n_rows_after_joint_state_filter = n_before_joint_filter - n_excluded_both_anomaly,
    n_excluded_both_anomaly = n_excluded_both_anomaly,
    n_ref_state_rows = ref_diagnostics$n_ref_state_rows,
    n_ref_state_cases = ref_diagnostics$n_ref_state_cases,
    n_ref_state_referents = ref_diagnostics$n_ref_state_referents,
    n_ref_state_sets = ref_diagnostics$n_ref_state_sets
  )
  save_dt(reference_diagnostics, file.path(DIR_TABLE, "table_00_reference_state_diagnostics.csv"))
  
  dt <- dt[!is.na(joint_state)]
  
  valid_id2 <- exclusion_details$valid_ids
  dt <- dt[id %in% valid_id2]
  rm(valid_id2, exclusion_details)
  gc()
  
  log_msg("Rows before removing within-period cold+heat records:", n_before_joint_filter)
  log_msg("Cases before removing within-period cold+heat records:", n_case_before_joint_filter)
  if (!nrow(dt) || !any(dt$case == 1L)) stop("No valid matched sets remain after joint-state exclusion.")
  save_dt(data.table(
    outcome = OUTCOME_NAME, source_reference = CURRENT_SHARED_DESIGN$source_name,
    n_rows_before_joint_filter = n_before_joint_filter,
    n_cases_before_joint_filter = n_case_before_joint_filter,
    n_rows_removed_for_cold_heat_overlap = n_excluded_both_anomaly,
    n_rows_final = nrow(dt), n_cases_final = sum(dt$case == 1L),
    n_cases_lost_total = n_case_before_joint_filter - sum(dt$case == 1L),
    case_retention_percent = 100 * sum(dt$case == 1L) / n_case_before_joint_filter,
    burden_denominator = "Final included cases in valid matched sets"
  ), file.path(DIR_TABLE, "table_00b_joint_analysis_sample_flow.csv"))
  log_msg("Rows after joint-state filter and matched-set refilter:", nrow(dt))
  log_msg("Cases after joint-state filter and matched-set refilter:", dt[, sum(case == 1)])
  
  dt[, joint_state := factor(joint_state, levels = STATE_LEVELS)]
  dt[, joint_state := relevel(joint_state, ref = REF_STATE)]
  
  threshold_summary <- data.table(
    outcome = OUTCOME_NAME,
    model_label = MODEL_LABEL,
    both_anomaly_policy = BOTH_ANOMALY_POLICY,
    state_threshold_strategy = STATE_THRESHOLD_STRATEGY,
    threshold_strategy_tag = state_thresholds$strategy_tag,
    reference_method = CURRENT_SHARED_DESIGN$settings$reference_method,
    reference_source_outcome = CURRENT_SHARED_DESIGN$source_name,
    reference_source_code = CURRENT_SHARED_DESIGN$source_code,
    calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
    rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance,
    day_cold_threshold = state_thresholds$day_cold_threshold,
    day_heat_threshold = state_thresholds$day_heat_threshold,
    night_cold_threshold = state_thresholds$night_cold_threshold,
    night_heat_threshold = state_thresholds$night_heat_threshold,
    threshold_distribution = state_thresholds$threshold_distribution,
    day_mmt = mmt_day,
    day_mmt_low = mmt_day_low,
    day_mmt_high = mmt_day_high,
    night_mmt = mmt_night,
    night_mmt_low = mmt_night_low,
    night_mmt_high = mmt_night_high,
    day_heat_lag_days_data_suffix = paste(DAY_HEAT_CLASS_LAG_DAYS, collapse = ","),
    day_cold_lag_days_data_suffix = paste(DAY_COLD_CLASS_LAG_DAYS, collapse = ","),
    night_heat_lag_days_data_suffix = paste(NIGHT_HEAT_CLASS_LAG_DAYS, collapse = ","),
    night_cold_lag_days_data_suffix = paste(NIGHT_COLD_CLASS_LAG_DAYS, collapse = ","),
    day_heat_lag_window = format_lag_window_from_suffix(DAY_HEAT_CLASS_LAG_DAYS),
    day_cold_lag_window = format_lag_window_from_suffix(DAY_COLD_CLASS_LAG_DAYS),
    night_heat_lag_window = format_lag_window_from_suffix(NIGHT_HEAT_CLASS_LAG_DAYS),
    night_cold_lag_window = format_lag_window_from_suffix(NIGHT_COLD_CLASS_LAG_DAYS),
    mmt_lag_knot_type = LAG_KNOT_TYPE,
    mmt_lag_n_internal_knots = LAG_N_INTERNAL_KNOTS,
    mmt_lag_intercept = LAG_INTERCEPT,
    mmt_lag_knots = paste(round(lag_spec$knots, 4), collapse = ","),
    mmt_lag_boundary_knots = paste(round(lag_spec$Boundary.knots, 4), collapse = ",")
  )
  save_dt(threshold_summary, file.path(DIR_TABLE, "table_02_joint_state_thresholds.csv"))
  save_obj(threshold_summary, file.path(DIR_MODEL, "joint_state_thresholds.qs"))
  
  state_counts <- dt[, .(
    n_rows = .N,
    n_cases = sum(case == 1),
    n_referents = sum(case == 0),
    n_sets = uniqueN(id)
  ), by = .(joint_state)]
  setnames(state_counts, "joint_state", "state")
  state_counts[, state := as.character(state)]
  state_counts <- merge(STATE_INFO, state_counts, by = "state", all.x = TRUE, sort = FALSE)
  state_counts[is.na(n_rows), `:=`(n_rows = 0L, n_cases = 0L, n_referents = 0L, n_sets = 0L)]
  
  state_counts[, support_status := fifelse(
    n_rows == 0L,
    "not_observed",
    fifelse(
      n_cases == 0L | n_referents == 0L,
      "observed_but_not_estimable",
      "observed_estimable"
    )
  )]
  setorder(state_counts, mask)
  save_dt(state_counts, file.path(DIR_TABLE, "table_03_joint_state_counts.csv"))
  
  matched_support <- summarize_matched_state_support(dt, STATE_LEVELS, REF_STATE)
  save_dt(matched_support$summary, file.path(DIR_TABLE, "table_03c_matched_state_comparison_support.csv"))
  save_dt(matched_support$comparisons, file.path(DIR_TABLE, "table_03d_case_referent_state_comparisons.csv"))
  log_msg("Matched sets with a within-set state contrast:",
          matched_support$summary$n_sets_with_any_state_contrast[1L], "of",
          matched_support$summary$n_matched_sets_total[1L])
  limited_support <- matched_support$summary[
    !comparison_support %in% c("not_observed", "both_directions_observed"), state]
  if (length(limited_support)) log_msg("Inspect within-set support diagnostics for:", paste(limited_support, collapse = ", "))
  rm(matched_support, limited_support)
  
  log_msg("Joint-state counts before model fitting:")
  for (ii in seq_len(nrow(state_counts))) {
    log_msg(
      "  ", state_counts$state[ii],
      "| rows:", state_counts$n_rows[ii],
      "| cases:", state_counts$n_cases[ii],
      "| referents:", state_counts$n_referents[ii],
      "| sets:", state_counts$n_sets[ii],
      "| support:", state_counts$support_status[ii]
    )
  }
  
  if (state_counts[state == REF_STATE, n_rows] == 0L) {
    stop("The reference state D0_N0 is not observed and the joint model cannot be fitted.")
  }
  
  unobserved_states <- state_counts[support_status == "not_observed", state]
  unsupported_observed_states <- state_counts[
    support_status == "observed_but_not_estimable",
    state
  ]
  
  if (length(unsupported_observed_states) > 0L) {
    stop(
      "The following states contain observations but have no cases or no referents and cannot be safely omitted: ",
      paste(unsupported_observed_states, collapse = ", "),
      ". Inspect table_03_joint_state_counts.csv."
    )
  }
  
  observed_state_levels <- STATE_LEVELS[!STATE_LEVELS %in% unobserved_states]
  estimable_non_ref_states <- setdiff(observed_state_levels, REF_STATE)
  
  if (length(unobserved_states) > 0L) {
    log_msg(
      "Structurally unobserved joint states will remain in reporting tables as NE but will not be included in the conditional logistic regression:",
      paste(unobserved_states, collapse = ", ")
    )
  }
  log_msg(
    "Joint model will estimate",
    length(estimable_non_ref_states),
    "non-reference state coefficients; global Wald df =",
    length(estimable_non_ref_states)
  )
  
  log_msg("Fitting joint categorical model:", MODEL_LABEL)
  joint_model_dt <- data.table(
    case = dt$case,
    id = dt$id,
    joint_state = factor(as.character(dt$joint_state), levels = observed_state_levels),
    holiday = dt$holiday,
    rh_lag01_03 = dt$rh_lag01_03
  )
  joint_model_dt[, joint_state := relevel(joint_state, ref = REF_STATE)]
  
  rh_basis <- as.data.table(ns(joint_model_dt$rh_lag01_03, df = RH_DF))
  setnames(rh_basis, paste0("rh_ns", seq_len(ncol(rh_basis))))
  joint_model_dt <- cbind(joint_model_dt, rh_basis)
  
  rhs_terms <- c("joint_state", names(rh_basis), "holiday")
  form_joint <- as.formula(paste0("case ~ ", paste(rhs_terms, collapse = " + "), " + strata(id)"))
  formula_joint_text <- paste(deparse(form_joint), collapse = " ")
  
  fit_joint <- clogit(
    form_joint,
    data = joint_model_dt,
    method = CLOGIT_METHOD,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  
  coef_joint <- coef(fit_joint)
  vcov_joint <- vcov(fit_joint)
  
  required_joint_terms <- paste0("joint_state", estimable_non_ref_states)
  names(required_joint_terms) <- estimable_non_ref_states
  missing_joint_terms <- setdiff(required_joint_terms, names(coef_joint))
  available_joint_terms <- intersect(required_joint_terms, names(coef_joint))
  nonfinite_joint_terms <- available_joint_terms[!is.finite(coef_joint[available_joint_terms])]
  
  if (length(missing_joint_terms) > 0L || length(nonfinite_joint_terms) > 0L) {
    failed_joint_compact <- list(
      formula_text = formula_joint_text,
      model_label = MODEL_LABEL,
      state_levels = STATE_LEVELS,
      observed_state_levels = observed_state_levels,
      estimable_non_reference_states = estimable_non_ref_states,
      unobserved_states = unobserved_states,
      ref_state = REF_STATE,
      coefficients = coef_joint,
      covariance = vcov_joint,
      missing_joint_terms = missing_joint_terms,
      nonfinite_joint_terms = nonfinite_joint_terms,
      state_counts = state_counts,
      threshold_summary = threshold_summary,
      storage_mode = "compact_failed_joint_model_diagnostics"
    )
    save_obj(
      failed_joint_compact,
      file.path(DIR_MODEL, paste0("model_joint_", MODEL_LABEL, "_FAILED_COMPACT.qs"))
    )
    
    rm(fit_joint, joint_model_dt, rh_basis, form_joint)
    invisible(gc(full = TRUE))
    
    stop(
      "At least one OBSERVED joint state was not estimable after removing structurally unobserved states. Missing terms: ",
      paste(missing_joint_terms, collapse = ", "),
      "; non-finite terms: ",
      paste(nonfinite_joint_terms, collapse = ", "),
      ". This indicates separation or insufficient within-stratum support rather than a structurally empty state."
    )
  }
  
  beta_info <- make_beta_state_vector(
    coefficients = coef_joint,
    covariance = vcov_joint,
    state_levels = STATE_LEVELS,
    ref_state = REF_STATE
  )
  beta_state <- beta_info$beta_state
  state_se <- beta_info$state_se
  
  covariance_state <- matrix(
    NA_real_,
    nrow = length(STATE_LEVELS),
    ncol = length(STATE_LEVELS),
    dimnames = list(STATE_LEVELS, STATE_LEVELS)
  )
  covariance_state[REF_STATE, REF_STATE] <- 0
  
  if (length(estimable_non_ref_states) > 0L) {
    estimable_terms <- required_joint_terms[estimable_non_ref_states]
    covariance_state[
      estimable_non_ref_states,
      estimable_non_ref_states
    ] <- vcov_joint[
      estimable_terms,
      estimable_terms,
      drop = FALSE
    ]
  }
  
  mu_non_reference <- setNames(
    as.numeric(coef_joint[required_joint_terms]),
    estimable_non_ref_states
  )
  covariance_non_reference <- vcov_joint[
    required_joint_terms,
    required_joint_terms,
    drop = FALSE
  ]
  dimnames(covariance_non_reference) <- list(
    estimable_non_ref_states,
    estimable_non_ref_states
  )
  
  joint_model_compact <- list(
    coefficients = coef_joint,
    covariance = vcov_joint,
    beta_state = beta_state,
    covariance_state = covariance_state,
    mu_non_reference = mu_non_reference,
    covariance_non_reference = covariance_non_reference,
    theoretical_state_levels = STATE_LEVELS,
    state_levels = STATE_LEVELS,
    observed_state_levels = observed_state_levels,
    estimable_non_reference_states = estimable_non_ref_states,
    unobserved_states = unobserved_states,
    n_estimable_non_reference = length(estimable_non_ref_states),
    state_terms = required_joint_terms,
    formula_text = formula_joint_text,
    model_label = MODEL_LABEL,
    both_anomaly_policy = BOTH_ANOMALY_POLICY,
    ref_state = REF_STATE,
    state_info = STATE_INFO,
    state_counts = state_counts,
    threshold_summary = threshold_summary,
    state_threshold_strategy = STATE_THRESHOLD_STRATEGY,
    threshold_strategy_tag = state_thresholds$strategy_tag,
    reference_method = CURRENT_SHARED_DESIGN$settings$reference_method,
    reference_source_outcome = CURRENT_SHARED_DESIGN$source_name,
    reference_source_code = CURRENT_SHARED_DESIGN$source_code,
    calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
    rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance,
    state_thresholds = state_thresholds,
    mmt_lag_spec = lag_spec,
    storage_mode = "compact_coefficients_and_covariance_only_dynamic_state_support"
  )
  
  save_obj(
    joint_model_compact,
    file.path(DIR_MODEL, paste0("model_joint_", MODEL_LABEL, ".qs"))
  )
  
  rm(fit_joint, joint_model_dt, rh_basis, form_joint)
  invisible(gc(full = TRUE))
  log_msg("Finished joint categorical model; retained compact coefficients/covariance only.")
  
  required_states <- unique(as.character(dt[case == 1, joint_state]))
  missing_required <- required_states[!is.finite(beta_state[required_states])]
  if (length(missing_required) > 0) {
    stop(
      "Some states required for burden calculation are not estimable in the joint model: ",
      paste(unique(missing_required), collapse = ", "),
      ". Check sparse states or the selected shared reference thresholds."
    )
  }
  
  rr_table <- make_rr_table(beta_state, state_se, STATE_INFO)
  rr_table <- merge(rr_table, state_counts[, .(state, n_rows, n_cases, n_referents, n_sets, support_status)], by = "state", all.x = TRUE, sort = FALSE)
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
    data.table(year = as.integer(y), state = STATE_LEVELS,
               n_cases = as.numeric(counts_by_year[[y]][STATE_LEVELS]),
               n_deaths = as.numeric(n_deaths_by_year[[y]]))
  }))
  save_dt(annual_counts, file.path(DIR_TABLE, "table_00c_annual_state_death_counts.csv"))
  export_covariance_csv(covariance_non_reference, file.path(DIR_TABLE, "table_04a_joint_logrr_covariance.csv"))
  analysis_meta <- list(
    outcome = OUTCOME_NAME,
    data_path = DATA_PATH,
    output_root = OUT_ROOT,
    model_label = MODEL_LABEL,
    both_anomaly_policy = BOTH_ANOMALY_POLICY,
    n_rows_model = nrow(dt),
    n_cases = dt[, sum(case == 1)],
    n_matched_sets = dt[, uniqueN(id)],
    mmt_summary = mmt_summary,
    threshold_summary = threshold_summary,
    state_levels = STATE_LEVELS,
    observed_state_levels = observed_state_levels,
    estimable_non_reference_states = estimable_non_ref_states,
    unobserved_states = unobserved_states,
    n_estimable_non_reference = length(estimable_non_ref_states),
    reference_state = REF_STATE,
    state_threshold_strategy = STATE_THRESHOLD_STRATEGY,
    threshold_strategy_tag = state_thresholds$strategy_tag,
    reference_method = CURRENT_SHARED_DESIGN$settings$reference_method,
    reference_source_outcome = CURRENT_SHARED_DESIGN$source_name,
    reference_source_code = CURRENT_SHARED_DESIGN$source_code,
    calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
    rr_tolerance = CURRENT_SHARED_DESIGN$thresholds$rr_tolerance,
    state_thresholds = state_thresholds,
    overlap_diagnostics = overlap_diagnostics,
    lag_response_setting_for_mmt = list(
      knot_type = LAG_KNOT_TYPE,
      n_internal_knots = LAG_N_INTERNAL_KNOTS,
      intercept = LAG_INTERCEPT,
      knots = lag_spec$knots,
      boundary_knots = lag_spec$Boundary.knots,
      note = "Natural cubic spline with intercept and three internal knots placed at equally spaced values on the log scale."
    ),
    day_heat_class_lag_days = DAY_HEAT_CLASS_LAG_DAYS,
    day_cold_class_lag_days = DAY_COLD_CLASS_LAG_DAYS,
    night_heat_class_lag_days = NIGHT_HEAT_CLASS_LAG_DAYS,
    night_cold_class_lag_days = NIGHT_COLD_CLASS_LAG_DAYS,
    auto_lag_selection = list(
      alpha = AUTO_LAG_ALPHA,
      selection_rule = AUTO_LAG_SELECTION_RULE,
      no_significant_lag_action = AUTO_LAG_NO_SIG_ACTION,
      fallback_heat_model_lags = AUTO_LAG_FALLBACK_HEAT,
      fallback_cold_model_lags = AUTO_LAG_FALLBACK_COLD,
      selection_summary = AUTO_LAG_SELECTION_SUMMARY
    ),
    n_sim = N_SIM,
    sim_seed = SIM_SEED,
    cap_negative_af = CAP_NEGATIVE_AF,
    shared_design = CURRENT_SHARED_DESIGN,
    uncertainty_scope = "Joint-model coefficient Monte Carlo; conditional on shared thresholds and lag selection; case counts treated as fixed",
    mean_annual_AF_definition = "Arithmetic mean of annual AF",
    burden_population = "Final included cases after complete-case, state-overlap and matched-set filtering",
    analysis_stage = "risk_estimation",
    n_years = length(years)
  )
  ready <- write_stage1_input_bundle(rr_table, annual_counts, case_dt, joint_model_compact, analysis_meta)
  log_msg("Stage-one risk estimation and input preparation completed:", OUT_ROOT)
  invisible(list(rr_table = rr_table, joint_model_compact = joint_model_compact,
                 annual_counts = annual_counts, threshold_summary = threshold_summary,
                 analysis_meta = analysis_meta, ready = ready))
}

# ============================================================
# Shared reference intervals and exposure windows
# ============================================================

validate_scalar_choice <- function(x, allowed, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% allowed)
    stop(label, " must be one of: ", paste(allowed, collapse = ", "))
}

validate_ma_days <- function(x, required_names, allow_partial = FALSE) {
  if (allow_partial && is.null(x)) return(invisible(TRUE))
  if (!is.numeric(x) || !length(x) || is.null(names(x)) ||
      anyNA(x) || any(!is.finite(x)) || any(x != floor(x)) ||
      any(x < 1L | x > length(LAG_INDEX)) || anyDuplicated(names(x)) ||
      any(!names(x) %in% required_names) ||
      (!allow_partial && !setequal(names(x), required_names))) {
    stop("Moving-average days must be named integers between 1 and ",
         length(LAG_INDEX), "; required names: ", paste(required_names, collapse = ", "))
  }
  invisible(TRUE)
}

validate_run_settings <- function(reference_method, rr_tolerance, lag_mode, manual_days,
                                  overrides, no_sig_action, fallback_days, allow_truncated) {
  validate_scalar_choice(reference_method, c("mmt_ci", "rr_tolerance"), "reference_method")
  validate_scalar_choice(lag_mode, c("auto", "manual"), "lag_mode")
  validate_scalar_choice(no_sig_action, c("fallback", "stop"), "no_sig_action")
  validate_scalar_choice(BOTH_ANOMALY_POLICY, c("exclude", "heat_priority", "cold_priority"), "BOTH_ANOMALY_POLICY")
  if (!is.numeric(rr_tolerance) || length(rr_tolerance) != 1L || !is.finite(rr_tolerance) || rr_tolerance <= 0)
    stop("rr_tolerance must be positive, e.g. 0.01.")
  if (!is.logical(allow_truncated) || length(allow_truncated) != 1L || is.na(allow_truncated))
    stop("allow_truncated_rr_band must be TRUE or FALSE.")
  keys <- c("day_cold", "day_heat", "night_cold", "night_heat")
  validate_ma_days(overrides, keys, TRUE)
  if (lag_mode == "manual") validate_ma_days(manual_days, keys)
  validate_ma_days(fallback_days, c("cold", "heat"))
  for (v in list(N_SIM, N_MMT_CI_SIM))
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 2 || v != floor(v))
      stop("N_SIM and N_MMT_CI_SIM must be integers >= 2.")
  if (!is.numeric(TEST_N) || length(TEST_N) != 1L || is.na(TEST_N) || TEST_N <= 0 ||
      (is.finite(TEST_N) && TEST_N != floor(TEST_N))) stop("TEST_N must be a positive integer or Inf.")
  invisible(TRUE)
}

cumulative_design_fast <- function(x, var_spec, lag_spec, prefix) {
  xb <- predict_ns_basis(x, var_spec)
  lb <- colSums(predict_ns_basis(LAG_INDEX, lag_spec))
  out <- do.call(cbind, lapply(seq_len(ncol(xb)), function(a) outer(xb[, a], lb)))
  colnames(out) <- unlist(lapply(seq_len(ncol(xb)), function(a)
    paste0(prefix, "_v", a, "_l", seq_along(lb))))
  out
}

refine_mmt <- function(temp_seq, eta_grid, eta_function) {
  n <- length(temp_seq)
  if (n < 3L || any(!is.finite(eta_grid))) stop("Invalid MMT search grid.")
  mids <- 2L:(n - 1L)
  candidates <- mids[eta_grid[mids] <= eta_grid[mids - 1L] &
                       eta_grid[mids] <= eta_grid[mids + 1L]]
  roots <- vapply(candidates, function(i)
    optimize(eta_function, c(temp_seq[i - 1L], temp_seq[i + 1L]), tol = 1e-7)$minimum,
    numeric(1))
  candidates_t <- unique(c(temp_seq[c(1L, n)], roots))
  values <- vapply(candidates_t, eta_function, numeric(1))
  candidates_t[which.min(values)]
}

reference_strategy_tag <- function() {
  if (is.null(CURRENT_SHARED_DESIGN)) stop("Shared design has not been initialized.")
  CURRENT_SHARED_DESIGN$reference_tag
}

compute_joint_state_thresholds <- function(dt = NULL) {
  if (is.null(CURRENT_SHARED_DESIGN)) stop("Shared design has not been initialized.")
  CURRENT_SHARED_DESIGN$thresholds
}

read_analysis_data <- function(paths) {
  if (!is.character(paths) || !length(paths) || anyNA(paths) ||
      any(!file.exists(paths))) stop("Input RDS file(s) missing: ", paste(paths[!file.exists(paths)], collapse = " | "))
  required <- c("id", "date", "case", "holiday",
                sprintf("temp_day_lag%02d", LAG_DAYS),
                sprintf("temp_night_lag%02d", LAG_DAYS),
                sprintf("rh_day_lag%02d", 1:3), sprintf("rh_night_lag%02d", 1:3))
  pieces <- vector("list", length(paths))
  next_id <- 0
  audit <- vector("list", length(paths))
  for (i in seq_along(paths)) {
    log_msg("Reading data:", paths[i])
    z <- as.data.table(readRDS(paths[i]))
    missing <- setdiff(required, names(z))
    if (length(missing)) stop("Missing columns in ", paths[i], ": ", paste(missing, collapse = ", "))
    geo_cols <- diagnostic_geography_columns(names(z))
    if (!is.null(DIAGNOSTIC_GEOGRAPHY_COLS)) {
      unavailable_geo <- setdiff(DIAGNOSTIC_GEOGRAPHY_COLS, names(z))
      if (length(unavailable_geo)) log_msg("Configured geography fields absent from this input:", paste(unavailable_geo, collapse = ", "))
    }
    remove_cols <- setdiff(names(z), c(required, geo_cols))
    if (length(remove_cols)) z[, (remove_cols) := NULL]
    for (field in geo_cols) set(z, j = field, value = as.character(z[[field]]))
    if (anyNA(z$id)) stop("Missing matched-set ids in ", paths[i])
    old_ids <- unique(z$id)
    z[, id := as.double(match(id, old_ids)) + next_id]
    next_id <- next_id + length(old_ids)
    z[, case := coerce_binary01(case, "case")]
    z[, holiday := coerce_binary01(holiday, "holiday")]
    z[, date := as.IDate(date)]
    if (anyNA(z$case) || anyNA(z$date)) stop("Unrecognized case code or missing/invalid date in ", paths[i])
    temp_rh <- setdiff(required, c("id", "date", "case", "holiday"))
    for (v in temp_rh) {
      if (!is.numeric(z[[v]])) stop("Expected numeric temperature/humidity column: ", v)
      if (any(is.infinite(z[[v]]))) stop("Infinite exposure values in ", v)
    }
    if (anyDuplicated(z, by = c("id", "date")))
      stop("Duplicate id/date rows in ", paths[i], "; check upstream matched-set construction.")
    audit[[i]] <- data.table(source_file = paths[i], n_rows = nrow(z),
                             n_cases = sum(z$case == 1L), n_sets = length(old_ids),
                             id_recode = "source-specific ids renumbered before row binding")
    pieces[[i]] <- z
    rm(z, old_ids)
    invisible(gc())
  }
  ans <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  audit_dt <- rbindlist(audit)
  stopifnot(nrow(ans) == sum(audit_dt$n_rows))
  save_dt(audit_dt, file.path(DIR_TABLE, "table_00_input_source_audit.csv"))
  rm(pieces)
  invisible(gc())
  ans
}

calibration_signature <- function(paths) {
  info <- file.info(paths)
  if (anyNA(info$size)) stop("Overall input files unavailable; cannot validate calibration cache.")
  list(version = "shared-design-20260925-v1", source_code = OVERALL_CODE,
       source_name = OVERALL_NAME,
       paths = unname(normalizePath(paths, winslash = "/", mustWork = TRUE)),
       sizes = unname(info$size), mtimes = unname(as.numeric(info$mtime)),
       settings = list(LAG_DAYS = LAG_DAYS, LAG_INDEX = LAG_INDEX, VAR_DF = VAR_DF,
                       RH_DF = RH_DF, LOG_LAG_N_INTERNAL_KNOTS = LOG_LAG_N_INTERNAL_KNOTS,
                       LOG_LAG_INTERCEPT = LOG_LAG_INTERCEPT, MMT_RANGE_PROBS = MMT_RANGE_PROBS,
                       CURVE_BY = CURVE_BY, N_MMT_CI_SIM = N_MMT_CI_SIM,
                       MMT_CI_SEED = MMT_CI_SEED, CLOGIT_METHOD = CLOGIT_METHOD, TEST_N = TEST_N))
}

get_shared_calibration <- function(overall_paths = unname(DATA_PATHS),
                                   shared_dir = SHARED_DIR, force = FALSE) {
  cache_file <- object_path(file.path(shared_dir, "overall_dlnm_calibration.qs"))
  signature <- calibration_signature(overall_paths)
  if (file.exists(cache_file) && !isTRUE(force)) {
    obj <- read_model_object(cache_file)
    if (identical(obj$signature, signature)) {
      obj$cache_file <- cache_file
      return(obj)
    }
  }
  runtime_names <- c("DATA_PATH", "OUT_BASE", "OUT_ROOT", "OUTCOME_CODE", "OUTCOME_NAME",
                     "DIR_TABLE", "DIR_MODEL", "DIR_PLOT_DATA", "DIR_FIG", "DIR_LOG", "DIR_CONTRIB", "LOG_FILE")
  env <- environment(get_shared_calibration)
  existed <- vapply(runtime_names, exists, logical(1), envir = env, inherits = FALSE)
  saved <- mget(runtime_names[existed], envir = env, inherits = FALSE)
  on.exit({
    for (nm in names(saved)) assign(nm, saved[[nm]], envir = env)
    for (nm in runtime_names[!existed]) if (exists(nm, envir = env, inherits = FALSE)) rm(list = nm, envir = env)
  }, add = TRUE)
  OUTCOME_CODE <<- OVERALL_CODE
  OUTCOME_NAME <<- OVERALL_NAME
  DATA_PATH <<- overall_paths
  OUT_BASE <<- shared_dir
  OUT_ROOT <<- file.path(shared_dir, "overall_single_dlnm")
  set_output_dirs(OUT_ROOT, "overall_calibration_log.txt")
  log_msg("Building overall DLNM calibration; thresholds/windows will be reused by all outcomes.")
  raw <- read_analysis_data(overall_paths)
  single <- run_single_temperature_models(raw)
  rm(raw)
  invisible(gc())
  obj <- list(signature = signature, source_code = OVERALL_CODE, source_name = OVERALL_NAME,
              created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
              calibration_id = paste0(OVERALL_CODE, "_", format(Sys.time(), "%Y%m%dT%H%M%OS6")),
              single_results = single, output_root = OUT_ROOT, cache_file = cache_file)
  temp_file <- tempfile("calibration_", tmpdir = shared_dir, fileext = paste0(".", OBJECT_FORMAT))
  on.exit(unlink(temp_file), add = TRUE)
  write_model_object(obj, temp_file)
  if (!file.copy(temp_file, cache_file, overwrite = TRUE)) stop("Could not save shared calibration: ", cache_file)
  log_msg("Shared calibration saved:", cache_file)
  obj
}

rr_tolerance_band <- function(curve, mmt, epsilon, allow_truncated = FALSE) {
  x <- copy(curve[, .(temperature, logrr)])
  setorder(x, temperature)
  if (anyDuplicated(x$temperature) || any(!is.finite(unlist(x)))) stop("Invalid cumulative curve.")
  k <- which.min(abs(x$temperature - mmt))
  if (abs(x$temperature[k] - mmt) > 1e-6) stop("MMT must be present in the prediction grid.")
  target <- log1p(epsilon)
  inside <- x$logrr <= target + 1e-12
  if (!inside[k]) stop("MMT is outside RR-tolerance band.")
  left <- right <- k
  while (left > 1L && inside[left - 1L]) left <- left - 1L
  while (right < nrow(x) && inside[right + 1L]) right <- right + 1L
  cross <- function(a, b) {
    x$temperature[a] + (target - x$logrr[a]) *
      (x$temperature[b] - x$temperature[a]) / (x$logrr[b] - x$logrr[a])
  }
  low <- if (left == 1L) x$temperature[1L] else cross(left - 1L, left)
  high <- if (right == nrow(x)) x$temperature[nrow(x)] else cross(right, right + 1L)
  truncated <- c(left = left == 1L, right = right == nrow(x))
  if (any(truncated) && !allow_truncated) {
    stop("RR-tolerance interval reaches the MMT search-domain boundary. Inspect the overall curve; ",
         "RUN_ALLOW_TRUNCATED_RR_BAND = TRUE explicitly permits a truncated band.")
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
    if (nrow(m) != 1L) stop("Expected one overall MMT per period.")
    band <- if (method == "mmt_ci") {
      list(low = m$mmt_low, high = m$mmt_high,
           truncated_left = m$mmt_low <= min(curve$temperature) + 1e-7,
           truncated_right = m$mmt_high >= max(curve$temperature) - 1e-7)
    } else rr_tolerance_band(curve, m$mmt, epsilon, allow_truncated)
    if (any(!is.finite(c(band$low, band$high))) || band$low >= band$high)
      stop("Empty or invalid overall reference band for ", period, "; inspect the overall MMT/curve.")
    data.table(source_outcome = calibration$source_name, source_code = calibration$source_code,
               calibration_id = calibration$calibration_id, exposure = period, method = method,
               rr_tolerance = if (method == "rr_tolerance") epsilon else NA_real_,
               mmt = m$mmt, mmt_low = m$mmt_low, mmt_high = m$mmt_high,
               reference_low = band$low, reference_high = band$high,
               contains_point_mmt = band$low <= m$mmt && band$high >= m$mmt,
               lower_search_limit = min(curve$temperature), upper_search_limit = max(curve$temperature),
               truncated_left = as.logical(band$truncated_left), truncated_right = as.logical(band$truncated_right),
               definition = if (method == "mmt_ci") "Operational band from overall MMT 95% CI" else
                 "Connected set containing overall MMT with fitted cumulative RR <= 1 + epsilon")
  })
  rbindlist(bands)
}

resolve_shared_lags <- function(single_results, lag_mode, manual_days, overrides, no_sig_action, fallback_days) {
  keys <- c("day_cold", "day_heat", "night_cold", "night_heat")
  validate_scalar_choice(lag_mode, c("auto", "manual"), "lag_mode")
  validate_scalar_choice(no_sig_action, c("fallback", "stop"), "no_sig_action")
  validate_ma_days(fallback_days, c("cold", "heat"))
  validate_ma_days(overrides, keys, allow_partial = TRUE)
  if (lag_mode == "manual") validate_ma_days(manual_days, keys)
  specs <- list(day_cold = c("daytime", "P2.5_vs_MMT", "Daytime cold", "cold"),
                day_heat = c("daytime", "P97.5_vs_MMT", "Daytime heat", "heat"),
                night_cold = c("nighttime", "P2.5_vs_MMT", "Nighttime cold", "cold"),
                night_heat = c("nighttime", "P97.5_vs_MMT", "Nighttime heat", "heat"))
  ans <- setNames(vector("list", length(keys)), keys)
  for (key in keys) {
    z <- specs[[key]]
    forced <- if (lag_mode == "manual") manual_days[[key]] else NULL
    if (key %in% names(overrides)) forced <- overrides[[key]]
    if (is.null(forced)) {
      a <- select_one_lag_window(single_results$lag_all_dt, z[1L], z[2L], z[3L], z[4L],
                                 alpha = AUTO_LAG_ALPHA, selection_rule = "zero_to_last",
                                 no_sig_action = no_sig_action,
                                 fallback_heat = seq_len(fallback_days[["heat"]]) - 1L,
                                 fallback_cold = seq_len(fallback_days[["cold"]]) - 1L)
    } else {
      d <- copy(single_results$lag_all_dt[exposure == z[1L] & contrast == z[2L]])
      setorder(d, lag_day)
      if (nrow(d) != length(LAG_INDEX)) stop("Missing overall lag diagnostics for ", key)
      d[, selection_ci_low := exp(logrr - qnorm(1 - AUTO_LAG_ALPHA / 2) * se)]
      d[, selection_ci_high := exp(logrr + qnorm(1 - AUTO_LAG_ALPHA / 2) * se)]
      d[, significant_harmful := selection_ci_low > 1]
      lags <- seq_len(forced) - 1L
      reason <- if (lag_mode == "manual") "manual_moving_average_days" else "manual_component_override"
      d[, `:=`(thermal = z[4L], component = z[3L], selected = lag_day %in% lags,
               selected_non_significant = (lag_day %in% lags) & !significant_harmful,
               fallback_used = FALSE, selection_reason = reason)]
      sig <- d[significant_harmful == TRUE, lag_day]
      s <- data.table(outcome = OVERALL_NAME, exposure = z[1L], thermal = z[4L], contrast = z[2L],
                      component = z[3L], significance_rule = "Not used for manual selection",
                      lag_selection_rule = "manual lag0 to lag(days - 1)", status = "manual",
                      n_significant_harmful_lags = length(sig), significant_harmful_lags = paste(sig, collapse = ","),
                      selected_lag_start = 0L, selected_lag_end = max(lags), selected_lag_days = paste(lags, collapse = ","),
                      selected_data_suffix = paste(lags + 1L, collapse = ","), selected_window = format_lag_window(lags),
                      n_selected_lags = length(lags), n_selected_non_significant = sum(d$selected_non_significant),
                      significant_at_lag20 = max(LAG_INDEX) %in% sig, fallback_used = FALSE,
                      selection_reason = reason, temperature = d$temperature[1L], reference_temperature = d$reference_temperature[1L])
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

make_shared_design <- function(calibration, reference_method, rr_tolerance, lag_mode,
                               manual_days, manual_overrides, no_sig_action, fallback_days,
                               allow_truncated) {
  bands <- make_reference_bands(calibration, reference_method, rr_tolerance, allow_truncated)
  lags <- resolve_shared_lags(calibration$single_results, lag_mode, manual_days,
                              manual_overrides, no_sig_action, fallback_days)
  ref_tag <- if (reference_method == "mmt_ci") "MMT_CI95" else
    paste0("RRtol_", gsub("[.]", "p", format(rr_tolerance, scientific = FALSE, trim = TRUE)))
  keys <- c("day_cold", "day_heat", "night_cold", "night_heat")
  lag_sizes <- vapply(lags[keys], function(x) length(x$lag_days), integer(1))
  selection_modes <- vapply(lags[keys], function(x) x$summary$status[1L], character(1))
  lag_tag <- paste0(lag_mode, "_", paste(paste0(c("DC", "DH", "NC", "NH"), lag_sizes), collapse = "_"))
  if (!is.null(manual_overrides)) lag_tag <- paste0(lag_tag, "_override_", paste(names(manual_overrides), collapse = "-"))
  if (any(selection_modes == "prespecified_fallback")) lag_tag <- paste0(lag_tag, "_fallback")
  day <- bands[exposure == "daytime"]
  night <- bands[exposure == "nighttime"]
  thresholds <- list(strategy = paste0("overall_", reference_method), strategy_tag = ref_tag,
                     reference_method = reference_method, rr_tolerance = if (reference_method == "rr_tolerance") rr_tolerance else NA_real_,
                     source_outcome = calibration$source_name, source_code = calibration$source_code,
                     calibration_id = calibration$calibration_id,
                     cold_quantile_prob = NA_real_, heat_quantile_prob = NA_real_,
                     day_cold_threshold = day$reference_low, day_heat_threshold = day$reference_high,
                     night_cold_threshold = night$reference_low, night_heat_threshold = night$reference_high,
                     threshold_distribution = "Overall-outcome full-lag cumulative DLNM; same thresholds applied to the four selected moving-average exposures")
  settings <- list(reference_method = reference_method, rr_tolerance = rr_tolerance,
                   lag_mode = lag_mode, manual_days = if (lag_mode == "manual") manual_days else NULL,
                   manual_overrides = manual_overrides, no_sig_action = no_sig_action,
                   fallback_days = fallback_days, alpha = AUTO_LAG_ALPHA,
                   allow_truncated = allow_truncated, overlap_policy = BOTH_ANOMALY_POLICY,
                   n_sim = N_SIM, sim_seed = SIM_SEED, cap_negative_af = CAP_NEGATIVE_AF)
  list(reference_tag = ref_tag, lag_tag = lag_tag,
       analysis_tag = paste(ref_tag, lag_tag, BOTH_ANOMALY_POLICY, sep = "__"),
       bands = bands, lags = lags, thresholds = thresholds, settings = settings,
       calibration_id = calibration$calibration_id, calibration_signature = calibration$signature,
       source_code = calibration$source_code, source_name = calibration$source_name,
       calibration_file = calibration$cache_file)
}

export_shared_design <- function(design) {
  save_dt(design$bands, file.path(DIR_TABLE, "table_02a_shared_reference_intervals.csv"))
  s <- copy(design$lags$summary)
  d <- copy(design$lags$diagnostics)
  s[, applied_outcome := OUTCOME_NAME]
  d[, applied_outcome := OUTCOME_NAME]
  save_dt(s, file.path(DIR_TABLE, "table_03a_auto_selected_lag_windows.csv"))
  save_dt(d, file.path(DIR_TABLE, "table_03b_auto_lag_selection_diagnostics.csv"))
  save_obj(s, file.path(DIR_MODEL, "auto_selected_lag_windows.qs"))
  save_obj(d, file.path(DIR_MODEL, "auto_lag_selection_diagnostics.qs"))
  save_obj(design, file.path(DIR_MODEL, "shared_analysis_design.qs"))
}

copy_calibration_outputs <- function(calibration) {
  root <- calibration$output_root
  for (subdir in c("tables", "models", "plot_data", "figures")) {
    from <- list.files(file.path(root, subdir), full.names = TRUE)
    if (!length(from)) stop("Overall calibration output files missing: ", root, "; rebuild calibration.")
    if (!all(file.copy(from, file.path(OUT_ROOT, subdir), overwrite = TRUE)))
      stop("Could not copy overall DLNM outputs.")
  }
}

# ============================================================
# Stage-one input bundle and pipeline registration
# ============================================================

pipeline_stage1_begin <- function(outcome_code, design) {
  active <- list(protocol = "daynight_two_stage_v1", calibration_id = design$calibration_id,
                 analysis_tag = design$analysis_tag, settings = design$settings,
                 object_format = OBJECT_FORMAT)
  index <- if (file.exists(PIPELINE_INDEX_FILE)) readRDS(PIPELINE_INDEX_FILE) else NULL
  if (is.null(index)) {
    index <- list(protocol = "daynight_two_stage_v1", active_design = active,
                  work_dir = WORK_DIR, results_root = normalizePath(RESULTS_ROOT, winslash = "/", mustWork = TRUE),
                  outcomes = list())
    log_msg("Starting a pipeline index for this reference/window/calibration configuration.")
  }
  index$active_design <- active
  index$outcomes[[outcome_code]] <- NULL
  saveRDS(index, PIPELINE_INDEX_FILE)
  owned_markers <- file.path(DIR_MODEL, c("stage1_completed.rds", "stage2_completed.rds",
                                          "run_completed.qs", "run_completed.rds"))
  for (p in owned_markers[file.exists(owned_markers)])
    if (!file.remove(p)) stop("Cannot invalidate completion marker: ", p)
  invisible(active)
}

pipeline_stage1_register <- function(ready) {
  index <- readRDS(PIPELINE_INDEX_FILE)
  index$outcomes[[ready$outcome_code]] <- ready
  index$updated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  saveRDS(index, PIPELINE_INDEX_FILE)
  log_msg("Stage-one outcome registered for stage two:", ready$outcome_code)
  log_msg("Pipeline index:", PIPELINE_INDEX_FILE)
  invisible(ready)
}

export_covariance_csv <- function(cv, file) {
  cols <- setNames(lapply(seq_len(ncol(cv)), function(j)
    formatC(cv[, j], digits = 17L, format = "g")), colnames(cv))
  save_dt(as.data.table(c(list(state = rownames(cv)), cols)), file)
}

write_stage1_input_bundle <- function(rr_table, annual_counts, case_dt, compact, metadata) {
  folder <- DIR_TABLE
  required_tables <- c("table_04_joint_9_state_rr.csv", "table_03_joint_state_counts.csv",
                       "table_02_joint_state_thresholds.csv", "table_02a_shared_reference_intervals.csv",
                       "table_03a_auto_selected_lag_windows.csv",
                       "table_00a_cold_heat_overlap_diagnostics.csv",
                       "table_00b_joint_analysis_sample_flow.csv",
                       "table_00d_overlap_exclusion_distribution.csv",
                       "table_03c_matched_state_comparison_support.csv",
                       "table_03d_case_referent_state_comparisons.csv")
  for (nm in required_tables) {
    if (!file.exists(file.path(folder, nm))) stop("Required stage-one table missing: ", nm)
  }
  save_dt(annual_counts, file.path(folder, "table_00c_annual_state_death_counts.csv"))
  cv <- compact$covariance_non_reference
  export_covariance_csv(cv, file.path(folder, "table_04a_joint_logrr_covariance.csv"))
  saveRDS(list(covariance_non_reference = cv, mu_non_reference = compact$mu_non_reference),
          file.path(DIR_MODEL, "joint_logrr_covariance.rds"))
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
  } else if (file.exists(case_path) && !file.remove(case_path)) {
    stop("Cannot clear an obsolete stage-one case-state file.")
  }
  writeLines(c(
    "Current analysis outputs and downstream calculation inputs.",
    "table_04 contains state-specific risk estimates; table_00c contains annual final-sample death counts.",
    "joint_logrr_covariance.rds stores the full-precision covariance; the CSV provides an inspection copy.",
    "source_parameters.rds stores reference definitions, lag windows, sample information and computational settings.",
    "table_00d describes overlap patterns, direct exclusions, subsequent matched-set losses and retention.",
    "Its percentages use pre-exclusion records within each group and record role; missing geography is retained.",
    "Calendar distributions use record dates; winter is December-February. Different partitions must not be added together.",
    "table_03c reports final-sample within-set state contrast support and reference connectivity.",
    "table_03d gives case-state by referent-state pair counts and distinct matched sets in each cell.",
    "Referents are records, not independent people; matched sets can contribute to several comparison cells.",
    "Comparison support and connectivity do not establish absence of separation, confounding or imprecision.",
    "Run script 02 from the same working directory to complete the analysis without raw data or model refitting.",
    "Both stages write tables, figures and model objects into this analysis directory.",
    "A specific analysis can be selected with run_stage2_outcome(rr_file='tables/table_04_joint_9_state_rr.csv').",
    "Downstream scripts read current files from tables, models and plot_data."
  ), file.path(DIR_MODEL, "README.txt"), useBytes = TRUE)
  ready <- list(protocol = "daynight_two_stage_v1", stage1_id = metadata$stage1_id,
                outcome_code = OUTCOME_CODE, outcome_name = OUTCOME_NAME,
                output_root = metadata$stage1_output_root, input_directory = metadata$stage1_input_directory,
                calibration_id = CURRENT_SHARED_DESIGN$calibration_id,
                analysis_tag = CURRENT_SHARED_DESIGN$analysis_tag,
                n_cases = sum(annual_counts$n_cases), completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  save_obj(metadata, file.path(DIR_MODEL, "analysis_metadata_stage1.qs"))
  saveRDS(ready, file.path(DIR_MODEL, "stage1_completed.rds"))
  ready
}

run_stage1_internal_checks <- function() {
  states <- make_state_levels()
  info <- state_dt_from_levels(states)
  if (length(states) != 9L || uniqueN(states) != 9L || REF_STATE != states[1L] ||
      uniqueN(info$mask) != 9L || uniqueN(info$component_mask) != 9L)
    stop("The nine-state definition failed the stage-one consistency check.")
  invisible(TRUE)
}

# ============================================================
# Outcome-specific risk analysis
# ============================================================

run_stage1_outcome <- function(data_path, outcome_code, outcome_name = outcome_code,
                               out_base = file.path(RESULTS_ROOT, outcome_code),
                               reference_method = RUN_REFERENCE_METHOD,
                               rr_tolerance = RUN_RR_TOLERANCE,
                               lag_mode = RUN_LAG_MODE,
                               manual_ma_days = RUN_MANUAL_MA_DAYS,
                               manual_overrides = RUN_MANUAL_OVERRIDES,
                               auto_lag_no_sig_action = RUN_NO_SIG_ACTION,
                               fallback_ma_days = RUN_FALLBACK_MA_DAYS,
                               allow_truncated_rr_band = RUN_ALLOW_TRUNCATED_RR_BAND,
                               shared_dir = SHARED_DIR,
                               overall_paths = unname(DATA_PATHS)) {
  validate_run_settings(reference_method, rr_tolerance, lag_mode, manual_ma_days,
                        manual_overrides, auto_lag_no_sig_action, fallback_ma_days,
                        allow_truncated_rr_band)
  if (!is.character(data_path) || !length(data_path) || anyNA(data_path) || any(!file.exists(data_path)))
    stop("Requested outcome data files are missing: ", paste(data_path, collapse = " | "))
  if (!OBJECT_FORMAT %in% c("qs", "rds")) stop("OBJECT_FORMAT must be qs or rds.")
  if (OBJECT_FORMAT == "qs" && !requireNamespace("qs", quietly = TRUE)) stop("Install qs or select OBJECT_FORMAT = 'rds'.")
  calibration <- get_shared_calibration(overall_paths = overall_paths, shared_dir = shared_dir)
  design <- make_shared_design(calibration, reference_method, rr_tolerance, lag_mode,
                               manual_ma_days, manual_overrides, auto_lag_no_sig_action,
                               fallback_ma_days, allow_truncated_rr_band)
  for (z in list(outcome_code, outcome_name, out_base))
    if (!is.character(z) || length(z) != 1L || is.na(z) || !nzchar(z)) stop("Invalid outcome/output name.")
  if (identical(outcome_code, OVERALL_CODE) &&
      !identical(unname(normalizePath(data_path, winslash = "/", mustWork = TRUE)), calibration$signature$paths))
    stop("Overall outcome call must use the same source files as the calibration.")
  CURRENT_SHARED_DESIGN <<- design
  STATE_THRESHOLD_STRATEGY <<- design$thresholds$strategy
  AUTO_LAG_NO_SIG_ACTION <<- auto_lag_no_sig_action
  AUTO_LAG_FALLBACK_HEAT <<- seq_len(fallback_ma_days[["heat"]]) - 1L
  AUTO_LAG_FALLBACK_COLD <<- seq_len(fallback_ma_days[["cold"]]) - 1L
  apply_auto_lag_windows(design$lags)
  DATA_PATH <<- data_path
  OUTCOME_CODE <<- outcome_code
  OUTCOME_NAME <<- outcome_name
  OUT_BASE <<- out_base
  OUT_ROOT <<- file.path(out_base, paste0("joint_analysis_", design$analysis_tag))
  set_output_dirs(OUT_ROOT, "stage1_risk_analysis_log.txt", SAVE_CASE_CONTRIB)
  manifest <- file.path(DIR_MODEL, "run_manifest.rds")
  identity <- list(design_id = design$calibration_id, design_settings = design$settings,
                   data = calibration_signature(data_path), outcome = outcome_code,
                   object_format = OBJECT_FORMAT)
  pipeline_stage1_begin(outcome_code, design)
  saveRDS(identity, manifest)
  log_msg("Starting:", OUTCOME_NAME, "| reference:", reference_method, "| shared calibration:", calibration$calibration_id)
  log_msg("Shared lag windows:", paste(design$lags$summary$selected_window, collapse = "; "))
  run_stage1_internal_checks()
  if (identical(outcome_code, OVERALL_CODE)) {
    copy_calibration_outputs(calibration)
    single <- calibration$single_results
  }
  raw <- read_analysis_data(data_path)
  if (!identical(outcome_code, OVERALL_CODE)) {
    single <- run_single_temperature_models(raw)
  }
  export_shared_design(design)
  risk <- run_joint_temperature_model(raw, calibration$single_results)
  rm(raw)
  invisible(gc())
  writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, "sessionInfo_stage1.txt"))
  pipeline_stage1_register(risk$ready)
  log_msg("Stage one completed. Script 02 will continue in:", OUT_ROOT)
  invisible(list(outcome_code = outcome_code, outcome_name = outcome_name,
                 output_root = OUT_ROOT, input_directory = risk$ready$input_directory,
                 shared_design = design, single_results = single,
                 auto_lag_results = design$lags, risk_results = risk))
}

# ============================================================
# Individual outcome calls
# ============================================================

if (!isTRUE(getOption("daynight.stage1.skip_calls", FALSE)) && RUN_OVERALL) {
  analysis_result <- run_stage1_outcome(
    data_path = unname(DATA_PATHS),
    outcome_code = OVERALL_CODE,
    outcome_name = OVERALL_NAME
  )
}

if (!isTRUE(getOption("daynight.stage1.skip_calls", FALSE)) && RUN_I10_I15) {
  analysis_result <- run_stage1_outcome(
    data_path = MATCHED_FILES[["I10_I15"]],
    outcome_code = "I10_I15", outcome_name = "I10-I15 Hypertensive diseases"
  )
}

if (!isTRUE(getOption("daynight.stage1.skip_calls", FALSE)) && RUN_I20_I25) {
  analysis_result <- run_stage1_outcome(
    data_path = MATCHED_FILES[["I20_I25"]],
    outcome_code = "I20_I25", outcome_name = "I20-I25 Ischemic heart disease"
  )
}

if (!isTRUE(getOption("daynight.stage1.skip_calls", FALSE)) && RUN_I50) {
  analysis_result <- run_stage1_outcome(
    data_path = MATCHED_FILES[["I50"]],
    outcome_code = "I50", outcome_name = "I50 Heart failure"
  )
}

if (!isTRUE(getOption("daynight.stage1.skip_calls", FALSE)) && RUN_I60_I62) {
  analysis_result <- run_stage1_outcome(
    data_path = MATCHED_FILES[["I60_I62"]],
    outcome_code = "I60_I62", outcome_name = "I60-I62 Hemorrhagic stroke"
  )
}

if (!isTRUE(getOption("daynight.stage1.skip_calls", FALSE)) && RUN_I63) {
  analysis_result <- run_stage1_outcome(
    data_path = MATCHED_FILES[["I63"]],
    outcome_code = "I63", outcome_name = "I63 Ischemic stroke"
  )
}

