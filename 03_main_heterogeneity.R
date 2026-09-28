#!/usr/bin/env Rscript
# =============================================================================
# Primary risk and burden-composition inference
#
# Purpose: test the joint association of temperature states, equality of four
# component ANs within each outcome, and heterogeneity of complete share profiles
# across the five mutually exclusive specific causes. The overall outcome is
# excluded from between-disease tests because it includes these causes.
#
# Inputs: current stage-one risks, full coefficient covariance, annual state counts,
# source_parameters.rds, and stage-two mean annual component burden tables in the
# selected primary model directories. Raw case-referent records are not required.
#
# Outputs: CSV tables and an Excel workbook in
# results/main/_core_heterogeneity_inference/<configuration>/tables; computational
# settings, diagnostics, and an RDS result in models. Existing outputs are replaced.
#
# Ten thousand correlated coefficient draws estimate component-share intervals
# and transformed covariance matrices. Holm adjustment applies separately to six
# global-risk tests, six within-outcome tests, and ten gated disease-pair tests.
# Cross-disease coefficient errors are treated as approximately independent.
# Reference intervals, windows, and counts are fixed during propagation. Profile
# inference requires stable nonzero total burden and adequate covariance rank.
# Functions only: options(daynight.stage3.skip_calls = TRUE); source(this_file).
# =============================================================================

# Local paths are configured through environment variables; defaults are relative.
PROJECT_ROOT <- normalizePath(Sys.getenv("DAYNIGHT_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
DATA_ROOT <- Sys.getenv("DAYNIGHT_DATA_ROOT", unset = file.path(PROJECT_ROOT, "data"))
RESULTS_BASE <- Sys.getenv("DAYNIGHT_RESULTS_ROOT", unset = file.path(PROJECT_ROOT, "results"))
WORK_DIR <- PROJECT_ROOT

if (.Platform$OS.type == "windows")
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)))

PIPELINE_INDEX_FILE <- file.path(RESULTS_BASE, "daynight_pipeline_index.rds")
ANALYSIS_DIRS <- NULL # Optional named vector: outcome code = analysis directory.
RESULTS_ROOT <- file.path(RESULTS_BASE, "main")
ANALYSIS_TAG <- NULL # Optional configuration folder name when several analyses exist.
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

script_directory <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    path <- frames[[i]]$ofile
    if (is.character(path) && length(path) == 1L && file.exists(path))
      return(dirname(normalizePath(path, winslash = "/", mustWork = TRUE)))
  }
  args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(args)) {
    path <- sub("^--file=", "", args[1L])
    if (file.exists(path)) return(dirname(normalizePath(path, winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
SCRIPT_DIR <- script_directory()

load_directory_index <- function(index_file, directories = ANALYSIS_DIRS) {
  codes <- c("I00_I52_I60_I69", "I10_I15", "I20_I25", "I50", "I60_I62", "I63")
  clean_tag <- function(x) sub("^joint_analysis_", "", x)
  scalar <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  risk_exists <- function(root) scalar(root) && dir.exists(root) &&
    any(file.exists(file.path(root, c("tables/table_04_joint_9_state_rr.csv", "table_04_joint_9_state_rr.csv"))))
  read_meta <- function(root) {
    path <- file.path(root, c("models/source_parameters.rds", "source_parameters.rds",
                              "models/analysis_metadata_stage1.rds", "models/analysis_metadata.rds"))
    for (f in path[file.exists(path)]) {
      value <- tryCatch(readRDS(f), error = function(e) NULL)
      if (is.list(value)) return(value)
    }
    list()
  }
  parent_chain <- function(path) {
    out <- path
    for (i in 1:3) out <- c(out, dirname(tail(out, 1L)))
    unique(out)
  }
  index_paths <- unique(c(index_file, file.path(SCRIPT_DIR, "daynight_pipeline_index.rds")))
  previous <- list()
  for (f in index_paths[file.exists(index_paths)]) {
    value <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.list(value) && length(value$outcomes)) { previous <- value; break }
  }
  selected <- list()
  records <- list()
  add_directory <- function(root, code = NULL) {
    if (!risk_exists(root) || basename(root) %in% c("tables", "models", "figures", "plot_data")) return(invisible(NULL))
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    meta <- read_meta(root)
    if (!scalar(code)) {
      code <- if (basename(dirname(root)) %in% codes) basename(dirname(root)) else if (basename(root) %in% codes) basename(root) else meta$outcome_code
    }
    if (!scalar(code) || !code %in% codes) return(invisible(NULL))
    tag <- if (startsWith(basename(root), "joint_analysis_")) clean_tag(basename(root)) else meta$shared_design$analysis_tag
    if (!scalar(tag)) tag <- basename(root)
    result_root <- if (basename(dirname(root)) == code) dirname(dirname(root)) else dirname(root)
    records[[root]] <<- list(code = code, root = root, tag = clean_tag(tag), results_root = result_root)
    invisible(NULL)
  }
  if (!is.null(directories)) {
    if (!is.character(directories) || is.null(names(directories)) || anyNA(names(directories)) ||
        any(!nzchar(names(directories))) || anyDuplicated(names(directories)))
      stop("ANALYSIS_DIRS must be a named character vector of outcome directories.")
    for (code in names(directories)) {
      path <- normalizePath(directories[[code]], winslash = "/", mustWork = TRUE)
      if (!risk_exists(path)) stop("Risk table not found in the configured outcome directory: ", path)
      add_directory(path, code)
      selected[[code]] <- records[[path]]
    }
  } else {
    if (!is.null(RESULTS_ROOT)) {
      if (!scalar(RESULTS_ROOT) || !dir.exists(RESULTS_ROOT)) stop("RESULTS_ROOT is not an existing directory: ", RESULTS_ROOT)
      bases <- normalizePath(RESULTS_ROOT, winslash = "/", mustWork = TRUE)
    } else {
      anchors <- unique(c(parent_chain(WORK_DIR), parent_chain(SCRIPT_DIR),
                          if (scalar(index_file)) dirname(index_file)))
      bases <- unique(c(file.path(anchors, "results", "main"),
                        file.path(anchors, "main"), anchors))
    }
    bases <- bases[dir.exists(bases)]
    for (base in bases) {
      add_directory(base)
      if (basename(base) %in% codes) {
        for (root in list.dirs(base, full.names = TRUE, recursive = FALSE)) add_directory(root, basename(base))
      }
      for (code in codes) {
        folder <- file.path(base, code)
        if (!dir.exists(folder)) next
        add_directory(folder, code)
        for (root in list.dirs(folder, full.names = TRUE, recursive = FALSE)) add_directory(root, code)
      }
    }
    discovered <- names(records)
    for (code in intersect(codes, names(previous$outcomes))) {
      path <- previous$outcomes[[code]]$output_root
      if (is.null(RESULTS_ROOT) && risk_exists(path)) add_directory(path, code)
    }
    if (length(records)) {
      tags <- unique(vapply(records, `[[`, character(1), "tag"))
      target <- if (scalar(ANALYSIS_TAG)) clean_tag(ANALYSIS_TAG) else NULL
      if (is.null(target)) {
        direct <- records[intersect(parent_chain(normalizePath(WORK_DIR, winslash = "/")), names(records))]
        if (length(direct)) target <- direct[[1L]]$tag
      }
      if (is.null(target) && scalar(previous$active_design$analysis_tag) &&
          clean_tag(previous$active_design$analysis_tag) %in% tags)
        target <- clean_tag(previous$active_design$analysis_tag)
      if (is.null(target) && length(tags) == 1L) target <- tags[1L]
      if (is.null(target)) stop("Several analysis configurations were found. Set ANALYSIS_TAG to one of: ",
                                paste(sort(tags), collapse = " | "), call. = FALSE)
      records <- records[vapply(records, function(x) identical(x$tag, target), logical(1))]
      for (code in codes) {
        candidates <- records[vapply(records, function(x) identical(x$code, code), logical(1))]
        local <- intersect(names(candidates), discovered)
        if (length(local)) candidates <- candidates[local]
        if (length(candidates) > 1L) {
          path <- previous$outcomes[[code]]$output_root
          if (scalar(path) && dir.exists(path)) {
            path <- normalizePath(path, winslash = "/", mustWork = TRUE)
            if (path %in% names(candidates)) candidates <- candidates[path]
          }
        }
        if (length(candidates) > 1L) stop("Several directories contain ", code,
                                          ". Set RESULTS_ROOT or ANALYSIS_DIRS to choose the result location: ", paste(names(candidates), collapse = " | "))
        if (length(candidates)) selected[[code]] <- candidates[[1L]]
      }
    }
  }
  if (!length(selected)) stop("No outcome risk tables were found. Expected <results root>/<outcome>/<analysis>/tables/table_04_joint_9_state_rr.csv. ",
                              "Set RESULTS_ROOT to the existing directory containing the outcome folders. Working directory: ", WORK_DIR,
                              "; script directory: ", SCRIPT_DIR, call. = FALSE)
  root <- selected[[1L]]$root
  meta <- readRDS(current_file(root, "source_parameters.rds"))
  active <- meta$shared_design
  if (!is.list(active)) active <- list()
  active$analysis_tag <- selected[[1L]]$tag
  if (is.null(active$settings)) active$settings <- list(reference_method = meta$reference_method)
  list(outcomes = lapply(selected, function(x) list(output_root = x$root)),
       results_root = if (scalar(RESULTS_ROOT)) normalizePath(RESULTS_ROOT, winslash = "/") else selected[[1L]]$results_root,
       active_design = active, work_dir = WORK_DIR)
}


OUTPUT_BASE <- NULL  # NULL uses the results root containing the selected outcome folders.
N_SIM <- 10000L
SIM_SEED <- 20260925L
ALPHA <- 0.05
WRITE_XLSX <- TRUE
SAVE_SIMULATION_DRAWS <- FALSE

# Prespecified numerical and Monte Carlo stability criteria.
EIGEN_REL_TOL <- 1e-10
RECONSTRUCTION_REL_TOL <- 1e-7
TOTAL_AN_REL_EPS <- 1e-8
MAX_TOTAL_SIGN_SWITCH_FRACTION <- 0.01
MAX_PROFILE_SPLIT_COV_REL_DIFF <- 0.50

OUTCOME_CONFIG <- data.frame(
  outcome_code = c("I00_I52_I60_I69", "I10_I15", "I20_I25", "I50", "I60_I62", "I63"),
  outcome = c("Overall cardiocerebrovascular mortality", "Hypertensive diseases",
              "Ischemic heart disease", "Heart failure", "Hemorrhagic stroke",
              "Ischemic stroke"),
  between_disease = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
RUN_OUTCOMES <- OUTCOME_CONFIG$outcome_code

REF_STATE <- "D0_N0"
STATE_LEVELS <- c("D0_N0", "DC_N0", "DH_N0", "D0_NC", "DC_NC", "DH_NC",
                  "D0_NH", "DC_NH", "DH_NH")
COMPONENTS <- c("D_C", "N_C", "D_H", "N_H")
COMPONENT_LABELS <- c(D_C = "Daytime cold", N_C = "Nighttime cold",
                      D_H = "Daytime heat", N_H = "Nighttime heat")
PROFILE_COORDINATES <- COMPONENTS[1:3]
PROTOCOL <- "daynight_two_stage_v1"

if (!requireNamespace("data.table", quietly = TRUE)) stop("Install package data.table.")
suppressPackageStartupMessages(library(data.table))

assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
  invisible(TRUE)
}

same_numeric <- function(a, b, tolerance = RECONSTRUCTION_REL_TOL) {
  a <- as.numeric(a); b <- as.numeric(b)
  length(a) == length(b) && all(is.finite(a)) && all(is.finite(b)) &&
    all(abs(a - b) <= tolerance * (1 + pmax(abs(a), abs(b))))
}

read_object <- function(path) {
  assert(file.exists(path), paste("Missing input:", path))
  extension <- tolower(tools::file_ext(path))
  if (extension == "rds") return(readRDS(path))
  if (extension == "qs") {
    assert(requireNamespace("qs", quietly = TRUE), "Install qs to read stage-two .qs metadata.")
    return(qs::qread(path))
  }
  stop("Unsupported object format: ", path)
}

describe_files <- function(paths) {
  paths <- unique(paths[file.exists(paths)])
  paths <- paths[!file.info(paths)$isdir]
  data.table(path = normalizePath(paths, winslash = "/", mustWork = TRUE),
             bytes = as.numeric(file.info(paths)$size))
}

check_psd <- function(V, label) {
  assert(is.matrix(V) && is.numeric(V) && nrow(V) == ncol(V) &&
           nrow(V) > 0L && all(is.finite(V)), paste("Invalid covariance:", label))
  scale <- max(abs(V), .Machine$double.eps)
  assert(max(abs(V - t(V))) <= 1e-8 * scale, paste("Asymmetric covariance:", label))
  V <- (V + t(V)) / 2
  eg <- eigen(V, symmetric = TRUE)
  threshold <- max(abs(eg$values), .Machine$double.eps) * EIGEN_REL_TOL
  assert(min(eg$values) >= -threshold, paste("Non-positive-semidefinite covariance:", label))
  list(matrix = V, eigen = eg, keep = eg$values > threshold)
}

wald_test <- function(estimate, covariance, required_rank = NULL) {
  estimate <- as.numeric(estimate)
  unavailable <- function(status, rank = NA_integer_) list(
    statistic = NA_real_, df = rank, p_value = NA_real_, status = status)
  if (any(!is.finite(estimate)) || any(!is.finite(covariance)))
    return(unavailable("not_estimable_nonfinite"))
  p <- check_psd(as.matrix(covariance), "Wald contrast")
  assert(length(estimate) == nrow(p$matrix), "Wald estimate/covariance size mismatch.")
  rank <- sum(p$keep)
  if (!is.null(required_rank) && rank != required_rank)
    return(unavailable("not_estimable_rank_deficient", rank))
  if (rank == 0L) return(unavailable("not_estimable_zero_variance", rank))
  if (any(!p$keep)) {
    unsupported <- crossprod(p$eigen$vectors[, !p$keep, drop = FALSE], estimate)
    if (max(abs(unsupported)) > 1e-7 * max(1, sqrt(sum(estimate^2))))
      return(unavailable("not_estimable_contrast_outside_covariance_range", rank))
  }
  coordinates <- crossprod(p$eigen$vectors[, p$keep, drop = FALSE], estimate)
  statistic <- sum(as.numeric(coordinates)^2 / p$eigen$values[p$keep])
  list(statistic = statistic, df = rank,
       p_value = pchisq(statistic, rank, lower.tail = FALSE),
       status = if (rank == length(estimate)) "estimated" else "estimated_reduced_rank")
}

inverse_full <- function(V) {
  p <- check_psd(V, "profile")
  assert(all(p$keep), "Profile covariance is rank deficient.")
  sweep(p$eigen$vectors, 2, 1 / p$eigen$values, "*") %*% t(p$eigen$vectors)
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
  out <- vapply(seq_len(ncol(draws)), function(j)
    quantile(draws[, j], c(.025, .975), type = 8, names = FALSE), numeric(2))
  dimnames(out) <- list(c("low", "high"), colnames(draws))
  out
}

# ============================================================
# Current analysis inputs
# ============================================================

design_key <- function(meta) {
  lag_names <- c(day_cold = "day_cold_class_lag_days", day_heat = "day_heat_class_lag_days",
                 night_cold = "night_cold_class_lag_days", night_heat = "night_heat_class_lag_days")
  thresholds <- c("day_cold_threshold", "day_heat_threshold", "night_cold_threshold", "night_heat_threshold")
  list(calibration_id = meta$calibration_id, analysis_tag = meta$shared_design$analysis_tag,
       reference_method = meta$reference_method, rr_tolerance = meta$rr_tolerance,
       thresholds = unlist(meta$state_thresholds[thresholds]),
       lag_column_indices = lapply(lag_names, function(name) as.integer(meta[[name]])),
       overlap_policy = meta$both_anomaly_policy, reference_source = meta$reference_source_code,
       settings = meta$shared_design$settings)
}

read_outcome_bundle <- function(code, index) {
  entry <- index$outcomes[[code]]
  assert(!is.null(entry), paste("No input directory configured for", code))
  root <- normalizePath(entry$output_root, winslash = "/", mustWork = TRUE)
  names <- c("source_parameters.rds", "table_04_joint_9_state_rr.csv",
             "table_00c_annual_state_death_counts.csv", "table_12_shapley_4_component_burden_mean_annual.csv")
  paths <- vapply(names, function(name) current_file(root, name), character(1))
  for (name in c("joint_logrr_covariance.rds", "table_04a_joint_logrr_covariance.csv"))
    current_file(root, name, required = FALSE)
  meta <- readRDS(paths[["source_parameters.rds"]])
  primary <- fread(paths[["table_12_shapley_4_component_burden_mean_annual.csv"]])
  assert(all(c("component", "AN", "AF_percent", "share_percent", "n_deaths") %in% names(primary)) &&
           nrow(primary) == 5L && !anyDuplicated(primary$component) && setequal(primary$component, c(COMPONENTS, "joint_total")),
         paste("Invalid component burden table:", code))
  ready <- list(stage1_id = meta$stage1_id,
                calibration_id = meta$calibration_id, analysis_tag = index$active_design$analysis_tag)
  list(code = code, root = root, folder = file.path(root, "tables"), ready = ready, metadata = meta,
       design = design_key(meta), primary = primary, input_files = describe_files(unname(paths)))
}

read_risk_counts <- function(bundle) {
  rr <- fread(file.path(bundle$folder, "table_04_joint_9_state_rr.csv"))
  need <- c("state", "beta", "se", "rr", "rr_low", "rr_high", "n_cases",
            "n_rows", "n_referents", "support_status")
  assert(all(need %in% names(rr)) && nrow(rr) == 9L && !anyDuplicated(rr$state) &&
           setequal(rr$state, STATE_LEVELS), "Invalid nine-state RR table.")
  rr <- rr[match(STATE_LEVELS, state)]
  for (v in c("n_cases", "n_rows", "n_referents"))
    assert(is.numeric(rr[[v]]) && all(is.finite(rr[[v]])) &&
             all(rr[[v]] >= 0 & rr[[v]] == floor(rr[[v]])), paste("Invalid counts:", v))
  assert(all(rr$n_rows == rr$n_cases + rr$n_referents), "Cases/referents do not sum to rows.")
  observed <- rr$n_rows > 0
  assert(all(rr$support_status == ifelse(observed, "observed_estimable", "not_observed")),
         "Invalid observed-state support flags.")
  for (v in c("beta", "se", "rr", "rr_low", "rr_high"))
    assert(is.numeric(rr[[v]]) && all(is.finite(rr[[v]][observed])), paste("Invalid risk column:", v))
  r <- rr[observed]
  assert(all(r$se >= 0 & r$rr > 0 & r$rr_low > 0 & r$rr_high > 0) &&
           same_numeric(log(r$rr), r$beta) &&
           same_numeric(log(r$rr_low), r$beta - 1.96 * r$se) &&
           same_numeric(log(r$rr_high), r$beta + 1.96 * r$se), "Inconsistent RR/beta/SE/CI values.")
  reference <- rr[state == REF_STATE]
  assert(reference$n_rows > 0 && reference$beta == 0 && reference$se == 0 && reference$rr == 1,
         "Invalid common reference state.")
  assert(all(rr[!observed, n_cases] == 0) && all(is.na(rr[!observed, beta])),
         "Unobserved states must retain non-estimable risks and zero case counts.")
  states <- rr[observed & state != REF_STATE, state]
  assert(length(states) > 0, "No estimable non-reference state.")
  mu <- setNames(rr$beta[match(states, rr$state)], states)
  csv <- file.path(bundle$folder, "table_04a_joint_logrr_covariance.csv")
  rds <- analysis_file(bundle$root, "joint_logrr_covariance.rds")
  candidates <- c(csv, rds); candidates <- candidates[file.exists(candidates)]
  assert(length(candidates) > 0, "A full coefficient covariance is required.")
  chosen <- candidates[which.max(as.numeric(file.info(candidates)$mtime))]
  if (tolower(tools::file_ext(chosen)) == "csv") {
    v <- fread(chosen); label <- intersect(c("state", "term"), names(v))
    assert(length(label) == 1L, "Covariance CSV requires one state/term column.")
    rn <- v[[label]]; v[, (label) := NULL]; V <- as.matrix(v); rownames(V) <- rn
  } else {
    saved <- readRDS(chosen)
    V <- if (is.matrix(saved)) saved else saved$covariance_non_reference
  }
  assert(is.matrix(V) && is.numeric(V) && !is.null(rownames(V)) && !is.null(colnames(V)), "Invalid covariance matrix.")
  rownames(V) <- sub("^joint_state", "", rownames(V)); colnames(V) <- sub("^joint_state", "", colnames(V))
  assert(!anyDuplicated(rownames(V)) && !anyDuplicated(colnames(V)) &&
           all(states %in% rownames(V)) && all(states %in% colnames(V)), "Covariance does not cover the current states.")
  V <- check_psd(V[states, states, drop = FALSE], "current covariance")$matrix
  supplied_se <- rr$se[match(states, rr$state)]; previous_se <- sqrt(diag(V))
  assert(!any(previous_se == 0 & supplied_se != 0), "A full covariance is needed for nonzero SEs with zero prior variance.")
  ratio <- ifelse(previous_se == 0, 1, supplied_se / previous_se)
  V <- V * outer(ratio, ratio)
  saveRDS(list(mu_non_reference = mu, covariance_non_reference = V), rds)
  fwrite(cbind(data.table(state = states), as.data.table(V)), csv)
  annual <- fread(file.path(bundle$folder, "table_00c_annual_state_death_counts.csv"))
  assert(all(c("year", "state", "n_cases", "n_deaths") %in% names(annual)),
         "Annual counts lack required columns.")
  for (v in c("year", "n_cases", "n_deaths"))
    assert(is.numeric(annual[[v]]) && all(is.finite(annual[[v]])) &&
             all(annual[[v]] == floor(annual[[v]])), paste("Invalid annual column:", v))
  assert(all(annual$n_cases >= 0 & annual$n_deaths > 0) &&
           !anyDuplicated(annual[, .(year, state)]), "Invalid or duplicate annual state counts.")
  years <- sort(unique(annual$year))
  assert(length(years) > 0, "No analysis year.")
  for (y in years) {
    z <- annual[year == y]
    assert(nrow(z) == 9L && setequal(z$state, STATE_LEVELS) && uniqueN(z$n_deaths) == 1L &&
             sum(z$n_cases) == z$n_deaths[1L], paste("Invalid state counts in year", y))
  }
  totals <- annual[, .(total = sum(n_cases)), by = state]
  assert(same_numeric(totals$total, rr$n_cases[match(totals$state, rr$state)], 0),
         "Annual counts do not match the current risk table.")
  beta <- setNames(rep(NA_real_, 9), STATE_LEVELS)
  beta[REF_STATE] <- 0; beta[states] <- mu
  list(rr = rr, states = states, mu = mu, covariance = V, beta = beta, annual = annual, years = years)
}

# ============================================================
# Shapley operators, point estimates, and conditional simulations
# ============================================================

state_components <- function(state) {
  tokens <- strsplit(state, "_", fixed = TRUE)[[1]]
  mapping <- c(DC = "D_C", DH = "D_H", NC = "N_C", NH = "N_H")
  unname(mapping[tokens[tokens %in% names(mapping)]])
}

SINGLETON <- c(D_C = "DC_N0", N_C = "D0_NC", D_H = "DH_N0", N_H = "D0_NH")

shapley_operator <- function(weights, estimable_states) {
  assert(setequal(names(weights), STATE_LEVELS) && all(is.finite(weights)) &&
           all(weights >= 0), "Invalid state weights.")
  weights <- weights[STATE_LEVELS]
  B <- matrix(0, 4, 9, dimnames = list(COMPONENTS, STATE_LEVELS))
  for (s in STATE_LEVELS) {
    active <- state_components(s)
    if (weights[[s]] == 0 || length(active) == 0L) next
    needed <- unique(c(s, unname(SINGLETON[active])))
    assert(all(needed %in% estimable_states),
           paste("Positive-weight state requires unavailable risks:", s,
                 paste(setdiff(needed, estimable_states), collapse = ", ")))
    if (length(active) == 1L) {
      B[active, s] <- B[active, s] + weights[[s]]
    } else {
      for (k in active) {
        other <- setdiff(active, k)
        B[k, SINGLETON[[k]]] <- B[k, SINGLETON[[k]]] + weights[[s]] / 2
        B[k, s] <- B[k, s] + weights[[s]] / 2
        B[k, SINGLETON[[other]]] <- B[k, SINGLETON[[other]]] - weights[[s]] / 2
      }
    }
  }
  expected <- weights; expected[REF_STATE] <- 0
  assert(same_numeric(colSums(B), expected), "Shapley additivity failed.")
  B[, estimable_states, drop = FALSE]
}

make_operators <- function(risk) {
  annual <- copy(risk$annual)
  summary <- annual[, .(mean_cases = mean(n_cases), mean_fraction = mean(n_cases / n_deaths)), by = state]
  counts <- setNames(summary$mean_cases, summary$state)
  fractions <- setNames(summary$mean_fraction, summary$state)
  assert(abs(sum(fractions) - 1) < 1e-10, "Mean annual state fractions do not sum to one.")
  list(AN = shapley_operator(counts, risk$states),
       AF_percent = 100 * shapley_operator(fractions, risk$states),
       mean_cases = counts, mean_fractions = fractions,
       mean_deaths = mean(unique(annual[, .(year, n_deaths)])$n_deaths))
}

validate_primary <- function(primary, an, af, mean_deaths) {
  p <- primary[match(c(COMPONENTS, "joint_total"), component)]
  target_an <- c(an, sum(an)); target_af <- c(af, sum(af))
  target_share <- if (sum(an) == 0) rep(NA_real_, 5) else c(an / sum(an) * 100, 100)
  data.table(metric = c("AN", "AF_percent", "share_percent"),
             max_absolute_difference = c(max(abs(p$AN - target_an)),
                                         max(abs(p$AF_percent - target_af)),
                                         if (all(is.na(target_share))) NA_real_ else max(abs(p$share_percent - target_share))))
}

profile_diagnostics <- function(an, draws) {
  draws <- as.matrix(draws)
  assert(ncol(draws) == 4L && nrow(draws) >= 4L && all(COMPONENTS %in% names(an)),
         "Profile diagnostics require four named components and at least four draws.")
  if (is.null(colnames(draws))) colnames(draws) <- COMPONENTS
  assert(setequal(colnames(draws), COMPONENTS), "Unexpected component names in profile draws.")
  draws <- draws[, COMPONENTS, drop = FALSE]
  an <- an[COMPONENTS]
  total <- sum(an)
  total_draws <- rowSums(draws)
  finite <- all(is.finite(draws)) && all(is.finite(total_draws))
  total_ci <- if (finite) quantile(total_draws, c(.025, .975), type = 8, names = FALSE) else c(NA_real_, NA_real_)
  epsilon <- TOTAL_AN_REL_EPS * max(1, sum(abs(an)))
  point_ok <- is.finite(total) && abs(total) > epsilon
  excludes_zero <- finite && (total_ci[1] > 0 || total_ci[2] < 0)
  sign_switch <- if (finite && point_ok) mean(sign(total_draws) != sign(total)) else NA_real_
  near_zero <- if (finite) mean(abs(total_draws) <= epsilon) else NA_real_
  ratio_draws <- matrix(NA_real_, nrow(draws), 4, dimnames = list(NULL, COMPONENTS))
  V <- matrix(NA_real_, 3, 3, dimnames = list(PROFILE_COORDINATES, PROFILE_COORDINATES))
  reasons <- character()
  if (!finite) reasons <- c(reasons, "nonfinite_AN_draws")
  if (!point_ok) reasons <- c(reasons, "total_AN_near_zero")
  if (!excludes_zero) reasons <- c(reasons, "total_AN_interval_includes_zero_or_unavailable")
  if (is.finite(sign_switch) && sign_switch > MAX_TOTAL_SIGN_SWITCH_FRACTION)
    reasons <- c(reasons, "excess_total_sign_switching")
  if (is.finite(near_zero) && near_zero > 0) reasons <- c(reasons, "near_zero_simulated_total")
  split_difference <- NA_real_
  rank <- NA_integer_
  max_abs_share <- NA_real_
  if (finite && all(total_draws != 0)) {
    ratio_draws <- sweep(draws, 1, total_draws, "/") * 100
    if (all(is.finite(ratio_draws))) {
      max_abs_share <- max(abs(ratio_draws))
      z <- ratio_draws[, PROFILE_COORDINATES, drop = FALSE]
      V <- cov(z)
      split <- floor(nrow(z) / 2)
      V1 <- cov(z[seq_len(split), , drop = FALSE])
      V2 <- cov(z[seq.int(split + 1L, nrow(z)), , drop = FALSE])
      split_difference <- norm(V1 - V2, "F") / max(norm(V, "F"), .Machine$double.eps)
      rank <- sum(check_psd(V, "simulated profile")$keep)
      if (rank != 3L) reasons <- c(reasons, "profile_covariance_rank_not_three")
      if (!is.finite(split_difference) || split_difference > MAX_PROFILE_SPLIT_COV_REL_DIFF)
        reasons <- c(reasons, "unstable_split_half_profile_covariance")
    } else reasons <- c(reasons, "nonfinite_share_draws")
  } else reasons <- c(reasons, "share_draws_unavailable")
  stable <- length(reasons) == 0L
  status <- if (stable) "eligible" else paste(unique(reasons), collapse = "; ")
  list(eligible = stable, status = status, covariance = V, share_draws = ratio_draws,
       row = data.table(
         total_mean_annual_AN = total, total_AN_CI_low = total_ci[1], total_AN_CI_high = total_ci[2],
         finite_AN_draw_fraction = mean(apply(is.finite(draws), 1, all)),
         total_sign_switch_fraction = sign_switch, near_zero_total_fraction = near_zero,
         max_absolute_simulated_share_percent = max_abs_share,
         profile_covariance_rank = rank, split_half_covariance_relative_difference = split_difference,
         profile_eligible = stable, profile_status = status,
         n_sim = nrow(draws), n_discarded_draws = 0L))
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
  covariance_an <- if (all(is.finite(an_draws))) cov(an_draws) else matrix(NA_real_, 4, 4)
  profile <- profile_diagnostics(an, an_draws)
  observed_share <- if (sum(an) == 0) setNames(rep(NA_real_, 4), COMPONENTS) else an / sum(an) * 100
  
  C <- matrix(c(1, -1, 0, 0, 1, 0, -1, 0, 1, 0, 0, -1), nrow = 3, byrow = TRUE,
              dimnames = list(c("D_C_minus_N_C", "D_C_minus_D_H", "D_C_minus_N_H"), COMPONENTS))
  within <- wald_test(C %*% an, C %*% covariance_an %*% t(C))
  global_risk <- wald_test(risk$mu, risk$covariance)
  peak <- risk$states[which.max(risk$mu)]
  peak_row <- risk$rr[state == peak]
  risk_row <- data.table(
    outcome_code = bundle$code, outcome = outcome_name,
    chi_square = global_risk$statistic, df = global_risk$df, p_value = global_risk$p_value,
    status = global_risk$status, n_estimable_nonreference_states = length(risk$states),
    unobserved_states = paste(setdiff(STATE_LEVELS, c(REF_STATE, risk$states)), collapse = "; "),
    highest_RR_state_descriptive = peak, highest_RR = peak_row$rr,
    highest_RR_CI_low = peak_row$rr_low, highest_RR_CI_high = peak_row$rr_high,
    null_hypothesis = "All estimable non-reference joint-state log relative risks equal zero",
    interpretation = "Global association with the reference included; the highest RR is descriptive")
  within_row <- data.table(
    outcome_code = bundle$code, outcome = outcome_name,
    chi_square = within$statistic, df = within$df, p_value = within$p_value, status = within$status,
    null_hypothesis = "Four mean annual component attributable numbers are equal",
    total_mean_annual_AN = sum(an), profile_eligible = profile$eligible, profile_status = profile$status)
  an_ci <- ci_columns(an_draws); af_ci <- ci_columns(af_draws)
  share_ci <- if (profile$eligible) ci_columns(profile$share_draws) else
    matrix(NA_real_, 2, 4, dimnames = list(c("low", "high"), COMPONENTS))
  components <- data.table(
    outcome_code = bundle$code, outcome = outcome_name, component = COMPONENTS,
    component_label = unname(COMPONENT_LABELS[COMPONENTS]),
    AN = as.numeric(an), AN_CI_low = an_ci[1, ], AN_CI_high = an_ci[2, ],
    AF_percent = as.numeric(af), AF_percent_CI_low = af_ci[1, ], AF_percent_CI_high = af_ci[2, ],
    share_percent = as.numeric(observed_share), share_CI_low = share_ci[1, ], share_CI_high = share_ci[2, ],
    profile_status = profile$status)
  diagnostics <- copy(profile$row)
  diagnostics[, `:=`(outcome_code = bundle$code, outcome = outcome_name, simulation_seed = seed,
                     negative_component_count = sum(an < 0), n_cases = sum(risk$annual$n_cases),
                     years = paste(risk$years, collapse = "; "),
                     max_AN_reconstruction_error = validation[metric == "AN", max_absolute_difference],
                     max_AF_reconstruction_error = validation[metric == "AF_percent", max_absolute_difference],
                     max_share_reconstruction_error = validation[metric == "share_percent", max_absolute_difference])]
  validation[, outcome_code := bundle$code]
  list(risk = risk_row, within = within_row, components = components, diagnostics = diagnostics,
       validation = validation, years = risk$years, AN = an, AF_percent = af,
       share = observed_share, profile = list(eligible = profile$eligible, status = profile$status,
                                              z = observed_share[PROFILE_COORDINATES], covariance = profile$covariance),
       covariance_AN = covariance_an, coefficient_mean = risk$mu, coefficient_covariance = risk$covariance,
       draws = if (SAVE_SIMULATION_DRAWS) list(AN = an_draws, AF_percent = af_draws,
                                               share_percent = profile$share_draws) else NULL)
}

# ============================================================
# Between-disease omnibus and hierarchical pairwise inference
# ============================================================

between_disease_tests <- function(results, alpha = ALPHA) {
  codes <- OUTCOME_CONFIG$outcome_code[OUTCOME_CONFIG$between_disease]
  available <- codes %in% names(results)
  reason <- character()
  if (!all(available)) reason <- c(reason, paste("missing_disease_outcomes", paste(codes[!available], collapse = ",")))
  if (all(available)) {
    unsuitable <- codes[!vapply(results[codes], function(x) x$profile$eligible, logical(1))]
    if (length(unsuitable)) reason <- c(reason, paste("ineligible_profiles", paste(unsuitable, collapse = ",")))
    if (!all(vapply(results[codes], function(x) identical(x$years, results[[codes[1]]]$years), logical(1))))
      reason <- c(reason, "different_analysis_years")
  }
  omnibus <- data.table(test = "Five-disease multivariate burden-profile Wald/Q test",
                        n_diseases_planned = 5L, profile_dimensions = 3L, chi_square = NA_real_, df = 12L,
                        p_value = NA_real_, status = if (length(reason)) paste(reason, collapse = "; ") else "estimated",
                        covariance_assumption = "Approximately independent disease-specific coefficient errors; shared calibration and counts fixed")
  if (!length(reason)) {
    precision <- lapply(results[codes], function(x) inverse_full(x$profile$covariance))
    pooled <- as.numeric(solve(Reduce(`+`, precision), Reduce(`+`, lapply(seq_along(codes), function(i)
      precision[[i]] %*% results[[codes[i]]]$profile$z))))
    Q <- sum(vapply(seq_along(codes), function(i) {
      difference <- results[[codes[i]]]$profile$z - pooled
      as.numeric(crossprod(difference, precision[[i]] %*% difference))
    }, numeric(1)))
    omnibus[, `:=`(chi_square = Q, p_value = pchisq(Q, 12L, lower.tail = FALSE))]
  }
  perform <- is.finite(omnibus$p_value) && omnibus$p_value < alpha
  pairs <- combn(codes, 2, simplify = FALSE)
  rows <- lapply(pairs, function(pair) {
    row <- data.table(disease_1 = pair[1], disease_2 = pair[2],
                      chi_square = NA_real_, df = 3L, p_value = NA_real_,
                      status = if (perform) "estimated" else if (length(reason)) "not_performed_omnibus_unavailable" else
                        "not_performed_omnibus_not_significant")
    for (k in COMPONENTS) row[, (paste0(k, "_share_difference_pp")) := NA_real_]
    if (perform) {
      x <- results[[pair[1]]]; y <- results[[pair[2]]]
      test <- wald_test(x$profile$z - y$profile$z,
                        x$profile$covariance + y$profile$covariance, required_rank = 3L)
      row[, `:=`(chi_square = test$statistic, df = test$df, p_value = test$p_value, status = test$status)]
      for (k in COMPONENTS) row[, (paste0(k, "_share_difference_pp")) := x$share[[k]] - y$share[[k]]]
    }
    row
  })
  pairwise <- rbindlist(rows)
  pairwise[, p_holm := p.adjust(p_value, method = "holm", n = 10L)]
  pairwise[, significant_after_holm := fifelse(is.finite(p_holm), p_holm < alpha, NA)]
  pairwise[, `:=`(disease_1_name = OUTCOME_CONFIG$outcome[match(disease_1, OUTCOME_CONFIG$outcome_code)],
                  disease_2_name = OUTCOME_CONFIG$outcome[match(disease_2, OUTCOME_CONFIG$outcome_code)])]
  list(omnibus = omnibus, pairwise = pairwise)
}

# ============================================================
# Reporting and execution
# ============================================================

add_holm_family <- function(x, alpha) {
  x <- copy(x)
  x[, p_holm_across_outcomes := p.adjust(p_value, method = "holm", n = nrow(x))]
  x[, significant_after_holm := fifelse(is.finite(p_holm_across_outcomes),
                                        p_holm_across_outcomes < alpha, NA)]
  x
}

format_p <- function(p) {
  out <- rep("NE", length(p))
  valid <- is.finite(p)
  out[valid & p < .001] <- "<0.001"
  out[valid & p >= .001] <- sprintf("%.3f", p[valid & p >= .001])
  out
}

write_workbook <- function(tables, file) {
  assert(requireNamespace("openxlsx", quietly = TRUE), "Install openxlsx or set WRITE_XLSX <- FALSE.")
  wb <- openxlsx::createWorkbook(creator = "Research analysis")
  title_style <- openxlsx::createStyle(fontSize = 14, textDecoration = "bold", fontColour = "#17365D")
  note_style <- openxlsx::createStyle(fontSize = 10, wrapText = TRUE, valign = "top")
  p_style <- openxlsx::createStyle(numFmt = "0.000E+00")
  setup_sheet <- function(sheet, title, note) {
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, title, startRow = 1, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = 1:10, rows = 1)
    openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = 1)
    openxlsx::writeData(wb, sheet, note, startRow = 2, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = 1:10, rows = 2)
    openxlsx::addStyle(wb, sheet, note_style, rows = 2, cols = 1)
    openxlsx::setRowHeights(wb, sheet, rows = 2, heights = 45)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 5, firstActiveCol = 3)
  }
  write_table <- function(sheet, x, row = 4L) {
    x <- as.data.frame(x)
    if (!nrow(x)) return(row)
    openxlsx::writeDataTable(wb, sheet, x, startRow = row, tableStyle = "TableStyleMedium2",
                             keepNA = TRUE, na.string = "NE")
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(x)), widths = 23)
    text_cols <- which(vapply(x, is.character, logical(1)))
    if (length(text_cols)) {
      openxlsx::setColWidths(wb, sheet, cols = text_cols, widths = 38)
      openxlsx::addStyle(wb, sheet, note_style, rows = seq.int(row + 1L, row + nrow(x)),
                         cols = text_cols, gridExpand = TRUE, stack = TRUE)
    }
    pcols <- which(grepl("^p_value$|^p_holm", names(x)))
    if (length(pcols)) openxlsx::addStyle(wb, sheet, p_style,
                                          rows = seq.int(row + 1L, row + nrow(x)), cols = pcols, gridExpand = TRUE, stack = TRUE)
    row + nrow(x) + 4L
  }
  setup_sheet("A_Risk_global", "Joint-state mortality risk: global association",
              "One Wald test per outcome. Holm correction across the selected outcomes. Highest-RR states are descriptive. NE denotes unavailable or unperformed estimates; see status columns.")
  write_table("A_Risk_global", tables$risk)
  setup_sheet("B_Burden_main", "Between-disease heterogeneity in four-component burden profiles",
              "Profiles are mean annual component AN / mean annual total AN. Overall mortality is excluded. Ten pairwise tests follow a significant omnibus test and use Holm correction. Covariances assume approximate cross-disease independence with calibration fixed.")
  row <- write_table("B_Burden_main", tables$between)
  write_table("B_Burden_main", tables$pairwise, row)
  setup_sheet("C_Burden_within", "Within-outcome equality of component attributable numbers",
              "Three linear contrasts test equality of four mean annual ANs. AF is the arithmetic mean of annual AF. Shares use AN and retain negative contributions. Unstable profiles have no share interval; AN/AF estimates remain available when finite.")
  row <- write_table("C_Burden_within", tables$within)
  write_table("C_Burden_within", tables$components, row)
  setup_sheet("D_Diagnostics", "Monte Carlo stability and reconstruction checks",
              "No simulation draws are trimmed. Profile eligibility requires a nonzero total with a 95% interval excluding zero, limited sign switching, no numerically near-zero denominators, and stable full-rank covariance. Maximum share magnitude is diagnostic only.")
  row <- write_table("D_Diagnostics", tables$diagnostics)
  write_table("D_Diagnostics", tables$validation, row)
  setup_sheet("E_Settings", "Analysis settings and inference scope",
              "All inference is conditional on the selected reference intervals, exposure windows and observed state death counts. Shared-calibration and between-disease covariance uncertainty are not propagated.")
  write_table("E_Settings", tables$settings)
  setup_sheet("F_Inputs", "Input provenance and file-integrity checks",
              "Stage-one bundles and stage-two completion markers are checked before analysis. Input files are hashed before and after execution. Source research results are not modified.")
  row <- write_table("F_Inputs", tables$sources)
  write_table("F_Inputs", tables$input_files, row)
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

run_stage3 <- function(pipeline_index_file = PIPELINE_INDEX_FILE,
                       outcome_codes = RUN_OUTCOMES, output_base = OUTPUT_BASE,
                       n_sim = N_SIM, sim_seed = SIM_SEED, alpha = ALPHA,
                       write_xlsx = WRITE_XLSX) {
  assert(length(n_sim) == 1 && is.finite(n_sim) && n_sim == as.integer(n_sim) && n_sim >= 1000L,
         "n_sim must be an integer of at least 1000.")
  assert(length(sim_seed) == 1 && is.finite(sim_seed) && sim_seed >= 0 &&
           sim_seed + 10000 < .Machine$integer.max && sim_seed == floor(sim_seed), "Invalid seed.")
  assert(length(alpha) == 1 && is.finite(alpha) && alpha > 0 && alpha < 1, "Invalid alpha.")
  assert(length(outcome_codes) > 0 && !anyDuplicated(outcome_codes) &&
           all(outcome_codes %in% OUTCOME_CONFIG$outcome_code), "Unknown or duplicate outcomes.")
  assert(is.finite(MAX_TOTAL_SIGN_SWITCH_FRACTION) && MAX_TOTAL_SIGN_SWITCH_FRACTION >= 0 &&
           MAX_TOTAL_SIGN_SWITCH_FRACTION < .5 && is.finite(MAX_PROFILE_SPLIT_COV_REL_DIFF) &&
           MAX_PROFILE_SPLIT_COV_REL_DIFF > 0 && is.finite(TOTAL_AN_REL_EPS) && TOTAL_AN_REL_EPS > 0,
         "Invalid profile stability settings.")
  if (write_xlsx) assert(requireNamespace("openxlsx", quietly = TRUE), "Install package openxlsx.")
  index <- load_directory_index(pipeline_index_file)
  pipeline_files <- if (file.exists(pipeline_index_file)) describe_files(pipeline_index_file) else data.table()
  if (is.null(output_base)) output_base <- file.path(index$results_root, "_core_heterogeneity_inference")
  tag <- index$active_design$analysis_tag
  assert(is.character(tag) && length(tag) == 1L && nzchar(tag) && !grepl("[/\\\\]", tag) &&
           !tag %in% c(".", ".."), "Invalid analysis tag.")
  base <- file.path(output_base, tag)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  output <- base
  assert(dir.exists(output), "Cannot create stage-three output directory.")
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  prepare_analysis_directory(output)
  log_file <- analysis_file(output, "core_heterogeneity_inference_log.txt")
  log_line <- function(...) {
    line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " "), sep = " | ")
    cat(line, "\n"); cat(line, "\n", file = log_file, append = TRUE)
  }
  log_line("Stage three started. Output:", output)
  success <- FALSE
  on.exit({if (!success) log_line("Stage three failed; no completion marker was written.")}, add = TRUE)
  bundles <- results <- list()
  common_design <- NULL
  for (code in outcome_codes) {
    log_line("Validating:", code)
    bundle <- read_outcome_bundle(code, index)
    if (is.null(common_design)) common_design <- bundle$design
    bundles[[code]] <- bundle
  }
  for (code in outcome_codes) {
    position <- match(code, OUTCOME_CONFIG$outcome_code)
    log_line("Analysing:", code, "with", n_sim, "coefficient draws")
    result <- analyse_outcome(bundles[[code]], OUTCOME_CONFIG$outcome[position], n_sim,
                              as.integer(sim_seed + position * 100L))
    results[[code]] <- result
    log_line("Profile status:", code, result$profile$status)
  }
  disease_codes <- setdiff(OUTCOME_CONFIG$outcome_code, "I00_I52_I60_I69")
  common <- TRUE
  if (all(disease_codes %in% names(bundles))) {
    design_fields <- c("reference_method", "rr_tolerance", "thresholds", "lag_column_indices", "overlap_policy")
    reference_design <- bundles[[disease_codes[1]]]$design[design_fields]
    common <- all(vapply(bundles[disease_codes], function(b)
      isTRUE(all.equal(b$design[design_fields], reference_design, tolerance = 1e-8)), logical(1)))
  }
  between <- between_disease_tests(if (common) results else list(), alpha)
  if (!common) {
    between$omnibus[, status := "Not performed: exposure definitions differ across outcomes"]
    between$pairwise[, status := "Not performed: exposure definitions differ across outcomes"]
  }
  sources <- rbindlist(lapply(outcome_codes, function(code) {
    b <- bundles[[code]]
    data.table(outcome_code = code, output_root = b$root, input_directory = b$folder,
               stage1_id = b$ready$stage1_id,
               calibration_id = b$ready$calibration_id, analysis_tag = b$ready$analysis_tag)
  }))
  inputs <- unique(rbindlist(c(list(pipeline_files), lapply(bundles, `[[`, "input_files"))))
  settings <- list(
    script_version = "daynight-stage3-1.0", completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    pipeline_index = normalizePath(pipeline_index_file, winslash = "/", mustWork = FALSE),
    outcomes = outcome_codes, n_sim = n_sim, seed_base = sim_seed, alpha = alpha,
    calibration_id = common_design$calibration_id, analysis_tag = common_design$analysis_tag,
    reference_method = common_design$reference_method, rr_tolerance = common_design$rr_tolerance,
    reference_thresholds = common_design$thresholds,
    model_lag_windows = lapply(common_design$lag_column_indices, function(x) x - 1L),
    overlap_policy = common_design$overlap_policy,
    profile_definition = "Mean annual component AN / mean annual total AN; equivalent to pooled study-period AN share",
    mean_AF_definition = "Arithmetic mean of annual AF; not used to normalize the burden profile",
    within_null = "Four mean annual attributable numbers equal; three independent linear contrasts",
    cross_disease_assumption = "Approximately independent coefficient errors across mutually exclusive underlying-cause outcomes",
    uncertainty_scope = "Coefficient uncertainty only; shared calibration, lag selection and death counts treated as fixed",
    negative_contributions = "Retained without truncation",
    profile_CI_rule = "Total AN percentile 95% interval must exclude zero",
    total_AN_relative_epsilon = TOTAL_AN_REL_EPS,
    max_total_sign_switch_fraction = MAX_TOTAL_SIGN_SWITCH_FRACTION,
    max_profile_split_half_covariance_relative_difference = MAX_PROFILE_SPLIT_COV_REL_DIFF,
    split_half_definition = "Frobenius norm of first-half minus second-half share covariance / norm of full share covariance",
    extreme_share_rule = "Maximum absolute simulated share is reported; no absolute share cap or draw trimming",
    multiplicity = "Separate Holm families for outcome risk tests, within-outcome tests, and 10 gated disease pairs",
    CI_definition = "2.5th and 97.5th empirical percentiles, R quantile type 8; point estimates match stage two",
    reconstruction_relative_tolerance = RECONSTRUCTION_REL_TOL,
    eigen_relative_tolerance = EIGEN_REL_TOL)
  settings_table <- rbindlist(lapply(names(settings), function(k)
    data.table(setting = k, value = paste(capture.output(dput(settings[[k]])), collapse = " "))))
  tables <- list(
    risk = add_holm_family(rbindlist(lapply(results, `[[`, "risk")), alpha),
    between = between$omnibus, pairwise = between$pairwise,
    within = add_holm_family(rbindlist(lapply(results, `[[`, "within")), alpha),
    components = rbindlist(lapply(results, `[[`, "components")),
    diagnostics = rbindlist(lapply(results, `[[`, "diagnostics")),
    validation = rbindlist(lapply(results, `[[`, "validation")),
    sources = sources, settings = settings_table, input_files = inputs)
  for (k in c("risk", "between", "pairwise", "within")) {
    tables[[k]][, p_text := format_p(p_value)]
    adjusted <- intersect(c("p_holm", "p_holm_across_outcomes"), names(tables[[k]]))
    if (length(adjusted)) tables[[k]][, p_holm_text := format_p(get(adjusted))]
  }
  filenames <- c(risk = "table_01_joint_state_risk_global.csv",
                 between = "table_02_between_disease_profile_global.csv",
                 pairwise = "table_03_between_disease_profile_pairwise.csv",
                 within = "table_04_within_outcome_equal_AN.csv",
                 components = "table_05_component_burden_and_profile.csv",
                 diagnostics = "table_06_monte_carlo_diagnostics.csv",
                 validation = "table_07_stage2_reconstruction_checks.csv",
                 sources = "table_08_source_outcomes.csv", settings = "table_09_analysis_settings.csv",
                 input_files = "table_10_input_files.csv")
  for (k in names(filenames)) fwrite(tables[[k]], analysis_file(output, filenames[[k]]), na = "NA")
  if (write_xlsx) write_workbook(tables, analysis_file(output, "table_core_heterogeneity_inference.xlsx"))
  saved <- list(settings = settings, tables = tables, outcomes = results, shared_design = common_design)
  saveRDS(saved, analysis_file(output, "core_heterogeneity_inference_results.rds"))
  writeLines(capture.output(sessionInfo()), analysis_file(output, "sessionInfo_stage3.txt"))
  output_files <- list.files(output, recursive = TRUE, full.names = TRUE)
  output_files <- output_files[basename(output_files) != basename(log_file) & !file.info(output_files)$isdir]
  output_manifest <- describe_files(output_files)
  fwrite(output_manifest, analysis_file(output, "output_manifest.csv"))
  saveRDS(list(protocol = "daynight_stage3_v1", calibration_id = common_design$calibration_id,
               analysis_tag = common_design$analysis_tag, input_input_files = inputs,
               completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
          analysis_file(output, "stage3_completed.rds"))
  success <- TRUE
  log_line("Between-disease test:", tables$between$status, "P =", tables$between$p_value)
  log_line("Stage three completed:", output)
  invisible(list(output_root = output, tables = tables, results = results))
}

if (!isTRUE(getOption("daynight.stage3.skip_calls", FALSE))) {
  run_stage3()
}
