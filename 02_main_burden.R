#!/usr/bin/env Rscript
# =============================================================================
# Primary attributable burden and Shapley allocation
#
# Purpose: transform nine-state log-RRs into signed attributable fractions and
# allocate each state to daytime cold, nighttime cold, daytime heat, and nighttime
# heat. No risk model is fitted by this script.
#
# Inputs: the selected primary model directory, including tables/table_04_joint_9_state_rr.csv,
# tables/table_04a_joint_logrr_covariance.csv (or models/joint_logrr_covariance.rds),
# tables/table_00c_annual_state_death_counts.csv, tables/table_02_joint_state_thresholds.csv,
# and models/source_parameters.rds. State rows and covariance labels must agree;
# beta, se, rr, rr_low, and rr_high must be internally consistent.
#
# Outputs: annual and mean annual nine-state and four-component AN, AF, and shares;
# Monte Carlo intervals; Shapley weights and additivity checks; component plots and
# plot data. Outputs replace the corresponding files in the same model directory.
#
# AN and AF intervals use 1,000 correlated coefficient draws by default. Negative
# contributions are retained. Mean annual AF is the arithmetic mean of annual AF;
# composition is normalized from mean annual AN. The covariance-source options
# are documented in README.md. Settings and optional outcome switches are below.
# Functions only: options(daynight.stage2.skip_calls = TRUE); source(this_file).
# =============================================================================

# Local paths are configured through environment variables; defaults are relative.
PROJECT_ROOT <- normalizePath(Sys.getenv("DAYNIGHT_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
DATA_ROOT <- Sys.getenv("DAYNIGHT_DATA_ROOT", unset = file.path(PROJECT_ROOT, "data"))
RESULTS_BASE <- Sys.getenv("DAYNIGHT_RESULTS_ROOT", unset = file.path(PROJECT_ROOT, "results"))
WORK_DIR <- PROJECT_ROOT

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


RUN_COVARIANCE_SOURCE <- "auto" # "auto", "rescale", "csv", or "rds".

# NULL uses the corresponding setting saved by stage one.
RUN_N_SIM <- NULL
RUN_SIM_SEED <- NULL
RUN_CAP_NEGATIVE_AF <- NULL

RUN_OVERALL <- TRUE
RUN_I10_I15 <- TRUE
RUN_I20_I25 <- TRUE
RUN_I50 <- TRUE
RUN_I60_I62 <- TRUE
RUN_I63 <- TRUE

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install packages first: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
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

# ============================================================
# Input/output helpers and publication figures
# ============================================================

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

save_obj <- function(obj, file) {
  if (OBJECT_FORMAT == "rds") {
    file <- sub("[.]qs$", ".rds", file)
    saveRDS(obj, file)
  } else if (OBJECT_FORMAT == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Saving .qs requires package qs.")
    qs::qsave(obj, file, preset = "fast")
  } else stop("OBJECT_FORMAT must be qs or rds.")
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

save_plot_device <- function(p, file, width, height, pdf = FALSE) {
  destination <- normalizePath(dirname(file), winslash = "/", mustWork = TRUE)
  previous <- getwd()
  on.exit(setwd(previous), add = TRUE)
  setwd(destination)
  if (pdf) grDevices::pdf(basename(file), width = width, height = height) else
    grDevices::png(basename(file), width = width, height = height, units = "in", res = FIG_DPI, bg = "white")
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

# ============================================================
# Risk uncertainty, exposure states and burden calculations
# ============================================================

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

safe_af_from_beta <- function(beta) {
  af <- 1 - exp(-beta)
  if (CAP_NEGATIVE_AF) af <- pmax(af, 0)
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
    stop(
      "At least one non-estimable joint state has observed cases and therefore cannot be assigned burden: ",
      paste(non_estimable_with_cases, collapse = ", ")
    )
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
    if (abs(sum(counts) - n_deaths) > 1e-8) {
      stop("State counts do not sum to the total number of deaths for year ", y, ".")
    }
    
    dt_y <- data.table(
      year = as.integer(y),
      state = state_order,
      n_cases = as.numeric(counts),
      n_deaths = n_deaths,
      state_af_value = af_state[state_order]
    )
    
    dt_y[, AN := fifelse(n_cases == 0, 0, n_cases * state_af_value)]
    if (any(!is.finite(dt_y$AN))) {
      bad <- dt_y[!is.finite(AN), state]
      stop("Non-finite attributable deaths were obtained for states: ", paste(bad, collapse = ", "))
    }
    
    joint_AN <- sum(dt_y$AN)
    dt_y[, AF := AN / n_deaths]
    dt_y[, AF_percent := AF * 100]
    if (!is.finite(joint_AN) || joint_AN == 0) {
      dt_y[, share_percent := NA_real_]
    } else {
      dt_y[, share_percent := AN / joint_AN * 100]
    }
    dt_y[]
  }))
  
  mean_annual <- annual[, .(
    n_cases = mean(n_cases, na.rm = TRUE),
    n_deaths = mean(n_deaths, na.rm = TRUE),
    AN = mean(AN, na.rm = TRUE),
    AF = mean(AF, na.rm = TRUE),
    AF_percent = mean(AF_percent, na.rm = TRUE)
  ), by = state]
  
  mean_joint_AN <- mean_annual[, sum(AN)]
  if (!is.finite(mean_joint_AN) || mean_joint_AN == 0) {
    mean_annual[, share_percent := NA_real_]
  } else {
    mean_annual[, share_percent := AN / mean_joint_AN * 100]
  }
  mean_annual[, period := "mean_annual"]
  
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
  
  phi <- matrix(
    NA_real_,
    nrow = nrow(state_info),
    ncol = length(COMPONENTS),
    dimnames = list(state_info$state, COMPONENTS)
  )
  phi[REF_STATE, ] <- 0
  
  for (ii in seq_len(nrow(state_info))) {
    st <- state_info$state[ii]
    if (st == REF_STATE) next
    
    if (!is.finite(af_state[st])) next
    
    active <- COMPONENTS[
      as.integer(unlist(state_info[ii, ..COMPONENTS], use.names = FALSE)) == 1L
    ]
    
    phi[st, ] <- 0
    
    if (length(active) == 0L) next
    
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
      stop(
        "An observed compound state requires a non-estimable single-component counterfactual for Shapley decomposition: ",
        st, ". Required single-component states: ", state_a, ", ", state_b, "."
      )
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
  if (length(counts_by_year) == 0L) stop("No annual death counts were available.")
  
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
    if (abs(sum(counts) - n_deaths) > 1e-8) {
      stop("State counts do not sum to the total number of deaths for year ", y, ".")
    }
    
    positive_states <- state_order[counts > 0]
    if (length(positive_states) > 0L) {
      finite_phi <- apply(
        phi_state[positive_states, COMPONENTS, drop = FALSE],
        1,
        function(z) all(is.finite(z))
      )
      if (any(!finite_phi)) {
        bad_states <- positive_states[!finite_phi]
        stop(
          "Observed case states lack estimable Shapley values in year ", y, ": ",
          paste(bad_states, collapse = ", ")
        )
      }
      
      component_an <- as.numeric(crossprod(
        counts[positive_states],
        phi_state[positive_states, COMPONENTS, drop = FALSE]
      ))
    } else {
      component_an <- rep(0, length(COMPONENTS))
    }
    
    names(component_an) <- COMPONENTS
    joint_an <- sum(component_an)
    
    out_y <- data.table(
      year = as.integer(y),
      component = COMPONENTS,
      n_deaths = n_deaths,
      AN = component_an
    )
    out_y[, AF := AN / n_deaths]
    out_y[, AF_percent := AF * 100]
    if (!is.finite(joint_an) || joint_an == 0) {
      out_y[, share_percent := NA_real_]
    } else {
      out_y[, share_percent := AN / joint_an * 100]
    }
    
    joint_y <- data.table(
      year = as.integer(y),
      component = "joint_total",
      n_deaths = n_deaths,
      AN = joint_an,
      AF = joint_an / n_deaths,
      AF_percent = joint_an / n_deaths * 100,
      share_percent = 100
    )
    
    rbindlist(list(out_y, joint_y), use.names = TRUE, fill = TRUE)
  }))
  
  mean_annual <- annual[, .(
    n_deaths = mean(n_deaths, na.rm = TRUE),
    AN = mean(AN, na.rm = TRUE),
    AF = mean(AF, na.rm = TRUE),
    AF_percent = mean(AF_percent, na.rm = TRUE)
  ), by = component]
  
  mean_joint_an <- mean_annual[component == "joint_total", AN]
  if (length(mean_joint_an) != 1L || !is.finite(mean_joint_an) || mean_joint_an == 0) {
    mean_annual[, share_percent := NA_real_]
  } else {
    mean_annual[component != "joint_total", share_percent := AN / mean_joint_an * 100]
    mean_annual[component == "joint_total", share_percent := 100]
  }
  mean_annual[, period := "mean_annual"]
  
  if (is.finite(mean_joint_an) && mean_joint_an != 0) {
    component_share_sum <- mean_annual[
      component %in% COMPONENTS,
      if (any(is.finite(share_percent))) sum(share_percent, na.rm = TRUE) else NA_real_
    ]
    if (!is.finite(component_share_sum) || abs(component_share_sum - 100) > 1e-8) {
      stop("Mean-annual Shapley component shares failed the 100% sum check.")
    }
  }
  
  annual <- merge(annual, COMPONENT_INFO, by = "component", all.x = TRUE, sort = FALSE)
  mean_annual <- merge(mean_annual, COMPONENT_INFO, by = "component", all.x = TRUE, sort = FALSE)
  setorder(annual, year, component_order)
  setorder(mean_annual, component_order)
  
  list(
    annual = annual,
    mean_annual = mean_annual,
    phi_by_state = phi_state
  )
}

add_ci <- function(obs, sims, by_cols, value_cols) {
  validate_unique_key(obs, by_cols, "Observed burden table")
  if (nrow(sims) == 0L) stop("The simulation table is empty.")
  
  missing_obs <- setdiff(c(by_cols, value_cols), names(obs))
  missing_sims <- setdiff(c(by_cols, value_cols), names(sims))
  if (length(missing_obs) > 0L) {
    stop("Observed burden table is missing columns: ", paste(missing_obs, collapse = ", "))
  }
  if (length(missing_sims) > 0L) {
    stop("Simulation burden table is missing columns: ", paste(missing_sims, collapse = ", "))
  }
  
  ci <- sims[
    ,
    c(
      setNames(
        lapply(.SD, safe_empirical_quantile, probability = 0.025),
        paste0(value_cols, "_low")
      ),
      setNames(
        lapply(.SD, safe_empirical_quantile, probability = 0.975),
        paste0(value_cols, "_high")
      )
    ),
    by = by_cols,
    .SDcols = value_cols
  ]
  
  validate_unique_key(ci, by_cols, "Simulation confidence-interval table")
  merge(obs, ci, by = by_cols, all.x = TRUE, sort = FALSE)
}

run_internal_checks <- function() {
  states <- make_state_levels()
  if (length(states) != 9L || uniqueN(states) != 9L || REF_STATE != states[1]) {
    stop("The joint-state definition must contain 9 unique states with D0_N0 as the reference.")
  }
  
  state_info <- state_dt_from_levels(states)
  if (nrow(state_info) != 9L ||
      uniqueN(state_info$mask) != 9L ||
      uniqueN(state_info$component_mask) != 9L) {
    stop("The joint-state metadata failed the uniqueness check.")
  }
  
  beta_test <- setNames(
    c(0, 0.04, 0.06, 0.03, 0.08, 0.10, 0.05, 0.09, 0.12),
    states
  )
  counts_test <- setNames(c(100, 20, 15, 25, 10, 8, 18, 7, 5), states)
  deaths_test <- sum(counts_test)
  counts_by_year_test <- list(`2013` = counts_test, `2014` = counts_test + 1)
  n_deaths_by_year_test <- list(
    `2013` = deaths_test,
    `2014` = sum(counts_test + 1)
  )
  
  shapley_test <- compute_shapley_burden_tables(
    beta_state = beta_test,
    counts_by_year = counts_by_year_test,
    n_deaths_by_year = n_deaths_by_year_test,
    state_info = state_info
  )
  
  if (nrow(shapley_test$annual) != 10L ||
      nrow(shapley_test$mean_annual) != 5L) {
    stop("The Shapley burden table dimensions failed the internal check.")
  }
  
  share_sum <- shapley_test$mean_annual[
    component %in% COMPONENTS,
    sum(share_percent)
  ]
  if (!is.finite(share_sum) || abs(share_sum - 100) > 1e-8) {
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

# ============================================================
# Risk, covariance and annual-count validation
# ============================================================

read_result_object <- function(path) {
  if (!file.exists(path)) stop("Result object not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(readRDS(path))
  if (ext == "qs") {
    if (!requireNamespace("qs", quietly = TRUE)) stop("Reading .qs requires package qs: ", path)
    return(qs::qread(path))
  }
  stop("Expected .qs or .rds object: ", path)
}

load_rr_input <- function(path) {
  rr <- if (is.data.frame(path)) copy(as.data.table(path)) else fread(path)
  need <- c("state", "rr", "rr_low", "rr_high", "n_cases", "n_rows", "n_referents", "n_sets", "support_status")
  missing <- setdiff(need, names(rr))
  if (length(missing)) stop("RR table missing columns from the original table_04: ", paste(missing, collapse = ", "))
  states <- make_state_levels()
  if (nrow(rr) != length(states) || anyDuplicated(rr$state) || !setequal(rr$state, states))
    stop("table_04 must contain exactly the nine unique original state names.")
  rr <- rr[match(states, state)]
  if (!"beta" %in% names(rr)) rr[, beta := log(rr)]
  if (!"se" %in% names(rr)) rr[, se := (log(rr_high) - log(rr_low)) / (2 * 1.96)]
  for (v in c("beta", "se", "rr", "rr_low", "rr_high", "n_cases", "n_rows", "n_referents", "n_sets"))
    if (!is.numeric(rr[[v]])) stop("Non-numeric RR table column: ", v)
  for (v in c("n_cases", "n_rows", "n_referents", "n_sets"))
    if (any(!is.finite(rr[[v]])) || any(rr[[v]] < 0) || any(rr[[v]] != floor(rr[[v]])))
      stop("Invalid counts in RR table column: ", v)
  if (any(rr$n_rows != rr$n_cases + rr$n_referents)) stop("RR n_rows differs from cases + referents.")
  valid_support <- c("observed_estimable", "not_observed")
  if (any(!rr$support_status %in% valid_support)) stop("Observed but non-estimable states cannot be resumed.")
  observed <- rr$n_rows > 0
  if (any((rr$support_status == "observed_estimable") != observed)) stop("Inconsistent RR state support flags/counts.")
  r <- rr[observed]
  if (any(!is.finite(r$beta)) || any(!is.finite(r$se)) || any(r$se < 0) ||
      any(!is.finite(r$rr)) || any(!is.finite(r$rr_low)) || any(!is.finite(r$rr_high)) ||
      any(r$rr <= 0 | r$rr_low <= 0 | r$rr_high <= 0)) stop("Invalid observed-state RR/beta/SE/CI.")
  close <- function(x, y) all(abs(x - y) <= 1e-6 * (1 + abs(y)))
  if (!close(log(r$rr), r$beta) || !close(log(r$rr_low), r$beta - 1.96 * r$se) ||
      !close(log(r$rr_high), r$beta + 1.96 * r$se)) {
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
  canonical <- merge(canonical, rr[, .(state, n_rows, n_cases, n_referents, n_sets, support_status)],
                     by = "state", all.x = TRUE, sort = FALSE)
  setorder(canonical, mask)
  list(table = canonical, original = rr,
       beta = setNames(rr$beta, rr$state), se = setNames(rr$se, rr$state),
       estimable = setdiff(rr[observed, state], REF_STATE), unobserved = rr[!observed, state])
}

normalise_annual_counts <- function(x, risks) {
  x <- copy(as.data.table(x))
  cols <- c("year", "state", "n_cases", "n_deaths")
  if (!all(cols %in% names(x))) stop("Annual counts require year, state, n_cases, n_deaths columns.")
  x <- x[, ..cols]
  validate_unique_key(x, c("year", "state"), "Annual state death counts")
  for (v in c("year", "n_cases", "n_deaths"))
    if (!is.numeric(x[[v]]) || any(!is.finite(x[[v]])) || any(x[[v]] != floor(x[[v]])))
      stop("Invalid integer annual count/year column: ", v)
  if (any(x$n_cases < 0) || any(x$n_deaths <= 0)) stop("Invalid annual deaths.")
  states <- make_state_levels()
  for (yr in unique(x$year)) {
    y <- x[year == yr]
    if (nrow(y) != 9L || !setequal(y$state, states) || uniqueN(y$n_deaths) != 1L ||
        sum(y$n_cases) != y$n_deaths[1L]) stop("Annual counts must contain all nine states and sum to n_deaths for year ", yr)
  }
  totals <- x[, .(total = sum(n_cases)), by = state]
  expected <- risks$original$n_cases[match(totals$state, risks$original$state)]
  if (any(totals$total != expected)) stop("Annual state counts do not match n_cases in table_04. Sample/threshold/window/outcome mismatch.")
  x[, state_order := match(state, states)]
  setorder(x, year, state_order)
  x[, state_order := NULL]
  x
}

# ============================================================
# Stage-two workflow
# ============================================================

read_stage1_package <- function(outcome_code = NULL, rr_file = NULL,
                                pipeline_index_file = PIPELINE_INDEX_FILE,
                                skip_if_missing = FALSE) {
  index <- NULL
  if (is.null(rr_file)) {
    index <- load_directory_index(pipeline_index_file)
    entry <- index$outcomes[[outcome_code]]
    if (is.null(entry)) {
      if (skip_if_missing) { message("Skipping unconfigured outcome: ", outcome_code); return(NULL) }
      stop("No input directory is configured for: ", outcome_code)
    }
    root <- normalizePath(entry$output_root, winslash = "/", mustWork = TRUE)
  } else {
    root <- normalizePath(dirname(rr_file), winslash = "/", mustWork = TRUE)
    if (basename(root) == "tables") root <- dirname(root)
  }
  required <- c("table_04_joint_9_state_rr.csv", "table_00c_annual_state_death_counts.csv",
                "source_parameters.rds", "table_02_joint_state_thresholds.csv")
  invisible(lapply(required, function(name) current_file(root, name)))
  for (name in c("joint_logrr_covariance.rds", "table_04a_joint_logrr_covariance.csv", "case_states.rds"))
    current_file(root, name, required = FALSE)
  params <- readRDS(current_file(root, "source_parameters.rds"))
  if (is.null(outcome_code)) outcome_code <- params$outcome_code
  ready <- list(protocol = "daynight_two_stage_v1", outcome_code = outcome_code,
                outcome_name = if (is.null(params$outcome)) outcome_code else params$outcome,
                stage1_id = if (is.null(params$stage1_id)) outcome_code else params$stage1_id,
                calibration_id = params$calibration_id,
                analysis_tag = if (is.null(params$shared_design$analysis_tag)) basename(root) else params$shared_design$analysis_tag,
                input_directory = file.path(root, "tables"), output_root = root)
  list(root = root, folder = file.path(root, "tables"), rr_file = current_file(root, "table_04_joint_9_state_rr.csv"),
       ready = ready, metadata = params, index = index,
       pipeline_index_file = pipeline_index_file)
}

read_working_covariance <- function(package, risks, source = RUN_COVARIANCE_SOURCE) {
  source <- match.arg(source, c("auto", "rescale", "csv", "rds"))
  states <- risks$estimable
  read_matrix <- function(file) {
    if (tolower(tools::file_ext(file)) == "csv") {
      x <- fread(file)
      labels <- intersect(c("state", "term"), names(x))
      if (length(labels) != 1L) stop("Covariance CSV requires one state or term column.")
      row_names <- x[[labels]]; x[, (labels) := NULL]
      if (!all(vapply(x, is.numeric, logical(1)))) stop("Covariance entries must be numeric.")
      V <- as.matrix(x); rownames(V) <- row_names
    } else {
      x <- readRDS(file)
      if (is.matrix(x)) V <- x else {
        V <- x$covariance_non_reference
        if (is.null(V)) V <- x$covariance
        if (is.null(V)) V <- x$covariance_state
      }
    }
    if (!is.matrix(V) || !is.numeric(V) || is.null(rownames(V)) || is.null(colnames(V)))
      stop("The covariance must be a named numeric matrix.")
    rownames(V) <- sub("^joint_state", "", rownames(V)); colnames(V) <- sub("^joint_state", "", colnames(V))
    if (anyDuplicated(rownames(V)) || anyDuplicated(colnames(V)) ||
        !all(states %in% rownames(V)) || !all(states %in% colnames(V)))
      stop("The covariance does not cover the estimable states.")
    V <- V[states, states, drop = FALSE]
    if (any(!is.finite(V))) stop("The covariance contains non-finite entries.")
    if (length(states)) {
      scale <- max(abs(V), .Machine$double.eps)
      if (max(abs(V - t(V))) > 1e-8 * scale) stop("The covariance is not symmetric.")
      V <- (V + t(V)) / 2
      if (min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < -1e-8 * scale)
        stop("The covariance is not positive semidefinite.")
    }
    V
  }
  rds_name <- "joint_logrr_covariance.rds"; csv_name <- "table_04a_joint_logrr_covariance.csv"
  can_rescale <- source %in% c("auto", "rescale")
  if (source == "auto") {
    candidates <- c(csv = analysis_file(package$root, csv_name), rds = analysis_file(package$root, rds_name))
    candidates <- candidates[file.exists(candidates)]
    if (!length(candidates)) stop("A full coefficient covariance matrix is required.")
    source <- names(candidates)[which.max(as.numeric(file.info(candidates)$mtime))]
  }
  V <- read_matrix(analysis_file(package$root, if (source == "csv") csv_name else rds_name))
  se <- risks$se[states]
  if (can_rescale && length(states)) {
    old_se <- sqrt(diag(V))
    if (any(old_se == 0 & se != 0)) stop("Cannot recover correlations for a zero-variance coefficient; provide a full covariance.")
    ratio <- ifelse(old_se == 0, 1, se / old_se)
    V <- V * outer(ratio, ratio)
  }
  if (length(states) && any(abs(diag(V) - se^2) > 1e-7 * pmax(1e-12, se^2)))
    stop("The supplied full covariance diagonal disagrees with SE squared in the risk table.")
  list(matrix = V, model = NULL, source = source)
}

save_working_inputs <- function(package, risks, covariance) {
  root <- package$root
  states <- risks$estimable
  saveRDS(list(mu_non_reference = risks$beta[states], covariance_non_reference = covariance$matrix),
          analysis_file(root, "joint_logrr_covariance.rds"))
  fwrite(cbind(data.table(state = states), as.data.table(covariance$matrix)),
         analysis_file(root, "table_04a_joint_logrr_covariance.csv"))
  package$metadata$stage1_input_directory <- file.path(root, "tables")
  package$metadata$output_root <- root
  saveRDS(package$metadata, analysis_file(root, "source_parameters.rds"))
  package
}

run_stage2_outcome <- function(outcome_code = NULL, rr_file = NULL,
                               pipeline_index_file = PIPELINE_INDEX_FILE,
                               skip_if_missing = FALSE,
                               n_sim = RUN_N_SIM, sim_seed = RUN_SIM_SEED,
                               cap_negative_af = RUN_CAP_NEGATIVE_AF,
                               covariance_source = RUN_COVARIANCE_SOURCE) {
  package <- read_stage1_package(outcome_code, rr_file, pipeline_index_file, skip_if_missing)
  if (is.null(package)) return(invisible(NULL))
  metadata <- package$metadata
  ready <- package$ready
  bundle <- list(directory = package$folder, rr_file = package$rr_file,
                 covariance_file = analysis_file(package$root, "joint_logrr_covariance.rds"),
                 annual_counts_file = file.path(package$folder, "table_00c_annual_state_death_counts.csv"),
                 metadata_file = analysis_file(package$root, "source_parameters.rds"),
                 threshold_file = file.path(package$folder, "table_02_joint_state_thresholds.csv"))
  risks <- load_rr_input(bundle$rr_file)
  covariance <- read_working_covariance(package, risks, covariance_source)
  annual <- normalise_annual_counts(fread(bundle$annual_counts_file), risks)
  n_sim <- if (is.null(n_sim)) metadata$n_sim else n_sim
  sim_seed <- if (is.null(sim_seed)) metadata$sim_seed else sim_seed
  cap_negative_af <- if (is.null(cap_negative_af)) metadata$cap_negative_af else cap_negative_af
  if (!is.numeric(n_sim) || length(n_sim) != 1L || !is.finite(n_sim) || n_sim < 2 || n_sim != floor(n_sim)) stop("n_sim must be an integer >= 2.")
  if (!is.numeric(sim_seed) || length(sim_seed) != 1L || !is.finite(sim_seed) || sim_seed < 0 || sim_seed > .Machine$integer.max || sim_seed != floor(sim_seed)) stop("Invalid Monte Carlo seed.")
  if (!is.logical(cap_negative_af) || length(cap_negative_af) != 1L || is.na(cap_negative_af)) stop("Invalid negative-AF setting.")
  if (!metadata$object_format %in% c("qs", "rds")) stop("Invalid stage-one object format.")
  if (metadata$object_format == "qs" && !requireNamespace("qs", quietly = TRUE)) stop("Install qs to retain stage one's object output format.")
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
  for (p in owned_markers[file.exists(owned_markers)])
    if (!file.remove(p)) stop("Cannot invalidate stage-two completion marker: ", p)
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
  compact <- if (file.exists(compact_file)) read_result_object(compact_file) else list()
  compact$coefficients <- coef_joint
  compact$covariance <- vcov_joint
  compact$beta_state <- beta_state
  compact$mu_non_reference <- risks$beta[estimable_non_ref_states]
  compact$covariance_non_reference <- covariance$matrix
  compact$covariance_state <- matrix(NA_real_, length(STATE_LEVELS), length(STATE_LEVELS),
                                     dimnames = list(STATE_LEVELS, STATE_LEVELS))
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
    if (!file.exists(analysis_file(package$root, "case_states.rds"))) stop("Stage one did not save required case-level input.")
    case_dt <- as.data.table(readRDS(analysis_file(package$root, "case_states.rds")))
    if (!all(c("id", "date", "year", "joint_state", "joint_mask") %in% names(case_dt)) ||
        nrow(case_dt) != sum(annual$n_cases) || anyNA(case_dt$year) ||
        any(!case_dt$joint_state %in% STATE_LEVELS) || anyDuplicated(case_dt, by = c("id", "date")))
      stop("Invalid saved case-state input.")
    actual <- case_dt[, .(n_cases = .N), by = .(year, state = as.character(joint_state))]
    compare <- merge(annual[, .(year, state, expected = n_cases)], actual, by = c("year", "state"), all.x = TRUE)
    compare[is.na(n_cases), n_cases := 0L]
    if (any(compare$expected != compare$n_cases)) stop("Case-level and annual counts disagree.")
    case_dt <- case_dt[, .(id, date, year, joint_state, joint_mask)]
  }
  
  plot_rr <- copy(rr_table)
  finite_rr <- plot_rr[is.finite(rr), rr]
  fill_limits <- range(finite_rr, na.rm = TRUE)
  if (!is.finite(fill_limits[1]) || !is.finite(fill_limits[2]) || fill_limits[1] == fill_limits[2]) {
    fill_limits <- c(0.9, 1.1)
  }
  
  p_heatmap <- ggplot(plot_rr, aes(x = day_axis, y = night_axis, fill = rr)) +
    geom_tile(colour = "white", size = 0.35) +
    geom_text(aes(label = rr_label), size = 2.6, family = FONT_FAMILY, lineheight = 0.9) +
    scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), limits = fill_limits, na.value = "white", name = "RR") +
    labs(x = "Daytime state", y = "Nighttime state") +
    coord_equal() +
    theme_nature(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), legend.position = "right")
  
  save_dt(plot_rr, file.path(DIR_PLOT_DATA, "plotdata_joint_3x3_rr_heatmap.csv"))
  save_obj(plot_rr, file.path(DIR_PLOT_DATA, "plotdata_joint_3x3_rr_heatmap.qs"))
  ggsave_out(p_heatmap, "fig_01_joint_3x3_rr_heatmap", width = 5.6, height = 4.8)
  
  log_msg("Computing observed exact nine-state attributable burden...")
  state_obs <- compute_state_burden_tables(beta_state, counts_by_year, n_deaths_by_year, STATE_INFO)
  
  rr_cols <- rr_table[, .(state, rr, rr_low, rr_high, support_status)]
  state_obs$annual <- merge(state_obs$annual, rr_cols, by = "state", all.x = TRUE, sort = FALSE)
  state_obs$mean_annual <- merge(state_obs$mean_annual, rr_cols, by = "state", all.x = TRUE, sort = FALSE)
  
  log_msg("Computing observed four-component Shapley burden from the nine-state model...")
  shapley_obs <- compute_shapley_burden_tables(
    beta_state,
    counts_by_year,
    n_deaths_by_year,
    STATE_INFO
  )
  
  phi_observed <- as.data.table(
    shapley_obs$phi_by_state,
    keep.rownames = "state"
  )
  
  missing_phi_columns <- setdiff(COMPONENTS, names(phi_observed))
  if (length(missing_phi_columns) > 0L) {
    stop(
      "The state-specific Shapley matrix is missing component columns: ",
      paste(missing_phi_columns, collapse = ", ")
    )
  }
  
  phi_columns <- paste0("phi_", COMPONENTS)
  setnames(phi_observed, COMPONENTS, phi_columns)
  
  phi_observed <- merge(
    phi_observed,
    STATE_INFO,
    by = "state",
    all.x = TRUE,
    sort = FALSE
  )
  phi_observed[, state_af_value := safe_af_from_beta(beta_state[state])]
  phi_observed[state == REF_STATE, state_af_value := 0]
  phi_observed[, shapley_sum := rowSums(.SD), .SDcols = phi_columns]
  phi_observed[, additivity_error := shapley_sum - state_af_value]
  setorder(phi_observed, mask)
  save_dt(
    phi_observed,
    file.path(DIR_TABLE, "table_07_state_specific_shapley_values.csv")
  )
  
  log_msg("Starting Monte Carlo uncertainty estimation for state and Shapley burden. N_SIM:", N_SIM)
  set.seed(SIM_SEED)
  
  non_ref_states <- setdiff(STATE_LEVELS, REF_STATE)
  term_names_all <- paste0("joint_state", non_ref_states)
  names(term_names_all) <- non_ref_states
  
  estimable_states <- estimable_non_ref_states
  estimable_terms <- term_names_all[estimable_states]
  
  if (length(unobserved_states) > 0L) {
    log_msg(
      "Monte Carlo simulation excludes structurally unobserved coefficient(s):",
      paste(unobserved_states, collapse = ", ")
    )
  }
  
  mu <- coef_joint[estimable_terms]
  Sigma <- vcov_joint[estimable_terms, estimable_terms, drop = FALSE]
  
  if (length(mu) > 0 && (any(!is.finite(mu)) || any(!is.finite(Sigma)))) {
    stop("Non-finite coefficient or covariance values were found in the joint model.")
  }
  
  if (length(mu) > 0) {
    draw_terms <- rmvnorm_eigen(N_SIM, mu, Sigma)
    colnames(draw_terms) <- estimable_states
  } else {
    draw_terms <- matrix(nrow = N_SIM, ncol = 0)
  }
  
  sim_state_annual_list <- vector("list", N_SIM)
  sim_state_mean_list <- vector("list", N_SIM)
  sim_shapley_annual_list <- vector("list", N_SIM)
  sim_shapley_mean_list <- vector("list", N_SIM)
  
  for (ss in seq_len(N_SIM)) {
    if (ss %% 100 == 0) log_msg("Monte Carlo simulation:", ss, "of", N_SIM)
    
    beta_draw <- setNames(rep(NA_real_, length(STATE_LEVELS)), STATE_LEVELS)
    beta_draw[REF_STATE] <- 0
    beta_draw[estimable_states] <- draw_terms[ss, estimable_states]
    
    state_sim <- compute_state_burden_tables(
      beta_draw,
      counts_by_year,
      n_deaths_by_year,
      STATE_INFO
    )
    shapley_sim <- compute_shapley_burden_tables(
      beta_draw,
      counts_by_year,
      n_deaths_by_year,
      STATE_INFO
    )
    
    sim_state_annual_list[[ss]] <- state_sim$annual[, .(
      sim = ss, year, state, AN, AF_percent, share_percent
    )]
    sim_state_mean_list[[ss]] <- state_sim$mean_annual[, .(
      sim = ss, period, state, AN, AF_percent, share_percent
    )]
    sim_shapley_annual_list[[ss]] <- shapley_sim$annual[, .(
      sim = ss, year, component, AN, AF_percent, share_percent
    )]
    sim_shapley_mean_list[[ss]] <- shapley_sim$mean_annual[, .(
      sim = ss, period, component, AN, AF_percent, share_percent
    )]
  }
  
  sim_state_annual <- rbindlist(sim_state_annual_list, use.names = TRUE, fill = TRUE)
  sim_state_mean <- rbindlist(sim_state_mean_list, use.names = TRUE, fill = TRUE)
  sim_shapley_annual <- rbindlist(sim_shapley_annual_list, use.names = TRUE, fill = TRUE)
  sim_shapley_mean <- rbindlist(sim_shapley_mean_list, use.names = TRUE, fill = TRUE)
  rm(
    sim_state_annual_list,
    sim_state_mean_list,
    sim_shapley_annual_list,
    sim_shapley_mean_list,
    draw_terms
  )
  gc()
  
  state_annual_ci <- add_ci(
    obs = state_obs$annual,
    sims = sim_state_annual,
    by_cols = c("year", "state"),
    value_cols = c("AN", "AF_percent", "share_percent")
  )
  
  state_mean_ci <- add_ci(
    obs = state_obs$mean_annual,
    sims = sim_state_mean,
    by_cols = c("period", "state"),
    value_cols = c("AN", "AF_percent", "share_percent")
  )
  
  shapley_annual_ci <- add_ci(
    obs = shapley_obs$annual,
    sims = sim_shapley_annual,
    by_cols = c("year", "component"),
    value_cols = c("AN", "AF_percent", "share_percent")
  )
  
  shapley_mean_ci <- add_ci(
    obs = shapley_obs$mean_annual,
    sims = sim_shapley_mean,
    by_cols = c("period", "component"),
    value_cols = c("AN", "AF_percent", "share_percent")
  )
  
  setorder(state_annual_ci, year, mask)
  setorder(state_mean_ci, mask)
  setorder(shapley_annual_ci, year, component_order)
  setorder(shapley_mean_ci, component_order)
  
  for (object_name in c("state_annual_ci", "state_mean_ci", "shapley_annual_ci", "shapley_mean_ci")) {
    object_value <- get(object_name)
    object_value[, `:=`(
      AN_CI = sprintf("%.3f (%.3f, %.3f)", AN, AN_low, AN_high),
      AF_percent_CI = sprintf("%.3f (%.3f, %.3f)", AF_percent, AF_percent_low, AF_percent_high),
      share_percent_CI = sprintf("%.3f (%.3f, %.3f)", share_percent, share_percent_low, share_percent_high)
    )]
    assign(object_name, object_value)
  }
  
  state_share_check_annual <- state_annual_ci[, .(
    state_share_sum = if (any(is.finite(share_percent))) sum(share_percent, na.rm = TRUE) else NA_real_
  ), by = year]
  state_share_check_mean <- state_mean_ci[, .(
    state_share_sum = if (any(is.finite(share_percent))) sum(share_percent, na.rm = TRUE) else NA_real_
  ), by = period]
  
  shapley_share_check_annual <- shapley_annual_ci[
    component %in% COMPONENTS,
    .(component_share_sum = if (any(is.finite(share_percent))) sum(share_percent, na.rm = TRUE) else NA_real_),
    by = year
  ]
  shapley_share_check_mean <- shapley_mean_ci[
    component %in% COMPONENTS,
    .(component_share_sum = if (any(is.finite(share_percent))) sum(share_percent, na.rm = TRUE) else NA_real_),
    by = period
  ]
  
  shapley_additivity_annual <- merge(
    state_annual_ci[, .(state_joint_AN = sum(AN, na.rm = TRUE)), by = year],
    shapley_annual_ci[component == "joint_total", .(year, shapley_joint_AN = AN)],
    by = "year",
    all = TRUE
  )
  shapley_additivity_annual[, difference := shapley_joint_AN - state_joint_AN]
  
  shapley_additivity_mean <- data.table(
    period = "mean_annual",
    state_joint_AN = state_mean_ci[, sum(AN, na.rm = TRUE)],
    shapley_joint_AN = shapley_mean_ci[component == "joint_total", AN]
  )
  shapley_additivity_mean[, difference := shapley_joint_AN - state_joint_AN]
  
  if (any(abs(state_share_check_annual$state_share_sum - 100) > 1e-6, na.rm = TRUE) ||
      any(abs(state_share_check_mean$state_share_sum - 100) > 1e-6, na.rm = TRUE)) {
    stop("The nine-state burden shares failed the 100% sum check.")
  }
  if (any(abs(shapley_share_check_annual$component_share_sum - 100) > 1e-6, na.rm = TRUE) ||
      any(abs(shapley_share_check_mean$component_share_sum - 100) > 1e-6, na.rm = TRUE)) {
    stop("The four-component Shapley burden shares failed the 100% sum check.")
  }
  if (any(abs(shapley_additivity_annual$difference) > 1e-8, na.rm = TRUE) ||
      any(abs(shapley_additivity_mean$difference) > 1e-8, na.rm = TRUE)) {
    stop("The Shapley burden totals failed the additivity check against the nine-state burden.")
  }
  
  save_dt(state_share_check_annual, file.path(DIR_TABLE, "table_05_state_share_sum_check_annual.csv"))
  save_dt(state_share_check_mean, file.path(DIR_TABLE, "table_06_state_share_sum_check_mean_annual.csv"))
  save_dt(shapley_share_check_annual, file.path(DIR_TABLE, "table_08_shapley_share_sum_check_annual.csv"))
  save_dt(shapley_share_check_mean, file.path(DIR_TABLE, "table_08_shapley_share_sum_check_mean_annual.csv"))
  save_dt(shapley_additivity_annual, file.path(DIR_TABLE, "table_08_shapley_additivity_check_annual.csv"))
  save_dt(shapley_additivity_mean, file.path(DIR_TABLE, "table_08_shapley_additivity_check_mean_annual.csv"))
  
  state_annual_out <- state_annual_ci[, .(
    year, state, mask, component_mask, D_C, D_H, N_C, N_H,
    day_code, night_code, day_status, night_status, support_status,
    n_cases, n_deaths, rr, rr_low, rr_high,
    AN, AN_low, AN_high, AN_CI,
    AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI,
    share_percent, share_percent_low, share_percent_high, share_percent_CI
  )]
  
  state_mean_out <- state_mean_ci[, .(
    period, state, mask, component_mask, D_C, D_H, N_C, N_H,
    day_code, night_code, day_status, night_status, support_status,
    n_cases, n_deaths, rr, rr_low, rr_high,
    AN, AN_low, AN_high, AN_CI,
    AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI,
    share_percent, share_percent_low, share_percent_high, share_percent_CI
  )]
  
  shapley_annual_out <- shapley_annual_ci[, .(
    year, component, component_order, component_label, component_group,
    n_deaths,
    AN, AN_low, AN_high, AN_CI,
    AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI,
    share_percent, share_percent_low, share_percent_high, share_percent_CI
  )]
  
  shapley_mean_out <- shapley_mean_ci[, .(
    period, component, component_order, component_label, component_group,
    n_deaths,
    AN, AN_low, AN_high, AN_CI,
    AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI,
    share_percent, share_percent_low, share_percent_high, share_percent_CI
  )]
  
  save_dt(state_annual_out, file.path(DIR_TABLE, "table_09_exact_9_state_burden_annual.csv"))
  save_dt(state_mean_out, file.path(DIR_TABLE, "table_10_exact_9_state_burden_mean_annual.csv"))
  save_dt(shapley_annual_out, file.path(DIR_TABLE, "table_11_shapley_4_component_burden_annual.csv"))
  save_dt(shapley_mean_out, file.path(DIR_TABLE, "table_12_shapley_4_component_burden_mean_annual.csv"))
  
  save_obj(state_annual_out, file.path(DIR_MODEL, "exact_9_state_burden_annual.qs"))
  save_obj(state_mean_out, file.path(DIR_MODEL, "exact_9_state_burden_mean_annual.qs"))
  save_obj(shapley_annual_out, file.path(DIR_MODEL, "shapley_4_component_burden_annual.qs"))
  save_obj(shapley_mean_out, file.path(DIR_MODEL, "shapley_4_component_burden_mean_annual.qs"))
  
  state_summary <- state_mean_out[, .(
    result_type = "exact_9_state",
    period,
    item = state,
    item_label = paste(day_status, "|", night_status),
    AN, AN_low, AN_high, AN_CI,
    AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI,
    share_percent, share_percent_low, share_percent_high, share_percent_CI
  )]
  
  shapley_summary <- shapley_mean_out[, .(
    result_type = "shapley_4_component",
    period,
    item = component,
    item_label = component_label,
    AN, AN_low, AN_high, AN_CI,
    AF_percent, AF_percent_low, AF_percent_high, AF_percent_CI,
    share_percent, share_percent_low, share_percent_high, share_percent_CI
  )]
  
  mean_summary_requested <- rbindlist(
    list(state_summary, shapley_summary),
    use.names = TRUE,
    fill = TRUE
  )
  save_dt(mean_summary_requested, file.path(DIR_TABLE, "table_13_requested_mean_annual_summary.csv"))
  save_obj(mean_summary_requested, file.path(DIR_MODEL, "requested_mean_annual_summary.qs"))
  
  shapley_plot_data <- copy(shapley_mean_out[component %in% COMPONENTS])
  shapley_plot_data[, component_label := factor(
    component_label,
    levels = COMPONENT_INFO[component %in% COMPONENTS, component_label]
  )]
  save_dt(shapley_plot_data, file.path(DIR_PLOT_DATA, "plotdata_shapley_4_component_burden_share.csv"))
  save_obj(shapley_plot_data, file.path(DIR_PLOT_DATA, "plotdata_shapley_4_component_burden_share.qs"))
  
  component_colours <- c(
    "Daytime cold" = "#2C7FB8",
    "Daytime heat" = "#F4A582",
    "Nighttime cold" = "#7FCDBB",
    "Nighttime heat" = "#C51B1D"
  )
  
  if (all(is.finite(shapley_plot_data$share_percent)) &&
      all(shapley_plot_data$share_percent >= 0) &&
      sum(shapley_plot_data$share_percent) > 0) {
    p_shapley <- ggplot(shapley_plot_data, aes(x = "", y = share_percent, fill = component_label)) +
      geom_col(width = 1, colour = "white", size = 0.35) +
      coord_polar(theta = "y") +
      geom_text(
        aes(label = sprintf("%.1f%%", share_percent)),
        position = position_stack(vjust = 0.5),
        size = 3.0,
        family = FONT_FAMILY
      ) +
      scale_fill_manual(values = component_colours, name = NULL) +
      labs(x = NULL, y = NULL) +
      theme_void(base_family = FONT_FAMILY) +
      theme(
        legend.position = "top",
        legend.text = element_text(size = 9, colour = "black"),
        plot.margin = margin(5, 5, 5, 5)
      )
    ggsave_out(
      p_shapley,
      "fig_02_shapley_4_component_burden_share",
      width = 5.2,
      height = 4.5
    )
  } else {
    for (extension in c("png", "pdf")) {
      old_figure <- file.path(DIR_FIG, paste0("fig_02_shapley_4_component_burden_share.", extension))
      if (file.exists(old_figure) && !file.remove(old_figure)) stop("Cannot remove an obsolete share figure: ", old_figure)
    }
    log_msg(
      "The Shapley share figure was not generated because one or more component shares were negative or non-finite."
    )
  }
  
  rm(
    sim_state_annual,
    sim_state_mean,
    sim_shapley_annual,
    sim_shapley_mean
  )
  gc()
  
  if (SAVE_CASE_CONTRIB) {
    log_msg("Saving individual case-level state AF contributions...")
    af_state_obs <- safe_af_from_beta(beta_state)
    af_state_obs[REF_STATE] <- 0
    contrib <- copy(case_dt)
    contrib[, state_af := af_state_obs[joint_state]]
    save_obj(contrib, file.path(DIR_MODEL, "case_level_state_af_contributions.qs"))
    save_dt(contrib, file.path(DIR_CONTRIB, "case_level_state_af_contributions.csv"))
  }
  
  analysis_meta <- metadata
  analysis_meta$output_root <- OUT_ROOT
  analysis_meta$n_sim <- N_SIM
  analysis_meta$sim_seed <- SIM_SEED
  analysis_meta$cap_negative_af <- CAP_NEGATIVE_AF
  analysis_meta$shapley_decomposition <- list(
    method = "State-specific Shapley allocation derived from the nine-state model",
    components = COMPONENTS, component_labels = COMPONENT_INFO, phi_by_state = phi_observed)
  analysis_meta$share_checks <- list(state_annual = state_share_check_annual, state_mean = state_share_check_mean,
                                     shapley_annual = shapley_share_check_annual, shapley_mean = shapley_share_check_mean,
                                     shapley_additivity_annual = shapley_additivity_annual, shapley_additivity_mean = shapley_additivity_mean)
  analysis_meta$stage2 <- list(stage1_id = ready$stage1_id, input_directory = bundle$directory,
                               no_model_refit = TRUE,
                               completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  save_obj(analysis_meta, file.path(DIR_MODEL, "analysis_metadata.qs"))
  writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, "sessionInfo_stage2.txt"))
  complete <- list(outcome = OUTCOME_CODE, stage1_id = ready$stage1_id,
                   completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  saveRDS(complete, file.path(DIR_MODEL, "stage2_completed.rds"))
  save_obj(complete, file.path(DIR_MODEL, "run_completed.qs"))
  log_msg("Both stages now form the complete analysis output:", OUT_ROOT)
  invisible(list(output_root = OUT_ROOT, rr_table = rr_table, annual_counts = annual,
                 state_annual_out = state_annual_out, state_mean_out = state_mean_out,
                 shapley_annual_out = shapley_annual_out, shapley_mean_out = shapley_mean_out,
                 phi_observed = phi_observed, analysis_meta = analysis_meta))
}

# if (!isTRUE(getOption("daynight.stage2.skip_calls", FALSE)) && RUN_OVERALL) {
#   analysis_result <- run_stage2_outcome("I00_I52_I60_I69", skip_if_missing = TRUE)
# }
# if (!isTRUE(getOption("daynight.stage2.skip_calls", FALSE)) && RUN_I10_I15) {
#   analysis_result <- run_stage2_outcome("I10_I15", skip_if_missing = TRUE)
# }
# if (!isTRUE(getOption("daynight.stage2.skip_calls", FALSE)) && RUN_I20_I25) {
#   analysis_result <- run_stage2_outcome("I20_I25", skip_if_missing = TRUE)
# }
# if (!isTRUE(getOption("daynight.stage2.skip_calls", FALSE)) && RUN_I50) {
#   analysis_result <- run_stage2_outcome("I50", skip_if_missing = TRUE)
# }
# if (!isTRUE(getOption("daynight.stage2.skip_calls", FALSE)) && RUN_I60_I62) {
#   analysis_result <- run_stage2_outcome("I60_I62", skip_if_missing = TRUE)
# }
if (!isTRUE(getOption("daynight.stage2.skip_calls", FALSE)) && RUN_I63) {
  analysis_result <- run_stage2_outcome("I63", skip_if_missing = TRUE)
}

