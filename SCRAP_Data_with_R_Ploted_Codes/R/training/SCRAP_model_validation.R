# SCRAP model validation and robustness diagnostics

rm(list = ls())
graphics.off()
set.seed(123)

required_packages <- c("tidyverse", "patchwork", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

# Explicit aliases prevent select()/filter()/rename() masking in long RStudio sessions.
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
rename <- dplyr::rename
left_join <- dplyr::left_join
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
relocate <- dplyr::relocate
transmute <- dplyr::transmute
pull <- dplyr::pull
count <- dplyr::count
if_else <- dplyr::if_else
case_when <- dplyr::case_when
all_of <- tidyselect::all_of

theme_nature <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "#1F1F1F"),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.08), color = "#111111"),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.88), color = "#333333", margin = ggplot2::margin(b = 4)),
      axis.title = ggplot2::element_text(color = "#222222"),
      axis.text = ggplot2::element_text(color = "#222222"),
      axis.line = ggplot2::element_line(color = "#222222", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(color = "#222222", linewidth = 0.28),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.90)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.84)),
      legend.key.size = grid::unit(0.35, "cm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", color = "#222222"),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

ggplot2::theme_set(theme_nature(base_size = 8.5))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && !is.null(ctx$path) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path, winslash = "/", mustWork = FALSE)))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
project_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

input_file <- file.path(project_root, "data", "training", "Training_segment_data.csv")
out_dir <- file.path(project_root, "results", "tables", "training", "SCRAP_model_validation")
fig_dir <- file.path(project_root, "results", "figures", "training", "SCRAP_model_validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

compatible_cutoff_nm <- 30
lambda_equal_weight <- 0.20

thresholds <- tribble(
  ~param,       ~term,       ~min_val, ~max_val, ~module,
  "Pro-Pol",   "Pro-Pol",   -38.1,    -18.0,    "Carrier integration",
  "Par-Pol",   "Par-Pol",   -19.0,      0.0,    "Carrier integration",
  "Mod-Pol",   "Mod-Pol",      NA,     -4.0,    "Carrier integration",
  "Pro-PCL",   "Pro-PCL",   -31.0,    -15.0,    "PCL anchoring",
  "Par-PCL",   "Par-PCL",   -15.6,      0.0,    "PCL anchoring",
  "Mod-PCL",   "Mod-PCL",      NA,     -3.5,    "PCL anchoring",
  "Pro-PEG",  "Pro-PEG",      NA,     -3.3,    "PEG-associated contacts",
  "Par-PEG",   "Par-PEG",    -6.5,      0.0,    "PEG-associated contacts",
  "Mod-PEG",   "Mod-PEG",    -7.0,      0.0,    "PEG-associated contacts",
  "Pro-Pro", "Pro-Pro",   -28.0,     -6.1,    "Prodrug self-association",
  "Par-Par",   "Par-Par",   -11.0,      0.0,    "Prodrug self-association",
  "Mod-Mod",   "Mod-Mod",      NA,     -0.42,   "Prodrug self-association",
  "Par-Mod",   "Par-Mod",    -9.0,     -1.0,    "Prodrug self-association"
) %>%
  dplyr::mutate(
    rule_type = case_when(
      is.na(min_val) ~ "max_only",
      !is.na(min_val) & !is.na(max_val) & max_val == 0 ~ "min_to_zero",
      TRUE ~ "interval"
    ),
    param = factor(param, levels = param)
  )

param_order <- as.character(thresholds$param)

score_pass <- function(x, min_val, max_val) {
  ok <- !is.na(x)
  if (!is.na(min_val)) ok <- ok & x >= min_val
  if (!is.na(max_val)) ok <- ok & x <= max_val
  as.integer(ok)
}

window_distance <- function(x, min_val, max_val) {
  d <- rep(0, length(x))
  if (!is.na(min_val)) d <- pmax(d, min_val - x, na.rm = FALSE)
  if (!is.na(max_val)) d <- pmax(d, x - max_val, na.rm = FALSE)
  d[is.na(x)] <- NA_real_
  width <- case_when(
    !is.na(min_val) & !is.na(max_val) & max_val != 0 ~ abs(max_val - min_val),
    !is.na(min_val) ~ max(abs(min_val), 1),
    !is.na(max_val) ~ max(abs(max_val), 1),
    TRUE ~ 1
  )
  d / width
}

calc_metrics <- function(actual, pred) {
  actual <- as.integer(actual)
  pred <- as.integer(pred)
  keep <- !is.na(actual) & !is.na(pred)
  actual <- actual[keep]
  pred <- pred[keep]
  tp <- sum(actual == 1 & pred == 1)
  tn <- sum(actual == 0 & pred == 0)
  fp <- sum(actual == 0 & pred == 1)
  fn <- sum(actual == 1 & pred == 0)
  n <- length(actual)
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  f1 <- ifelse(!is.na(precision) & !is.na(sensitivity) & (precision + sensitivity) > 0, 2 * precision * sensitivity / (precision + sensitivity), NA_real_)
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(denom > 0, ((tp * tn) - (fp * fn)) / denom, 0)
  tibble(n = n, TP = tp, TN = tn, FP = fp, FN = fn, accuracy = (tp + tn) / n, balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE), sensitivity = sensitivity, specificity = specificity, precision = precision, F1 = f1, MCC = mcc)
}

entropy_bin <- function(y) {
  y <- y[!is.na(y)]
  if (length(y) == 0) return(NA_real_)
  p <- mean(y == 1)
  if (p <= 0 || p >= 1) return(0)
  -p * log2(p) - (1 - p) * log2(1 - p)
}

information_gain <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  if (length(unique(x)) <= 1 || length(unique(y)) <= 1) return(0)
  h0 <- entropy_bin(y)
  h_cond <- 0
  for (lv in unique(x)) {
    idx <- x == lv
    h_cond <- h_cond + mean(idx) * entropy_bin(y[idx])
  }
  h0 - h_cond
}

norm01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || abs(rng[2] - rng[1]) < .Machine$double.eps) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

calc_auc <- function(score, actual) {
  keep <- !is.na(score) & !is.na(actual)
  score <- as.numeric(score[keep])
  actual <- as.integer(actual[keep])
  n_pos <- sum(actual == 1)
  n_neg <- sum(actual == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[actual == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

score_with_thresholds <- function(data, threshold_tbl) {
  long <- data %>%
    dplyr::select(.row_id, Lib, Mol, Sizes, outcome, outcome_bin, dplyr::all_of(param_order)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(param_order),
      names_to = "param",
      values_to = "energy"
    ) %>%
    dplyr::mutate(param = factor(param, levels = param_order)) %>%
    dplyr::left_join(
      threshold_tbl %>% dplyr::mutate(param = factor(param, levels = param_order)),
      by = "param"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      pass = score_pass(energy, min_val, max_val),
      distance_to_window = window_distance(energy, min_val, max_val)
    ) %>%
    dplyr::ungroup()
  
  score <- long %>%
    dplyr::group_by(.row_id, Lib, Mol, Sizes, outcome, outcome_bin) %>%
    dplyr::summarise(
      score_13 = sum(pass, na.rm = TRUE),
      failed_terms = paste(term[pass == 0], collapse = ";"),
      mean_distance_to_window = if (all(is.na(distance_to_window))) NA_real_ else mean(distance_to_window, na.rm = TRUE),
      max_distance_to_window = if (all(is.na(distance_to_window))) NA_real_ else max(distance_to_window, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      score_unweighted_100 = 100 * score_13 / length(param_order),
      pred_strict_13 = score_13 == length(param_order)
    )
  
  duplicate_keys <- long %>%
    dplyr::count(.row_id, param, name = "n") %>%
    dplyr::filter(n != 1)
  
  if (nrow(duplicate_keys) > 0) {
    stop(
      "Internal scoring error: duplicate row/parameter combinations were detected.",
      call. = FALSE
    )
  }
  
  wide_pass <- long %>%
    dplyr::select(.row_id, param, pass) %>%
    tidyr::pivot_wider(
      names_from = param,
      values_from = pass
    ) %>%
    dplyr::arrange(.row_id)
  
  list(long = long, score = score, wide_pass = wide_pass)
}

raw_data <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "na", "NaN")
) %>%
  dplyr::rename_with(~ gsub("^\\ufeff", "", .x))

if (!"Lib"%in% names(raw_data)) {
  raw_data$Lib <- "Training"
}

required_cols <- c("Lib", "Mol", "Sizes", param_order)
missing_cols <- setdiff(required_cols, names(raw_data))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# Remove trailing/blank CSV rows before any scoring.
raw_data <- raw_data %>%
  dplyr::mutate(
    Lib = as.character(Lib),
    Mol = trimws(as.character(Mol))
  ) %>%
  dplyr::filter(
    !is.na(Mol),
    Mol != "",
    toupper(Mol) != "NA"
  ) %>%
  dplyr::mutate(
    Lib = dplyr::if_else(is.na(Lib) | trimws(Lib) == "", "Training", Lib)
  )

if (nrow(raw_data) == 0) {
  stop("No valid formulation rows remain after removing blank rows.", call. = FALSE)
}

duplicate_mol <- raw_data %>%
  dplyr::count(Mol, name = "n") %>%
  dplyr::filter(n > 1)

if (nrow(duplicate_mol) > 0) {
  stop(
    "Duplicate molecule names were found in Training_segment_data.csv: ",
    paste(duplicate_mol$Mol, collapse = ", "),
    ". Each formulation must have one unique Mol identifier.",
    call. = FALSE
  )
}

raw_data <- raw_data %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(c("Sizes", param_order)),
      ~ suppressWarnings(as.numeric(.x))
    )
  )

# Every valid formulation must contain all 13 interaction descriptors.
missing_descriptor_rows <- raw_data %>%
  dplyr::filter(dplyr::if_any(dplyr::all_of(param_order), is.na)) %>%
  dplyr::pull(Mol)

if (length(missing_descriptor_rows) > 0) {
  stop(
    "Missing/non-numeric interaction descriptor(s) for: ",
    paste(unique(missing_descriptor_rows), collapse = ", "),
    call. = FALSE
  )
}

raw_data <- raw_data %>%
  dplyr::mutate(
    .row_id = dplyr::row_number(),
    outcome = dplyr::if_else(
      !is.na(Sizes) & Sizes > 0 & Sizes < compatible_cutoff_nm,
      "Compatible",
      "Incompatible"
    ),
    outcome = factor(outcome, levels = c("Incompatible", "Compatible")),
    outcome_bin = as.integer(outcome == "Compatible")
  )

message("Valid training formulations: ", nrow(raw_data))

scored <- score_with_thresholds(raw_data, thresholds)
score_long <- scored$long
score_table <- scored$score %>% dplyr::arrange(.row_id)

pass_by_row <- score_table %>%
  dplyr::select(.row_id, Mol) %>%
  dplyr::left_join(scored$wide_pass, by = ".row_id") %>%
  dplyr::arrange(.row_id)

for (p in param_order) {
  pass_by_row[[p]] <- as.integer(pass_by_row[[p]])
}
full_metrics <- calc_metrics(score_table$outcome_bin, score_table$pred_strict_13) %>% dplyr::mutate(model = "Strict 13/13 rule")

importance_tbl <- map_dfr(param_order, function(p) {
  tmp <- score_long %>% dplyr::filter(as.character(param) == p)
  pass_vec <- tmp$pass
  y <- tmp$outcome_bin
  compat_pass_rate <- mean(pass_vec[y == 1], na.rm = TRUE)
  incompat_pass_rate <- mean(pass_vec[y == 0], na.rm = TRUE)
  score_without <- score_table$score_13 - pass_by_row[[p]]
  pred_without <- score_without == (length(param_order) - 1)
  without_metrics <- calc_metrics(score_table$outcome_bin, pred_without)
  dist_auc <- calc_auc(-tmp$distance_to_window, y)
  tibble(param = p, compatible_pass_rate = compat_pass_rate, incompatible_pass_rate = incompat_pass_rate, failure_gap = compat_pass_rate - incompat_pass_rate, information_gain = information_gain(pass_vec, y), phi_abs = abs(calc_metrics(y, pass_vec)$MCC[[1]]), distance_auc = dist_auc, distance_auc_gain = pmax(dist_auc - 0.5, 0), ablation_MCC_drop = full_metrics$MCC[[1]] - without_metrics$MCC[[1]])
}) %>%
  dplyr::left_join(thresholds %>% dplyr::mutate(param = as.character(param)) %>% dplyr::select(param, term, module, rule_type, min_val, max_val), by = "param") %>%
  dplyr::mutate(norm_failure_gap = norm01(failure_gap), norm_information_gain = norm01(information_gain), norm_phi_abs = norm01(phi_abs), norm_distance_auc_gain = norm01(distance_auc_gain), norm_ablation_MCC_drop = norm01(pmax(ablation_MCC_drop, 0)), importance_index = rowMeans(cbind(norm_failure_gap, norm_information_gain, norm_phi_abs, norm_distance_auc_gain, norm_ablation_MCC_drop), na.rm = TRUE))

if (sum(importance_tbl$importance_index, na.rm = TRUE) <= 0) importance_tbl <- importance_tbl %>% dplyr::mutate(importance_index = 1)

importance_tbl <- importance_tbl %>%
  dplyr::mutate(weight_fraction = lambda_equal_weight / length(param_order) + (1 - lambda_equal_weight) * importance_index / sum(importance_index, na.rm = TRUE), weight_percent = 100 * weight_fraction) %>%
  dplyr::arrange(desc(importance_index), desc(failure_gap), desc(ablation_MCC_drop), param) %>%
  dplyr::mutate(importance_rank = row_number()) %>%
  dplyr::relocate(importance_rank, param, term, module, rule_type, min_val, max_val, importance_index, weight_fraction, weight_percent)

weight_vec <- importance_tbl %>% dplyr::select(param, weight_fraction) %>% deframe()
weight_vec <- weight_vec[param_order]
pass_matrix <- pass_by_row %>% dplyr::select(dplyr::all_of(param_order)) %>% as.matrix()
storage.mode(pass_matrix) <- "numeric"
weighted_score_table <- score_table %>%
  dplyr::mutate(score_weighted_100 = as.numeric(pass_matrix %*% as.numeric(weight_vec)) * 100, pred_weighted_strict = score_weighted_100 >= (100 - 1e-8))
weighted_metrics <- calc_metrics(weighted_score_table$outcome_bin, weighted_score_table$pred_weighted_strict) %>% dplyr::mutate(model = "Importance-weighted strict 100% rule")
auc_tbl <- tibble(score_type = c("13-point score", "Unweighted percent score", "Importance-weighted 100% score"), AUC = c(calc_auc(score_table$score_13, score_table$outcome_bin), calc_auc(score_table$score_unweighted_100, score_table$outcome_bin), calc_auc(weighted_score_table$score_weighted_100, weighted_score_table$outcome_bin)))

module_importance_tbl <- importance_tbl %>%
  group_by(module) %>%
  dplyr::summarise(n_terms = n(), total_weight_percent = sum(weight_percent, na.rm = TRUE), total_importance_index = sum(importance_index, na.rm = TRUE), mean_importance_index = mean(importance_index, na.rm = TRUE), top_terms = paste(term[order(-importance_index)][seq_len(min(3, n()))], collapse = ";"), .groups = "drop") %>%
  dplyr::arrange(desc(total_weight_percent))

module_terms <- split(thresholds$param %>% as.character(), thresholds$module)
module_ablation_tbl <- map_dfr(names(module_terms), function(m) {
  removed <- module_terms[[m]]
  kept <- setdiff(param_order, removed)
  score_without <- rowSums(pass_matrix[, kept, drop = FALSE], na.rm = TRUE)
  met <- calc_metrics(score_table$outcome_bin, score_without == length(kept))
  met %>% dplyr::mutate(module = m, removed_terms = paste(thresholds$term[as.character(thresholds$param) %in% removed], collapse = ";"), n_removed = length(removed), MCC_drop = full_metrics$MCC[[1]] - MCC, balanced_accuracy_drop = full_metrics$balanced_accuracy[[1]] - balanced_accuracy)
}) %>%
  dplyr::relocate(module, removed_terms, n_removed, MCC_drop, balanced_accuracy_drop)

module_failure_tbl <- score_long %>%
  dplyr::mutate(param = as.character(param), failed = pass == 0) %>%
  dplyr::left_join(importance_tbl %>% dplyr::select(param, weight_percent), by = "param") %>%
  group_by(Lib, Mol, Sizes, outcome, module) %>%
  dplyr::summarise(module_score = sum(pass, na.rm = TRUE), n_failed = sum(failed, na.rm = TRUE), module_weight_lost_percent = sum(if_else(failed, weight_percent, 0), na.rm = TRUE), failed_terms = paste(term[failed], collapse = ";"), .groups = "drop")

recalculate_thresholds <- function(train_data, base_thresholds) {
  train_compatible <- train_data %>% dplyr::filter(outcome_bin == 1)
  base_thresholds %>%
    rowwise() %>%
    dplyr::mutate(train_min = suppressWarnings(min(train_compatible[[as.character(param)]], na.rm = TRUE)), train_max = suppressWarnings(max(train_compatible[[as.character(param)]], na.rm = TRUE)), min_val = case_when(rule_type %in% c("interval", "min_to_zero") ~ train_min, TRUE ~ NA_real_), max_val = case_when(rule_type == "interval"~ train_max, rule_type == "max_only"~ train_max, rule_type == "min_to_zero"~ 0, TRUE ~ max_val)) %>%
    ungroup() %>%
    dplyr::select(param, term, min_val, max_val, module, rule_type)
}

loocv_results <- map_dfr(seq_len(nrow(raw_data)), function(i) {
  score_with_thresholds(raw_data[i, ], recalculate_thresholds(raw_data[-i, ], thresholds))$score %>%
    dplyr::mutate(fold = i, left_out = raw_data$Mol[[i]], pred_loocv_strict = score_13 == length(param_order)) %>%
    dplyr::select(fold, left_out, Lib, Mol, Sizes, outcome, outcome_bin, score_13, score_unweighted_100, failed_terms, pred_loocv_strict)
})
loocv_metrics <- calc_metrics(loocv_results$outcome_bin, loocv_results$pred_loocv_strict) %>% dplyr::mutate(model = "LOOCV with recalculated compatible-window thresholds")
loocv_thresholds <- map_dfr(seq_len(nrow(raw_data)), function(i) recalculate_thresholds(raw_data[-i, ], thresholds) %>% dplyr::mutate(fold = i, left_out = raw_data$Mol[[i]]))
threshold_stability <- loocv_thresholds %>%
  pivot_longer(cols = c(min_val, max_val), names_to = "boundary", values_to = "value") %>%
  dplyr::filter(!is.na(value)) %>%
  group_by(param, term, boundary) %>%
  dplyr::summarise(mean_value = mean(value, na.rm = TRUE), sd_value = sd(value, na.rm = TRUE), cv_percent = ifelse(abs(mean_value) > 1e-8, 100 * sd_value / abs(mean_value), NA_real_), .groups = "drop") %>%
  dplyr::left_join(thresholds %>% dplyr::mutate(param = as.character(param)) %>% pivot_longer(cols = c(min_val, max_val), names_to = "boundary", values_to = "full_threshold") %>% dplyr::filter(!is.na(full_threshold)) %>% dplyr::select(param, boundary, full_threshold), by = c("param", "boundary"))

perturb_thresholds <- function(base_thresholds, perturb_fraction) {
  base_thresholds %>%
    rowwise() %>%
    dplyr::mutate(center = ifelse(!is.na(min_val) & !is.na(max_val), (min_val + max_val) / 2, NA_real_), half_width = ifelse(!is.na(min_val) & !is.na(max_val), abs(max_val - min_val) / 2, NA_real_), min_new = case_when(rule_type == "interval"~ center - half_width * (1 + perturb_fraction), rule_type == "min_to_zero"~ min_val - abs(min_val) * perturb_fraction, TRUE ~ NA_real_), max_new = case_when(rule_type == "interval"~ center + half_width * (1 + perturb_fraction), rule_type == "max_only"~ max_val + abs(max_val) * perturb_fraction, rule_type == "min_to_zero"~ 0, TRUE ~ max_val)) %>%
    ungroup() %>%
    transmute(param, term, min_val = min_new, max_val = max_new, module, rule_type)
}

threshold_perturbation_tbl <- map_dfr(seq(-0.20, 0.20, by = 0.05), function(pct) {
  sc <- score_with_thresholds(raw_data, perturb_thresholds(thresholds, pct))$score %>% dplyr::mutate(pred = score_13 == length(param_order))
  calc_metrics(sc$outcome_bin, sc$pred) %>% dplyr::mutate(perturbation_fraction = pct, perturbation_percent = 100 * pct)
})

correlation_tbl <- suppressWarnings(cor(raw_data %>% dplyr::select(dplyr::all_of(param_order)), method = "spearman", use = "pairwise.complete.obs")) %>%
  as.data.frame() %>%
  rownames_to_column("param_x") %>%
  pivot_longer(-param_x, names_to = "param_y", values_to = "spearman_rho")

pca_input <- raw_data %>% dplyr::select(dplyr::all_of(param_order)) %>% as.data.frame()
pca_complete <- complete.cases(pca_input)
pca_fit <- prcomp(pca_input[pca_complete, ], center = TRUE, scale. = TRUE)
pca_scores <- as_tibble(pca_fit$x[, 1:2], .name_repair = "minimal") %>% bind_cols(raw_data[pca_complete, ] %>% dplyr::select(Lib, Mol, Sizes, outcome, outcome_bin))
pca_loadings <- as_tibble(pca_fit$rotation[, 1:2], rownames = "param") %>% dplyr::left_join(thresholds %>% dplyr::mutate(param = as.character(param)) %>% dplyr::select(param, term, module), by = "param")
pca_variance <- tibble(PC = names(summary(pca_fit)$importance[2, ]), variance_percent = as.numeric(summary(pca_fit)$importance[2, ] * 100))

readr::write_csv(thresholds, file.path(out_dir, "SCRAP_thresholds.csv"))
readr::write_csv(score_long %>% dplyr::select(-.row_id), file.path(out_dir, "SCRAP_score_long.csv"))
readr::write_csv(score_table %>% dplyr::select(-.row_id), file.path(out_dir, "SCRAP_score_summary.csv"))
readr::write_csv(full_metrics, file.path(out_dir, "SCRAP_strict_13_metrics.csv"))
readr::write_csv(importance_tbl, file.path(out_dir, "SCRAP_term_importance.csv"))
readr::write_csv(weighted_score_table %>% dplyr::select(-.row_id), file.path(out_dir, "SCRAP_weighted_score_summary.csv"))
readr::write_csv(weighted_metrics, file.path(out_dir, "SCRAP_weighted_score_metrics.csv"))
readr::write_csv(auc_tbl, file.path(out_dir, "SCRAP_AUC_summary.csv"))
readr::write_csv(module_importance_tbl, file.path(out_dir, "SCRAP_module_importance.csv"))
readr::write_csv(module_ablation_tbl, file.path(out_dir, "SCRAP_module_ablation.csv"))
readr::write_csv(module_failure_tbl, file.path(out_dir, "SCRAP_module_failure_decomposition.csv"))
readr::write_csv(loocv_results, file.path(out_dir, "SCRAP_LOOCV_results.csv"))
readr::write_csv(loocv_metrics, file.path(out_dir, "SCRAP_LOOCV_metrics.csv"))
readr::write_csv(threshold_stability, file.path(out_dir, "SCRAP_threshold_stability.csv"))
readr::write_csv(threshold_perturbation_tbl, file.path(out_dir, "SCRAP_threshold_perturbation.csv"))
readr::write_csv(correlation_tbl, file.path(out_dir, "SCRAP_spearman_correlation.csv"))
readr::write_csv(pca_scores, file.path(out_dir, "SCRAP_PCA_scores.csv"))
readr::write_csv(pca_loadings, file.path(out_dir, "SCRAP_PCA_loadings.csv"))
readr::write_csv(pca_variance, file.path(out_dir, "SCRAP_PCA_variance.csv"))


# -----------------------------
# SI diagnostic figures
# -----------------------------

is_taxane_related_mol <- function(mol) {
  mol <- as.character(mol)
  grepl("^(PTX|DTX|Cbz|CBZ|P2|P7|D2|C2)", mol)
}

add_compatibility_sort_keys <- function(df) {
  if (!"Mol"%in% names(df)) return(df)
  
  mol_chr <- as.character(df$Mol)
  df$mol_taxane_rank <- ifelse(is_taxane_related_mol(mol_chr), 1L, 2L)
  
  outcome_chr <- if ("outcome"%in% names(df)) as.character(df$outcome) else rep(NA_character_, nrow(df))
  compatibility_chr <- if ("compatibility"%in% names(df)) as.character(df$compatibility) else rep(NA_character_, nrow(df))
  
  df$mol_outcome_rank <- dplyr::case_when(
    outcome_chr %in% c("Compatible", "Sub-30 nm", "Sub-30 nm particles") ~ 1L,
    compatibility_chr == "Compatible"~ 1L,
    outcome_chr %in% c("Incompatible", "Oversized", "Precipitated", "Oversized particles") ~ 2L,
    compatibility_chr == "Incompatible"~ 2L,
    TRUE ~ 3L
  )
  
  status_chr <- rep(NA_character_, nrow(df))
  if ("assembly_status"%in% names(df)) status_chr <- as.character(df$assembly_status)
  if ("assembly_type"%in% names(df)) status_chr <- as.character(df$assembly_type)
  if ("outcome"%in% names(df)) {
    status_chr <- ifelse(is.na(status_chr) | status_chr == "", outcome_chr, status_chr)
  }
  
  size_val <- rep(NA_real_, nrow(df))
  if ("Sizes_num"%in% names(df)) size_val <- suppressWarnings(as.numeric(df$Sizes_num))
  if ("Sizes"%in% names(df)) size_val <- suppressWarnings(as.numeric(df$Sizes))
  if ("Size"%in% names(df)) size_val <- suppressWarnings(as.numeric(df$Size))
  
  precip_flag <- rep(FALSE, nrow(df))
  if ("is_precipitated"%in% names(df)) precip_flag <- as.logical(df$is_precipitated)
  precip_flag <- precip_flag | is.na(size_val) | size_val <= 0 | grepl("Precip", status_chr, ignore.case = TRUE)
  
  df$mol_assembly_rank <- dplyr::case_when(
    df$mol_outcome_rank == 1L ~ 1L,
    precip_flag ~ 3L,
    grepl("Oversized", status_chr, ignore.case = TRUE) ~ 2L,
    df$mol_outcome_rank == 2L & !precip_flag ~ 2L,
    TRUE ~ 4L
  )
  
  df$mol_size_sort <- ifelse(is.na(size_val), Inf, size_val)
  df
}

module_order <- c(
  "Prodrug self-association",
  "PCL anchoring",
  "Carrier integration",
  "PEG-associated contacts"
)

heatmap_module_order <- c(
  "Carrier integration",
  "PCL anchoring",
  "PEG-associated contacts",
  "Prodrug self-association"
)

term_display <- c(
  "Pro-Pol"= "Pro\u2013Pol",
  "Par-Pol"= "Par\u2013Pol",
  "Mod-Pol"= "Mod\u2013Pol",
  "Pro-PCL"= "Pro\u2013PCL",
  "Par-PCL"= "Par\u2013PCL",
  "Mod-PCL"= "Mod\u2013PCL",
  "Pro-PEG"= "Pro\u2013PEG",
  "Par-PEG"= "Par\u2013PEG",
  "Mod-PEG"= "Mod\u2013PEG",
  "Pro-Pro"= "Pro\u2013Pro",
  "Par-Par"= "Par\u2013Par",
  "Mod-Mod"= "Mod\u2013Mod",
  "Par-Mod"= "Par\u2013Mod"
)

# Recreate the triangular Spearman table used by the original plotted panel.
energy_matrix <- raw_data %>%
  dplyr::select(dplyr::all_of(param_order))

cor_mat <- suppressWarnings(
  stats::cor(energy_matrix, method = "spearman", use = "pairwise.complete.obs")
)

cor_long_full <- as.data.frame(as.table(cor_mat)) %>%
  tibble::as_tibble() %>%
  dplyr::rename(param_x = Var1, param_y = Var2, spearman_rho = Freq) %>%
  dplyr::mutate(
    param_x = as.character(param_x),
    param_y = as.character(param_y),
    param_x_index = match(param_x, param_order),
    param_y_index = match(param_y, param_order)
  )

cor_long <- cor_long_full %>%
  dplyr::filter(param_x_index <= param_y_index) %>%
  dplyr::mutate(
    param_x_display = factor(
      term_display[param_x],
      levels = term_display[param_order]
    ),
    param_y_display = factor(
      term_display[param_y],
      levels = rev(term_display[param_order])
    )
  )

pca_var <- summary(pca_fit)$importance[2, 1:2] * 100
perturbation_results <- threshold_perturbation_tbl
failure_module_tbl <- module_failure_tbl

module_palette <- c(
  "Carrier integration"= "#6F94B9",
  "PCL anchoring"= "#D89C45",
  "PEG-associated contacts"= "#7FAE77",
  "Prodrug self-association"= "#C96A70"
)

outcome_palette <- c(
  "Incompatible"= "#C96A70",
  "Compatible"= "#66AAA0"
)

window_palette <- c(
  "Outside window"= "#EBEBEB",
  "Within window"= "#66AAA0"
)

metric_palette <- c(
  "Accuracy"= "#5E8C84",
  "Balanced accuracy"= "#6F94B9",
  "Sensitivity"= "#D89C45",
  "Specificity"= "#9D7BA5"
)

corr_palette <- c(
  low = "#D9A0A4",
  mid = "#F7F7F7",
  high = "#7EA6C8"
)

axis_joined_theme <- ggplot2::theme(
  axis.line.x = ggplot2::element_line(color = "#222222", linewidth = 0.35),
  axis.line.y = ggplot2::element_line(color = "#222222", linewidth = 0.35),
  axis.ticks = ggplot2::element_line(color = "#222222", linewidth = 0.28),
  panel.border = ggplot2::element_blank()
)

module_importance_plot_tbl <- module_importance_tbl %>%
  dplyr::mutate(module = factor(as.character(module), levels = module_order))

module_ablation_plot_tbl <- module_ablation_tbl %>%
  dplyr::mutate(module = factor(as.character(module), levels = module_order))

failure_plot_data <- failure_module_tbl %>%
  dplyr::filter(outcome == "Incompatible") %>%
  dplyr::group_by(Mol) %>%
  dplyr::mutate(total_weight_lost_percent = sum(module_weight_lost_percent, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(module_weight_lost_percent > 0) %>%
  dplyr::mutate(
    Mol = forcats::fct_reorder(Mol, total_weight_lost_percent),
    module = factor(as.character(module), levels = module_order)
  )

p_module_importance <- module_importance_plot_tbl %>%
  dplyr::mutate(module = forcats::fct_reorder(module, total_weight_percent)) %>%
  ggplot2::ggplot(ggplot2::aes(x = total_weight_percent, y = module, fill = module)) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.20) +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(round(total_weight_percent, 1), "%")),
    hjust = -0.10, size = 2.4
  ) +
  ggplot2::scale_fill_manual(values = module_palette) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::labs(
    x = "Total weighted contribution (%)",
    y = NULL,
    fill = "Module",
    title = "Module-level diagnostic contribution",
    subtitle = "Aggregated from parameter-level importance weights"
  ) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8),
    plot.margin = ggplot2::margin(5, 6, 5, 6)
  ) +
  axis_joined_theme

p_module_ablation <- module_ablation_plot_tbl %>%
  dplyr::mutate(module = forcats::fct_reorder(module, MCC_drop)) %>%
  ggplot2::ggplot(ggplot2::aes(x = MCC_drop, y = module, fill = module)) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.20) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.3f", MCC_drop)),
    hjust = -0.10, size = 2.4
  ) +
  ggplot2::scale_fill_manual(values = module_palette) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.22))) +
  ggplot2::labs(
    x = "MCC decrease after removing module",
    y = NULL,
    fill = "Module",
    title = "Module-level ablation test",
    subtitle = "Whole-module removal evaluates process-level necessity"
  ) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8)
  )

p_failure_modules <- ggplot2::ggplot(
  failure_plot_data,
  ggplot2::aes(x = Mol, y = module_weight_lost_percent, fill = module)
) +
  ggplot2::geom_col(width = 0.70, color = "white", linewidth = 0.18) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = module_palette) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.06))) +
  ggplot2::labs(
    x = NULL,
    y = "Lost weighted score caused by failed windows (%)",
    fill = "Failed module",
    title = "Failure-mode decomposition of incompatible formulations",
    subtitle = "Higher lost weight indicates more influential failed interaction modules"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE, title.position = "top")
  ) +
  ggplot2::theme(
    legend.position = "right",
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 7),
    legend.key.size = grid::unit(0.28, "cm"),
    legend.text = ggplot2::element_text(size = 6.5),
    legend.title = ggplot2::element_text(size = 7),
    plot.margin = ggplot2::margin(5, 6, 5, 6)
  ) +
  axis_joined_theme

importance_plot_tbl <- importance_tbl %>%
  dplyr::mutate(
    param_display = factor(
      term_display[as.character(param)],
      levels = term_display[param_order]
    ),
    module = factor(as.character(module), levels = module_order)
  )

p_importance <- importance_plot_tbl %>%
  dplyr::mutate(param_display = forcats::fct_reorder(param_display, importance_index)) %>%
  ggplot2::ggplot(ggplot2::aes(x = importance_index, y = param_display, fill = module)) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.20) +
  ggplot2::scale_fill_manual(values = module_palette) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.03))) +
  ggplot2::labs(
    x = "Integrated importance index",
    y = NULL,
    fill = "Module",
    title = "Parameter importance by integrated diagnostic evidence",
    subtitle = "Weights combine equal baseline contribution with data-derived importance"
  ) +
  ggplot2::theme(
    legend.position = "right",
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8),
    plot.margin = ggplot2::margin(5, 6, 5, 6)
  ) +
  axis_joined_theme

score_plot_data <- weighted_score_table %>%
  dplyr::mutate(
    Mol = forcats::fct_reorder(Mol, score_weighted_100),
    compatibility_label = dplyr::if_else(
      outcome == "Compatible", "Compatible", "Incompatible"
    )
  )

incompatible_mol_order_for_heatmap <- weighted_score_table %>%
  dplyr::filter(outcome == "Incompatible") %>%
  add_compatibility_sort_keys() %>%
  dplyr::arrange(
    mol_assembly_rank,
    mol_taxane_rank,
    dplyr::desc(score_13),
    dplyr::desc(score_weighted_100),
    mol_size_sort,
    Mol
  ) %>%
  dplyr::pull(Mol) %>%
  as.character()

incompatible_mol_labels <- weighted_score_table %>%
  dplyr::filter(as.character(Mol) %in% incompatible_mol_order_for_heatmap) %>%
  dplyr::mutate(
    Mol = as.character(Mol),
    size_label = dplyr::if_else(
      is.na(Sizes) | Sizes <= 0,
      "precip.",
      paste0(round(Sizes, 1), "nm")
    ),
    mol_label = paste0(Mol, "\n", size_label, "\n", score_13, "/13")
  ) %>%
  dplyr::select(Mol, mol_label) %>%
  tibble::deframe()

incompatible_mol_label_levels <- incompatible_mol_labels[incompatible_mol_order_for_heatmap]

pass_heatmap_data <- score_long %>%
  dplyr::filter(outcome == "Incompatible") %>%
  dplyr::left_join(
    weighted_score_table %>%
      dplyr::select(Mol, score_13, score_weighted_100),
    by = "Mol"
  ) %>%
  dplyr::mutate(
    param_chr = as.character(param),
    param_display = factor(
      term_display[param_chr],
      levels = rev(term_display[param_order])
    ),
    module = factor(as.character(module), levels = heatmap_module_order),
    Mol = factor(
      as.character(Mol),
      levels = incompatible_mol_order_for_heatmap,
      labels = incompatible_mol_label_levels
    ),
    pass_label = dplyr::if_else(pass == 1, "Within window", "Outside window"),
    pass_label = factor(
      pass_label,
      levels = c("Outside window", "Within window")
    )
  )

n_incompatible_for_heatmap <- length(incompatible_mol_order_for_heatmap)
heatmap_width <- max(8.2, 2.4 + 0.42 * n_incompatible_for_heatmap)

p_pass_heatmap <- ggplot2::ggplot(
  pass_heatmap_data,
  ggplot2::aes(x = Mol, y = param_display, fill = pass_label)
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.62,
    width = 0.94,
    height = 0.86
  ) +
  ggplot2::facet_grid(
    module ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::scale_fill_manual(values = window_palette) +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title = "Segmented window map of incompatible formulations",
    subtitle = "Grey marks failed interaction windows; teal marks terms within the favourable window"
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = ggplot2::element_text(face = "bold"),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    strip.text.y.left = ggplot2::element_text(
      angle = 0, hjust = 1, face = "plain", size = 8.2
    ),
    axis.text.x.top = ggplot2::element_text(
      size = 6.7, lineheight = 0.90, margin = ggplot2::margin(b = 4)
    ),
    axis.text.y = ggplot2::element_text(size = 7.2),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.spacing.y = grid::unit(0.12, "lines"),
    panel.border = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(6, 8, 6, 6)
  )

p_weighted_score <- ggplot2::ggplot(
  score_plot_data,
  ggplot2::aes(x = Mol, y = score_weighted_100, fill = outcome)
) +
  ggplot2::geom_col(width = 0.70, color = "white", linewidth = 0.20) +
  ggplot2::geom_hline(
    yintercept = 100,
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = outcome_palette) +
  ggplot2::scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Importance-weighted compatibility score (%)",
    fill = "Observed outcome",
    title = "Importance-weighted 100% compatibility score",
    subtitle = "Strict compatible call requires complete satisfaction of all 13 windows"
  ) +
  ggplot2::theme(
    legend.position = "right",
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 7)
  )

p_ablation <- importance_plot_tbl %>%
  dplyr::mutate(
    param_display = forcats::fct_reorder(param_display, ablation_MCC_drop)
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(x = ablation_MCC_drop, y = param_display, fill = module)
  ) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.20) +
  ggplot2::scale_fill_manual(values = module_palette) +
  ggplot2::labs(
    x = "MCC decrease after removing a single parameter",
    y = NULL,
    fill = "Module",
    title = "Ablation-based non-redundant necessity test",
    subtitle = "Single-parameter removal identifies terms that cannot be compensated by correlated descriptors"
  ) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8)
  ) +
  axis_joined_theme

p_perturb <- perturbation_results %>%
  dplyr::select(
    perturbation_percent,
    accuracy,
    sensitivity,
    specificity,
    balanced_accuracy
  ) %>%
  tidyr::pivot_longer(
    cols = c(accuracy, sensitivity, specificity, balanced_accuracy),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      accuracy = "Accuracy",
      sensitivity = "Sensitivity",
      specificity = "Specificity",
      balanced_accuracy = "Balanced accuracy"
    ),
    metric = factor(
      metric,
      levels = c(
        "Accuracy",
        "Balanced accuracy",
        "Sensitivity",
        "Specificity"
      )
    )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = perturbation_percent,
      y = 100 * value,
      color = metric,
      shape = metric,
      linetype = metric
    )
  ) +
  ggplot2::geom_line(linewidth = 0.45) +
  ggplot2::geom_point(size = 1.7, stroke = 0.35) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.30,
    color = "#6B6B6B"
  ) +
  ggplot2::scale_color_manual(values = metric_palette) +
  ggplot2::scale_linetype_manual(values = c("solid", "solid", "22", "22")) +
  ggplot2::scale_shape_manual(values = c(16, 17, 15, 3)) +
  ggplot2::scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::scale_x_continuous(breaks = seq(-20, 20, 10)) +
  ggplot2::labs(
    x = "Uniform threshold perturbation (%)",
    y = "Performance (%)",
    color = "Metric",
    linetype = "Metric",
    shape = "Metric",
    title = "Threshold perturbation analysis"
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  ) +
  ggplot2::theme(
    legend.position = "right",
    legend.key.width = grid::unit(0.55, "cm"),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(5, 6, 5, 6)
  ) +
  axis_joined_theme

p_corr <- ggplot2::ggplot(
  cor_long,
  ggplot2::aes(
    x = param_x_display,
    y = param_y_display,
    fill = spearman_rho
  )
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.18) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", spearman_rho)),
    size = 1.55,
    color = "#222222"
  ) +
  ggplot2::scale_fill_gradient2(
    low = corr_palette["low"],
    mid = corr_palette["mid"],
    high = corr_palette["high"],
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish
  ) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = NULL,
    y = NULL,
    fill = "Spearman rho",
    title = "Correlation among segment-specific interaction energies",
    subtitle = "One triangular half is shown because the correlation matrix is symmetric"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1, size = 7
    ),
    axis.text.y = ggplot2::element_text(size = 7),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.margin = ggplot2::margin(5, 6, 5, 6)
  ) +
  axis_joined_theme

p_pca <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, fill = outcome, size = Sizes)
) +
  geom_point(
    shape = 21,
    colour = "white",
    alpha = 0.90,
    stroke = 0.30
  ) +
  scale_fill_manual(values = outcome_palette) +
  scale_size_continuous(
    range = c(2.2, 6.5),
    breaks = c(0, 50, 100)
  ) +
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        fill = "grey70",
        colour = "grey35",
        alpha = 1,
        stroke = 0.4
      )
    ),
    fill = guide_legend(order = 2)
  ) +
  labs(
    x = paste0("PC1 (", round(pca_var[[1]], 1), "% variance)"),
    y = paste0("PC2 (", round(pca_var[[2]], 1), "% variance)"),
    fill = "Observed outcome",
    size = "DLS size (nm)",
    title = "Low-dimensional structure of the 13-parameter space"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(5, 6, 5, 6)
  ) +
  axis_joined_theme

# Individual figures: filenames, dimensions and 600-dpi PNG settings retained.
ggplot2::ggsave(
  file.path(fig_dir, "incompatible_13window_pass_heatmap.pdf"),
  p_pass_heatmap, width = heatmap_width, height = 5.8
)
ggplot2::ggsave(
  file.path(fig_dir, "incompatible_13window_pass_heatmap.png"),
  p_pass_heatmap, width = heatmap_width, height = 5.8, dpi = 600
)
ggplot2::ggsave(
  file.path(fig_dir, "13window_pass_heatmap.pdf"),
  p_pass_heatmap, width = heatmap_width, height = 5.8
)
ggplot2::ggsave(
  file.path(fig_dir, "13window_pass_heatmap.png"),
  p_pass_heatmap, width = heatmap_width, height = 5.8, dpi = 600
)

ggplot2::ggsave(
  file.path(fig_dir, "module_importance.pdf"),
  p_module_importance, width = 6.2, height = 3.2
)
ggplot2::ggsave(
  file.path(fig_dir, "module_importance.png"),
  p_module_importance, width = 6.2, height = 3.2, dpi = 600
)

ggplot2::ggsave(
  file.path(fig_dir, "module_ablation.pdf"),
  p_module_ablation, width = 6.2, height = 3.2
)
ggplot2::ggsave(
  file.path(fig_dir, "module_ablation.png"),
  p_module_ablation, width = 6.2, height = 3.2, dpi = 600
)

ggplot2::ggsave(
  file.path(fig_dir, "incompatible_failure_decomposition.pdf"),
  p_failure_modules, width = 7.4, height = 4.8
)
ggplot2::ggsave(
  file.path(fig_dir, "incompatible_failure_decomposition.png"),
  p_failure_modules, width = 7.4, height = 4.8, dpi = 600
)

ggplot2::ggsave(
  file.path(fig_dir, "parameter_importance.pdf"),
  p_importance, width = 7.2, height = 4.6
)
ggplot2::ggsave(
  file.path(fig_dir, "parameter_importance.png"),
  p_importance, width = 7.2, height = 4.6, dpi = 600
)

ggplot2::ggsave(
  file.path(fig_dir, "weighted_100_score.pdf"),
  p_weighted_score, width = 7.4, height = 6.4
)
ggplot2::ggsave(
  file.path(fig_dir, "weighted_100_score.png"),
  p_weighted_score, width = 7.4, height = 6.4, dpi = 600
)

ggplot2::ggsave(
  file.path(fig_dir, "ablation.pdf"),
  p_ablation, width = 6.2, height = 4.4
)
ggplot2::ggsave(
  file.path(fig_dir, "ablation.png"),
  p_ablation, width = 6.2, height = 4.4, dpi = 600
)

ggplot2::ggsave(
  file.path(fig_dir, "threshold_perturbation.pdf"),
  p_perturb, width = 6.4, height = 4.2
)
ggplot2::ggsave(
  file.path(fig_dir, "threshold_perturbation.png"),
  p_perturb, width = 6.4, height = 4.2, dpi = 600
)

ggplot2::ggsave(
  file.path(fig_dir, "spearman_correlation.pdf"),
  p_corr, width = 6.4, height = 5.6
)
ggplot2::ggsave(
  file.path(fig_dir, "spearman_correlation.png"),
  p_corr, width = 6.4, height = 5.6, dpi = 600
)

ggplot2::ggsave(
  file.path(fig_dir, "pca.pdf"),
  p_pca, width = 5.8, height = 4.6
)
ggplot2::ggsave(
  file.path(fig_dir, "pca.png"),
  p_pca, width = 5.8, height = 4.6, dpi = 600
)

combined_main <- (p_importance / p_weighted_score) +
  patchwork::plot_layout(heights = c(0.85, 1.15))

ggplot2::ggsave(
  file.path(fig_dir, "combined_importance_and_weighted_score.pdf"),
  combined_main, width = 8.0, height = 10.0
)
ggplot2::ggsave(
  file.path(fig_dir, "combined_importance_and_weighted_score.png"),
  combined_main, width = 8.0, height = 10.0, dpi = 600
)

combined_validation <- (p_ablation | p_perturb) /
  (p_corr | p_pca) +
  patchwork::plot_layout(heights = c(0.8, 1.2))

ggplot2::ggsave(
  file.path(fig_dir, "combined_model_validation.pdf"),
  combined_validation, width = 11.0, height = 9.0
)
ggplot2::ggsave(
  file.path(fig_dir, "combined_model_validation.png"),
  combined_validation, width = 11.0, height = 9.0, dpi = 600
)

combined_module_diagnostics <- (
  p_module_importance | p_module_ablation
) / p_failure_modules +
  patchwork::plot_layout(heights = c(0.72, 1.28))

ggplot2::ggsave(
  file.path(fig_dir, "combined_module_diagnostics.pdf"),
  combined_module_diagnostics, width = 11.0, height = 8.2
)
ggplot2::ggsave(
  file.path(fig_dir, "combined_module_diagnostics.png"),
  combined_module_diagnostics, width = 11.0, height = 8.2, dpi = 600
)

extended_row_importance <- (
  (p_importance + ggplot2::theme(legend.position = "none")) |
    (p_module_importance + ggplot2::theme(legend.position = "none"))
) +
  patchwork::plot_layout(widths = c(1.02, 0.98))

extended_row_failure <- (
  (p_failure_modules + ggplot2::theme(legend.position = "bottom")) |
    (p_corr + ggplot2::theme(legend.position = "right"))
) +
  patchwork::plot_layout(widths = c(0.82, 1.18))

extended_row_validation <- (
  (p_perturb + ggplot2::theme(legend.position = "bottom")) |
    (p_pca + ggplot2::theme(legend.position = "right"))
) +
  patchwork::plot_layout(widths = c(1, 1))

extended_model_diagnostics <- (
  p_pass_heatmap /
    extended_row_importance /
    extended_row_failure /
    extended_row_validation
) +
  patchwork::plot_layout(
    heights = c(1.12, 0.95, 1.25, 1.00)
  ) +
  patchwork::plot_annotation(tag_levels = "a") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 11),
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  )

ggplot2::ggsave(
  file.path(fig_dir, "extended_segmented_model_diagnostics.pdf"),
  extended_model_diagnostics,
  width = 13.0, height = 12.4, limitsize = FALSE
)
ggplot2::ggsave(
  file.path(fig_dir, "extended_segmented_model_diagnostics.png"),
  extended_model_diagnostics,
  width = 13.0, height = 12.4, dpi = 600, limitsize = FALSE
)

message(
  "Analysis complete. Tables written to: ",
  normalizePath(out_dir, winslash = "/", mustWork = FALSE)
)
message(
  "Figures written to: ",
  normalizePath(fig_dir, winslash = "/", mustWork = FALSE)
)
