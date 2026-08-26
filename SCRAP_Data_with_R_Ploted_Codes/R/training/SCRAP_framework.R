# SCRAP framework analysis
# Input: data/training/Training_segment_data.csv


rm(list = ls())
graphics.off()
set.seed(123)

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_arg_match <- grep(file_arg, cmd_args)
  if (length(file_arg_match) > 0) {
    script_path <- sub(file_arg, "", cmd_args[file_arg_match[1]])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
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
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

required_packages <- c("tidyverse", "patchwork", "ggrepel", "scales", "grid")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
})

input_file <- file.path(project_root, "data", "training", "Training_segment_data.csv")
figure_dir <- file.path(project_root, "results", "figures", "training", "SCRAP_framework")
table_dir <- file.path(project_root, "results", "tables", "training", "SCRAP_framework")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Unified theme and palettes
# -----------------------------
theme_nature <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(family = "sans", colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      axis.title = element_text(colour = "#222222"),
      axis.line = element_line(linewidth = 0.30, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#222222"),
      axis.ticks.length = unit(1.2, "mm"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.8),
      plot.title = element_text(face = "bold", size = base_size + 0.9, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.6, hjust = 0, colour = "#555555"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
}

col_compatible   <- "#2F5D8C"
col_incompatible <- "#C76F5A"
col_pass         <- "#2A9D8F"
col_fail         <- "#ECECEC"

outcome_cols <- c(
  "Compatible"= col_compatible,
  "Incompatible"= col_incompatible
)

module_order <- c(
  "Carrier integration",
  "PCL anchoring",
  "PEG-associated contacts",
  "Prodrug self-association"
)

module_cols <- c(
  "Carrier integration"= "#6F94B9",
  "PCL anchoring"= "#D9A441",
  "PEG-associated contacts"= "#8BBE84",
  "Prodrug self-association"= "#C96A70"
)

# -----------------------------
# 2. Parameters, labels and scoring windows
# -----------------------------
params <- c(
  "Pro-Pol", "Par-Pol", "Mod-Pol",
  "Pro-PCL", "Par-PCL", "Mod-PCL",
  "Pro-PEG", "Par-PEG", "Mod-PEG",
  "Pro-Pro", "Par-Par", "Mod-Mod", "Par-Mod"
)

param_labels <- c(
  "Pro-Pol"= "Pro–Pol",
  "Par-Pol"= "Par–Pol",
  "Mod-Pol"= "Mod–Pol",
  "Pro-PCL"= "Pro–PCL",
  "Par-PCL"= "Par–PCL",
  "Mod-PCL"= "Mod–PCL",
  "Pro-PEG"= "Pro–PEG",
  "Par-PEG"= "Par–PEG",
  "Mod-PEG"= "Mod–PEG",
  "Pro-Pro"= "Pro–Pro",
  "Par-Par"= "Par–Par",
  "Mod-Mod"= "Mod–Mod",
  "Par-Mod"= "Par–Mod"
)

thresholds <- tribble(
  ~param,       ~min_val, ~max_val, ~module,
  "Pro-Pol",   -38.1,    -18,      "Carrier integration",
  "Par-Pol",   -19,       0,       "Carrier integration",
  "Mod-Pol",   NA,       -4,       "Carrier integration",
  "Pro-PCL",   -31,      -15,      "PCL anchoring",
  "Par-PCL",   -15.6,     0,       "PCL anchoring",
  "Mod-PCL",   NA,       -3.5,     "PCL anchoring",
  "Pro-PEG",  NA,       -3.3,     "PEG-associated contacts",
  "Par-PEG",   -6.5,      0,       "PEG-associated contacts",
  "Mod-PEG",   -7,        0,       "PEG-associated contacts",
  "Pro-Pro", -28,      -6.1,     "Prodrug self-association",
  "Par-Par",   -11,       0,       "Prodrug self-association",
  "Mod-Mod",   NA,       -0.42,    "Prodrug self-association",
  "Par-Mod",   -9,       -1,       "Prodrug self-association"
) %>%
  mutate(
    param = factor(param, levels = params),
    module = factor(module, levels = module_order)
  )

in_window <- function(x, min_val, max_val) {
  ok_min <- ifelse(is.na(min_val), TRUE, x >= min_val)
  ok_max <- ifelse(is.na(max_val), TRUE, x <= max_val)
  ok_min & ok_max
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
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(denom > 0, ((tp * tn) - (fp * fn)) / denom, 0)
  tibble(
    TP = tp, TN = tn, FP = fp, FN = fn,
    accuracy = (tp + tn) / length(actual),
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    sensitivity = sensitivity,
    specificity = specificity,
    MCC = mcc
  )
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
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || abs(rng[2] - rng[1]) < .Machine$double.eps) {
    return(rep(0, length(x)))
  }
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

# -----------------------------
# 3. Read data and calculate 13/13 score
# -----------------------------
raw_data <- readr::read_csv(input_file, show_col_types = FALSE) %>%
  rename_with(~ gsub("^\\ufeff", "", .x)) %>%
  filter(
    !is.na(Mol),
    trimws(as.character(Mol)) != "",
    toupper(trimws(as.character(Mol))) != "NA"
  )

missing_cols <- setdiff(c("Lib", "Mol", "Sizes", params), colnames(raw_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

df <- raw_data %>%
  mutate(
    Mol = as.character(Mol),
    across(all_of(c("Sizes", params)), as.numeric),
    size_label = case_when(
      is.na(Sizes) | Sizes <= 0 ~ "precip.",
      TRUE ~ paste0(round(Sizes, 1), "nm")
    ),
    size_plot = if_else(is.na(Sizes) | Sizes <= 0, 1, Sizes),
    outcome = if_else(!is.na(Sizes) & Sizes > 0 & Sizes < 30, "Compatible", "Incompatible"),
    outcome = factor(outcome, levels = c("Compatible", "Incompatible")),
    outcome_bin = as.integer(outcome == "Compatible")
  )

score_long <- df %>%
  pivot_longer(cols = all_of(params), names_to = "param", values_to = "energy") %>%
  left_join(thresholds, by = "param") %>%
  rowwise() %>%
  mutate(
    pass = as.integer(in_window(energy, min_val, max_val)),
    distance_to_window = window_distance(energy, min_val, max_val)
  ) %>%
  ungroup() %>%
  mutate(
    param = factor(as.character(param), levels = params),
    param_label = factor(param_labels[as.character(param)], levels = rev(param_labels[params])),
    module = factor(module, levels = module_order),
    pass_label = factor(if_else(pass == 1, "Within window", "Outside window"),
                        levels = c("Outside window", "Within window"))
  )

score_summary <- score_long %>%
  group_by(Lib, Mol, Sizes, size_label, size_plot, outcome, outcome_bin) %>%
  summarise(
    score_13 = sum(pass, na.rm = TRUE),
    n_failed = 13 - score_13,
    .groups = "drop"
  ) %>%
  arrange(outcome, Sizes, desc(score_13), Mol)

mol_order <- score_summary %>%
  arrange(outcome, desc(score_13), size_plot, Mol) %>%
  pull(Mol) %>%
  unique()

score_summary <- score_summary %>%
  mutate(Mol = factor(Mol, levels = mol_order)) %>%
  filter(!is.na(Mol))

score_long <- score_long %>%
  left_join(
    score_summary %>%
      mutate(Mol = as.character(Mol)) %>%
      select(Lib, Mol, Sizes, outcome, outcome_bin, score_13, n_failed),
    by = c("Lib", "Mol", "Sizes", "outcome", "outcome_bin")
  ) %>%
  mutate(Mol = factor(as.character(Mol), levels = mol_order)) %>%
  filter(!is.na(Mol))

readr::write_csv(score_summary, file.path(table_dir, "score_summary.csv"))
readr::write_csv(score_long, file.path(table_dir, "score_long.csv"))

# -----------------------------
# 4. Integrated importance index
# -----------------------------
wide_pass <- score_long %>%
  select(Mol, param, pass) %>%
  distinct() %>%
  pivot_wider(names_from = param, values_from = pass) %>%
  mutate(Mol = factor(as.character(Mol), levels = mol_order)) %>%
  arrange(Mol)

score_summary_for_metrics <- score_summary %>%
  mutate(Mol = factor(as.character(Mol), levels = mol_order)) %>%
  arrange(Mol)

full_metrics <- calc_metrics(score_summary_for_metrics$outcome_bin, score_summary_for_metrics$score_13 == 13)

importance_tbl <- map_dfr(params, function(p) {
  tmp <- score_long %>% filter(as.character(param) == p)
  pass_vec <- tmp$pass
  y <- tmp$outcome_bin
  compat_pass_rate <- mean(pass_vec[y == 1], na.rm = TRUE)
  incompat_pass_rate <- mean(pass_vec[y == 0], na.rm = TRUE)
  failure_gap <- compat_pass_rate - incompat_pass_rate
  ig <- information_gain(pass_vec, y)
  mcc_single <- abs(calc_metrics(y, pass_vec)$MCC[[1]])
  dist_auc <- calc_auc(-tmp$distance_to_window, y)
  dist_auc_gain <- pmax(dist_auc - 0.5, 0)
  score_without <- score_summary_for_metrics$score_13 - wide_pass[[p]]
  pred_without <- score_without == 12
  without_metrics <- calc_metrics(score_summary_for_metrics$outcome_bin, pred_without)
  tibble(
    param = p,
    compatible_pass_rate = compat_pass_rate,
    incompatible_pass_rate = incompat_pass_rate,
    failure_gap = failure_gap,
    information_gain = ig,
    phi_abs = mcc_single,
    distance_auc_gain = dist_auc_gain,
    ablation_MCC_drop = full_metrics$MCC[[1]] - without_metrics$MCC[[1]]
  )
}) %>%
  left_join(thresholds %>% mutate(param = as.character(param)) %>% select(param, module), by = "param") %>%
  mutate(
    norm_failure_gap = norm01(failure_gap),
    norm_information_gain = norm01(information_gain),
    norm_phi_abs = norm01(phi_abs),
    norm_distance_auc_gain = norm01(distance_auc_gain),
    norm_ablation_MCC_drop = norm01(pmax(ablation_MCC_drop, 0)),
    importance_index = rowMeans(
      cbind(norm_failure_gap, norm_information_gain, norm_phi_abs,
            norm_distance_auc_gain, norm_ablation_MCC_drop),
      na.rm = TRUE
    )
  )

if (sum(importance_tbl$importance_index, na.rm = TRUE) <= 0) {
  importance_tbl <- importance_tbl %>% mutate(importance_index = 1)
}

importance_tbl <- importance_tbl %>%
  mutate(
    module = factor(module, levels = module_order),
    weight_fraction = 0.20 / 13 + 0.80 * importance_index / sum(importance_index, na.rm = TRUE),
    weight_percent = 100 * weight_fraction,
    param_label = factor(param, levels = params, labels = param_labels[params])
  ) %>%
  arrange(desc(importance_index), param)

module_importance_tbl <- importance_tbl %>%
  group_by(module) %>%
  summarise(
    total_weight_percent = sum(weight_percent, na.rm = TRUE),
    total_importance_index = sum(importance_index, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(module = factor(module, levels = module_order))

readr::write_csv(importance_tbl, file.path(table_dir, "term_importance.csv"))
readr::write_csv(module_importance_tbl, file.path(table_dir, "module_importance.csv"))

# -----------------------------
# 5. Panel c: interaction-energy heatmap
# -----------------------------
p_c_heatmap <- ggplot(score_long, aes(x = param_label, y = Mol, fill = energy)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  facet_grid(. ~ module, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(
    colours = c("#1F4E79", "#7FA6C7", "#F5F5F5", "#E3B6AD", "#B85648"),
    values = scales::rescale(c(-40, -25, -15, -5, 2)),
    limits = c(-40, 2),
    oob = squish,
    name = expression(Delta * "G\n(kcal mol"^-1 * ")")
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, title = "Segment interaction-energy landscape") +
  theme_nature(7.0) +
  theme(
    axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 5.0),
    axis.text.y = element_text(size = 5.3),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    strip.text = element_text(size = 6.0),
    panel.spacing.x = unit(0.35, "lines"),
    legend.position = "bottom",
    legend.key.width = unit(1.0, "cm"),
    legend.key.height = unit(0.18, "cm")
  )

# -----------------------------
# 6. Panel d: favourable-window strip
# -----------------------------
strip_df <- score_long %>%
  mutate(
    row_id = as.numeric(param_label),
    jitter_y = row_id + runif(n(), -0.12, 0.12)
  )

strip_windows <- thresholds %>%
  mutate(
    param_label = factor(param_labels[as.character(param)], levels = rev(param_labels[params])),
    row_id = as.numeric(param_label),
    plot_min = if_else(is.na(min_val), -42, min_val),
    plot_max = if_else(is.na(max_val), 2, max_val),
    ymin = row_id - 0.22,
    ymax = row_id + 0.22
  )

p_d_windows <- ggplot() +
  geom_rect(
    data = strip_windows,
    aes(xmin = plot_min, xmax = plot_max, ymin = ymin, ymax = ymax),
    fill = "#B7DCCB",
    colour = "#7AAE92",
    linewidth = 0.20,
    alpha = 0.75
  ) +
  geom_segment(
    data = strip_windows,
    aes(x = -42, xend = 2, y = row_id, yend = row_id),
    colour = "#A0A0A0",
    linewidth = 0.18
  ) +
  geom_point(
    data = strip_df,
    aes(x = energy, y = jitter_y, fill = outcome),
    shape = 21,
    size = 1.25,
    colour = "white",
    stroke = 0.12,
    alpha = 0.78
  ) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = outcome_cols, name = NULL) +
  scale_x_continuous(limits = c(-42, 2), breaks = c(-40, -30, -20, -10, 0)) +
  scale_y_continuous(
    breaks = sort(unique(strip_windows$row_id)),
    labels = levels(score_long$param_label),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(
    x = expression(Delta * "G (kcal mol"^-1 * ")"),
    y = NULL,
    title = "Favourable windows and raw values"
  ) +
  theme_nature(7.0) +
  theme(
    axis.text.y = element_text(size = 5.4),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 5.8, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    panel.spacing.y = unit(0.30, "lines")
  )

# -----------------------------
# 7. Panel e: 13-point score
# -----------------------------
p_e_score <- ggplot(score_summary, aes(x = Mol, y = score_13)) +
  geom_hline(yintercept = 13, linetype = "dashed", linewidth = 0.32, colour = "#555555") +
  geom_segment(aes(xend = Mol, y = 0, yend = score_13, colour = outcome), linewidth = 0.38, alpha = 0.65) +
  geom_point(aes(fill = outcome), shape = 21, colour = "#222222", size = 2.15, stroke = 0.25) +
  geom_text(
    data = score_summary %>% filter(outcome == "Incompatible"),
    aes(label = score_13),
    colour = "white",
    size = 1.7,
    fontface = "bold"
  ) +
  scale_fill_manual(values = outcome_cols, name = NULL) +
  scale_colour_manual(values = outcome_cols, guide = "none") +
  scale_y_continuous(breaks = seq(1, 13, 2), limits = c(0, 13.4), expand = expansion(mult = c(0, 0.01))) +
  labs(
    x = NULL,
    y = "Segmented compatibility\nscore",
    title = "13-point score separates DLS outcomes"
  ) +
  theme_nature(7.0) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.0),
    legend.position = "right"
  )

# -----------------------------
# 8. Panel f: parameter importance
# -----------------------------
p_f_importance <- importance_tbl %>%
  mutate(param_label = fct_reorder(param_label, importance_index)) %>%
  ggplot(aes(x = importance_index, y = param_label, fill = module)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.20) +
  scale_fill_manual(values = module_cols, name = NULL) +
  scale_x_continuous(limits = c(0, max(importance_tbl$importance_index, na.rm = TRUE) * 1.05),
                     expand = expansion(mult = c(0, 0.01))) +
  labs(
    x = "Integrated importance index",
    y = NULL,
    title = "Parameter-level diagnostic importance"
  ) +
  theme_nature(7.0) +
  theme(
    axis.text.y = element_text(size = 5.5),
    legend.position = "none"
  )

# -----------------------------
# 9. Panel g: module importance
# -----------------------------
p_g_module <- module_importance_tbl %>%
  mutate(module = fct_reorder(module, total_weight_percent)) %>%
  ggplot(aes(x = total_weight_percent, y = module, fill = module)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.20) +
  geom_text(aes(label = paste0(round(total_weight_percent, 1), "%")),
            hjust = -0.10, size = 2.3) +
  scale_fill_manual(values = module_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    x = "Total weighted contribution (%)",
    y = NULL,
    title = "Module-level contribution"
  ) +
  theme_nature(7.0) +
  theme(
    axis.text.y = element_text(size = 5.8),
    legend.position = "none"
  )

# -----------------------------
# 10. Panel h: failed-window burden heatmap
# -----------------------------
# Summarize failed-window burden by incompatibility severity.
severity_levels <- c(
  "Near-miss\n(11–12)",
  "Intermediate\n(8–10)",
  "Severe\n(≤7/precip.)"
)

incompat_severity <- score_summary %>%
  filter(outcome == "Incompatible") %>%
  mutate(
    Mol_chr = as.character(Mol),
    severity = case_when(
      is.na(Sizes) | Sizes <= 0 | score_13 <= 7 ~ "Severe\n(≤7/precip.)",
      score_13 >= 8 & score_13 <= 10 ~ "Intermediate\n(8–10)",
      score_13 >= 11 & score_13 <= 12 ~ "Near-miss\n(11–12)",
      TRUE ~ "Intermediate\n(8–10)"
    ),
    severity = factor(severity, levels = severity_levels)
  ) %>%
  select(Mol_chr, severity)

severity_n <- incompat_severity %>%
  count(severity, name = "n_formulations") %>%
  complete(severity = factor(severity_levels, levels = severity_levels),
           fill = list(n_formulations = 0))

severity_label_map <- severity_n %>%
  mutate(label = paste0(as.character(severity), "\nn=", n_formulations)) %>%
  select(severity, label) %>%
  deframe()

module_n_terms <- thresholds %>%
  count(module, name = "n_terms") %>%
  mutate(module = factor(module, levels = module_order))

failure_heatmap_data <- score_long %>%
  filter(outcome == "Incompatible") %>%
  mutate(Mol_chr = as.character(Mol)) %>%
  left_join(incompat_severity, by = "Mol_chr") %>%
  group_by(severity, module) %>%
  summarise(failed_windows = sum(pass == 0, na.rm = TRUE), .groups = "drop") %>%
  right_join(
    tidyr::crossing(
      severity = factor(severity_levels, levels = severity_levels),
      module = factor(module_order, levels = module_order)
    ),
    by = c("severity", "module")
  ) %>%
  mutate(failed_windows = replace_na(failed_windows, 0L)) %>%
  left_join(severity_n, by = "severity") %>%
  left_join(module_n_terms, by = "module") %>%
  mutate(
    total_possible_windows = n_formulations * n_terms,
    failed_fraction = if_else(total_possible_windows > 0,
                              100 * failed_windows / total_possible_windows,
                              NA_real_),
    mean_failed_per_formulation = if_else(n_formulations > 0,
                                          failed_windows / n_formulations,
                                          NA_real_),
    tile_label = case_when(
      n_formulations > 0 ~ paste0(round(failed_fraction), "%"),
      TRUE ~ "–"
    ),
    label_colour = if_else(!is.na(failed_fraction) & failed_fraction >= 70, "white", "#222222"),
    module = factor(module, levels = rev(module_order)),
    severity = factor(severity, levels = severity_levels)
  )

p_h_failure <- ggplot(failure_heatmap_data, aes(x = severity, y = module, fill = failed_fraction)) +
  geom_tile(colour = "white", linewidth = 0.70, width = 0.94, height = 0.86) +
  geom_text(aes(label = tile_label, colour = label_colour), size = 2.35, fontface = "bold") +
  scale_colour_identity() +
  scale_x_discrete(labels = severity_label_map, position = "top") +
  scale_y_discrete(labels = c(
    "Carrier integration"= "Carrier\nintegration",
    "PCL anchoring"= "PCL\nanchoring",
    "PEG-associated contacts"= "PEG-associated\ncontacts",
    "Prodrug self-association"= "Self-\nassociation"
  )) +
  scale_fill_gradientn(
    colours = c("#DDEAF6", "#F7F7F7", "#E3B1B5", "#C6656D"),
    values = scales::rescale(c(0, 35, 70, 100)),
    limits = c(0, 100),
    oob = scales::squish,
    na.value = "#F1F1F1",
    name = "Failed\nwindows (%)"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Failed-window burden by severity"
  ) +
  theme_nature(7.0) +
  theme(
    axis.text.x.top = element_text(size = 5.4, lineheight = 0.86, margin = margin(b = 2.5)),
    axis.text.y = element_text(size = 5.6, lineheight = 0.86),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.56, "cm"),
    legend.key.width = unit(0.16, "cm"),
    legend.title = element_text(size = 5.8, lineheight = 0.85),
    legend.text = element_text(size = 5.2),
    plot.title = element_text(size = 7.8),
    plot.margin = margin(3, 3, 3, 3)
  )

readr::write_csv(incompat_severity, file.path(table_dir, "incompatible_severity.csv"))
readr::write_csv(failure_heatmap_data, file.path(table_dir, "failed_window_burden.csv"))

# -----------------------------
# 11. Optional placeholders for panel a and b
# -----------------------------
p_a_placeholder <- ggplot() +
  annotate("text", x = 0.5, y = 0.56, label = "Panel a\n13 segment-specific\ninteraction terms",
           size = 3.5, fontface = "bold", colour = "#0B7F78", lineheight = 0.95) +
  annotate("text", x = 0.5, y = 0.34, label = "Replace with schematic artwork",
           size = 2.5, colour = "#555555") +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, colour = "#D0D0D0", linewidth = 0.35))

p_b_placeholder <- ggplot() +
  annotate("text", x = 0.5, y = 0.56, label = "Panel b\nMD/MM-PBSA\nworkflow",
           size = 3.5, fontface = "bold", colour = "#222222", lineheight = 0.95) +
  annotate("text", x = 0.5, y = 0.34, label = "Replace with schematic artwork",
           size = 2.5, colour = "#555555") +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, colour = "#D0D0D0", linewidth = 0.35))

# -----------------------------
# 12. Export panels
# -----------------------------
ggsave(file.path(figure_dir, "interaction_energy_matrix.pdf"), p_c_heatmap, width = 5.6, height = 4.2, device = cairo_pdf)
ggsave(file.path(figure_dir, "interaction_energy_matrix.png"), p_c_heatmap, width = 5.6, height = 4.2, dpi = 900, bg = "white")

ggsave(file.path(figure_dir, "favourable_energy_windows.pdf"), p_d_windows, width = 5.6, height = 4.2, device = cairo_pdf)
ggsave(file.path(figure_dir, "favourable_energy_windows.png"), p_d_windows, width = 5.6, height = 4.2, dpi = 900, bg = "white")

ggsave(file.path(figure_dir, "compatibility_score.pdf"), p_e_score, width = 7.0, height = 3.0, device = cairo_pdf)
ggsave(file.path(figure_dir, "compatibility_score.png"), p_e_score, width = 7.0, height = 3.0, dpi = 900, bg = "white")

ggsave(file.path(figure_dir, "interaction_term_importance.pdf"), p_f_importance, width = 3.6, height = 3.0, device = cairo_pdf)
ggsave(file.path(figure_dir, "interaction_term_importance.png"), p_f_importance, width = 3.6, height = 3.0, dpi = 900, bg = "white")

ggsave(file.path(figure_dir, "interaction_module_importance.pdf"), p_g_module, width = 3.3, height = 3.0, device = cairo_pdf)
ggsave(file.path(figure_dir, "interaction_module_importance.png"), p_g_module, width = 3.3, height = 3.0, dpi = 900, bg = "white")

ggsave(file.path(figure_dir, "incompatibility_failure_modes.pdf"), p_h_failure, width = 4.25, height = 3.0, device = cairo_pdf)
ggsave(file.path(figure_dir, "incompatibility_failure_modes.png"), p_h_failure, width = 4.25, height = 3.0, dpi = 900, bg = "white")

# -----------------------------
# 13. Export combined layouts
# -----------------------------
combined_c_to_h <- (p_c_heatmap | p_d_windows) /
  p_e_score /
  (p_f_importance | p_g_module | p_h_failure) +
  plot_layout(heights = c(1.18, 0.78, 0.92), widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = list(c("c", "d", "e", "f", "g", "h"))) &
  theme(
    plot.tag = element_text(face = "bold", size = 11, family = "sans"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(file.path(figure_dir, "SCRAP_analysis_panels.pdf"), combined_c_to_h,
       width = 11.8, height = 10.0, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(figure_dir, "SCRAP_analysis_panels.png"), combined_c_to_h,
       width = 11.8, height = 10.0, dpi = 900, bg = "white", limitsize = FALSE)

full_with_ab <- (p_a_placeholder | p_b_placeholder | p_c_heatmap) /
  (p_d_windows | p_e_score) /
  (p_f_importance | p_g_module | p_h_failure) +
  plot_layout(heights = c(0.95, 1.15, 0.90)) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e", "f", "g", "h"))) &
  theme(
    plot.tag = element_text(face = "bold", size = 11, family = "sans"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(file.path(figure_dir, "SCRAP_framework_with_placeholders.pdf"), full_with_ab,
       width = 12.2, height = 12.0, device = cairo_pdf, limitsize = FALSE)
ggsave(file.path(figure_dir, "SCRAP_framework_with_placeholders.png"), full_with_ab,
       width = 12.2, height = 12.0, dpi = 900, bg = "white", limitsize = FALSE)

message("Figure 3 analysis completed.")
message("Figure outputs: ", normalizePath(figure_dir, winslash = "/", mustWork = FALSE))
message("Table outputs: ", normalizePath(table_dir, winslash = "/", mustWork = FALSE))
