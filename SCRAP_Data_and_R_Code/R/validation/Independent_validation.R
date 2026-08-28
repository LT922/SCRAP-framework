# independent validation
# Input: data/validation/test_lib.xlsx


rm(list = ls())
graphics.off()
set.seed(123)

# -----------------------------
# 0. Script directory
# -----------------------------
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

figure_dir <- file.path(project_root, "results", "figures", "validation", "independent_validation")
table_dir <- file.path(project_root, "results", "tables", "validation", "independent_validation")

main_figure_dir <- figure_dir
si_figure_dir <- figure_dir
main_table_dir <- table_dir
si_table_dir <- table_dir
panel_dir_main <- figure_dir
panel_dir_si <- figure_dir

for (d in c(figure_dir, table_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# -----------------------------
# 1. Packages
# -----------------------------
install_missing <- TRUE
packages <- c(
  "readxl", "tidyverse", "patchwork", "scales",
  "ggrepel", "grid"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (install_missing) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    } else {
      stop("Package not installed: ", pkg, call. = FALSE)
    }
  }
}

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(ggrepel)
  library(grid)
})

# -----------------------------
# 2. Theme and colours
# -----------------------------
theme_nature <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(family = "sans", colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      axis.title = element_text(colour = "#222222"),
      axis.line = element_line(linewidth = 0.32, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.28, colour = "#222222"),
      axis.ticks.length = unit(1.4, "mm"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.title = element_text(size = base_size, colour = "#222222"),
      legend.text = element_text(size = base_size - 0.7, colour = "#222222"),
      plot.title = element_text(face = "bold", size = base_size + 1.0, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.25, hjust = 0, colour = "#555555"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 0.1, colour = "#222222"),
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

col_compatible   <- "#2F5D8C"
col_incompatible <- "#C76F5A"
col_pass         <- "#2A9D8F"
col_fail         <- "#ECECEC"
col_window       <- "#B7DCCB"
col_window_edge  <- "#7AAE92"
col_text         <- "#222222"
col_grey         <- "#777777"

outcome_cols <- c(
  "Compatible"  = col_compatible,
  "Incompatible"= col_incompatible
)

assembly_cols <- c(
  "Compatible"  = col_compatible,
  "Oversized"   = col_incompatible,
  "Precipitated"= "#9E5A4D"
)

module_cols <- c(
  "Carrier integration"     = "#6F94B9",
  "PCL anchoring"           = "#D89C45",
  "PEG-associated contacts" = "#7FAE77",
  "Prodrug self-association"= "#C96A70"
)

# -----------------------------
# 3. Segment thresholds from the trained 13-point model
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
    module = factor(module, levels = names(module_cols)),
    param_label = factor(param_labels[as.character(param)],
                         levels = rev(param_labels[params])),
    window_label = case_when(
      is.na(min_val) & !is.na(max_val) ~ paste0("\u0394G \u2264 ", max_val),
      !is.na(min_val) & !is.na(max_val) ~ paste0(min_val, "\u2264 \u0394G \u2264 ", max_val),
      TRUE ~ ""
    )
  )

in_window <- function(x, min_val, max_val) {
  ok_min <- ifelse(is.na(min_val), TRUE, x >= min_val)
  ok_max <- ifelse(is.na(max_val), TRUE, x <= max_val)
  ok_min & ok_max
}

# -----------------------------
# 4. Read test library
# -----------------------------
input_file <- file.path(project_root, "data", "validation", "test_lib.xlsx")
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

df_raw <- readxl::read_excel(input_file, sheet = "test") %>%
  rename_with(~ gsub("^\\ufeff", "", .x)) %>%
  mutate(Mol = as.character(Mol), Lib = as.character(Lib))

required_cols <- c("Lib", "Mol", "Sizes", params, "Log P", "HLB", "χ")
missing_cols <- setdiff(required_cols, colnames(df_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

df <- df_raw %>%
  mutate(
    across(all_of(c("Sizes", params, "Log P", "HLB", "χ")), as.numeric),
    is_precipitated = is.na(Sizes) | Sizes <= 0,
    assembly_type = case_when(
      is_precipitated ~ "Precipitated",
      Sizes < 30 ~ "Compatible",
      TRUE ~ "Oversized"
    ),
    outcome = if_else(assembly_type == "Compatible", "Compatible", "Incompatible"),
    outcome = factor(outcome, levels = c("Compatible", "Incompatible")),
    assembly_type = factor(assembly_type, levels = c("Compatible", "Oversized", "Precipitated")),
    size_plot = if_else(is_precipitated, max(Sizes, na.rm = TRUE) * 1.15, Sizes),
    size_label = if_else(is_precipitated, "precip.", paste0(round(Sizes, 1), "nm"))
  )

# -----------------------------
# 5. Calculate 13-point score
# -----------------------------
score_long <- df %>%
  select(Lib, Mol, Sizes, size_plot, size_label, outcome, assembly_type, all_of(params)) %>%
  pivot_longer(cols = all_of(params), names_to = "param", values_to = "energy") %>%
  mutate(param = factor(param, levels = params)) %>%
  left_join(thresholds, by = "param") %>%
  mutate(
    pass = as.integer(in_window(energy, min_val, max_val)),
    pass_label = factor(if_else(pass == 1, "Within window", "Outside window"),
                        levels = c("Outside window", "Within window")),
    module = factor(module, levels = names(module_cols)),
    param_label = factor(param_labels[as.character(param)],
                         levels = rev(param_labels[params]))
  )

score_summary <- score_long %>%
  group_by(Lib, Mol, Sizes, size_plot, size_label, outcome, assembly_type) %>%
  summarise(
    score = sum(pass, na.rm = TRUE),
    failed_terms = paste(param_labels[as.character(param[pass == 0])], collapse = ";"),
    n_failed = sum(pass == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pred_outcome = if_else(score == 13, "Compatible", "Incompatible"),
    pred_outcome = factor(pred_outcome, levels = c("Compatible", "Incompatible")),
    correct = as.character(pred_outcome) == as.character(outcome)
  )

score_long <- score_long %>%
  left_join(score_summary %>% select(Mol, score, pred_outcome, n_failed), by = "Mol")

mol_order <- score_summary %>%
  mutate(
    outcome_rank = if_else(outcome == "Compatible", 1L, 2L),
    assembly_rank = case_when(
      assembly_type == "Compatible"~ 1L,
      assembly_type == "Oversized"~ 2L,
      TRUE ~ 3L
    )
  ) %>%
  arrange(outcome_rank, assembly_rank, desc(score), Sizes, Mol) %>%
  pull(Mol)

score_summary <- score_summary %>%
  mutate(Mol = factor(Mol, levels = mol_order))

score_long <- score_long %>%
  mutate(
    Mol = factor(Mol, levels = mol_order),
    mol_display = paste0(as.character(Mol), "\n", size_label, "\n", score, "/13"),
    mol_display = factor(mol_display, levels = unique(paste0(mol_order, "\n",
                                                             score_summary$size_label[match(mol_order, as.character(score_summary$Mol))],
                                                             "\n",
                                                             score_summary$score[match(mol_order, as.character(score_summary$Mol))],
                                                             "/13")))
  )

# -----------------------------
# 6. Performance metrics
# -----------------------------
calc_metrics <- function(actual, pred) {
  actual <- factor(actual, levels = c("Compatible", "Incompatible"))
  pred <- factor(pred, levels = c("Compatible", "Incompatible"))
  
  tp <- sum(actual == "Compatible"& pred == "Compatible", na.rm = TRUE)
  tn <- sum(actual == "Incompatible"& pred == "Incompatible", na.rm = TRUE)
  fp <- sum(actual == "Incompatible"& pred == "Compatible", na.rm = TRUE)
  fn <- sum(actual == "Compatible"& pred == "Incompatible", na.rm = TRUE)
  n <- tp + tn + fp + fn
  
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
  ppv <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  npv <- ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_)
  accuracy <- (tp + tn) / n
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(denom > 0, ((tp * tn) - (fp * fn)) / denom, NA_real_)
  
  tibble(
    n = n, TP = tp, TN = tn, FP = fp, FN = fn,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    PPV = ppv,
    NPV = npv,
    MCC = mcc
  )
}

metrics_segmented <- calc_metrics(score_summary$outcome, score_summary$pred_outcome) %>%
  mutate(model = "Segmented 13/13 score") %>%
  relocate(model)

confusion_segmented <- score_summary %>%
  count(outcome, pred_outcome, name = "n") %>%
  complete(outcome = factor(c("Compatible", "Incompatible"), levels = c("Compatible", "Incompatible")),
           pred_outcome = factor(c("Compatible", "Incompatible"), levels = c("Compatible", "Incompatible")),
           fill = list(n = 0))

# -----------------------------
# 7. Module-level failed windows
# -----------------------------
failure_counts <- score_long %>%
  filter(pass == 0) %>%
  mutate(
    Mol_chr = as.character(Mol),
    module_chr = as.character(module)
  ) %>%
  count(Mol_chr, module_chr, name = "failed_terms")

failure_mol_order <- score_summary %>%
  filter(outcome == "Incompatible") %>%
  arrange(assembly_type, Sizes, desc(score), Mol) %>%
  pull(Mol) %>%
  as.character()

failure_module <- tidyr::expand_grid(
  Mol_chr = as.character(mol_order),
  module_chr = names(module_cols)
) %>%
  left_join(failure_counts, by = c("Mol_chr", "module_chr")) %>%
  mutate(failed_terms = replace_na(failed_terms, 0L)) %>%
  left_join(
    score_summary %>%
      mutate(Mol_chr = as.character(Mol)) %>%
      select(Mol_chr, outcome, assembly_type, score, Sizes),
    by = "Mol_chr"
  ) %>%
  filter(outcome == "Incompatible") %>%
  transmute(
    Mol = factor(Mol_chr, levels = rev(failure_mol_order)),
    outcome,
    assembly_type,
    score,
    Sizes,
    module = factor(module_chr, levels = names(module_cols)),
    failed_terms
  )

module_score_summary <- score_long %>%
  group_by(Mol, outcome, assembly_type, module) %>%
  summarise(module_pass_fraction = mean(pass, na.rm = TRUE), .groups = "drop")

module_validation <- module_score_summary %>%
  group_by(outcome, module) %>%
  summarise(
    mean_pass_fraction = mean(module_pass_fraction, na.rm = TRUE),
    sd_pass_fraction = sd(module_pass_fraction, na.rm = TRUE),
    n_molecules = n_distinct(Mol),
    .groups = "drop"
  ) %>%
  mutate(
    outcome = factor(outcome, levels = c("Compatible", "Incompatible")),
    module = factor(module, levels = names(module_cols))
  )

# -----------------------------
# 8. Export processed data
# -----------------------------
readr::write_csv(thresholds, file.path(main_table_dir, "segment_thresholds.csv"))
readr::write_csv(score_long, file.path(main_table_dir, "test_segment_score_long.csv"))
readr::write_csv(score_summary, file.path(main_table_dir, "test_segment_score_summary.csv"))
readr::write_csv(metrics_segmented, file.path(main_table_dir, "test_segment_score_metrics.csv"))
readr::write_csv(confusion_segmented, file.path(main_table_dir, "test_segment_confusion_matrix.csv"))
readr::write_csv(failure_module, file.path(main_table_dir, "test_failed_windows_by_module.csv"))
readr::write_csv(module_score_summary, file.path(main_table_dir, "test_module_pass_fraction.csv"))
readr::write_csv(module_validation, file.path(main_table_dir, "test_module_window_retention_by_outcome.csv"))

# -----------------------------
# 9. Helper function for saving figures
# -----------------------------
save_pdf_png_tiff <- function(plot, filename_prefix, width, height, dpi = 600) {
  ggsave(paste0(filename_prefix, ".pdf"), plot = plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
  ggsave(paste0(filename_prefix, ".png"), plot = plot, width = width, height = height,
         dpi = dpi, bg = "white")
  ggsave(paste0(filename_prefix, ".tiff"), plot = plot, width = width, height = height,
         dpi = dpi, compression = "lzw", bg = "white")
}

# -----------------------------
# 10. Main-text panels
# -----------------------------

# a: placeholder for structures
p_A <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           fill = "white", colour = "#BBBBBB", linewidth = 0.35) +
  annotate("text", x = 0.5, y = 0.56,
           label = "Test-library\nchemical structures",
           size = 3.2, fontface = "bold", lineheight = 0.92, colour = col_text) +
  annotate("text", x = 0.5, y = 0.36,
           label = "placeholder for panel a",
           size = 2.3, colour = col_grey) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(3, 4, 3, 4))

# b: score gate versus DLS size
p_B <- ggplot(score_summary, aes(x = score, y = size_plot)) +
  annotate("rect", xmin = 12.65, xmax = 13.35, ymin = 1, ymax = 30,
           fill = col_compatible, alpha = 0.08) +
  geom_hline(yintercept = 30, linetype = "dashed", linewidth = 0.32, colour = "#555555") +
  geom_vline(xintercept = 13, linetype = "dashed", linewidth = 0.32, colour = "#555555") +
  geom_point(aes(fill = assembly_type), shape = 21, size = 2.8,
             colour = "#222222", stroke = 0.25, alpha = 0.90,
             position = position_jitter(width = 0.06, height = 0)) +
  geom_text_repel(
    data = score_summary %>% filter(outcome == "Incompatible"| score < 13),
    aes(label = as.character(Mol), colour = outcome),
    size = 2.1, family = "sans", box.padding = 0.20,
    point.padding = 0.12, segment.size = 0.18,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_fill_manual(values = assembly_cols, name = NULL) +
  scale_colour_manual(values = outcome_cols, guide = "none") +
  scale_x_continuous(breaks = seq(5, 13, 1), limits = c(5, 13.35)) +
  scale_y_log10(
    breaks = c(10, 20, 30, 50, 100),
    labels = c("10", "20", "30", "50", "100"),
    limits = c(max(8, min(score_summary$size_plot, na.rm = TRUE) * 0.8),
               max(130, max(score_summary$size_plot, na.rm = TRUE) * 1.15))
  ) +
  labs(
    x = "Segmented compatibility score",
    y = "Hydrodynamic diameter (nm)",
    title = "13/13 score validates DLS outcome"
  ) +
  theme_nature(7.2) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 6.2))

# c: 13-window heatmap
p_C <- ggplot(score_long, aes(x = mol_display, y = param_label, fill = pass_label)) +
  geom_tile(colour = "white", linewidth = 0.22, width = 0.93, height = 0.90) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("Outside window"= col_fail, "Within window"= col_pass),
                    name = NULL) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL,
    title = "Test-library 13-window map"
  ) +
  theme_nature(6.8) +
  theme(
    axis.text.x = element_text(size = 5.0, lineheight = 0.86, vjust = 0.5),
    axis.text.y = element_text(size = 5.6),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 5.6),
    panel.spacing.y = unit(0.20, "lines")
  )

# d: anchoring / self-association diagnostic map
balance_df <- df %>%
  left_join(score_summary %>% select(Mol, score, outcome, assembly_type), by = c("Mol", "outcome", "assembly_type"))

p_D <- ggplot(balance_df, aes(x = `Pro-PCL`, y = `Pro-Pro`)) +
  annotate("rect", xmin = -31, xmax = -15, ymin = -28, ymax = -6.1,
           fill = col_window, colour = col_window_edge, alpha = 0.35, linewidth = 0.25) +
  geom_vline(xintercept = c(-31, -15), linetype = "dashed", linewidth = 0.25, colour = col_window_edge) +
  geom_hline(yintercept = c(-28, -6.1), linetype = "dashed", linewidth = 0.25, colour = col_window_edge) +
  geom_point(aes(fill = assembly_type), shape = 21, size = 2.7,
             colour = "#222222", stroke = 0.25, alpha = 0.88) +
  geom_text_repel(
    data = balance_df %>% filter(outcome == "Incompatible"),
    aes(label = paste0(Mol, "\n", score, "/13"), colour = outcome),
    size = 2.0, family = "sans", lineheight = 0.82,
    label.size = 0.12, label.padding = unit(0.07, "lines"),
    fill = alpha("white", 0.82),
    box.padding = 0.18, point.padding = 0.10,
    segment.size = 0.15, min.segment.length = 0,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_fill_manual(values = assembly_cols, guide = "none") +
  scale_colour_manual(values = outcome_cols, guide = "none") +
  labs(
    x = expression("drug-PCL "* Delta * "G (kcal mol"^-1 * ")"),
    y = expression("drug-drug "* Delta * "G (kcal mol"^-1 * ")"),
    title = "Pro–PCL/Pro–Pro interaction map"
  ) +
  theme_nature(7.0)

# e: failed windows by module in incompatible formulations
if (nrow(failure_module) > 0) {
  p_E <- ggplot(failure_module, aes(x = failed_terms, y = Mol, fill = module)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.18) +
    scale_fill_manual(values = module_cols, name = NULL) +
    scale_x_continuous(breaks = 0:6, expand = expansion(mult = c(0, 0.08))) +
    labs(
      x = "Number of failed windows",
      y = NULL,
      title = "Failure modes in oversized candidates"
    ) +
    theme_nature(7.0) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(size = 5.8))
} else {
  p_E <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No failed windows", size = 3) +
    theme_void()
}

# f: module window retention by observed outcome
p_F <- ggplot(module_validation,
              aes(x = module, y = 100 * mean_pass_fraction, fill = outcome)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.60,
           colour = "white", linewidth = 0.18, alpha = 0.92) +
  geom_hline(yintercept = 100, linetype = "dashed", linewidth = 0.28, colour = "#555555") +
  scale_fill_manual(values = outcome_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL,
    y = "Window retention (%)",
    title = "Module-level validation"
  ) +
  theme_nature(7.0) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.text.x = element_text(angle = 25, hjust = 1, size = 5.8)
  )

# g: confusion matrix / accuracy
acc_label <- paste0("Accuracy = ", round(100 * metrics_segmented$accuracy[[1]], 1), "%")
p_G <- ggplot(confusion_segmented, aes(x = pred_outcome, y = outcome, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = n), size = 4.6, fontface = "bold", colour = "#222222") +
  scale_fill_gradient(low = "#F2F2F2", high = "#7FAE77", name = "n") +
  labs(
    x = "Predicted outcome",
    y = "Observed DLS outcome",
    title = acc_label,
    subtitle = "Strict 13/13 rule"
  ) +
  theme_nature(7.0) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  )

# Save main panels
main_panels <- list(
  structure_placeholder = p_A,
  window_satisfaction_map = p_C,
  score_vs_DLS_size = p_B,
  module_window_retention = p_F,
  confusion_matrix = p_G,
  failed_window_decomposition = p_E,
  ProPCL_ProPro_interaction_map = p_D
)

panel_sizes <- list(
  structure_placeholder = c(3.8, 2.8),
  window_satisfaction_map = c(7.8, 3.4),
  score_vs_DLS_size = c(3.2, 2.45),
  module_window_retention = c(3.2, 2.45),
  confusion_matrix = c(2.7, 2.45),
  failed_window_decomposition = c(5.4, 2.5),
  ProPCL_ProPro_interaction_map = c(3.5, 3.5)
)

for (nm in names(main_panels)) {
  sz <- panel_sizes[[nm]]
  save_pdf_png_tiff(main_panels[[nm]],
                    file.path(panel_dir_main, nm),
                    width = sz[1], height = sz[2], dpi = 600)
}

main_design <- "
AAACCE
AAADDE
BBBBGG
FFFFGG
"

combined_main <- wrap_plots(
  A = p_A,
  B = p_C,
  C = p_B,
  D = p_F,
  E = p_G,
  F = p_E,
  G = p_D,
  design = main_design
) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 10, family = "sans"),
      plot.margin = margin(3, 3, 3, 3)
    )
  )

save_pdf_png_tiff(combined_main,
                  file.path(main_figure_dir, "independent_validation"),
                  width = 11.0, height = 8.6, dpi = 600)

# -----------------------------
# 11. SI figure: whole-molecule descriptor benchmark
# -----------------------------
descriptor_thresholds <- tribble(
  ~descriptor, ~column,    ~direction, ~boundary, ~display_label,
  "Log P",     "Log P",    "greater",     4.8,   "Log P > 4.8",
  "HLB",       "HLB",      "less",        9.8,   "HLB < 9.8",
  "\u03c7",    "\u03c7",   "less",        7.5,   "\u03c7 < 7.5",
  "Pro–Pol ΔG", "Pro-Pol", "less", -18.7, "Pro–Pol ΔG < -18.7"
)

descriptor_long <- descriptor_thresholds %>%
  mutate(data = purrr::map(column, function(col_nm) {
    df %>%
      select(Lib, Mol, Sizes, size_plot, size_label, outcome, assembly_type, all_of(col_nm)) %>%
      rename(value = all_of(col_nm))
  })) %>%
  select(-column) %>%
  tidyr::unnest(data) %>%
  mutate(
    pass = case_when(
      direction == "greater"~ value > boundary,
      direction == "less"~ value < boundary,
      TRUE ~ FALSE
    ),
    descriptor = factor(descriptor, levels = descriptor_thresholds$descriptor),
    pass_label = factor(if_else(pass, "Within descriptor range", "Outside range"),
                        levels = c("Outside range", "Within descriptor range"))
  )

descriptor_summary <- descriptor_long %>%
  group_by(Mol, outcome, assembly_type, Sizes, size_plot, size_label) %>%
  summarise(
    descriptor_score = sum(pass, na.rm = TRUE),
    failed_descriptors = paste(as.character(descriptor[!pass]), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    pred_descriptor = if_else(descriptor_score == 4, "Compatible", "Incompatible"),
    pred_descriptor = factor(pred_descriptor, levels = c("Compatible", "Incompatible")),
    Mol = factor(Mol, levels = mol_order)
  )

metrics_descriptor <- calc_metrics(descriptor_summary$outcome, descriptor_summary$pred_descriptor) %>%
  mutate(model = "Whole-molecule descriptors, 4/4 rule") %>%
  relocate(model)

confusion_descriptor <- descriptor_summary %>%
  count(outcome, pred_descriptor, name = "n") %>%
  complete(outcome = factor(c("Compatible", "Incompatible"), levels = c("Compatible", "Incompatible")),
           pred_descriptor = factor(c("Compatible", "Incompatible"), levels = c("Compatible", "Incompatible")),
           fill = list(n = 0))

readr::write_csv(descriptor_thresholds, file.path(si_table_dir, "whole_descriptor_thresholds.csv"))
readr::write_csv(descriptor_long, file.path(si_table_dir, "whole_descriptor_long.csv"))
readr::write_csv(descriptor_summary, file.path(si_table_dir, "whole_descriptor_summary.csv"))
readr::write_csv(metrics_descriptor, file.path(si_table_dir, "whole_descriptor_metrics.csv"))
readr::write_csv(confusion_descriptor, file.path(si_table_dir, "whole_descriptor_confusion_matrix.csv"))

# SI-a: descriptor scatter
p_S1 <- ggplot(descriptor_long, aes(x = value, y = size_plot)) +
  geom_hline(yintercept = 30, linetype = "dashed", linewidth = 0.28, colour = "#555555") +
  geom_vline(aes(xintercept = boundary), data = descriptor_thresholds,
             linetype = "dashed", linewidth = 0.28, colour = col_grey) +
  geom_point(aes(fill = assembly_type), shape = 21, colour = "#222222",
             stroke = 0.22, size = 2.1, alpha = 0.88) +
  facet_wrap(~ display_label, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = assembly_cols, name = NULL) +
  scale_y_log10(
    breaks = c(10, 20, 30, 50, 100),
    labels = c("10", "20", "30", "50", "100")
  ) +
  labs(
    x = "Whole-molecule descriptor value",
    y = "Hydrodynamic diameter (nm)",
    title = "Whole-molecule descriptor ranges"
  ) +
  theme_nature(7.0) +
  theme(legend.position = "top",
        legend.direction = "horizontal")

# SI-b: descriptor pass matrix
descriptor_mol_labels <- descriptor_summary %>%
  mutate(label = paste0(as.character(Mol), "\n", size_label, "\n", descriptor_score, "/4")) %>%
  select(Mol, label) %>%
  deframe()

p_S2 <- descriptor_long %>%
  mutate(
    Mol = factor(Mol, levels = mol_order),
    mol_display = factor(descriptor_mol_labels[as.character(Mol)],
                         levels = descriptor_mol_labels[mol_order])
  ) %>%
  ggplot(aes(x = mol_display, y = descriptor, fill = pass_label)) +
  geom_tile(colour = "white", linewidth = 0.24, width = 0.92, height = 0.86) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(
    values = c("Outside range"= col_fail, "Within descriptor range"= col_pass),
    name = NULL
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Whole-descriptor range map"
  ) +
  theme_nature(6.8) +
  theme(
    axis.text.x = element_text(size = 5.0, lineheight = 0.86),
    axis.text.y = element_text(size = 6.0),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# SI-c: descriptor score vs size
p_S3 <- ggplot(descriptor_summary, aes(x = descriptor_score, y = size_plot)) +
  annotate("rect", xmin = 3.65, xmax = 4.35, ymin = 1, ymax = 30,
           fill = col_compatible, alpha = 0.06) +
  geom_hline(yintercept = 30, linetype = "dashed", linewidth = 0.30, colour = "#555555") +
  geom_vline(xintercept = 4, linetype = "dashed", linewidth = 0.30, colour = "#555555") +
  geom_point(aes(fill = assembly_type), shape = 21, colour = "#222222",
             size = 2.6, stroke = 0.24, alpha = 0.90,
             position = position_jitter(width = 0.05, height = 0)) +
  geom_text_repel(
    data = descriptor_summary %>% filter(pred_descriptor != outcome),
    aes(label = as.character(Mol), colour = outcome),
    size = 2.0, family = "sans", box.padding = 0.18,
    point.padding = 0.10, segment.size = 0.15,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_fill_manual(values = assembly_cols, guide = "none") +
  scale_colour_manual(values = outcome_cols, guide = "none") +
  scale_x_continuous(breaks = 0:4, limits = c(0, 4.25)) +
  scale_y_log10(breaks = c(10, 20, 30, 50, 100),
                labels = c("10", "20", "30", "50", "100")) +
  labs(
    x = "Whole-descriptor criteria satisfied",
    y = "Hydrodynamic diameter (nm)",
    title = "4/4 descriptor rule gives false calls"
  ) +
  theme_nature(7.0)

# SI-d: descriptor confusion matrix
desc_acc_label <- paste0("Accuracy = ", round(100 * metrics_descriptor$accuracy[[1]], 1), "%")
p_S4 <- ggplot(confusion_descriptor, aes(x = pred_descriptor, y = outcome, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = n), size = 4.4, fontface = "bold", colour = "#222222") +
  scale_fill_gradient(low = "#F2F2F2", high = "#C96A70", name = "n") +
  labs(
    x = "Predicted by 4/4 descriptors",
    y = "Observed DLS outcome",
    title = desc_acc_label
  ) +
  theme_nature(7.0) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "none"
  )

# SI-e: individual descriptor pass rate by outcome
p_S5 <- descriptor_long %>%
  group_by(descriptor, outcome) %>%
  summarise(pass_rate = mean(pass, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = descriptor, y = 100 * pass_rate, fill = outcome)) +
  geom_col(position = position_dodge(width = 0.70), width = 0.62,
           colour = "white", linewidth = 0.18) +
  scale_fill_manual(values = outcome_cols, name = NULL) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL,
    y = "Pass rate (%)",
    title = "Descriptor criteria are not outcome-specific"
  ) +
  theme_nature(7.0) +
  theme(legend.position = "top")

# SI-f: segmented versus whole-descriptor score comparison
compare_scores <- score_summary %>%
  select(Mol, outcome, assembly_type, segmented_score = score) %>%
  left_join(descriptor_summary %>% select(Mol, descriptor_score), by = "Mol") %>%
  mutate(Mol = factor(Mol, levels = mol_order))

p_S6 <- ggplot(compare_scores, aes(x = descriptor_score, y = segmented_score)) +
  geom_point(aes(fill = assembly_type), shape = 21, colour = "#222222",
             size = 2.8, stroke = 0.25, alpha = 0.90) +
  geom_vline(xintercept = 4, linetype = "dashed", linewidth = 0.30, colour = "#555555") +
  geom_hline(yintercept = 13, linetype = "dashed", linewidth = 0.30, colour = "#555555") +
  geom_text_repel(
    data = compare_scores %>% filter(outcome == "Incompatible"| descriptor_score < 4),
    aes(label = as.character(Mol), colour = outcome),
    size = 1.9, family = "sans", box.padding = 0.16,
    point.padding = 0.10, segment.size = 0.15,
    min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_fill_manual(values = assembly_cols, name = NULL) +
  scale_colour_manual(values = outcome_cols, guide = "none") +
  scale_x_continuous(breaks = 0:4, limits = c(0, 4.25)) +
  scale_y_continuous(breaks = seq(5, 13, 2), limits = c(4.7, 13.35)) +
  labs(
    x = "Whole-descriptor criteria satisfied",
    y = "Segmented score",
    title = "Segmented score resolves descriptor false calls"
  ) +
  theme_nature(7.0) +
  theme(legend.position = "top")

si_panels <- list(
  descriptor_scatter = p_S1,
  descriptor_pass_matrix = p_S2,
  descriptor_score_vs_size = p_S3,
  descriptor_confusion_matrix = p_S4,
  descriptor_pass_rate = p_S5,
  segmented_vs_descriptor_score = p_S6
)

si_panel_sizes <- list(
  descriptor_scatter = c(5.4, 4.0),
  descriptor_pass_matrix = c(7.2, 2.8),
  descriptor_score_vs_size = c(3.2, 2.7),
  descriptor_confusion_matrix = c(3.0, 2.6),
  descriptor_pass_rate = c(3.5, 2.6),
  segmented_vs_descriptor_score = c(3.5, 2.8)
)

for (nm in names(si_panels)) {
  sz <- si_panel_sizes[[nm]]
  save_pdf_png_tiff(si_panels[[nm]],
                    file.path(panel_dir_si, nm),
                    width = sz[1], height = sz[2], dpi = 600)
}

combined_si <- (p_S1 | p_S2) / (p_S3 | p_S4 | p_S5) / p_S6 +
  plot_layout(heights = c(1.05, 0.82, 0.78)) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 10, family = "sans"),
      plot.margin = margin(3, 3, 3, 3)
    )
  )

save_pdf_png_tiff(combined_si,
                  file.path(si_figure_dir, "whole_molecule_descriptor_benchmark"),
                  width = 11.0, height = 8.6, dpi = 600)

# -----------------------------
# 12. Short console summary
# -----------------------------
message("Figure 4 outputs: ", normalizePath(main_figure_dir, winslash = "/", mustWork = FALSE))
message("Figure 4 tables: ", normalizePath(main_table_dir, winslash = "/", mustWork = FALSE))
message("Descriptor benchmark outputs: ", normalizePath(si_figure_dir, winslash = "/", mustWork = FALSE))
message("Main figure: independent_validation.pdf/.png/.tiff")
message("SI figure: whole_molecule_descriptor_benchmark.pdf/.png/.tiff")
message("Tables: processed score and descriptor CSV files")
