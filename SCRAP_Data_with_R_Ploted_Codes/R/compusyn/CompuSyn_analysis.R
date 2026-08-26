# CompuSyn CI provenance and supporting plots
# rev221: preserve HTML table-cell separators so CI rows parse correctly.

rm(list = ls())
graphics.off()
set.seed(123)

required_packages <- c("tidyverse", "stringr", "readr", "rvest", "xml2", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(readr)
  library(rvest)
  library(xml2)
  library(scales)
})

# -----------------------------
# 0. Paths and input
# -----------------------------
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

input_dir <- file.path(project_root, "data", "compusyn")
fig_dir <- file.path(project_root, "results", "figures", "compusyn")
table_dir <- file.path(project_root, "results", "tables", "compusyn")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

candidate_files <- c("report.html", "report.htm", "report.txt", "Pasted text.txt")
input_file <- candidate_files[file.exists(file.path(input_dir, candidate_files))][1]
if (is.na(input_file)) {
  stop("No CompuSyn report found in: ", input_dir, call. = FALSE)
}
input_path <- file.path(input_dir, input_file)
message("Reading CompuSyn report: ", input_path)

read_report_text <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("html", "htm")) {
    txt <- tryCatch({
      page <- xml2::read_html(path)
      rvest::html_text2(page)
    }, error = function(e) {
      paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    })
  } else {
    txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }
  txt %>%
    str_replace_all("\r\n|\r", "\n") %>%
    str_replace_all("\u00A0", "")
}

report_text <- read_report_text(input_path)
lines <- str_split(report_text, "\n", simplify = FALSE)[[1]] %>%
  str_replace_all("\t", " ") %>%
  str_replace_all("\u00A0", "") %>%
  str_squish()
lines <- lines[nzchar(lines)]

# -----------------------------
# 1. Basic parsing helpers
# -----------------------------
num_pat <- "[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[Ee][-+]?\\d+)?"

extract_numbers <- function(x) {
  vals <- str_extract_all(x, num_pat)[[1]]
  if (length(vals) == 0) return(numeric(0))
  as.numeric(vals)
}

first_num <- function(x) {
  if (length(x) == 0 || is.na(x)) return(NA_real_)
  vals <- extract_numbers(x)
  if (length(vals) == 0) NA_real_ else vals[1]
}

parse_ratio <- function(x) {
  m <- str_match(x, "\\[(\\d+(?:\\.\\d+)?):(\\d+(?:\\.\\d+)?)\\]")
  if (all(is.na(m))) return(c(NA_real_, NA_real_))
  c(as.numeric(m[2]), as.numeric(m[3]))
}

parse_combo_label <- function(x) {
  # Examples:
  # Drug Combo: 1/20 (1/20) (dtx+Cri [1:20])
  # Data for Drug Combo: 1/20 (dtx+Cri [1:20])
  lab <- str_match(x, "Drug Combo:\\s*([^\\(]+?)\\s*\\(")[, 2]
  if (is.na(lab)) lab <- str_match(x, "Data for Drug Combo:\\s*([^\\(]+?)\\s*\\(")[, 2]
  str_trim(lab)
}

parse_drug_label <- function(x) {
  lab <- str_match(x, "Data for Drug:\\s*([^\\[]+?)\\s*\\[")[, 2]
  str_trim(lab)
}

next_start <- function(idx, start_patterns) {
  if (idx >= length(lines)) return(length(lines) + 1)
  pat <- paste(start_patterns, collapse = "|")
  j <- which(seq_along(lines) > idx & str_detect(lines, pat))
  if (length(j) == 0) length(lines) + 1 else j[1]
}

# Ratio order: low DTX fraction to high DTX fraction.
get_ratio_levels <- function(df) {
  df %>%
    filter(type == "Combination", !is.na(ratio_a), !is.na(ratio_b)) %>%
    distinct(group, ratio_a, ratio_b) %>%
    mutate(dtx_fraction = ratio_a / (ratio_a + ratio_b)) %>%
    arrange(dtx_fraction) %>%
    pull(group)
}

# -----------------------------
# 2. Dose-effect data and median-effect parameters
# -----------------------------
block_starts <- which(str_detect(lines, "^Data for Drug:|^Data for Drug Combo:"))

dose_effect_list <- list()
model_list <- list()

for (idx in block_starts) {
  header <- lines[idx]
  is_combo <- str_detect(header, "^Data for Drug Combo:")
  end_idx <- next_start(
    idx,
    c("^Data for Drug:", "^Data for Drug Combo:", "^Dose-Effect Curve", "^Median-Effect Plot",
      "^CI Data for Drug Combo:", "^DRI Data for Drug Combo:", "^Summary Table")
  )
  block <- lines[idx:(end_idx - 1)]
  
  if (is_combo) {
    group <- parse_combo_label(header)
    ratio_vals <- parse_ratio(header)
    ratio_a <- ratio_vals[1]
    ratio_b <- ratio_vals[2]
    type <- "Combination"
  } else {
    group <- parse_drug_label(header)
    ratio_a <- NA_real_
    ratio_b <- NA_real_
    type <- "Single drug"
  }
  
  # Numeric rows in the dose-effect section.
  dat_rows <- map_dfr(block, function(l) {
    if (!str_detect(l, paste0("^", num_pat))) return(tibble())
    if (str_detect(l, "data points entered")) return(tibble())
    clean <- str_replace_all(l, "\\+", "")
    vals <- extract_numbers(clean)
    if (length(vals) < 2) return(tibble())
    tibble(dose_a = vals[1], effect = vals[2])
  })
  
  if (nrow(dat_rows) > 0) {
    dat_rows <- dat_rows %>%
      mutate(
        group = group,
        type = type,
        ratio_a = ratio_a,
        ratio_b = ratio_b,
        total_dose = if_else(
          type == "Combination"& !is.na(ratio_a) & ratio_a > 0,
          dose_a * (ratio_a + ratio_b) / ratio_a,
          dose_a
        ),
        viability_pct = (1 - effect) * 100
      ) %>%
      select(group, type, ratio_a, ratio_b, dose_a, total_dose, effect, viability_pct)
    dose_effect_list[[length(dose_effect_list) + 1]] <- dat_rows
  }
  
  model_list[[length(model_list) + 1]] <- tibble(
    group = group,
    type = type,
    ratio_a = ratio_a,
    ratio_b = ratio_b,
    m = first_num(block[str_detect(block, "^m:")][1]),
    Dm = first_num(block[str_detect(block, "^Dm:")][1]),
    r = first_num(block[str_detect(block, "^r:")][1])
  )
}

dose_effect <- bind_rows(dose_effect_list)
model_summary <- bind_rows(model_list) %>%
  filter(!is.na(group), !is.na(Dm), !is.na(m))

if (nrow(dose_effect) == 0 || nrow(model_summary) == 0) {
  stop("Failed to parse dose-effect data or median-effect parameters from the report.")
}

ratio_levels <- get_ratio_levels(model_summary)
drug_levels <- model_summary %>% filter(type == "Single drug") %>% pull(group)
all_levels <- c(drug_levels, ratio_levels)

# -----------------------------
# 3. CI curves and actual-point CI values
# -----------------------------
ci_starts <- which(str_detect(lines, "^CI Data for Drug Combo:"))
ci_curve_list <- list()
ci_actual_list <- list()

for (idx in ci_starts) {
  header <- lines[idx]
  group <- parse_combo_label(header)
  ratio_vals <- parse_ratio(header)
  ratio_a <- ratio_vals[1]
  ratio_b <- ratio_vals[2]
  
  end_idx <- next_start(idx, c("^CI Data for Drug Combo:", "^DRI Data for Drug Combo:",
                               "^Combination Index Plot", "^DRI Plot", "^Isobologram", "^Summary Table"))
  block <- lines[idx:(end_idx - 1)]
  
  state <- "none"
  for (l in block[-1]) {
    if (str_detect(l, "^Fa\\s+CI Value\\s+Total Dose")) {
      state <- "curve"
      next
    }
    if (str_detect(l, "^CI values for actual experimental points")) {
      state <- "actual_wait_header"
      next
    }
    if (state == "actual_wait_header"&& str_detect(l, "^Total Dose\\s+Fa\\s+CI Value")) {
      state <- "actual"
      next
    }
    if (!str_detect(l, paste0("^", num_pat))) next
    vals <- extract_numbers(l)
    
    if (state == "curve"&& length(vals) >= 3) {
      ci_curve_list[[length(ci_curve_list) + 1]] <- tibble(
        group = group,
        ratio_a = ratio_a,
        ratio_b = ratio_b,
        fa = vals[1],
        ci = vals[2],
        total_dose = vals[3]
      )
    }
    if (state == "actual"&& length(vals) >= 3) {
      ci_actual_list[[length(ci_actual_list) + 1]] <- tibble(
        group = group,
        ratio_a = ratio_a,
        ratio_b = ratio_b,
        total_dose = vals[1],
        fa = vals[2],
        ci = vals[3]
      )
    }
  }
}

fa_ci <- bind_rows(ci_curve_list)
actual_ci <- bind_rows(ci_actual_list)

if (nrow(fa_ci) == 0 || !"fa" %in% names(fa_ci)) {
  stop(
    "Failed to parse Fa-CI data blocks. ",
    "The CompuSyn HTML table structure was not preserved during text extraction.",
    call. = FALSE
  )
}

fa_ci <- fa_ci %>%
  mutate(group = factor(group, levels = ratio_levels)) %>%
  arrange(group, fa)

if (nrow(actual_ci) > 0 && all(c("group", "total_dose") %in% names(actual_ci))) {
  actual_ci <- actual_ci %>%
    mutate(group = factor(group, levels = ratio_levels)) %>%
    arrange(group, total_dose)
} else {
  actual_ci <- tibble(
    group = factor(character(), levels = ratio_levels),
    ratio_a = numeric(),
    ratio_b = numeric(),
    total_dose = numeric(),
    fa = numeric(),
    ci = numeric()
  )
}

# Correct CI50 extraction: one ED50 row per ratio.
ci50_summary <- fa_ci %>%
  filter(abs(fa - 0.5) < 1e-8) %>%
  transmute(
    ratio = as.character(group),
    ratio_a,
    ratio_b,
    dtx_fraction = ratio_a / (ratio_a + ratio_b),
    CI50_report = ci,
    combo_total_Dm = total_dose,
    dtx_dose_at_ED50 = total_dose * ratio_a / (ratio_a + ratio_b),
    cri_dose_at_ED50 = total_dose * ratio_b / (ratio_a + ratio_b)
  ) %>%
  arrange(dtx_fraction)

if (nrow(ci50_summary) == 0) {
  stop("No Fa = 0.5 rows were found in the CI curves. Cannot build CI50 summary.")
}

# ED-level CI summary for manuscript-level comparison.
# These values are taken directly from each ratio-specific CI curve in the CompuSyn report.
ed_levels <- c(0.50, 0.75, 0.90, 0.95)
ci_ed_summary <- fa_ci %>%
  filter(fa %in% ed_levels) %>%
  transmute(
    ratio = as.character(group),
    ratio_a,
    ratio_b,
    dtx_fraction = ratio_a / (ratio_a + ratio_b),
    ED = case_when(
      abs(fa - 0.50) < 1e-8 ~ "ED50",
      abs(fa - 0.75) < 1e-8 ~ "ED75",
      abs(fa - 0.90) < 1e-8 ~ "ED90",
      abs(fa - 0.95) < 1e-8 ~ "ED95",
      TRUE ~ paste0("Fa", fa)
    ),
    fa,
    CI = ci,
    total_dose = total_dose,
    dtx_dose = total_dose * ratio_a / (ratio_a + ratio_b),
    cri_dose = total_dose * ratio_b / (ratio_a + ratio_b)
  ) %>%
  mutate(
    ratio = factor(ratio, levels = ratio_levels),
    ED = factor(ED, levels = c("ED50", "ED75", "ED90", "ED95"))
  ) %>%
  arrange(ratio, ED)

# Wide table for direct checking and Supplementary Data export.
ci_ed_wide <- ci_ed_summary %>%
  select(ratio, ED, CI) %>%
  pivot_wider(names_from = ED, values_from = CI) %>%
  arrange(factor(ratio, levels = ratio_levels))

# Optional non-exclusive CI50 calculation for reference only.
# The manuscript figure can use CI50_report if strictly reproducing CompuSyn output.
dtx_ic50 <- model_summary %>% filter(type == "Single drug") %>% slice(1) %>% pull(Dm)
cri_ic50 <- model_summary %>% filter(type == "Single drug") %>% slice(2) %>% pull(Dm)

if (length(dtx_ic50) == 1 && length(cri_ic50) == 1 && !is.na(dtx_ic50) && !is.na(cri_ic50)) {
  ci50_summary <- ci50_summary %>%
    mutate(
      CI50_nonexclusive_calculated =
        dtx_dose_at_ED50 / dtx_ic50 +
        cri_dose_at_ED50 / cri_ic50 +
        (dtx_dose_at_ED50 * cri_dose_at_ED50) / (dtx_ic50 * cri_ic50)
    )
}

best_ratio <- ci50_summary %>% slice_min(CI50_report, n = 1, with_ties = FALSE)
message("Lowest CompuSyn-reported CI50: ", best_ratio$ratio, "= ", round(best_ratio$CI50_report, 4))

# -----------------------------
# 4. Fitted curves from median-effect parameters
# -----------------------------
fit_curve <- map_dfr(seq_len(nrow(model_summary)), function(i) {
  row <- model_summary[i, ]
  obs <- dose_effect %>% filter(group == row$group, total_dose > 0)
  if (nrow(obs) == 0 || is.na(row$Dm) || is.na(row$m) || row$Dm <= 0) return(tibble())
  min_dose <- min(obs$total_dose, na.rm = TRUE)
  max_dose <- max(obs$total_dose, na.rm = TRUE)
  if (!is.finite(min_dose) || !is.finite(max_dose) || min_dose <= 0 || max_dose <= min_dose) return(tibble())
  dose_seq <- 10 ^ seq(log10(min_dose), log10(max_dose), length.out = 200)
  fa <- 1 / (1 + (row$Dm / dose_seq)^row$m)
  tibble(
    group = row$group,
    type = row$type,
    dose = dose_seq,
    effect_fit = fa,
    viability_fit = (1 - fa) * 100
  )
})

median_data <- dose_effect %>%
  filter(effect > 0, effect < 1, total_dose > 0) %>%
  mutate(
    log_dose = log10(total_dose),
    log_fa_fu = log10(effect / (1 - effect))
  )

median_fit <- map_dfr(seq_len(nrow(model_summary)), function(i) {
  row <- model_summary[i, ]
  obs <- median_data %>% filter(group == row$group)
  if (nrow(obs) == 0 || is.na(row$m) || is.na(row$Dm) || row$Dm <= 0) return(tibble())
  xseq <- seq(min(obs$log_dose, na.rm = TRUE), max(obs$log_dose, na.rm = TRUE), length.out = 100)
  tibble(
    group = row$group,
    type = row$type,
    log_dose = xseq,
    log_fa_fu = row$m * (xseq - log10(row$Dm))
  )
})

# -----------------------------
# 5. Plot settings and save helper
# -----------------------------
ratio_palette <- c(
  "1/20"= "#C84F55",
  "1/10"= "#4F78B8",
  "1/4" = "#55A868",
  "1/1" = "#8172B2",
  "4/1" = "#CCB974",
  "10/1"= "#64B5CD"
)
# In case additional ratios exist, assign remaining colours automatically.
extra_ratios <- setdiff(ratio_levels, names(ratio_palette))
if (length(extra_ratios) > 0) {
  extra_cols <- scales::hue_pal()(length(extra_ratios))
  names(extra_cols) <- extra_ratios
  ratio_palette <- c(ratio_palette, extra_cols)
}

base_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", colour = "black"),
      plot.subtitle = element_text(colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(colour = "black"),
      panel.border = element_blank()
    )
}

save_plot <- function(p, name, width = 5.4, height = 4.0) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height, useDingbats = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 600)
}

# -----------------------------
# 6. Dose-effect plots
# -----------------------------
plot_dose_effect <- function(data_filter, fit_filter, title, subtitle) {
  ggplot() +
    geom_line(
      data = fit_curve %>% filter({{ fit_filter }}),
      aes(x = dose, y = effect_fit, colour = group),
      linewidth = 0.8
    ) +
    geom_point(
      data = dose_effect %>% filter({{ data_filter }}),
      aes(x = total_dose, y = effect, colour = group),
      size = 2.1
    ) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.1)) +
    scale_colour_manual(values = c("dtx"= "#333333", "Cri"= "#777777", ratio_palette), drop = FALSE) +
    labs(x = "Dose (µM; total dose for combinations)", y = "Fraction affected (Fa)",
         title = title, subtitle = subtitle) +
    base_theme()
}

p_dose_drugs <- plot_dose_effect(type == "Single drug", type == "Single drug",
                                 "Dose-effect curves for single-prodrug NPs",
                                 "Median-effect model fits from the CompuSyn report")
save_plot(p_dose_drugs, "CompuSyn_dose_effect_single", 5.0, 3.8)

p_dose_combos <- plot_dose_effect(type == "Combination", type == "Combination",
                                  "Dose-effect curves for fixed-ratio combinations",
                                  "Dose is total DTX-SI-C18 + Cri-SI-C18 concentration")
save_plot(p_dose_combos, "CompuSyn_dose_effect_combinations", 5.7, 4.1)

p_dose_all <- plot_dose_effect(TRUE, TRUE,
                               "Dose-effect curves for drugs and fixed-ratio combinations",
                               "Points are report values; lines are median-effect fits")
save_plot(p_dose_all, "CompuSyn_dose_effect_all", 6.2, 4.4)

# -----------------------------
# 7. Median-effect plots
# -----------------------------
plot_median <- function(data_type, title, subtitle) {
  ggplot() +
    geom_line(
      data = median_fit %>% filter(type == data_type),
      aes(x = log_dose, y = log_fa_fu, colour = group),
      linewidth = 0.8
    ) +
    geom_point(
      data = median_data %>% filter(type == data_type),
      aes(x = log_dose, y = log_fa_fu, colour = group),
      size = 2.0
    ) +
    scale_colour_manual(values = c("dtx"= "#333333", "Cri"= "#777777", ratio_palette), drop = FALSE) +
    labs(
      x = expression(log[10]~"dose"),
      y = expression(log[10]~(F[a]/F[u])),
      title = title,
      subtitle = subtitle
    ) +
    base_theme()
}

p_median_drugs <- plot_median("Single drug", "Median-effect plots for single-prodrug NPs",
                              "Linearized dose-effect relationship used for Dm and m estimation")
save_plot(p_median_drugs, "CompuSyn_median_effect_single", 4.8, 3.8)

p_median_combos <- plot_median("Combination", "Median-effect plots for fixed-ratio combinations",
                               "Each ratio was fitted independently in the CompuSyn report")
save_plot(p_median_combos, "CompuSyn_median_effect_combinations", 5.8, 4.1)

# -----------------------------
# 8. Fa-CI and log(Fa-CI) profiles
# -----------------------------
p_faci <- ggplot(fa_ci, aes(x = fa, y = ci, colour = group)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.45, colour = "grey35") +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.35, alpha = 0.85) +
  geom_point(data = fa_ci %>% filter(fa %in% ed_levels), size = 2.2) +
  scale_colour_manual(values = ratio_palette, drop = FALSE) +
  scale_x_continuous(breaks = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95), limits = c(0.03, 0.98)) +
  labs(
    x = "Fraction affected (Fa)",
    y = "Combination index (CI)",
    title = "Fa-CI profiles for tested molar ratios",
    subtitle = "All CompuSyn-reported CI values are shown; larger points mark ED50, ED75, ED90 and ED95"
  ) +
  base_theme(11)
save_plot(p_faci, "CompuSyn_Fa_CI", 6.4, 4.5)

p_falogci <- ggplot(fa_ci, aes(x = fa, y = log10(ci), colour = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "grey35") +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.35, alpha = 0.85) +
  geom_point(data = fa_ci %>% filter(fa %in% ed_levels), size = 2.2) +
  scale_colour_manual(values = ratio_palette, drop = FALSE) +
  scale_x_continuous(breaks = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95), limits = c(0.03, 0.98)) +
  labs(
    x = "Fraction affected (Fa)",
    y = expression(log[10]~"CI"),
    title = "Log-transformed Fa-CI profiles",
    subtitle = "The dashed line marks CI = 1"
  ) +
  base_theme(11)
save_plot(p_falogci, "CompuSyn_Fa_logCI", 6.4, 4.5)

# -----------------------------
# 9. CI50 ratio optimization
# -----------------------------
ci50_plot_data <- ci50_summary %>%
  mutate(
    ratio = factor(ratio, levels = ratio_levels),
    label = sprintf("%.2f", CI50_report),
    is_best = ratio == best_ratio$ratio
  )

p_ci50 <- ggplot(ci50_plot_data, aes(x = ratio, y = CI50_report, fill = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.45, colour = "grey35") +
  geom_col(width = 0.68, colour = "black", linewidth = 0.35) +
  geom_text(aes(label = label), vjust = -0.35, size = 3.5) +
  scale_fill_manual(values = ratio_palette, drop = FALSE) +
  scale_y_continuous(limits = c(0, max(1.05, max(ci50_plot_data$CI50_report, na.rm = TRUE) * 1.18))) +
  labs(
    x = "DTX-SI-C18:Cri-SI-C18 molar ratio",
    y = expression(CI[50]),
    title = expression("Ratio optimization by "* CI[50]),
    subtitle = paste0("CompuSyn-reported CI50 values; ", best_ratio$ratio, "is the lowest ratio in this report")
  ) +
  guides(fill = "none") +
  base_theme(12)
save_plot(p_ci50, "CompuSyn_CI50_ratio_optimization", 5.2, 4.0)

# -----------------------------
# 9b. All ED-level CI summaries
# -----------------------------
ci_ed_plot_data <- ci_ed_summary %>%
  mutate(label = sprintf("%.2f", CI))

p_ci_ed_heatmap <- ggplot(ci_ed_plot_data, aes(x = ED, y = ratio, fill = CI)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = label), size = 3.2, colour = "black") +
  scale_fill_gradient2(low = "#3B86C4", mid = "white", high = "#D45A5A", midpoint = 1,
                       limits = c(0, max(1.25, max(ci_ed_plot_data$CI, na.rm = TRUE))),
                       oob = scales::squish) +
  labs(
    x = "Effect level",
    y = "DTX-SI-C18:Cri-SI-C18 molar ratio",
    fill = "CI",
    title = "CI values at ED50, ED75, ED90 and ED95",
    subtitle = "Values are extracted from the ratio-specific Fa-CI curves in the CompuSyn report"
  ) +
  base_theme(11) +
  theme(panel.grid = element_blank())
save_plot(p_ci_ed_heatmap, "CompuSyn_CI_ED_heatmap", 5.6, 4.0)

p_ci_ed_bar <- ggplot(ci_ed_plot_data, aes(x = ratio, y = CI, fill = ED)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.45, colour = "grey35") +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, colour = "black", linewidth = 0.25) +
  labs(
    x = "DTX-SI-C18:Cri-SI-C18 molar ratio",
    y = "Combination index (CI)",
    fill = "Effect level",
    title = "Ratio comparison using all reported CI levels",
    subtitle = "Lower CI indicates stronger combination effect under the Chou-Talalay definition"
  ) +
  scale_fill_manual(values = c("ED50"= "#C84E52", "ED75"= "#4E76B8", "ED90"= "#55A868", "ED95"= "#8A75B6"), drop = FALSE) +
  scale_y_continuous(limits = c(0, max(1.12, max(ci_ed_plot_data$CI, na.rm = TRUE) * 1.18))) +
  base_theme(11)
save_plot(p_ci_ed_bar, "CompuSyn_CI_ED_bar", 6.2, 4.2)

# -----------------------------
# 10. Isobologram-style ED50 plot
# -----------------------------
if (length(dtx_ic50) == 1 && length(cri_ic50) == 1 && !is.na(dtx_ic50) && !is.na(cri_ic50)) {
  additive_line <- tibble(
    dtx = c(0, dtx_ic50),
    cri = c(cri_ic50, 0)
  )
  
  p_iso <- ggplot() +
    geom_line(data = additive_line, aes(x = dtx, y = cri), linetype = "dashed", colour = "grey35", linewidth = 0.65) +
    geom_point(data = ci50_plot_data, aes(x = dtx_dose_at_ED50, y = cri_dose_at_ED50, colour = ratio), size = 3.0) +
    geom_text(data = ci50_plot_data, aes(x = dtx_dose_at_ED50, y = cri_dose_at_ED50, label = ratio, colour = ratio),
              nudge_y = cri_ic50 * 0.035, size = 3.2, show.legend = FALSE) +
    annotate("point", x = dtx_ic50, y = 0, shape = 21, fill = "white", colour = "black", size = 2.8) +
    annotate("point", x = 0, y = cri_ic50, shape = 21, fill = "white", colour = "black", size = 2.8) +
    annotate("text", x = dtx_ic50, y = cri_ic50 * 0.045, label = "DTX IC50", size = 3.0, hjust = 1) +
    annotate("text", x = dtx_ic50 * 0.04, y = cri_ic50, label = "Cri IC50", size = 3.0, hjust = 0) +
    scale_colour_manual(values = ratio_palette, drop = FALSE) +
    coord_cartesian(xlim = c(0, dtx_ic50 * 1.08), ylim = c(0, cri_ic50 * 1.08), expand = FALSE) +
    labs(
      x = expression("DTX-SI-C18 dose at ED"[50]~"(µM)"),
      y = expression("Cri-SI-C18 dose at ED"[50]~"(µM)"),
      title = expression("Isobologram-style summary at ED"[50]),
      subtitle = "Points below the dashed additivity line indicate CI50 < 1"
    ) +
    base_theme(11)
  save_plot(p_iso, "CompuSyn_isobologram_ED50", 5.2, 4.4)
}

# -----------------------------
# 11. Polygonogram-style CI50 summary
# -----------------------------
poly_data <- ci50_plot_data %>%
  mutate(
    id = row_number(),
    angle = 2 * pi * (id - 1) / n(),
    x = CI50_report * sin(angle),
    y = CI50_report * cos(angle),
    x_outer = sin(angle),
    y_outer = cos(angle),
    x_lab = 1.12 * sin(angle),
    y_lab = 1.12 * cos(angle)
  )

circle_data <- tibble(theta = seq(0, 2 * pi, length.out = 361), x = sin(theta), y = cos(theta))

p_poly <- ggplot() +
  geom_path(data = circle_data, aes(x = x, y = y), linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_segment(data = poly_data, aes(x = 0, y = 0, xend = x, yend = y, colour = ratio), linewidth = 0.75) +
  geom_polygon(data = poly_data, aes(x = x, y = y, group = 1), fill = NA, colour = "grey35", linewidth = 0.5) +
  geom_point(data = poly_data, aes(x = x, y = y, colour = ratio), size = 3.2) +
  geom_text(data = poly_data, aes(x = x_lab, y = y_lab, label = ratio), size = 3.3, colour = "black") +
  geom_text(data = poly_data, aes(x = x * 1.12, y = y * 1.12, label = sprintf("%.2f", CI50_report)), size = 3.0, colour = "black") +
  scale_colour_manual(values = ratio_palette, drop = FALSE) +
  coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25), expand = FALSE) +
  labs(
    title = expression("Polygonogram-style summary of "* CI[50]),
    subtitle = "Radial distance represents CI50; lower values indicate stronger combination effect"
  ) +
  guides(colour = "none") +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", colour = "black"),
    plot.subtitle = element_text(colour = "black")
  )
save_plot(p_poly, "CompuSyn_polygonogram_CI50", 5.0, 5.0)

# -----------------------------
# 12. Export parsed data and summary
# -----------------------------
readr::write_csv(dose_effect, file.path(table_dir, "CompuSyn_dose_effect_data.csv"))
readr::write_csv(model_summary, file.path(table_dir, "CompuSyn_median_effect_parameters.csv"))
readr::write_csv(fa_ci, file.path(table_dir, "CompuSyn_Fa_CI_curve.csv"))
readr::write_csv(actual_ci, file.path(table_dir, "CompuSyn_actual_point_CI.csv"))
readr::write_csv(ci50_summary, file.path(table_dir, "CompuSyn_CI50_summary.csv"))
readr::write_csv(ci_ed_summary, file.path(table_dir, "CompuSyn_CI_ED_summary_long.csv"))
readr::write_csv(ci_ed_wide, file.path(table_dir, "CompuSyn_CI_ED_summary_wide.csv"))

summary_lines <- c(
  "CompuSyn report parsing summary",
  paste0("Input file: ", basename(input_path)),
  paste0("Single-drug Dm values: ", paste(model_summary %>% filter(type == "Single drug") %>% transmute(txt = paste0(group, "= ", signif(Dm, 5), "µM")) %>% pull(txt), collapse = "; ")),
  paste0("Lowest CompuSyn-reported CI50: ", best_ratio$ratio, "= ", signif(best_ratio$CI50_report, 5)),
  "CI50 values were extracted from the Fa = 0.5 row of each independent CI Data block.",
  "ED50/ED75/ED90/ED95 CI values are exported in ci_ed_summary_long.csv and ci_ed_summary_wide.csv.",
  "Use CI50_report to reproduce CompuSyn output. CI50_nonexclusive_calculated is provided only when single-drug Dm values are available."
)
writeLines(summary_lines, con = file.path(table_dir, "CompuSyn_summary.txt"))

message("Done.")
message("Figures saved to: ", fig_dir)
message("Tables saved to: ", table_dir)
message("CI50 summary:")
print(ci50_summary %>% select(ratio, CI50_report, combo_total_Dm, dtx_dose_at_ED50, cri_dose_at_ED50))
message("ED-level CI summary:")
print(ci_ed_wide)
