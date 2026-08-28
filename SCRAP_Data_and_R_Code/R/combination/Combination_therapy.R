# Combination screening and combination therapy
# Inputs: data/combination/combination_index.xlsx and data/combination/combination_therapy.xlsx



rm(list = ls())

# -----------------------------
# 0. Packages and paths
# -----------------------------
required_pkgs <- c("readxl", "tidyverse", "scales", "drc", "multcomp")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(scales)
  library(drc)
})

# Defensive aliases: some packages can mask dplyr verbs, especially select().
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate
summarise <- dplyr::summarise

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_path <- sub(file_arg, "", cmd_args[grepl(file_arg, cmd_args)])
  if (length(file_path) > 0 && nzchar(file_path[1])) {
    return(dirname(normalizePath(file_path[1], winslash = "/", mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && !is.null(ctx$path) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path, winslash = "/", mustWork = FALSE)))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

input_file <- file.path(project_root, "data", "combination", "combination_index.xlsx")
main_figure_dir <- file.path(project_root, "results", "figures", "combination")
out_dir <- file.path(main_figure_dir, "screening")
table_dir <- file.path(project_root, "results", "tables", "combination", "screening")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

message("Input workbook: ", input_file)
message("Screening figures: ", out_dir)
message("Screening tables: ", table_dir)

# -----------------------------
# 1. Global style and colours
# -----------------------------
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.45, colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0, size = base_size, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
}

partner_order <- c("PTX-SI-C18", "DOX-SI-C18", "Abi-SI-C18", "Cri-SI-C18", "PPT-SI-C18", "IDB-SI-C18")
partner_short <- c(
  "PTX-SI-C18"= "PTX",
  "DOX-SI-C18"= "DOX",
  "Abi-SI-C18"= "Abi",
  "Cri-SI-C18"= "Cri",
  "PPT-SI-C18"= "PPT",
  "IDB-SI-C18"= "IDB"
)

partner_cols <- c(
  "PTX-SI-C18"= "#AFC6E9",
  "DOX-SI-C18"= "#C7B6E8",
  "Abi-SI-C18"= "#E9B6E8",
  "Cri-SI-C18"= "#D5E9B6",
  "PPT-SI-C18"= "#92C6D8",
  "IDB-SI-C18"= "#A99BC7"
)

ratio_order <- c("1/20", "1/10", "1/4", "1/1", "4/1", "10/1")
ratio_cols <- c(
  "1/20"= "#C94E52",
  "1/10"= "#4C76B8",
  "1/4" = "#56AA68",
  "1/1" = "#8372B6",
  "4/1" = "#CCBA71",
  "10/1"= "#67B7CF"
)

# Optional manual correction for the ratio-screening curve sheet.
# Leave as NULL for automatic parsing. If the workbook has ambiguous duplicated
# ratio headers, set this vector to the intended block names in left-to-right order.
# Example: RATIO_CURVE_BLOCK_NAME_OVERRIDE <- c("1/20", "1/10", "1/1", "1/4", "10/1")
RATIO_CURVE_BLOCK_NAME_OVERRIDE <- NULL
# For ratio-screening dose-response sheets, each ratio should be represented by one block of three replicate columns.
# If the sheet contains duplicated ratio headers, later duplicate blocks are ignored by default to avoid plotting block 1/block 2/block 3 as separate curves.
DROP_DUPLICATE_RATIO_CURVE_BLOCKS <- TRUE


save_plot <- function(p, name, width, height, dpi = 600) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height, useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height, dpi = dpi)
}

mean_sd_label <- function(mean_value, sd_value, digits = 1) {
  ifelse(
    is.na(mean_value),
    "NA",
    paste0(format(round(mean_value, digits), nsmall = digits), "± ", format(round(sd_value, digits), nsmall = digits))
  )
}

# -----------------------------
# 2. Helpers for workbook sheets
# -----------------------------
read_replicate_sheet <- function(sheet_name, value_name) {
  raw <- readxl::read_excel(input_file, sheet = sheet_name, col_names = FALSE)
  raw <- raw %>% filter(!if_all(everything(), ~ is.na(.x)))
  if (ncol(raw) < 2) stop("Sheet ", sheet_name, "must contain at least two columns.")
  
  long <- raw %>%
    rename(item = 1) %>%
    mutate(item = as.character(item)) %>%
    pivot_longer(
      cols = -item,
      names_to = "replicate_raw",
      values_to = "value_raw"
    ) %>%
    mutate(
      replicate = as.integer(str_extract(replicate_raw, "\\d+")),
      value = suppressWarnings(as.numeric(value_raw)),
      value_text = as.character(value_raw),
      value_type = if_else(is.na(value) & !is.na(value_text), "text", "numeric"),
      sheet = sheet_name
    ) %>%
    filter(!is.na(item), item != "", !(is.na(value) & (is.na(value_text) | value_text == "NA"))) %>%
    dplyr::select(sheet, item, replicate, value, value_text, value_type) %>%
    rename(!!value_name := value)
  
  long
}

summarise_replicates <- function(dat, value_col) {
  value_col <- rlang::ensym(value_col)
  dat %>%
    group_by(item) %>%
    summarise(
      n = sum(!is.na(!!value_col)),
      mean = mean(!!value_col, na.rm = TRUE),
      sd = sd(!!value_col, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mean = if_else(is.nan(mean), NA_real_, mean),
      sd = if_else(is.na(sd), 0, sd)
    )
}

# -----------------------------
# 3. Read IC50, single viability and combination viability
# -----------------------------
ic50_long <- read_replicate_sheet("IC50", "ic50_uM") %>%
  mutate(
    item_clean = case_when(
      str_detect(item, regex("PTX", ignore_case = TRUE)) ~ "PTX-SI-C18 NPs",
      str_detect(item, regex("DTX", ignore_case = TRUE)) ~ "DTX-SI-C18 NPs",
      str_detect(item, regex("Abi", ignore_case = TRUE)) ~ "Abi-SI-C18 NPs",
      str_detect(item, regex("Dox", ignore_case = TRUE)) ~ "Dox-SI-C18 NPs",
      str_detect(item, regex("Cri", ignore_case = TRUE)) ~ "Cri-SI-C18 NPs",
      str_detect(item, regex("PPT", ignore_case = TRUE)) ~ "PPT-SI-C18 NPs",
      str_detect(item, regex("IDB", ignore_case = TRUE)) ~ "IDB-SI-C18 NPs",
      str_detect(item, regex("Dex", ignore_case = TRUE)) ~ "Dex-SI-C18 NPs",
      TRUE ~ item
    )
  )

ic50_order <- c("PTX-SI-C18 NPs", "DTX-SI-C18 NPs", "Abi-SI-C18 NPs", "Dox-SI-C18 NPs",
                "Cri-SI-C18 NPs", "PPT-SI-C18 NPs", "IDB-SI-C18 NPs", "Dex-SI-C18 NPs")

ic50_summary <- ic50_long %>%
  group_by(item_clean) %>%
  summarise(
    n = sum(!is.na(ic50_uM)),
    mean = mean(ic50_uM, na.rm = TRUE),
    sd = sd(ic50_uM, na.rm = TRUE),
    text_note = paste(na.omit(unique(value_text[value_type == "text"])), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    mean = if_else(is.nan(mean), NA_real_, mean),
    sd = if_else(is.na(sd), 0, sd),
    item_clean = factor(item_clean, levels = ic50_order),
    label = if_else(is.na(mean) & nzchar(text_note), text_note, sprintf("%.1f", mean))
  ) %>%
  arrange(item_clean)

single_long <- read_replicate_sheet("single", "viability_pct") %>%
  mutate(
    item = factor(item, levels = partner_order),
    condition = "Single NPs"
  )

combination_long <- read_replicate_sheet("combination", "viability_pct") %>%
  mutate(
    item = factor(item, levels = partner_order),
    condition = "DTX-SI-C18 co-loaded NPs"
  )

viability_long <- bind_rows(single_long, combination_long) %>%
  mutate(
    condition = factor(condition, levels = c("Single NPs", "DTX-SI-C18 co-loaded NPs")),
    short = recode(as.character(item), !!!partner_short)
  )

viability_summary <- viability_long %>%
  group_by(item, condition) %>%
  summarise(
    n = sum(!is.na(viability_pct)),
    mean = mean(viability_pct, na.rm = TRUE),
    sd = sd(viability_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = mean_sd_label(mean, sd, digits = 1),
    item = factor(item, levels = partner_order)
  )

# Single/combination viability ratio
# Calculation used for Figure 6:
#   1) For each partner prodrug, calculate the mean viability of the single-prodrug NP group.
#   2) Use this single-group mean as the common numerator.
#   3) Divide it by EACH replicate value of the corresponding DTX-SI-C18 co-loaded group.
#   4) Plot the resulting replicate-level ratios and summarize them as mean ± SD.
#
# This is intentionally NOT a replicate-paired single_i / combination_i calculation.

single_reference <- single_long %>%
  group_by(item) %>%
  summarise(
    single_mean = mean(viability_pct, na.rm = TRUE),
    n_single = sum(!is.na(viability_pct)),
    .groups = "drop"
  )

ratio_long <- combination_long %>%
  dplyr::select(
    item,
    replicate,
    combination_viability = viability_pct
  ) %>%
  left_join(single_reference, by = "item") %>%
  mutate(
    ratio = single_mean / combination_viability,
    item = factor(as.character(item), levels = partner_order),
    short = recode(as.character(item), !!!partner_short)
  )

ratio_summary <- ratio_long %>%
  group_by(item, short) %>%
  summarise(
    n = sum(!is.na(ratio)),
    single_mean_reference = first(single_mean),
    combination_mean = mean(combination_viability, na.rm = TRUE),
    mean = mean(ratio, na.rm = TRUE),
    sd = sd(ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%.2f", mean))

# -----------------------------
# 4. Read CI sheet
# -----------------------------
ci_raw <- readxl::read_excel(input_file, sheet = "CI", col_names = FALSE)
ci_raw <- ci_raw %>% filter(!if_all(everything(), ~ is.na(.x)))

ci_ratio_header <- as.character(unlist(ci_raw[1, -1]))
ci_values <- ci_raw[-1, ]

ci_long <- ci_values %>%
  rename(effect_level = 1) %>%
  mutate(effect_level = as.character(effect_level)) %>%
  pivot_longer(cols = -effect_level, names_to = "col", values_to = "CI") %>%
  mutate(
    col_index = as.integer(str_extract(col, "\\d+")),
    ratio = ci_ratio_header[col_index - 1],
    CI = suppressWarnings(as.numeric(CI))
  ) %>%
  filter(!is.na(effect_level), effect_level != "", !is.na(CI), !is.na(ratio), ratio != "") %>%
  group_by(effect_level, ratio) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  mutate(
    ratio = factor(ratio, levels = ratio_order),
    effect_level = factor(effect_level, levels = c("ED50", "ED75", "ED90", "ED95"))
  )

ci_summary <- ci_long %>%
  group_by(ratio, effect_level) %>%
  summarise(
    n = sum(!is.na(CI)),
    mean = mean(CI, na.rm = TRUE),
    sd = sd(CI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sd = if_else(is.na(sd), 0, sd),
    label = sprintf("%.2f", mean),
    label_mean_sd = sprintf("%.2f ± %.2f", mean, sd),
    ratio = factor(as.character(ratio), levels = ratio_order),
    effect_level = factor(as.character(effect_level), levels = c("ED50", "ED75", "ED90", "ED95"))
  )

ci_ed50 <- ci_long %>% filter(effect_level == "ED50")
ci_ed50_summary <- ci_summary %>% filter(effect_level == "ED50")

# -----------------------------
# 5. Export parsed data
# -----------------------------
readr::write_csv(ic50_long, file.path(table_dir, "parsed_ic50_replicates.csv"))
readr::write_csv(ic50_summary, file.path(table_dir, "parsed_ic50_summary.csv"))
readr::write_csv(viability_long, file.path(table_dir, "parsed_viability_replicates.csv"))
readr::write_csv(viability_summary, file.path(table_dir, "parsed_viability_summary.csv"))
readr::write_csv(ratio_long, file.path(table_dir, "single_combination_ratio_replicates.csv"))
readr::write_csv(ratio_summary, file.path(table_dir, "single_combination_ratio_summary.csv"))
readr::write_csv(ci_long, file.path(table_dir, "ci_replicates_long.csv"))
readr::write_csv(ci_summary, file.path(table_dir, "ci_summary_long.csv"))

# -----------------------------
# 5b. Read and plot representative dose-response curves
#     Newly added workbook sheets:
#       singledrugic50curve      : representative single-prodrug NP dose-response data
#       ratioscreeningic50curve  : representative fixed-ratio DTX-SI-C18/Cri-SI-C18 dose-response data
# -----------------------------
normalize_curve_name <- function(x) {
  x <- as.character(x)
  case_when(
    str_detect(x, regex("DTX", ignore_case = TRUE)) ~ "DTX-SI-C18",
    str_detect(x, regex("PTX", ignore_case = TRUE)) ~ "PTX-SI-C18",
    str_detect(x, regex("DOX|Dox", ignore_case = TRUE)) ~ "DOX-SI-C18",
    str_detect(x, regex("Abi", ignore_case = TRUE)) ~ "Abi-SI-C18",
    str_detect(x, regex("Cri", ignore_case = TRUE)) ~ "Cri-SI-C18",
    str_detect(x, regex("PPT", ignore_case = TRUE)) ~ "PPT-SI-C18",
    str_detect(x, regex("IDB", ignore_case = TRUE)) ~ "IDB-SI-C18",
    str_detect(x, regex("Dex", ignore_case = TRUE)) ~ "Dex-SI-C18",
    TRUE ~ x
  )
}

single_curve_order <- c("DTX-SI-C18", "PTX-SI-C18", "DOX-SI-C18", "Abi-SI-C18",
                        "Cri-SI-C18", "PPT-SI-C18", "IDB-SI-C18", "Dex-SI-C18")

single_curve_cols <- c(
  "DTX-SI-C18"= "#D9A441",
  "PTX-SI-C18"= "#AFC6E9",
  "DOX-SI-C18"= "#C7B6E8",
  "Abi-SI-C18"= "#E9B6E8",
  "Cri-SI-C18"= "#D5E9B6",
  "PPT-SI-C18"= "#92C6D8",
  "IDB-SI-C18"= "#A99BC7",
  "Dex-SI-C18"= "#D5D5D5"
)

read_triplicate_curve_sheet <- function(sheet_name, type = c("single", "ratio")) {
  type <- match.arg(type)
  raw <- readxl::read_excel(input_file, sheet = sheet_name, col_names = FALSE, col_types = "text")
  raw <- raw %>% dplyr::filter(!if_all(everything(), ~ is.na(.x) | trimws(as.character(.x)) == ""))
  if (nrow(raw) < 3 || ncol(raw) < 4) {
    stop("Sheet ", sheet_name, "does not contain the expected triplicate dose-response structure.")
  }
  
  dose_vec <- suppressWarnings(as.numeric(unlist(raw[-1, 1])))
  
  # Use non-empty first-row cells as the start of each merged triplicate block.
  # This is more robust than assuming 2,5,8,... when a header is shifted or merged.
  header_all <- trimws(as.character(unlist(raw[1, ])))
  block_starts <- which(!is.na(header_all) & header_all != "")
  block_starts <- block_starts[block_starts != 1]
  block_starts <- block_starts[block_starts + 2 <= ncol(raw)]
  if (length(block_starts) == 0) {
    stop("No curve blocks were detected in sheet ", sheet_name, ".")
  }
  
  block_names_raw <- header_all[block_starts]
  block_names <- block_names_raw
  
  if (type == "single") {
    block_names <- normalize_curve_name(block_names)
  }
  
  if (type == "ratio") {
    # If an explicit override is supplied, use it. The override should contain one
    # label per detected three-column block.
    if (!is.null(RATIO_CURVE_BLOCK_NAME_OVERRIDE)) {
      if (length(RATIO_CURVE_BLOCK_NAME_OVERRIDE) != length(block_names)) {
        stop("RATIO_CURVE_BLOCK_NAME_OVERRIDE has ", length(RATIO_CURVE_BLOCK_NAME_OVERRIDE),
             "entries, but ", length(block_names), "curve blocks were detected in ", sheet_name, ".")
      }
      block_names <- RATIO_CURVE_BLOCK_NAME_OVERRIDE
    } else if (anyDuplicated(block_names) > 0) {
      duplicated_labels <- unique(block_names[duplicated(block_names)])
      if (isTRUE(DROP_DUPLICATE_RATIO_CURVE_BLOCKS)) {
        keep_blocks <- !duplicated(block_names)
        message("Warning: duplicated ratio headers were detected in ", sheet_name,
                ": ", paste(duplicated_labels, collapse = ", "),
                ". Later duplicate blocks were ignored so each ratio is plotted as one n = 3 curve. ",
                "Correct the Excel headers or set RATIO_CURVE_BLOCK_NAME_OVERRIDE if these blocks are distinct ratios.")
        block_starts <- block_starts[keep_blocks]
        block_names_raw <- block_names_raw[keep_blocks]
        block_names <- block_names[keep_blocks]
      } else {
        block_names <- ave(block_names, block_names, FUN = function(z) {
          if (length(z) == 1) z else paste0(z, "block ", seq_along(z))
        })
      }
    }
  }
  
  diagnostics <- tibble(
    sheet = sheet_name,
    block_index = seq_along(block_starts),
    start_column = block_starts,
    raw_header = block_names_raw,
    parsed_group = block_names,
    replicate_columns = "three columns per plotted ratio"
  )
  readr::write_csv(diagnostics, file.path(table_dir, paste0(sheet_name, "_block_diagnostics.csv")))
  
  block_list <- purrr::map2(seq_along(block_starts), block_names, function(i, group_name) {
    start_col <- block_starts[[i]]
    cols <- start_col + 0:2
    value_mat <- raw[-1, cols, drop = FALSE]
    names(value_mat) <- paste0("rep", seq_len(ncol(value_mat)))
    
    value_mat %>%
      dplyr::mutate(
        concentration_uM = dose_vec,
        dplyr::across(starts_with("rep"), ~ as.character(.x))
      ) %>%
      tidyr::pivot_longer(
        cols = starts_with("rep"),
        names_to = "replicate_raw",
        values_to = "viability_raw",
        values_transform = list(viability_raw = as.character)
      ) %>%
      dplyr::mutate(
        replicate = as.integer(stringr::str_extract(replicate_raw, "\\d+")),
        group = group_name,
        viability_pct = suppressWarnings(as.numeric(viability_raw)),
        sheet = sheet_name
      ) %>%
      dplyr::filter(!is.na(concentration_uM), !is.na(viability_pct)) %>%
      dplyr::select(sheet, group, concentration_uM, replicate, viability_pct)
  })
  
  dplyr::bind_rows(block_list)
}

curve_summary_by_dose <- function(dat) {
  dat %>%
    group_by(group, concentration_uM) %>%
    summarise(
      n = sum(!is.na(viability_pct)),
      mean = mean(viability_pct, na.rm = TRUE),
      sd = sd(viability_pct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(sd = if_else(is.na(sd), 0, sd))
}

fit_dose_response_curve <- function(dat, group_name) {
  dat_fit <- dat %>%
    filter(group == group_name, concentration_uM > 0, !is.na(viability_pct))
  
  if (nrow(dat_fit) < 5 || length(unique(dat_fit$concentration_uM)) < 5) {
    return(list(pred = tibble(), summary = tibble(group = group_name, fitted_IC50_uM = NA_real_, fit_status = "insufficient data")))
  }
  
  fit <- tryCatch(
    drc::drm(
      viability_pct ~ concentration_uM,
      data = dat_fit,
      fct = drc::LL.4(names = c("slope", "lower", "upper", "IC50"))
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(list(pred = tibble(), summary = tibble(group = group_name, fitted_IC50_uM = NA_real_, fit_status = fit$message)))
  }
  
  dose_grid <- exp(seq(log(min(dat_fit$concentration_uM, na.rm = TRUE)),
                       log(max(dat_fit$concentration_uM, na.rm = TRUE)),
                       length.out = 250))
  pred <- tryCatch(
    tibble(
      group = group_name,
      concentration_uM = dose_grid,
      fitted_viability = as.numeric(predict(fit, newdata = data.frame(concentration_uM = dose_grid)))
    ),
    error = function(e) tibble()
  )
  
  ic50_val <- tryCatch(as.numeric(drc::ED(fit, 50, type = "absolute", display = FALSE)[1, 1]), error = function(e) NA_real_)
  
  list(
    pred = pred,
    summary = tibble(group = group_name, fitted_IC50_uM = ic50_val, fit_status = "ok")
  )
}

format_ic50_annotation <- function(drug_name) {
  lookup <- case_when(
    drug_name == "DTX-SI-C18"~ "DTX-SI-C18 NPs",
    drug_name == "PTX-SI-C18"~ "PTX-SI-C18 NPs",
    drug_name == "DOX-SI-C18"~ "Dox-SI-C18 NPs",
    drug_name == "Abi-SI-C18"~ "Abi-SI-C18 NPs",
    drug_name == "Cri-SI-C18"~ "Cri-SI-C18 NPs",
    drug_name == "PPT-SI-C18"~ "PPT-SI-C18 NPs",
    drug_name == "IDB-SI-C18"~ "IDB-SI-C18 NPs",
    drug_name == "Dex-SI-C18"~ "Dex-SI-C18 NPs",
    TRUE ~ drug_name
  )
  row <- ic50_summary %>% filter(as.character(item_clean) == lookup) %>% slice(1)
  if (nrow(row) == 0) return("IC50 = NA")
  if (is.na(row$mean[1]) && nzchar(row$text_note[1])) {
    return(paste0("IC50 ", row$text_note[1], "μM"))
  }
  paste0("IC50 = ", sprintf("%.1f ± %.1f μM", row$mean[1], row$sd[1]))
}

# Read new sheets if present. The rest of the script continues if either optional sheet is absent.
wb_sheets <- readxl::excel_sheets(input_file)

if ("singledrugic50curve"%in% wb_sheets) {
  single_curve_long <- read_triplicate_curve_sheet("singledrugic50curve", type = "single") %>%
    mutate(group = factor(group, levels = single_curve_order)) %>%
    filter(!is.na(group))
  
  single_curve_summary <- curve_summary_by_dose(single_curve_long) %>%
    mutate(group = factor(group, levels = single_curve_order))
  
  single_curve_fit <- map(single_curve_order, ~ fit_dose_response_curve(single_curve_long, .x))
  single_curve_pred <- map_dfr(single_curve_fit, "pred") %>%
    mutate(group = factor(group, levels = single_curve_order))
  single_curve_fit_summary <- map_dfr(single_curve_fit, "summary") %>%
    mutate(group = factor(group, levels = single_curve_order))
  
  single_curve_ann <- tibble(
    group = factor(single_curve_order, levels = single_curve_order),
    label = map_chr(single_curve_order, format_ic50_annotation),
    x = 0.06,
    y = 117
  )
  
  readr::write_csv(single_curve_long, file.path(table_dir, "single_drug_ic50_curve_replicates.csv"))
  readr::write_csv(single_curve_summary, file.path(table_dir, "single_drug_ic50_curve_summary.csv"))
  readr::write_csv(single_curve_fit_summary, file.path(table_dir, "single_drug_ic50_curve_fit_summary.csv"))
  
  p_single_curves <- ggplot() +
    geom_point(
      data = single_curve_long,
      aes(x = concentration_uM, y = viability_pct, colour = group),
      size = 1.3, alpha = 0.45
    ) +
    geom_errorbar(
      data = single_curve_summary,
      aes(x = concentration_uM, ymin = mean - sd, ymax = mean + sd, colour = group),
      width = 0, linewidth = 0.25, alpha = 0.85
    ) +
    geom_point(
      data = single_curve_summary,
      aes(x = concentration_uM, y = mean, fill = group),
      shape = 21, size = 2.0, colour = "black", stroke = 0.25
    ) +
    geom_line(
      data = single_curve_pred,
      aes(x = concentration_uM, y = fitted_viability, colour = group),
      linewidth = 0.75, na.rm = TRUE
    ) +
    geom_text(
      data = single_curve_ann,
      aes(x = x, y = y, label = label),
      hjust = 0, vjust = 1, size = 2.65, colour = "black"
    ) +
    scale_x_log10(
      breaks = c(0.05, 0.5, 5, 50, 500, 1000),
      labels = c("0.05", "0.5", "5", "50", "500", "1000"),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30), expand = expansion(mult = c(0.01, 0.03))) +
    scale_colour_manual(values = single_curve_cols, guide = "none") +
    scale_fill_manual(values = single_curve_cols, guide = "none") +
    facet_wrap(~ group, ncol = 4) +
    labs(
      x = "Concentration (μM)",
      y = "C4-2 cell viability (%)",
      title = "Representative single-prodrug NP dose-response curves"
    ) +
    theme_pub(base_size = 9) +
    theme(strip.text = element_text(face = "bold", size = 9))
  
  save_plot(p_single_curves, "single_drug_IC50_curves_combined", 8.6, 5.8)
  
  # Individual panels are also exported for flexible figure assembly.
  for (g in single_curve_order) {
    p_one <- ggplot() +
      geom_point(
        data = single_curve_long %>% filter(as.character(group) == g),
        aes(x = concentration_uM, y = viability_pct),
        size = 1.5, alpha = 0.55, colour = single_curve_cols[[g]]
      ) +
      geom_errorbar(
        data = single_curve_summary %>% filter(as.character(group) == g),
        aes(x = concentration_uM, ymin = mean - sd, ymax = mean + sd),
        width = 0, linewidth = 0.28, colour = single_curve_cols[[g]]
      ) +
      geom_point(
        data = single_curve_summary %>% filter(as.character(group) == g),
        aes(x = concentration_uM, y = mean),
        shape = 21, fill = single_curve_cols[[g]], colour = "black", stroke = 0.25, size = 2.2
      ) +
      geom_line(
        data = single_curve_pred %>% filter(as.character(group) == g),
        aes(x = concentration_uM, y = fitted_viability),
        linewidth = 0.8, colour = single_curve_cols[[g]], na.rm = TRUE
      ) +
      annotate("text", x = 0.06, y = 116, label = format_ic50_annotation(g), hjust = 0, vjust = 1, size = 3.0) +
      scale_x_log10(breaks = c(0.05, 0.5, 5, 50, 500, 1000), labels = c("0.05", "0.5", "5", "50", "500", "1000")) +
      scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30), expand = expansion(mult = c(0.01, 0.03))) +
      labs(x = "Concentration (μM)", y = "C4-2 cell viability (%)", title = g) +
      theme_pub(base_size = 10)
    
    safe_name <- str_replace_all(g, "[^A-Za-z0-9]+", "_")
    save_plot(p_one, paste0("single_drug_IC50_curve_", safe_name), 3.4, 3.0)
  }
}

if ("ratioscreeningic50curve"%in% wb_sheets) {
  ratio_curve_long <- read_triplicate_curve_sheet("ratioscreeningic50curve", type = "ratio") %>%
    dplyr::filter(!is.na(group), !is.na(concentration_uM), !is.na(viability_pct))
  
  # Keep only groups actually present in the curve sheet. This prevents empty legend
  # entries such as 4/1 when that ratio is not present in the representative curve sheet.
  ratio_groups_present <- unique(as.character(ratio_curve_long$group))
  ratio_curve_long <- ratio_curve_long %>%
    dplyr::mutate(group = factor(group, levels = ratio_groups_present))
  
  ratio_curve_summary <- curve_summary_by_dose(ratio_curve_long) %>%
    dplyr::mutate(group = factor(group, levels = ratio_groups_present))
  
  # Build a palette for present groups. For corrected labels such as "10/1 block 2",
  # use the colour of the base ratio before "block".
  base_ratio_name <- function(x) stringr::str_replace(as.character(x), "block \\d+$", "")
  ratio_curve_cols <- setNames(
    vapply(ratio_groups_present, function(g) {
      base <- base_ratio_name(g)
      if (base %in% names(ratio_cols)) ratio_cols[[base]] else "#808080"
    }, character(1)),
    ratio_groups_present
  )
  
  # Fit one representative curve per ratio using the dose-level means. Raw
  # triplicates are still shown as light points; the displayed fitted curve and
  # mean ± SD symbols correspond to one n = 3 curve per ratio.
  ratio_curve_fit_input <- ratio_curve_summary %>%
    dplyr::transmute(group, concentration_uM, viability_pct = mean)
  
  ratio_curve_fit <- purrr::map(ratio_groups_present, ~ fit_dose_response_curve(ratio_curve_fit_input, .x))
  ratio_curve_pred <- purrr::map_dfr(ratio_curve_fit, "pred") %>%
    dplyr::mutate(group = factor(group, levels = ratio_groups_present))
  ratio_curve_fit_summary <- purrr::map_dfr(ratio_curve_fit, "summary") %>%
    dplyr::mutate(group = factor(group, levels = ratio_groups_present))
  
  readr::write_csv(ratio_curve_long, file.path(table_dir, "ratio_screening_ic50_curve_replicates.csv"))
  readr::write_csv(ratio_curve_summary, file.path(table_dir, "ratio_screening_ic50_curve_summary.csv"))
  readr::write_csv(ratio_curve_fit_input, file.path(table_dir, "ratio_screening_ic50_curve_fit_input_mean.csv"))
  readr::write_csv(ratio_curve_fit_summary, file.path(table_dir, "ratio_screening_ic50_curve_fit_summary.csv"))
  
  # Cleaner representative curve plot: each ratio is shown as one fitted curve
  # from the n = 3 dose-level means. Light points indicate individual repeats,
  # and prominent points with error bars show mean ± SD at each concentration.
  p_ratio_curves <- ggplot() +
    geom_point(
      data = ratio_curve_long,
      aes(x = concentration_uM, y = viability_pct, colour = group),
      size = 0.75, alpha = 0.20,
      position = position_jitter(width = 0.02, height = 0, seed = 11),
      show.legend = FALSE
    ) +
    geom_errorbar(
      data = ratio_curve_summary,
      aes(x = concentration_uM, ymin = mean - sd, ymax = mean + sd, colour = group),
      width = 0, linewidth = 0.32, alpha = 0.90,
      show.legend = FALSE
    ) +
    geom_line(
      data = ratio_curve_pred,
      aes(x = concentration_uM, y = fitted_viability, colour = group),
      linewidth = 1.0, na.rm = TRUE
    ) +
    geom_point(
      data = ratio_curve_summary,
      aes(x = concentration_uM, y = mean, fill = group),
      shape = 21, size = 2.2, colour = "black", stroke = 0.25
    ) +
    scale_x_log10(
      breaks = c(0.05, 0.15, 0.5, 1.5, 5, 15, 50, 150, 300, 500, 750, 1000),
      labels = c("0.05", "0.15", "0.5", "1.5", "5", "15", "50", "150", "300", "500", "750", "1000"),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30), expand = expansion(mult = c(0.01, 0.03))) +
    scale_colour_manual(values = ratio_curve_cols, drop = TRUE) +
    scale_fill_manual(values = ratio_curve_cols, drop = TRUE) +
    labs(
      x = "Total prodrug concentration (μM)",
      y = "C4-2 cell viability (%)",
      title = "Representative fixed-ratio DTX-SI-C18/Cri-SI-C18 dose-response curves"
    ) +
    theme_pub(base_size = 10) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  save_plot(p_ratio_curves, "ratio_screening_IC50_curves", 6.8, 4.2)
}



# -----------------------------
# 6. Figure: single-prodrug IC50 values
# -----------------------------
ic50_plot_data <- ic50_summary %>% filter(!is.na(item_clean))
ic50_point_data <- ic50_long %>%
  mutate(item_clean = factor(item_clean, levels = ic50_order)) %>%
  filter(!is.na(item_clean), !is.na(ic50_uM))

ic50_max <- max(ic50_plot_data$mean + ic50_plot_data$sd, na.rm = TRUE)
if (!is.finite(ic50_max)) ic50_max <- 500

p_ic50 <- ggplot(ic50_plot_data, aes(x = item_clean, y = mean)) +
  geom_col(data = ic50_plot_data %>% filter(!is.na(mean)), fill = "#D6ECF4", colour = "black", width = 0.68, linewidth = 0.35) +
  geom_errorbar(data = ic50_plot_data %>% filter(!is.na(mean)), aes(ymin = mean - sd, ymax = mean + sd), width = 0.18, linewidth = 0.35) +
  geom_point(data = ic50_point_data, aes(y = ic50_uM), shape = 21, size = 2.2, stroke = 0.35, fill = "white", colour = "#4A4A4A", position = position_jitter(width = 0.08, height = 0)) +
  geom_text(aes(y = if_else(is.na(mean), ic50_max * 0.08, mean + sd + ic50_max * 0.045), label = label), size = 3.2, angle = 90, hjust = 0, na.rm = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16)), limits = c(0, ic50_max * 1.22)) +
  labs(x = NULL, y = expression(IC[50]~"("*mu*"M)"), title = "Single-prodrug NP activity") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1))

save_plot(p_ic50, "single_np_ic50", 4.6, 3.8)

# -----------------------------
# 7. Figure: viability heatmap
# -----------------------------
p_viability_heatmap <- ggplot(viability_summary, aes(x = condition, y = item, fill = mean)) +
  geom_tile(colour = "white", linewidth = 1.0) +
  geom_text(aes(label = label), size = 3.5, colour = "black") +
  scale_fill_gradient2(
    low = "#CFEAF3",
    mid = "#F7F4ED",
    high = "#CC6B54",
    midpoint = 50,
    limits = c(0, 100),
    name = "Cell viability\n(%)"
  ) +
  labs(x = NULL, y = "Partner prodrug", title = "IC50-based combination screening") +
  theme_pub(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 9, colour = "black"),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position = "right"
  )

save_plot(p_viability_heatmap, "viability_heatmap", 4.6, 3.4)

# -----------------------------
# 8. Figure: single/combination viability ratio
# -----------------------------
p_ratio <- ggplot(ratio_summary, aes(x = item, y = mean, fill = item)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#8A8A8A", linewidth = 0.35) +
  geom_col(width = 0.58, colour = "black", linewidth = 0.35, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.16, linewidth = 0.35) +
  geom_point(data = ratio_long, aes(x = item, y = ratio, fill = item), shape = 21, size = 2.2, stroke = 0.35, colour = "#4A4A4A", position = position_jitter(width = 0.08, height = 0)) +
  geom_text(aes(y = mean + sd + 0.10, label = label), size = 3.1) +
  scale_fill_manual(values = partner_cols, guide = "none") +
  scale_x_discrete(labels = partner_short) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, max(ratio_summary$mean + ratio_summary$sd, na.rm = TRUE) * 1.22)) +
  labs(x = "DTX-SI-C18/partner combination", y = "Single/combination viability ratio", title = "Partner-dependent viability reduction") +
  theme_pub(base_size = 10)

save_plot(p_ratio, "single_combination_viability_ratio", 4.4, 3.4)

# -----------------------------
# Dunnett multiple-comparisons helper
# -----------------------------
run_dunnett_vs_1_20 <- function(dat, effect_label) {
  dat_test <- dat %>%
    dplyr::filter(!is.na(CI), as.character(ratio) %in% ratio_order) %>%
    dplyr::mutate(
      ratio = factor(as.character(ratio), levels = ratio_order)
    )
  
  if (nrow(dat_test) == 0 || dplyr::n_distinct(dat_test$ratio) < 2) {
    return(tibble::tibble())
  }
  
  fit <- stats::aov(CI ~ ratio, data = dat_test)
  overall_p <- tryCatch(
    summary(fit)[[1]][["Pr(>F)"]][1],
    error = function(e) NA_real_
  )
  
  glht_fit <- multcomp::glht(
    fit,
    linfct = multcomp::mcp(ratio = "Dunnett")
  )
  glht_summary <- summary(
    glht_fit,
    test = multcomp::adjusted("single-step")
  )
  
  comparison <- names(glht_summary$test$coefficients)
  group2 <- stringr::str_replace(
    comparison,
    stringr::fixed(" - 1/20"),
    ""
  )
  
  tibble::tibble(
    effect_level = effect_label,
    group1 = "1/20",
    group2 = group2,
    estimate = as.numeric(glht_summary$test$coefficients),
    p_value = as.numeric(glht_summary$test$pvalues),
    overall_ANOVA_p = overall_p,
    method = "Ordinary one-way ANOVA followed by Dunnett's multiple-comparisons test"
  ) %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        is.na(p_value) ~ "",
        p_value < 0.0001 ~ "****",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

# -----------------------------
# 9. Figure: ED50 CI ratio optimization
# -----------------------------
# Ordinary one-way ANOVA followed by Dunnett's multiple-comparisons test.
# The 1/20 ratio is the control group; all other ratios are compared only with 1/20.
pairwise_ed50 <- run_dunnett_vs_1_20(
  ci_ed50,
  effect_label = "ED50"
)

readr::write_csv(
  pairwise_ed50,
  file.path(table_dir, "ci_ed50_dunnett_vs_1_20.csv")
)

p_ci50 <- ggplot(ci_ed50_summary, aes(x = ratio, y = mean, fill = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#6A6A6A", linewidth = 0.4) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.35, alpha = 0.90) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15, linewidth = 0.35) +
  geom_point(data = ci_ed50, aes(x = ratio, y = CI, fill = ratio), shape = 21, size = 2.4, colour = "#4A4A4A", stroke = 0.35, position = position_jitter(width = 0.08, height = 0)) +
  geom_text(aes(y = mean + sd + 0.055, label = sprintf("%.2f", mean)), size = 3.2) +
  scale_fill_manual(values = ratio_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, max(1.05, max(ci_ed50_summary$mean + ci_ed50_summary$sd, na.rm = TRUE) * 1.22))) +
  labs(x = "DTX-SI-C18:Cri-SI-C18 molar ratio", y = expression(CI[50]), title = expression("Ratio optimization by "*CI[50])) +
  theme_pub(base_size = 10)

save_plot(p_ci50, "CI50_ratio_optimization", 4.6, 3.5)

# -----------------------------
# 9b. Figure: CI50 statistics with pairwise significance
# -----------------------------
# This separate panel is intended for the ratio-optimization result.
# It shows all individual CI50 replicate values, mean ± SD and pairwise
# comparisons of the selected 1/20 ratio versus the other tested ratios.
ci50_y_base <- max(1.05, max(ci_ed50_summary$mean + ci_ed50_summary$sd, na.rm = TRUE))

sig_ci50 <- pairwise_ed50 %>%
  filter(!is.na(p_value)) %>%
  arrange(match(group2, ratio_order)) %>%
  mutate(
    group1 = factor(group1, levels = ratio_order),
    group2 = factor(group2, levels = ratio_order),
    x_mid = (as.numeric(group1) + as.numeric(group2)) / 2,
    y = ci50_y_base + 0.08 * row_number(),
    y_tip = y - 0.025,
    label = significance
  )
readr::write_csv(sig_ci50, file.path(table_dir, "ci_ed50_dunnett_for_plot.csv"))

ci50_plot_stats <- ci_ed50_summary %>%
  mutate(
    ratio_chr = as.character(ratio),
    ratio_index = as.numeric(factor(ratio_chr, levels = ratio_order)),
    label_y = mean + sd + 0.045
  )

ci50_point_stats <- ci_ed50 %>%
  mutate(
    ratio_chr = as.character(ratio),
    ratio_index = as.numeric(factor(ratio_chr, levels = ratio_order))
  )

ci50_y_limit <- max(ci50_y_base + 0.08 * (nrow(sig_ci50) + 1),
                    max(ci50_plot_stats$label_y, na.rm = TRUE) + 0.05,
                    1.15)

p_ci50_stats <- ggplot(ci50_plot_stats, aes(x = ratio_index, y = mean, fill = ratio_chr)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#6A6A6A", linewidth = 0.4) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.35, alpha = 0.90) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.15, linewidth = 0.35) +
  geom_point(
    data = ci50_point_stats,
    aes(x = ratio_index, y = CI, fill = ratio_chr),
    shape = 21, size = 2.5, colour = "#4A4A4A", stroke = 0.35,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  geom_text(aes(y = label_y, label = sprintf("%.2f", mean)), size = 3.2) +
  geom_segment(
    data = sig_ci50,
    aes(x = as.numeric(group1), xend = as.numeric(group2), y = y, yend = y),
    inherit.aes = FALSE, linewidth = 0.35
  ) +
  geom_segment(
    data = sig_ci50,
    aes(x = as.numeric(group1), xend = as.numeric(group1), y = y, yend = y_tip),
    inherit.aes = FALSE, linewidth = 0.35
  ) +
  geom_segment(
    data = sig_ci50,
    aes(x = as.numeric(group2), xend = as.numeric(group2), y = y, yend = y_tip),
    inherit.aes = FALSE, linewidth = 0.35
  ) +
  geom_text(
    data = sig_ci50,
    aes(x = x_mid, y = y + 0.025, label = label),
    inherit.aes = FALSE, size = 3.2
  ) +
  scale_x_continuous(breaks = seq_along(ratio_order), labels = ratio_order) +
  scale_fill_manual(values = ratio_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)), limits = c(0, ci50_y_limit)) +
  labs(
    x = "DTX-SI-C18:Cri-SI-C18 molar ratio",
    y = expression(CI[50]),
    title = expression(CI[50]~"statistics for ratio optimization"),
    subtitle = "Bars show mean ± SD; points show individual values. Statistics: one-way ANOVA with Dunnett's test versus 1/20."
  ) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_ci50_stats, "CI50_statistics_with_significance", 5.4, 4.0)

# -----------------------------
# 10. Figure: CI values at ED50, ED75, ED90 and ED95 heatmap
# -----------------------------
p_ci_heatmap <- ggplot(ci_summary, aes(x = effect_level, y = ratio, fill = mean)) +
  geom_tile(colour = "white", linewidth = 1.0) +
  geom_text(aes(label = label_mean_sd), size = 3.1, colour = "black") +
  scale_fill_gradient2(
    low = "#8FB6DD",
    mid = "#F7F7F7",
    high = "#EBA6A1",
    midpoint = 1,
    limits = c(0, max(1.45, max(ci_summary$mean, na.rm = TRUE))),
    name = "CI"
  ) +
  labs(x = "Effect level", y = "DTX-SI-C18:Cri-SI-C18 molar ratio", title = "CI values across effect levels") +
  theme_pub(base_size = 10) +
  theme(panel.background = element_rect(fill = "white", colour = NA))

save_plot(p_ci_heatmap, "CI_ED_heatmap", 5.2, 3.6)

# -----------------------------
# 11. Figure: CI values at ED levels, showing replicate points
#     plus pairwise significance versus 1/20 within each ED level
# -----------------------------
# For each effect level, perform an ordinary one-way ANOVA across all ratios,
# followed by Dunnett's multiple-comparisons test using 1/20 as the control.
pairwise_all_ed <- purrr::map_dfr(
  c("ED50", "ED75", "ED90", "ED95"),
  function(ed) {
    run_dunnett_vs_1_20(
      ci_long %>%
        dplyr::filter(as.character(effect_level) == ed),
      effect_label = ed
    )
  }
) %>%
  dplyr::mutate(
    effect_level = factor(
      effect_level,
      levels = c("ED50", "ED75", "ED90", "ED95")
    ),
    group1 = factor(group1, levels = ratio_order),
    group2 = factor(group2, levels = ratio_order)
  )

readr::write_csv(
  pairwise_all_ed,
  file.path(table_dir, "ci_all_ed_dunnett_vs_1_20.csv")
)

# Annotation positions are calculated separately for each ED facet.
ci_ed_y_base <- ci_summary %>%
  group_by(effect_level) %>%
  summarise(
    y_base = max(mean + sd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_base = if_else(is.finite(y_base), pmax(y_base, 1.02), 1.02)
  )

sig_all_ed <- pairwise_all_ed %>%
  filter(!is.na(p_value)) %>%
  left_join(ci_ed_y_base, by = "effect_level") %>%
  group_by(effect_level) %>%
  arrange(match(as.character(group2), ratio_order), .by_group = TRUE) %>%
  mutate(
    row_id = row_number(),
    y = y_base + 0.085 * row_id,
    y_tip = y - 0.026,
    x1 = as.numeric(group1),
    x2 = as.numeric(group2),
    x_mid = (x1 + x2) / 2,
    label = significance
  ) %>%
  ungroup()

readr::write_csv(sig_all_ed, file.path(table_dir, "ci_all_ed_dunnett_for_plot.csv"))

ci_summary_bar <- ci_summary %>%
  mutate(
    ratio_chr = as.character(ratio),
    ratio_index = as.numeric(factor(ratio_chr, levels = ratio_order)),
    label_y = mean + sd + 0.045
  )

ci_long_bar <- ci_long %>%
  mutate(
    ratio_chr = as.character(ratio),
    ratio_index = as.numeric(factor(ratio_chr, levels = ratio_order))
  )

ci_ed_y_limit <- max(
  max(ci_summary_bar$label_y, na.rm = TRUE) + 0.05,
  max(sig_all_ed$y, na.rm = TRUE) + 0.08,
  1.15,
  na.rm = TRUE
)
if (!is.finite(ci_ed_y_limit)) ci_ed_y_limit <- 1.5

p_ci_ed_bar <- ggplot(ci_summary_bar, aes(x = ratio_index, y = mean, fill = ratio_chr)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#6A6A6A", linewidth = 0.35) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.25, alpha = 0.88) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.14, linewidth = 0.30) +
  geom_point(
    data = ci_long_bar,
    aes(x = ratio_index, y = CI, fill = ratio_chr),
    shape = 21, size = 1.9, colour = "#4A4A4A", stroke = 0.25,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  geom_text(aes(y = label_y, label = sprintf("%.2f", mean)), size = 2.5, vjust = 0) +
  geom_segment(
    data = sig_all_ed,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE, linewidth = 0.25
  ) +
  geom_segment(
    data = sig_all_ed,
    aes(x = x1, xend = x1, y = y, yend = y_tip),
    inherit.aes = FALSE, linewidth = 0.25
  ) +
  geom_segment(
    data = sig_all_ed,
    aes(x = x2, xend = x2, y = y, yend = y_tip),
    inherit.aes = FALSE, linewidth = 0.25
  ) +
  geom_text(
    data = sig_all_ed,
    aes(x = x_mid, y = y + 0.020, label = label),
    inherit.aes = FALSE, size = 2.35
  ) +
  scale_x_continuous(breaks = seq_along(ratio_order), labels = ratio_order) +
  scale_fill_manual(values = ratio_cols, guide = "none") +
  facet_wrap(~ effect_level, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)), limits = c(0, ci_ed_y_limit)) +
  labs(
    x = "DTX-SI-C18:Cri-SI-C18 molar ratio",
    y = "Combination index (CI)",
    title = "CI profiles for ratio optimization",
    subtitle = "Bars show mean ± SD; points show individual values. Within each ED level: one-way ANOVA with Dunnett's test versus 1/20."
  ) +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_ci_ed_bar, "CI_ED_bar_with_replicates", 9.6, 3.8)

# -----------------------------
# 12. Figure: ED profiles by ratio, mean ± SD
# -----------------------------
p_ci_profile <- ggplot(ci_summary, aes(x = effect_level, y = mean, group = ratio, colour = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#6A6A6A", linewidth = 0.35) +
  geom_line(linewidth = 0.75) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.08, linewidth = 0.35) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = ratio_cols) +
  labs(x = "Effect level", y = "Combination index (CI)", title = "Effect-level CI profiles") +
  theme_pub(base_size = 10) +
  theme(legend.position = "right")

save_plot(p_ci_profile, "CI_ED_profile", 5.2, 3.6)

# -----------------------------
# 13. Write compact summary text
# -----------------------------
best_ci50 <- ci_ed50_summary %>% arrange(mean) %>% slice(1)
best_ratio <- as.character(best_ci50$ratio[1])
best_value <- best_ci50$mean[1]

summary_lines <- c(
  paste0("Input workbook: ", input_file),
  paste0("Output folder: ", out_dir),
  "",
  "Key outputs:",
  "single_np_ic50.pdf/png",
  "viability_heatmap.pdf/png",
  "single_combination_viability_ratio.pdf/png",
  "CI50_ratio_optimization.pdf/png",
  "CI50_statistics_with_significance.pdf/png (one-way ANOVA + Dunnett vs 1/20)",
  "CI_ED_heatmap.pdf/png (mean ± SD labels)",
  "CI_ED_bar_with_replicates.pdf/png (one-way ANOVA + Dunnett vs 1/20 within each ED level)",
  "CI_ED_profile.pdf/png",
  "single_drug_IC50_curves_combined.pdf/png",
  "single_drug_IC50_curve_<drug>.pdf/png",
  "ratio_screening_IC50_curves.pdf/png",
  "",
  paste0("Lowest mean CI50 ratio: ", best_ratio, "(mean CI50 = ", sprintf("%.3f", best_value), ")"),
  "CI < 1 was interpreted according to the Chou-Talalay definition."
)
writeLines(summary_lines, con = file.path(table_dir, "analysis_summary.txt"))

message("Done. Lowest mean CI50 ratio: ", best_ratio, "(", sprintf("%.3f", best_value), ")")
message("Screening figures saved to: ", out_dir)
message("Screening tables saved to: ", table_dir)


# ============================================================================
# END OF MODULE 1 / START OF MODULE 2
# The following cleanup reproduces running the second original script in a
# fresh R session; it does not alter any plotting/statistical code below.
# ============================================================================
rm(list = ls())
graphics.off()

# Combination-therapy module
# Combination therapy figure generation with raw data points overlaid on statistical panels.
# Input: data/combination/combination_therapy.xlsx
# Outputs are written to results/figures/main/Fig6/therapy and results/tables/Fig6/therapy.
#

options(stringsAsFactors = FALSE)

# -----------------------------
# 0. Packages and paths
# -----------------------------
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "stringr",
  "forcats", "patchwork", "readr", "purrr", "scales"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grepl(file_arg, cmd_args)])
  if (length(script_path) > 0 && nzchar(script_path[1])) {
    return(normalizePath(dirname(script_path[1]), winslash = "/", mustWork = TRUE))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(normalizePath(dirname(ctx$path), winslash = "/", mustWork = TRUE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

base_dir <- get_script_dir()
project_root <- normalizePath(file.path(base_dir, "..", ".."), winslash = "/", mustWork = TRUE)

input_file <- file.path(project_root, "data", "combination", "combination_therapy.xlsx")
screening_input_file <- file.path(project_root, "data", "combination", "combination_index.xlsx")
ic50_input_file <- screening_input_file
main_figure_dir <- file.path(project_root, "results", "figures", "combination")
out_dir <- file.path(main_figure_dir, "therapy")
table_dir <- file.path(project_root, "results", "tables", "combination", "therapy")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}
if (!file.exists(screening_input_file)) {
  stop("Screening input file not found: ", screening_input_file)
}

required_screening_sheets <- c(
  "IC50",
  "singledrugic50curve",
  "ratioscreeningic50curve",
  "single",
  "combination",
  "CI"
)
missing_screening_sheets <- setdiff(
  required_screening_sheets,
  readxl::excel_sheets(screening_input_file)
)
if (length(missing_screening_sheets) > 0) {
  stop(
    "Missing sheet(s) in combination_index.xlsx: ",
    paste(missing_screening_sheets, collapse = ", "),
    call. = FALSE
  )
}

required_therapy_sheets <- c("tumor volume", "body weight")
missing_therapy_sheets <- setdiff(
  required_therapy_sheets,
  readxl::excel_sheets(input_file)
)
if (length(missing_therapy_sheets) > 0) {
  stop(
    "Missing sheet(s) in combination_therapy.xlsx: ",
    paste(missing_therapy_sheets, collapse = ", "),
    call. = FALSE
  )
}

dir.create(main_figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

bar_error_type <- "SD"
curve_error_type <- "SEM"
plot_base_family <- "sans"

# -----------------------------
# 1. Labels, colors, and theme
# -----------------------------
clean_file_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

clean_drug <- function(x) {
  x2 <- stringr::str_to_upper(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x2, "PTX") ~ "PTX",
    stringr::str_detect(x2, "DTX|D2") ~ "DTX",
    stringr::str_detect(x2, "DOX") ~ "Dox",
    stringr::str_detect(x2, "ABI") ~ "Abi",
    stringr::str_detect(x2, "CRI") ~ "Cri",
    stringr::str_detect(x2, "PPT") ~ "PPT",
    stringr::str_detect(x2, "IDB") ~ "IDB",
    stringr::str_detect(x2, "DEX") ~ "Dex",
    TRUE ~ as.character(x)
  )
}

clean_treatment <- function(x) {
  x2 <- stringr::str_to_upper(as.character(x))
  dplyr::case_when(
    stringr::str_detect(x2, "^PBS") ~ "PBS",
    stringr::str_detect(x2, "NPS-SEQ|SEQ") ~ "DTX-SI-C18 NPs + Cri-SI-C18 NPs",
    stringr::str_detect(x2, "NPS-COM|CO-LOADED|COLOADED|COM") ~ "DTX-SI-C18/Cri-SI-C18 NPs",
    stringr::str_detect(x2, "CRISIC18|CRI-SI") ~ "Cri-SI-C18 NPs",
    stringr::str_detect(x2, "D2SIC18|DTXSIC18|DTX-SI") ~ "DTX-SI-C18 NPs",
    TRUE ~ as.character(x)
  )
}

drug_levels_ic50 <- c("PTX", "DTX", "Abi", "Dox", "Cri", "PPT", "IDB", "Dex")
drug_levels_screen <- c("PTX", "Dox", "Abi", "Cri", "PPT", "IDB")
ratio_levels <- c("1/20", "1/10", "1/4", "1/1", "4/1", "10/1")
treatment_levels <- c(
  "PBS",
  "DTX-SI-C18 NPs",
  "Cri-SI-C18 NPs",
  "DTX-SI-C18 NPs + Cri-SI-C18 NPs",
  "DTX-SI-C18/Cri-SI-C18 NPs"
)

col_axis <- "#222222"
col_grey_0 <- "#FFFFFF"
col_grey_1 <- "#FAFAF8"
col_grey_2 <- "#E8E3D8"
col_grey_3 <- "#8C867A"
col_panel_bg <- "#FCFBF8"
col_line <- "#D8D2C6"
col_blue <- "#4E79A7"
col_blue_light <- "#DCE7F4"
col_teal <- "#2F9C95"
col_teal_dark <- "#217C76"
col_orange <- "#D99537"
col_green <- "#67A567"
col_red <- "#C95C57"
col_beige <- "#E6DED1"
col_selected <- col_teal
col_heat_low <- "#D9EEF7"
col_heat_mid <- "#F4EFE6"
col_heat_high <- "#A8422D"
col_heat_high_mid <- "#D99077"

condition_cols <- c(
  "Single prodrug NP"= col_blue,
  "DTX-SI-C18 combination"= col_teal
)

treatment_cols <- c(
  "PBS"= "#3F3F3F",
  "DTX-SI-C18 NPs"= col_blue,
  "Cri-SI-C18 NPs"= col_orange,
  "DTX-SI-C18 NPs + Cri-SI-C18 NPs"= col_green,
  "DTX-SI-C18/Cri-SI-C18 NPs"= col_red
)

raw_point_style <- list(
  shape = 21,
  size = 1.45,
  stroke = 0.28,
  fill = "white",
  color = col_axis,
  alpha = 0.90
)

make_jitter <- function(width = 0.12, height = 0, seed = 1) {
  ggplot2::position_jitter(width = width, height = height, seed = seed)
}

theme_nature <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size, base_family = plot_base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = col_axis),
      axis.title = ggplot2::element_text(size = base_size + 0.5),
      axis.text = ggplot2::element_text(size = base_size, color = col_axis),
      axis.line = ggplot2::element_line(linewidth = 0.35, color = col_axis),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, color = col_axis),
      axis.ticks.length = grid::unit(1.4, "mm"),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.18, color = "#E8E6DE"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = base_size, color = col_axis),
      legend.text = ggplot2::element_text(size = base_size, color = col_axis),
      legend.key.size = grid::unit(3.0, "mm"),
      legend.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = ggplot2::margin(4, 5, 4, 5)
    )
}

theme_heatmap <- function(base_size = 7) {
  ggplot2::theme_minimal(base_size = base_size, base_family = plot_base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = col_axis),
      axis.title = ggplot2::element_text(size = base_size + 0.5),
      axis.text = ggplot2::element_text(size = base_size, color = col_axis),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size),
      legend.key.height = grid::unit(8, "mm"),
      legend.key.width = grid::unit(2.4, "mm"),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.margin = ggplot2::margin(4, 5, 4, 5)
    )
}

mean_summary <- function(df, group_cols, value_col) {
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      N = sum(!is.na(.data[[value_col]])),
      Mean = mean(.data[[value_col]], na.rm = TRUE),
      SD = stats::sd(.data[[value_col]], na.rm = TRUE),
      SEM = SD / sqrt(N),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      SD = dplyr::if_else(is.na(SD), 0, SD),
      SEM = dplyr::if_else(is.na(SEM), 0, SEM)
    )
}

get_error <- function(df, error_type = "SD") {
  if (toupper(error_type) == "SEM") df$SEM else df$SD
}

save_plot <- function(plot_obj, file_stub, width, height) {
  pdf_file <- file.path(out_dir, paste0(file_stub, ".pdf"))
  png_file <- file.path(out_dir, paste0(file_stub, ".png"))
  tryCatch(
    ggplot2::ggsave(pdf_file, plot_obj, width = width, height = height, units = "in", device = grDevices::cairo_pdf, bg = "white"),
    error = function(e) ggplot2::ggsave(pdf_file, plot_obj, width = width, height = height, units = "in", device = "pdf", bg = "white")
  )
  ggplot2::ggsave(png_file, plot_obj, width = width, height = height, units = "in", dpi = 900, bg = "white")
}

# -----------------------------
# 2. Data import helpers
# -----------------------------
read_screening_triplicates <- function(sheet_name) {
  raw <- readxl::read_excel(
    screening_input_file,
    sheet = sheet_name,
    col_names = FALSE,
    col_types = "text"
  ) %>%
    dplyr::filter(
      !dplyr::if_all(
        dplyr::everything(),
        ~ is.na(.x) | trimws(as.character(.x)) == ""
      )
    )
  
  if (ncol(raw) < 4) {
    stop(
      "Sheet '", sheet_name,
      "' in combination_index.xlsx must contain a formulation column and three replicate columns.",
      call. = FALSE
    )
  }
  
  raw <- raw[, 1:4, drop = FALSE]
  names(raw) <- c("ItemRaw", "rep1", "rep2", "rep3")
  
  raw %>%
    dplyr::mutate(ItemRaw = stringr::str_squish(as.character(ItemRaw))) %>%
    tidyr::pivot_longer(
      cols = c(rep1, rep2, rep3),
      names_to = "ReplicateRaw",
      values_to = "ValueRaw"
    ) %>%
    dplyr::mutate(
      Replicate = as.integer(stringr::str_extract(ReplicateRaw, "\\d+")),
      Value = suppressWarnings(as.numeric(ValueRaw))
    ) %>%
    dplyr::filter(
      !is.na(ItemRaw),
      ItemRaw != "",
      !is.na(Value),
      is.finite(Value)
    ) %>%
    dplyr::select(ItemRaw, Replicate, Value)
}

read_ci_sheet <- function() {
  raw <- readxl::read_excel(
    screening_input_file,
    sheet = "CI",
    col_names = FALSE,
    col_types = "text"
  ) %>%
    dplyr::filter(
      !dplyr::if_all(
        dplyr::everything(),
        ~ is.na(.x) | trimws(as.character(.x)) == ""
      )
    )
  
  if (nrow(raw) < 2 || ncol(raw) < 4) {
    stop("Sheet 'CI' in combination_index.xlsx has an unexpected structure.", call. = FALSE)
  }
  
  ratio_header <- stringr::str_squish(
    as.character(unlist(raw[1, -1, drop = TRUE], use.names = FALSE))
  )
  effect_levels <- stringr::str_squish(as.character(raw[[1]][-1]))
  values <- raw[-1, -1, drop = FALSE]
  
  rep_index <- ave(
    seq_along(ratio_header),
    ratio_header,
    FUN = seq_along
  )
  
  purrr::map_dfr(seq_along(effect_levels), function(i) {
    vals <- suppressWarnings(
      as.numeric(unlist(values[i, ], use.names = FALSE))
    )
    
    tibble::tibble(
      Effect_level = effect_levels[i],
      DTX_Cri_Ratio = ratio_header,
      Replicate = rep_index,
      CI = vals
    )
  }) %>%
    dplyr::filter(
      !is.na(Effect_level),
      Effect_level != "",
      DTX_Cri_Ratio %in% ratio_levels,
      !is.na(CI),
      is.finite(CI)
    ) %>%
    dplyr::mutate(
      Effect_level = factor(
        Effect_level,
        levels = c("ED50", "ED75", "ED90", "ED95")
      ),
      DTX_Cri_Ratio = factor(
        DTX_Cri_Ratio,
        levels = ratio_levels
      )
    )
}

parse_longitudinal_sheet <- function(sheet_name, value_name = "Value") {
  raw <- readxl::read_excel(input_file, sheet = sheet_name, col_names = FALSE)
  days <- suppressWarnings(as.numeric(raw[[1]][-1]))
  groups_raw <- as.character(unlist(raw[1, -1], use.names = FALSE))
  value_mat <- raw[-1, -1, drop = FALSE]
  rep_df <- tibble::tibble(Column = seq_along(groups_raw), GroupRaw = groups_raw) %>%
    dplyr::filter(!is.na(GroupRaw), GroupRaw != "") %>%
    dplyr::group_by(GroupRaw) %>%
    dplyr::mutate(Replicate = dplyr::row_number()) %>%
    dplyr::ungroup()
  purrr::map_dfr(rep_df$Column, function(j) {
    values <- suppressWarnings(as.numeric(unlist(value_mat[, j], use.names = FALSE)))
    group_clean <- clean_treatment(groups_raw[j])
    tibble::tibble(
      Day = days,
      GroupRaw = groups_raw[j],
      Group = factor(group_clean, levels = treatment_levels),
      Replicate = rep_df$Replicate[rep_df$Column == j],
      Mouse = paste0(clean_file_name(group_clean), "_", rep_df$Replicate[rep_df$Column == j]),
      "{value_name}" := values
    )
  }) %>%
    dplyr::filter(!is.na(Day), !is.na(.data[[value_name]]), !is.na(Group))
}

# -----------------------------
# 3. Import and process data
# -----------------------------
# b. IC50.
# Source: data/combination/combination_index.xlsx, sheet "IC50".
# Sheet structure: column A = NP name; columns B-D = three IC50 replicates.
# Text entries such as ">1000" are treated as not determined within the tested range.
ic50_raw <- readxl::read_excel(
  ic50_input_file,
  sheet = "IC50",
  col_names = FALSE,
  col_types = "text"
)

ic50_raw <- ic50_raw %>%
  dplyr::filter(!dplyr::if_all(dplyr::everything(), ~ is.na(.x) | trimws(as.character(.x)) == ""))

if (nrow(ic50_raw) < 2 || ncol(ic50_raw) < 4) {
  stop(
    "Sheet 'IC50' in combination_index.xlsx must contain column A with NP names and columns B-D with three IC50 replicates.",
    call. = FALSE
  )
}

ic50_body <- ic50_raw[-1, 1:4, drop = FALSE]
names(ic50_body) <- c("NP", "rep1", "rep2", "rep3")

ic50_long_module2 <- ic50_body %>%
  dplyr::mutate(NP = stringr::str_squish(as.character(NP))) %>%
  tidyr::pivot_longer(
    cols = c(rep1, rep2, rep3),
    names_to = "Replicate_raw",
    values_to = "IC50_raw"
  ) %>%
  dplyr::mutate(
    Replicate = as.integer(stringr::str_extract(Replicate_raw, "\\d+")),
    IC50_raw = stringr::str_squish(as.character(IC50_raw)),
    IC50_value = suppressWarnings(as.numeric(IC50_raw)),
    Is_censored = !is.na(IC50_raw) & stringr::str_detect(IC50_raw, "^>"),
    Censor_threshold = suppressWarnings(
      as.numeric(stringr::str_extract(IC50_raw, "\\d+\\.?\\d*"))
    ),
    Drug = factor(clean_drug(NP), levels = drug_levels_ic50)
  ) %>%
  dplyr::filter(!is.na(NP), NP != "", !is.na(Drug))

ic50_points <- ic50_long_module2 %>%
  dplyr::filter(!is.na(IC50_value), is.finite(IC50_value)) %>%
  dplyr::select(
    NP,
    Drug,
    Replicate,
    IC50_value
  )

ic50_data <- ic50_long_module2 %>%
  dplyr::group_by(NP, Drug) %>%
  dplyr::summarise(
    N_numeric = sum(!is.na(IC50_value)),
    IC50_mean_numeric = ifelse(
      N_numeric > 0,
      mean(IC50_value, na.rm = TRUE),
      NA_real_
    ),
    IC50_sd_numeric = ifelse(
      N_numeric > 1,
      stats::sd(IC50_value, na.rm = TRUE),
      NA_real_
    ),
    Not_determined = any(Is_censored, na.rm = TRUE) && N_numeric == 0,
    Censor_threshold = ifelse(
      any(Is_censored, na.rm = TRUE),
      max(Censor_threshold[Is_censored], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    IC50_mean = dplyr::if_else(
      Not_determined,
      Censor_threshold,
      IC50_mean_numeric
    ),
    IC50_sd = IC50_sd_numeric,
    IC50_plot = dplyr::if_else(Not_determined, NA_real_, IC50_mean),
    IC50_label = dplyr::if_else(
      Not_determined,
      paste0("ND\n(>", format(Censor_threshold, trim = TRUE, scientific = FALSE), ")"),
      sprintf("%.1f", IC50_mean)
    )
  ) %>%
  dplyr::filter(!is.na(Drug), !is.na(IC50_mean))

single_screening_long <- read_screening_triplicates("single") %>%
  dplyr::transmute(
    DrugRaw = ItemRaw,
    Drug = factor(clean_drug(ItemRaw), levels = drug_levels_screen),
    Condition = factor(
      "Single prodrug NP",
      levels = c("Single prodrug NP", "DTX-SI-C18 combination")
    ),
    Replicate = Replicate,
    Cell_viability = Value
  ) %>%
  dplyr::filter(!is.na(Drug), !is.na(Cell_viability))

combination_screening_long <- read_screening_triplicates("combination") %>%
  dplyr::transmute(
    DrugRaw = ItemRaw,
    Drug = factor(clean_drug(ItemRaw), levels = drug_levels_screen),
    Condition = factor(
      "DTX-SI-C18 combination",
      levels = c("Single prodrug NP", "DTX-SI-C18 combination")
    ),
    Replicate = Replicate,
    Cell_viability = Value
  ) %>%
  dplyr::filter(!is.na(Drug), !is.na(Cell_viability))

combination_long <- dplyr::bind_rows(
  single_screening_long,
  combination_screening_long
)

combination_summary <- mean_summary(
  combination_long,
  c("Drug", "Condition"),
  "Cell_viability"
) %>%
  dplyr::mutate(
    Error = get_error(., bar_error_type),
    ymin = pmax(Mean - Error, 0),
    ymax = Mean + Error,
    Heatmap_label = sprintf("%.1f ± %.1f", Mean, SD),
    Heatmap_text_color = col_axis,
    Condition_short = dplyr::recode(
      as.character(Condition),
      "Single prodrug NP" = "Single NP",
      "DTX-SI-C18 combination" = "+ DTX-SI-C18"
    ),
    Condition_short = factor(
      Condition_short,
      levels = c("Single NP", "+ DTX-SI-C18")
    )
  )

# Activity-enhancement ratio:
# mean viability of the single-prodrug group / each replicate of the corresponding combination group.
single_reference_module2 <- single_screening_long %>%
  dplyr::group_by(Drug) %>%
  dplyr::summarise(
    Single_mean_reference = mean(Cell_viability, na.rm = TRUE),
    .groups = "drop"
  )

ratio_long <- combination_screening_long %>%
  dplyr::select(
    DrugRaw,
    Drug,
    Replicate,
    Combination_viability = Cell_viability
  ) %>%
  dplyr::left_join(single_reference_module2, by = "Drug") %>%
  dplyr::mutate(
    Single_to_combination_ratio =
      Single_mean_reference / Combination_viability,
    Highlight = ifelse(as.character(Drug) == "Cri", "Selected", "Other")
  ) %>%
  dplyr::filter(!is.na(Single_to_combination_ratio))

ratio_summary <- mean_summary(
  ratio_long,
  c("Drug"),
  "Single_to_combination_ratio"
) %>%
  dplyr::mutate(
    Error = get_error(., bar_error_type),
    ymin = pmax(Mean - Error, 0),
    ymax = Mean + Error,
    Highlight = ifelse(as.character(Drug) == "Cri", "Selected", "Other"),
    Mean_label = sprintf("%.2f", Mean)
  )

# DTX-SI-C18/Cri-SI-C18 combination-index data from the CI sheet.
ci_long <- read_ci_sheet()

ci_summary <- mean_summary(
  ci_long,
  c("DTX_Cri_Ratio", "Effect_level"),
  "CI"
) %>%
  dplyr::mutate(
    Error = get_error(., bar_error_type),
    ymin = pmax(Mean - Error, 0),
    ymax = Mean + Error,
    Mean_label = sprintf("%.2f", Mean),
    Heatmap_label = sprintf("%.2f ± %.2f", Mean, SD)
  )

# ED50 is used for the compact DTX/Cri ratio-optimization panel.
dtxcri_long <- ci_long %>%
  dplyr::filter(Effect_level == "ED50") %>%
  dplyr::transmute(
    DTX_Cri_Ratio,
    Replicate,
    CI
  )

dtxcri_summary <- mean_summary(
  dtxcri_long,
  c("DTX_Cri_Ratio"),
  "CI"
) %>%
  dplyr::mutate(
    Error = get_error(., bar_error_type),
    ymin = pmax(Mean - Error, 0),
    ymax = Mean + Error
  )

selected_ratio <- dtxcri_summary %>%
  dplyr::filter(Mean == min(Mean, na.rm = TRUE)) %>%
  dplyr::slice(1) %>%
  dplyr::pull(DTX_Cri_Ratio) %>%
  as.character()

dtxcri_summary <- dtxcri_summary %>%
  dplyr::mutate(
    Highlight = ifelse(
      as.character(DTX_Cri_Ratio) == selected_ratio,
      "Selected",
      "Other"
    ),
    Mean_label = sprintf("%.2f", Mean)
  )

dtxcri_long <- dtxcri_long %>%
  dplyr::mutate(
    Highlight = ifelse(
      as.character(DTX_Cri_Ratio) == selected_ratio,
      "Selected",
      "Other"
    )
  )

tumor_long <- parse_longitudinal_sheet("tumor volume", value_name = "Tumor_volume")
tumor_summary <- mean_summary(tumor_long, c("Group", "Day"), "Tumor_volume") %>%
  dplyr::mutate(Error = get_error(., curve_error_type), ymin = pmax(Mean - Error, 0), ymax = Mean + Error)

body_weight_long <- parse_longitudinal_sheet("body weight", value_name = "Body_weight")
body_weight_change <- body_weight_long %>%
  dplyr::group_by(Mouse) %>%
  dplyr::arrange(Day, .by_group = TRUE) %>%
  dplyr::mutate(Baseline_weight = Body_weight[which.min(Day)], Body_weight_change = (Body_weight / Baseline_weight - 1) * 100) %>%
  dplyr::ungroup()

body_weight_summary <- mean_summary(body_weight_change, c("Group", "Day"), "Body_weight_change") %>%
  dplyr::mutate(Error = get_error(., curve_error_type), ymin = Mean - Error, ymax = Mean + Error)

endpoint_day <- max(tumor_long$Day, na.rm = TRUE)
tumor_delta <- tumor_long %>%
  dplyr::group_by(Group, Mouse) %>%
  dplyr::arrange(Day, .by_group = TRUE) %>%
  dplyr::summarise(
    Baseline_volume = Tumor_volume[which.min(Day)],
    Endpoint_volume = Tumor_volume[which.max(Day)],
    Endpoint_day = Day[which.max(Day)],
    Delta_volume = Endpoint_volume - Baseline_volume,
    .groups = "drop"
  ) %>%
  dplyr::filter(Endpoint_day == endpoint_day)

control_delta_mean <- tumor_delta %>%
  dplyr::filter(Group == "PBS") %>%
  dplyr::summarise(Control_delta_mean = mean(Delta_volume, na.rm = TRUE)) %>%
  dplyr::pull(Control_delta_mean)

tgi_long <- tumor_delta %>%
  dplyr::mutate(TGI_raw = 100 * (1 - Delta_volume / control_delta_mean), TGI = dplyr::if_else(Group == "PBS", 0, TGI_raw))

tgi_summary <- mean_summary(tgi_long, c("Group"), "TGI") %>%
  dplyr::mutate(Error = get_error(., curve_error_type), xmin = Mean - Error, xmax = Mean + Error)

# -----------------------------
# 4. Export processed data
# -----------------------------
readr::write_csv(ic50_data, file.path(table_dir, "ic50_processed.csv"))
readr::write_csv(ic50_points, file.path(table_dir, "ic50_raw_points_if_available.csv"))
readr::write_csv(combination_long, file.path(table_dir, "combination_viability_long.csv"))
readr::write_csv(combination_summary, file.path(table_dir, "combination_viability_summary.csv"))
readr::write_csv(ratio_long, file.path(table_dir, "single_to_combination_ratio_long.csv"))
readr::write_csv(ratio_summary, file.path(table_dir, "single_to_combination_ratio_summary.csv"))
readr::write_csv(ci_long, file.path(table_dir, "combination_index_long.csv"))
readr::write_csv(ci_summary, file.path(table_dir, "combination_index_summary.csv"))
readr::write_csv(dtxcri_long, file.path(table_dir, "dtx_cri_ratio_ci_long.csv"))
readr::write_csv(dtxcri_summary, file.path(table_dir, "dtx_cri_ratio_ci_summary.csv"))
readr::write_csv(tumor_long, file.path(table_dir, "tumor_volume_long.csv"))
readr::write_csv(tumor_summary, file.path(table_dir, "tumor_volume_summary.csv"))
readr::write_csv(body_weight_change, file.path(table_dir, "body_weight_change_long.csv"))
readr::write_csv(body_weight_summary, file.path(table_dir, "body_weight_change_summary.csv"))
readr::write_csv(tumor_delta, file.path(table_dir, "tumor_delta_endpoint.csv"))
readr::write_csv(tgi_long, file.path(table_dir, "tumor_inhibition_long.csv"))
readr::write_csv(tgi_summary, file.path(table_dir, "tumor_inhibition_summary.csv"))

readr::write_lines(
  c(
    "Figure 6 combination-therapy output notes",
    "Raw replicate points are overlaid on panels d, e, f and i. Panel b overlays raw IC50 points only if replicate-level IC50 columns are present; otherwise only the summarized mean point can be shown.",
    "Panel b: censored IC50 values such as '>1000' are treated as not determined within the tested range and are not plotted as numeric bars.",
    "Panel c: 'Single prodrug NP' indicates the single candidate prodrug nanomedicine; '+ DTX-SI-C18' indicates DTX-SI-C18 NPs fixed at IC50 and combined with the indicated partner prodrug nanomedicine below IC50.",
    "Panel d: ratio is single-prodrug NP cell viability divided by DTX-SI-C18 combination cell viability.",
    "Panel e shows CI across ED50, ED75, ED90 and ED95; panel f summarizes ED50 for DTX-SI-C18/Cri-SI-C18 ratio optimization.",
    paste0("Panel f: selected DTX-SI-C18/Cri-SI-C18 ratio based on the lowest mean ED50 CI = ", selected_ratio, "."),
    paste0("Panel i: endpoint TGI (%) = [1 - DeltaV_treated / mean(DeltaV_PBS)] × 100, using DeltaV = V_endpoint - V_baseline at day ", endpoint_day, "."),
    paste0("Bar/dot panels: mean ± ", bar_error_type, "."),
    paste0("Tumor/body-weight/TGI panels: mean ± ", curve_error_type, ".")
  ),
  file.path(table_dir, "README_data_outputs.txt")
)

# -----------------------------
# 5. Plot panels
# -----------------------------
# a. Schematic placeholder.
p_a <- ggplot2::ggplot() +
  ggplot2::annotate("rect", xmin = 0.03, xmax = 0.97, ymin = 0.14, ymax = 0.86, fill = col_panel_bg, color = col_line, linewidth = 0.35, linetype = "dashed") +
  ggplot2::annotate("text", x = 0.50, y = 0.56, label = "Schematic placeholder", size = 3.0, fontface = "bold", color = col_axis) +
  ggplot2::annotate("text", x = 0.50, y = 0.43, label = "compatibility-gated DTX-SI-C18 combination screening", size = 2.3, color = col_grey_3) +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  ggplot2::theme_void(base_family = plot_base_family) +
  ggplot2::theme(plot.margin = ggplot2::margin(4, 5, 4, 5))

# b. IC50.
ic50_measured <- ic50_data %>%
  dplyr::filter(!Not_determined) %>%
  dplyr::mutate(
    ymin = pmax(IC50_plot - dplyr::coalesce(IC50_sd, 0), 0),
    ymax = IC50_plot + dplyr::coalesce(IC50_sd, 0),
    Label_y = ymax + max(ymax, na.rm = TRUE) * 0.045,
    IC50_bar_label = dplyr::if_else(abs(IC50_mean - round(IC50_mean)) < 0.05, sprintf("%.0f", IC50_mean), sprintf("%.1f", IC50_mean))
  )

ic50_not_determined <- ic50_data %>%
  dplyr::filter(Not_determined) %>%
  dplyr::mutate(label_y = max(ic50_measured$ymax, na.rm = TRUE) * 0.08)

p_b <- ggplot2::ggplot(ic50_measured, ggplot2::aes(x = Drug, y = IC50_plot)) +
  ggplot2::geom_col(width = 0.62, fill = "#D9EEF7", color = "black", linewidth = 0.38) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = 0.16, linewidth = 0.32, color = "black") +
  {
    if (nrow(ic50_points) > 0) {
      ggplot2::geom_point(data = ic50_points, ggplot2::aes(x = Drug, y = IC50_value), inherit.aes = FALSE, position = make_jitter(0.12, 0, 11), shape = 21, size = 1.35, stroke = 0.25, fill = "white", color = col_axis, alpha = 0.90)
    } else {
      ggplot2::geom_point(ggplot2::aes(y = IC50_plot), shape = 21, size = 1.45, stroke = 0.25, fill = "white", color = col_axis)
    }
  } +
  ggplot2::geom_text(ggplot2::aes(y = Label_y, label = IC50_bar_label), size = 2.05, color = col_axis, vjust = 0) +
  ggplot2::geom_text(data = ic50_not_determined, ggplot2::aes(x = Drug, y = label_y, label = IC50_label), inherit.aes = FALSE, size = 2.0, lineheight = 0.82, color = col_grey_3) +
  ggplot2::scale_x_discrete(drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
  ggplot2::labs(x = NULL, y = expression(IC[50]~"("*mu*"M)"), title = "Single-prodrug potency") +
  theme_nature() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1))

# c. Cell viability heatmap.
p_c <- ggplot2::ggplot(combination_summary, ggplot2::aes(x = Condition_short, y = Drug, fill = Mean)) +
  ggplot2::geom_tile(color = "#F7F5EF", linewidth = 0.50, width = 0.92, height = 0.86) +
  ggplot2::geom_text(ggplot2::aes(label = Heatmap_label, color = Heatmap_text_color), size = 2.05) +
  ggplot2::scale_color_identity() +
  ggplot2::geom_tile(data = dplyr::filter(combination_summary, Drug == "Cri", Condition == "DTX-SI-C18 combination"), fill = NA, color = "#B84E4A", linewidth = 0.58, width = 0.92, height = 0.86) +
  ggplot2::scale_y_discrete(limits = rev(drug_levels_screen), drop = FALSE) +
  ggplot2::scale_fill_gradientn(
    colours = c(col_heat_low, "#B9D9E8", col_heat_mid, col_heat_high_mid, col_heat_high),
    values = scales::rescale(c(0, 20, 40, 60, 100)),
    limits = c(0, 100),
    breaks = c(0, 20, 40, 60, 80, 100),
    name = "Cell viability\n(%)"
  ) +
  ggplot2::labs(x = NULL, y = "Partner prodrug", title = "C4-2 cell viability (%)") +
  theme_heatmap() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5), legend.position = "right")

# d. Single/combination viability ratio.
p_d <- ggplot2::ggplot(ratio_summary, ggplot2::aes(x = Drug, y = Mean)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.28, color = "#9C9A92") +
  ggplot2::geom_point(data = ratio_long, ggplot2::aes(x = Drug, y = Single_to_combination_ratio), inherit.aes = FALSE, position = make_jitter(0.13, 0, 21), shape = 21, size = 1.35, stroke = 0.25, fill = "white", color = col_axis, alpha = 0.9) +
  ggplot2::geom_linerange(ggplot2::aes(ymin = ymin, ymax = ymax, color = Highlight), linewidth = 0.42) +
  ggplot2::geom_point(ggplot2::aes(fill = Highlight), shape = 21, size = 2.10, stroke = 0.30, color = col_axis) +
  ggplot2::geom_text(ggplot2::aes(y = ymax + 0.14, label = Mean_label, color = Highlight), size = 1.85, vjust = 0, show.legend = FALSE) +
  ggplot2::scale_fill_manual(values = c("Other"= col_beige, "Selected"= col_selected), guide = "none") +
  ggplot2::scale_color_manual(values = c("Other"= col_grey_3, "Selected"= col_selected), guide = "none") +
  ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.20))) +
  ggplot2::labs(x = "Partner prodrug", y = "Single/combination\nviability ratio", title = "Activity-enhancement ratio") +
  theme_nature() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1), plot.margin = ggplot2::margin(4, 4, 4, 4))

# e. CI across effect levels.
p_e <- ggplot2::ggplot(
  ci_summary,
  ggplot2::aes(x = Effect_level, y = DTX_Cri_Ratio, fill = Mean)
) +
  ggplot2::geom_tile(
    color = "#F7F5EF",
    linewidth = 0.50,
    width = 0.92,
    height = 0.86
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = Heatmap_label),
    size = 1.75,
    color = col_axis
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#8FB6DD",
    mid = "#F7F7F7",
    high = "#EBA6A1",
    midpoint = 1,
    limits = c(0, max(1.45, max(ci_summary$Mean, na.rm = TRUE))),
    name = "CI"
  ) +
  ggplot2::labs(
    x = "Effect level",
    y = "DTX-SI-C18:Cri-SI-C18 ratio",
    title = "Combination index across effect levels"
  ) +
  theme_heatmap() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
    axis.text.y = ggplot2::element_text(size = 6.0),
    legend.position = "right"
  )

# f. DTX/Cri ratio optimization.
p_f <- ggplot2::ggplot(dtxcri_summary, ggplot2::aes(x = DTX_Cri_Ratio, y = Mean)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.30, color = "#7E7C74") +
  ggplot2::geom_col(ggplot2::aes(fill = Highlight), width = 0.62, color = NA) +
  ggplot2::geom_point(data = dtxcri_long, ggplot2::aes(x = DTX_Cri_Ratio, y = CI), inherit.aes = FALSE, position = make_jitter(0.12, 0, 23), shape = 21, size = 1.35, stroke = 0.25, fill = "white", color = col_axis, alpha = 0.95) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = 0.16, linewidth = 0.30, color = col_axis) +
  ggplot2::geom_text(ggplot2::aes(y = ymax + 0.035, label = Mean_label), size = 1.85, vjust = 0, color = col_axis) +
  ggplot2::scale_fill_manual(values = c("Other"= col_beige, "Selected"= col_selected), guide = "none") +
  ggplot2::scale_y_continuous(limits = c(0, 1.30), breaks = seq(0, 1.0, 0.25), expand = ggplot2::expansion(mult = c(0, 0.02))) +
  ggplot2::labs(x = "DTX-SI-C18/Cri-SI-C18\nmolar ratio", y = "Combination index (CI)", title = "DTX/Cri ratio optimization") +
  theme_nature() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1))

# g. Tumor growth curve.
p_g <- ggplot2::ggplot(tumor_summary, ggplot2::aes(x = Day, y = Mean, color = Group, group = Group)) +
  ggplot2::geom_line(linewidth = 0.66) +
  ggplot2::geom_point(size = 1.35) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = 0.42, linewidth = 0.26) +
  ggplot2::scale_color_manual(values = treatment_cols, drop = FALSE, name = NULL) +
  ggplot2::scale_x_continuous(breaks = sort(unique(tumor_summary$Day)), expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
  ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::labs(x = "Time after treatment\ninitiation (days)", y = expression("Tumor volume (mm"^3*")"), title = "Tumor growth") +
  theme_nature() +
  ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 6.2), legend.key.width = grid::unit(4.2, "mm"), legend.margin = ggplot2::margin(t = -2, b = -4), legend.box.margin = ggplot2::margin(t = -4, b = -4))

# h. Body-weight change.
bw_lim <- max(abs(body_weight_summary$ymin), abs(body_weight_summary$ymax), 10, na.rm = TRUE)
bw_lim <- ceiling(bw_lim / 5) * 5
p_h <- ggplot2::ggplot(body_weight_summary, ggplot2::aes(x = Day, y = Mean, color = Group, group = Group)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "#A5A39B") +
  ggplot2::geom_line(linewidth = 0.64) +
  ggplot2::geom_point(size = 1.28) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = 0.42, linewidth = 0.26) +
  ggplot2::scale_color_manual(values = treatment_cols, drop = FALSE, name = NULL) +
  ggplot2::scale_x_continuous(breaks = sort(unique(body_weight_summary$Day)), expand = ggplot2::expansion(mult = c(0.02, 0.04))) +
  ggplot2::scale_y_continuous(limits = c(-bw_lim, bw_lim), breaks = scales::pretty_breaks(n = 5)) +
  ggplot2::labs(x = "Time after treatment\ninitiation (days)", y = "Body-weight change (%)", title = "Body-weight change") +
  theme_nature() +
  ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 6.4), legend.key.width = grid::unit(4.2, "mm"), legend.margin = ggplot2::margin(t = -2, b = -4), legend.box.margin = ggplot2::margin(t = -4, b = -4))

# i. Tumor-growth inhibition rate with raw mouse-level points.
x_min_tgi <- min(0, floor(min(tgi_summary$xmin, na.rm = TRUE) / 10) * 10)
x_max_tgi <- ceiling(max(tgi_summary$xmax, na.rm = TRUE) / 10) * 10
p_i <- ggplot2::ggplot(tgi_summary, ggplot2::aes(y = Group, x = Mean, fill = Group)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, color = "#A5A39B") +
  ggplot2::geom_col(width = 0.62, color = NA, alpha = 0.95) +
  ggplot2::geom_point(data = tgi_long, ggplot2::aes(x = TGI, y = Group), inherit.aes = FALSE, position = ggplot2::position_jitter(width = 0, height = 0.13, seed = 24), shape = 21, size = 1.55, stroke = 0.30, fill = "white", color = col_axis, alpha = 0.95) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = xmin, xmax = xmax), height = 0.18, linewidth = 0.30, color = col_axis) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", Mean), x = Mean + 3), size = 2.1, hjust = 0, color = col_axis) +
  ggplot2::scale_fill_manual(values = treatment_cols, guide = "none", drop = FALSE) +
  ggplot2::scale_y_discrete(limits = rev(treatment_levels), drop = FALSE) +
  ggplot2::coord_cartesian(xlim = c(x_min_tgi, x_max_tgi + 12), clip = "off") +
  ggplot2::labs(x = "Tumor-growth inhibition (%)", y = NULL, title = "Endpoint tumor inhibition") +
  theme_nature() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_line(linewidth = 0.18, color = "#E8E6DE"), panel.grid.major.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 6.4), legend.position = "none")

# -----------------------------
# 6. Save single panels
# -----------------------------
save_plot(p_a, "workflow", 7.4, 1.35)
save_plot(p_b, "single_prodrug_IC50", 3.1, 2.45)
save_plot(p_c, "cell_viability_heatmap", 5.5, 2.45)
save_plot(p_d, "single_to_combination_ratio", 3.0, 2.45)
save_plot(p_e, "combination_index_screen", 2.85, 2.35)
save_plot(p_f, "DTX_Cri_ratio_optimization", 3.0, 2.45)
save_plot(p_g, "tumor_growth", 4.9, 2.75)
save_plot(p_h, "body_weight", 4.9, 2.75)
save_plot(p_i, "tumor_growth_inhibition", 3.9, 2.75)

# -----------------------------
# 7. Combined figure
# -----------------------------
combined_design <- "
AAAAAA
BBBCCC
DDEEFF
GGGIII
HHHIII
"

p_c_comb <- p_c + ggplot2::theme(legend.position = "right")
p_g_comb <- p_g
p_h_comb <- p_h + ggplot2::theme(legend.position = "none")

combined_figure <- p_a + p_b + p_c_comb + p_d + p_e + p_f + p_g_comb + p_h_comb + p_i +
  patchwork::plot_layout(design = combined_design, heights = c(0.60, 1.05, 1.02, 1.10, 1.10), widths = rep(1, 6), guides = "collect") +
  patchwork::plot_annotation(tag_levels = "a") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(face = "bold", size = 8.5, family = plot_base_family),
    plot.tag.position = c(0.005, 0.985),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 6.2),
    legend.title = ggplot2::element_blank(),
    legend.key.width = grid::unit(4.2, "mm"),
    legend.margin = ggplot2::margin(t = -2, b = -3),
    legend.box.margin = ggplot2::margin(t = -3, b = -3)
  )

save_plot(combined_figure, "combination_therapy", 8.9, 9.4)

# Canonical main-figure copy at results/figures/main/Fig6/.
ggplot2::ggsave(
  file.path(main_figure_dir, "combination_therapy.pdf"),
  combined_figure, width = 8.9, height = 9.4, units = "in",
  device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
)
ggplot2::ggsave(
  file.path(main_figure_dir, "combination_therapy.png"),
  combined_figure, width = 8.9, height = 9.4, units = "in",
  dpi = 600, bg = "white", limitsize = FALSE
)

save(
  ic50_data, ic50_points,
  combination_long, combination_summary,
  ratio_long, ratio_summary,
  ci_long, ci_summary,
  dtxcri_long, dtxcri_summary,
  tumor_long, tumor_summary,
  body_weight_change, body_weight_summary,
  tumor_delta, tgi_long, tgi_summary,
  selected_ratio, endpoint_day, control_delta_mean,
  file = file.path(table_dir, "combination_therapy_processed_data.RData")
)

message("Screening source: ", screening_input_file, " [IC50, single, combination, CI]")
message("Therapy source: ", input_file, " [tumor volume, body weight]")
message("Figure 6 therapy figures: ", out_dir)
message("Figure 6 therapy tables: ", table_dir)
message("Canonical main figure: ", file.path(main_figure_dir, "combination_therapy.pdf"))
