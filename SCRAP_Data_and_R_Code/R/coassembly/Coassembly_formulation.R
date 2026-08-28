# Co-assembly formulation analysis
# Co-delivery size data are read from data/coassembly/Co_delivery_sizes.xlsx.
rm(list = ls())
graphics.off()
set.seed(123)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(scales)
  library(grid)
  library(readxl)
})

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required. Please install it with install.packages('patchwork').")
}
library(patchwork)

# -----------------------------
# 0. Project paths
# -----------------------------
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grepl(file_arg, cmd_args)])
  if (length(script_path) > 0 && nzchar(script_path[1])) {
    return(dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)))
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

source_workbook <- file.path(project_root, "data", "coassembly", "Co_delivery_sizes.xlsx")
figure_dir <- file.path(project_root, "results", "figures", "coassembly")
table_dir <- file.path(project_root, "results", "tables", "coassembly")

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Raw input data
# -----------------------------
if (!file.exists(source_workbook)) {
  stop("Input workbook not found: ", source_workbook, call. = FALSE)
}

read_replicate_sheet <- function(path, sheet_name) {
  x <- readxl::read_excel(
    path,
    sheet = sheet_name,
    col_names = FALSE,
    .name_repair = "minimal"
  )
  
  x <- x[rowSums(!is.na(x)) > 0, , drop = FALSE]
  
  if (nrow(x) < 2 || ncol(x) < 2) {
    stop("Sheet '", sheet_name, "' does not contain a valid replicate table.", call. = FALSE)
  }
  
  x
}

# Binary co-assembly with increasing DTX-SI-C18 fraction.
# withDTX layout:
#   column 1    = candidate
#   columns 2-4 = 0% DTX triplicates
#   columns 5-7 = 25% DTX triplicates
#   columns 8-10 = 50% DTX triplicates
#   columns 11-13 = 75% DTX triplicates
withdtx_sheet <- read_replicate_sheet(source_workbook, "withDTX")
withdtx_body <- as.data.frame(withdtx_sheet[-1, , drop = FALSE], stringsAsFactors = FALSE)

fraction_columns <- list(
  "0% DTX" = 2:4,
  "25% DTX" = 5:7,
  "50% DTX" = 8:10,
  "75% DTX" = 11:13
)

binary_replicates <- dplyr::bind_rows(
  lapply(seq_len(nrow(withdtx_body)), function(i) {
    candidate <- trimws(as.character(withdtx_body[i, 1][[1]]))
    if (is.na(candidate) || candidate == "") return(tibble())
    
    dplyr::bind_rows(
      lapply(names(fraction_columns), function(frac) {
        vals <- suppressWarnings(
          as.numeric(
            unlist(
              withdtx_body[i, fraction_columns[[frac]], drop = FALSE],
              use.names = FALSE
            )
          )
        )
        
        tibble(
          Candidate = candidate,
          DTX_fraction = frac,
          Replicate = seq_along(vals),
          Size = vals
        )
      })
    )
  })
) %>%
  filter(!is.na(Size))

binary_n_check <- binary_replicates %>%
  count(Candidate, DTX_fraction, name = "n")

if (
  nrow(binary_n_check) == 0 ||
  any(binary_n_check$n != 3) ||
  n_distinct(binary_replicates$DTX_fraction) != 4
) {
  stop(
    "Sheet 'withDTX' must contain exactly three replicates for each of four DTX fractions.",
    call. = FALSE
  )
}

binary_summary_raw <- binary_replicates %>%
  group_by(Candidate, DTX_fraction) %>%
  summarise(
    Mean = mean(Size),
    SD = stats::sd(Size),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  rename(Size = Mean)

if (any(is.na(binary_summary_raw$SD))) {
  stop(
    "SD calculation failed for at least one withDTX condition. Check that all triplicate cells are numeric.",
    call. = FALSE
  )
}

binary_wide <- binary_summary_raw %>%
  select(Candidate, DTX_fraction, Size) %>%
  pivot_wider(names_from = DTX_fraction, values_from = Size)

# Multi-prodrug co-formulations.
# multi layout:
#   column 1 = manuscript-facing formulation name
#   columns 2-4 = three DLS replicates
multi_sheet <- read_replicate_sheet(source_workbook, "multi")
multi_body <- as.data.frame(multi_sheet[-1, , drop = FALSE], stringsAsFactors = FALSE)

formulation_order <- trimws(as.character(multi_body[[1]]))
formulation_order <- formulation_order[!is.na(formulation_order) & formulation_order != ""]

multi_replicates <- dplyr::bind_rows(
  lapply(seq_len(nrow(multi_body)), function(i) {
    formulation <- trimws(as.character(multi_body[i, 1][[1]]))
    if (is.na(formulation) || formulation == "") return(tibble())
    
    vals <- suppressWarnings(
      as.numeric(
        unlist(multi_body[i, 2:4, drop = FALSE], use.names = FALSE)
      )
    )
    
    tibble(
      Formulation = formulation,
      Replicate = seq_along(vals),
      Size = vals
    )
  })
) %>%
  filter(!is.na(Size))

multi_n_check <- multi_replicates %>%
  count(Formulation, name = "n")

if (nrow(multi_n_check) == 0 || any(multi_n_check$n != 3)) {
  stop(
    "Sheet 'multi' must contain exactly three size replicates for every formulation.",
    call. = FALSE
  )
}

multi_df <- multi_replicates %>%
  group_by(Formulation) %>%
  summarise(
    Mean = mean(Size),
    SD = stats::sd(Size),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  rename(Size = Mean) %>%
  mutate(.row_order = match(Formulation, formulation_order)) %>%
  arrange(.row_order) %>%
  select(-.row_order)

if (any(is.na(multi_df$SD))) {
  stop(
    "SD calculation failed for at least one multi-prodrug formulation. Check that all triplicate cells are numeric.",
    call. = FALSE
  )
}

# Single-prodrug DLS size and current 13-point scores.
# Candidate single-prodrug sizes are taken from the 0% DTX means in the
# withDTX sheet. DTX-SI-C18 is not present in that sheet, so its existing
# 18.8-nm anchor is retained.
single_scores <- tibble::tribble(
  ~Candidate,       ~Score13,
  "DTX-SI-C18",        13,
  "Dex-SI-C18",        13,
  "Cri-SI-C18",        13,
  "PPT-SI-C18",        13,
  "Dox-SI-C18",        13,
  "Ful-SI-C18",         7,
  "IDB-SI-C18",        11,
  "Abi-SI-C18",        10
)

candidate_single_size <- binary_summary_raw %>%
  filter(DTX_fraction == "0% DTX") %>%
  select(Candidate, Single_size = Size)

single_si <- bind_rows(
  tibble(Candidate = "DTX-SI-C18", Single_size = 18.8),
  candidate_single_size
) %>%
  left_join(single_scores, by = "Candidate")

# Optional MD placeholder table for future DTX/Abi or DTX/IDB multi-molecule simulations.
# Replace NA values with real MD results after analysis. The plotting code automatically
# exports an empty template if values are not available.
md_template <- tibble::tribble(
  ~System,                    ~DTX_fraction, ~System_Rg_nm, ~Drug_cluster_Rg_nm, ~Drug_polymer_contact_nm2, ~Drug_drug_contact_index,
  "Abi-SI-C18",                         0,             NA,                  NA,                       NA,                     NA,
  "DTX/Abi 25/75",                     25,             NA,                  NA,                       NA,                     NA,
  "DTX/Abi 50/50",                     50,             NA,                  NA,                       NA,                     NA,
  "DTX/Abi 75/25",                     75,             NA,                  NA,                       NA,                     NA,
  "DTX-SI-C18",                       100,             NA,                  NA,                       NA,                     NA
)

# -----------------------------
# 2. Labels and derived variables
# -----------------------------
full_payload <- c(
  "DTX-SI-C18"= "Docetaxel",
  "Abi-SI-C18"= "Abiraterone",
  "IDB-SI-C18"= "Idebenone",
  "Ful-SI-C18"= "Fulvestrant",
  "Dox-SI-C18"= "Doxorubicin",
  "Cri-SI-C18"= "Crizotinib",
  "PPT-SI-C18"= "Podophyllotoxin",
  "Dex-SI-C18"= "Dexamethasone"
)

candidate_label_map <- c(
  "Abi-SI-C18"= "Abiraterone\n(Abi-SI-C18)",
  "IDB-SI-C18"= "Idebenone\n(IDB-SI-C18)",
  "Ful-SI-C18"= "Fulvestrant\n(Ful-SI-C18)",
  "Dox-SI-C18"= "Doxorubicin\n(Dox-SI-C18)",
  "Cri-SI-C18"= "Crizotinib\n(Cri-SI-C18)",
  "PPT-SI-C18"= "Podophyllotoxin\n(PPT-SI-C18)",
  "Dex-SI-C18"= "Dexamethasone\n(Dex-SI-C18)"
)

candidate_group_map <- c(
  "Abi-SI-C18"= "Initially oversized",
  "IDB-SI-C18"= "Initially oversized",
  "Ful-SI-C18"= "Initially oversized",
  "Dox-SI-C18"= "Initially sub-30 nm",
  "Cri-SI-C18"= "Initially sub-30 nm",
  "PPT-SI-C18"= "Initially sub-30 nm",
  "Dex-SI-C18"= "Initially sub-30 nm"
)

candidate_order <- c(
  "Abi-SI-C18", "IDB-SI-C18", "Ful-SI-C18",
  "Dox-SI-C18", "Cri-SI-C18", "PPT-SI-C18", "Dex-SI-C18"
)

binary_long <- binary_summary_raw %>%
  mutate(
    DTX_fraction = factor(DTX_fraction, levels = c("0% DTX", "25% DTX", "50% DTX", "75% DTX")),
    DTX_numeric = recode(as.character(DTX_fraction),
                         "0% DTX"= 0, "25% DTX"= 25,
                         "50% DTX"= 50, "75% DTX"= 75) %>% as.numeric(),
    Candidate = factor(Candidate, levels = candidate_order),
    Candidate_label = factor(candidate_label_map[as.character(Candidate)],
                             levels = candidate_label_map[candidate_order]),
    Initial_group = factor(candidate_group_map[as.character(Candidate)],
                           levels = c("Initially oversized", "Initially sub-30 nm")),
    Outcome = ifelse(Size < 30, "Sub-30 nm", "Oversized"),
    Outcome = factor(Outcome, levels = c("Sub-30 nm", "Oversized")),
    Size_label = sprintf("%.1f ± %.1f", Size, SD)
  )

summary_df <- binary_long %>%
  group_by(Candidate, Candidate_label, Initial_group) %>%
  summarise(
    baseline_size = Size[DTX_numeric == 0],
    min_size = min(Size, na.rm = TRUE),
    min_size_fraction = DTX_numeric[which.min(Size)][1],
    size_75 = Size[DTX_numeric == 75],
    rescue_fraction = ifelse(any(Size < 30), min(DTX_numeric[Size < 30]), NA_real_),
    abs_reduction_min = baseline_size - min_size,
    percent_reduction_min = 100 * (baseline_size - min_size) / baseline_size,
    .groups = "drop"
  ) %>%
  mutate(
    Rescue_label = ifelse(is.na(rescue_fraction), "Not reached", paste0(rescue_fraction, "%")),
    Candidate_label = factor(Candidate_label, levels = candidate_label_map[candidate_order])
  )

oversized_summary <- summary_df %>%
  filter(Initial_group == "Initially oversized") %>%
  arrange(desc(baseline_size)) %>%
  mutate(Candidate_label = factor(Candidate_label, levels = Candidate_label))

single_si_plot <- single_si %>%
  mutate(
    Candidate = factor(Candidate, levels = c("DTX-SI-C18", candidate_order)),
    Payload = full_payload[as.character(Candidate)],
    Label = paste0(Payload, "\n(", Candidate, ")"),
    Outcome = ifelse(Single_size < 30, "Sub-30 nm", "Oversized"),
    Outcome = factor(Outcome, levels = c("Sub-30 nm", "Oversized")),
    Score_group = case_when(
      Score13 == 13 ~ "13/13",
      Score13 >= 9 ~ "9-12/13",
      TRUE ~ "≤8/13"
    ),
    Score_group = factor(Score_group, levels = c("13/13", "9-12/13", "≤8/13"))
  )

multi_order <- formulation_order[formulation_order %in% multi_df$Formulation]

multi_df <- multi_df %>%
  mutate(
    No_prodrugs = stringr::str_count(Formulation, stringr::fixed(" and ")) + 1L,
    Outcome = ifelse(Size < 30, "Sub-30 nm", "Oversized"),
    Formulation_short = str_wrap(Formulation, width = 32),
    Formulation_short = factor(Formulation_short, levels = rev(str_wrap(multi_order, width = 32)))
  )

multi_replicates_plot <- multi_replicates %>%
  left_join(
    multi_df %>% select(Formulation, No_prodrugs, Formulation_short),
    by = "Formulation"
  )


message("withDTX triplicate summary (mean ± SD):")
print(
  binary_summary_raw %>%
    mutate(mean_SD = sprintf("%.1f ± %.1f", Size, SD)) %>%
    select(Candidate, DTX_fraction, mean_SD)
)

message("Multi-prodrug triplicate summary (mean ± SD):")
print(
  multi_df %>%
    mutate(mean_SD = sprintf("%.1f ± %.1f", Size, SD)) %>%
    select(Formulation, mean_SD)
)

# Exploratory composition-weighted single-prodrug score.
# This is NOT the rigorous mixed 13-point score. It only indicates how the single-prodrug
# score distribution changes when DTX-SI-C18 is added. Keep the subtitle/label in the figure.
score_mix <- binary_long %>%
  left_join(single_si %>% select(Candidate, Candidate_score = Score13), by = c("Candidate"= "Candidate")) %>%
  mutate(
    DTX_score = 13,
    Exploratory_weighted_score = (DTX_numeric / 100) * DTX_score + (1 - DTX_numeric / 100) * Candidate_score,
    Score_label = sprintf("%.1f", Exploratory_weighted_score)
  )

# -----------------------------
# 3. Plot style
# -----------------------------
# Restrained, colorblind-friendly, Nature-like palette.
col_blue <- "#2F6F9F"
col_blue2 <- "#6EA9C8"
col_blue_light <- "#EAF4F8"
col_teal <- "#078C82"
col_teal_light <- "#CFE9E5"
col_orange <- "#D77A61"
col_orange_dark <- "#A94935"
col_cream <- "#F5E6D3"
col_grey <- "#D9D9D9"
col_grey2 <- "#8B8B8B"
col_dark <- "#222222"
col_purple <- "#8D7AB8"

pal_group <- c("Initially oversized"= col_orange, "Initially sub-30 nm"= col_blue)
pal_outcome <- c("Sub-30 nm"= col_teal, "Oversized"= col_orange)
pal_score <- c("13/13"= col_teal, "9-12/13"= col_blue2, "≤8/13"= col_orange)
pal_multidrug <- c("3"= "#9CC9C3", "4"= "#8FB3D9", "5"= "#8D7AB8")

base_theme <- theme_classic(base_size = 8) +
  theme(
    text = element_text(colour = col_dark, family = "sans"),
    axis.text = element_text(colour = col_dark, size = 7),
    axis.title = element_text(colour = col_dark, size = 8),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.title = element_text(face = "bold", hjust = 0, size = 10, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 7.2, colour = "grey30", margin = margin(b = 3)),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 7.2),
    legend.title = element_text(size = 7.3),
    legend.text = element_text(size = 6.8),
    legend.key.size = unit(0.35, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

# Small helper to save panels.
save_panel <- function(plot, name, width = 4.6, height = 3.3) {
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot, width = width, height = height, useDingbats = FALSE)
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 600)
  ggsave(file.path(figure_dir, paste0(name, ".tiff")), plot, width = width, height = height, dpi = 600, compression = "lzw")
}

# -----------------------------
# 4. Panel a: R-drawn workflow schematic
# -----------------------------
# A lightweight vector schematic generated directly in R. It is intended as an editable
# guide for AI/Illustrator, not as a replacement for a hand-polished graphical abstract panel.
p_a <- ggplot() +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 4), expand = FALSE) +
  theme_void(base_size = 8) +
  labs(title = "a") +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0, vjust = 1, margin = margin(b = 0))) +
  # vial
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = 1.15, ymax = 2.65, fill = "#F4FAFD", colour = "grey35", linewidth = 0.35) +
  annotate("rect", xmin = 0.50, xmax = 1.50, ymin = 2.65, ymax = 2.95, fill = "#D9D9D9", colour = "grey35", linewidth = 0.35) +
  annotate("segment", x = 0.55, xend = 1.45, y = 2.1, yend = 2.1, colour = "#9FCBE2", linewidth = 0.45) +
  annotate("point", x = c(0.78, 0.98, 1.18, 1.25, 0.88), y = c(1.55, 1.85, 1.55, 2.0, 2.25),
           size = 3.0, shape = 21, stroke = 0.25,
           fill = c("#078C82", "#6EA9C8", "#D77A61", "#8D7AB8", "#E08B3E"), colour = "grey30") +
  annotate("text", x = 1.0, y = 0.55, label = "Mixed SI\nprodrugs", size = 2.8, lineheight = 0.95) +
  # arrow to beaker
  annotate("segment", x = 1.75, xend = 2.75, y = 2.0, yend = 2.0,
           arrow = arrow(length = unit(0.16, "cm")), linewidth = 0.45, colour = col_dark) +
  # beaker and polymer
  annotate("rect", xmin = 3.05, xmax = 5.35, ymin = 1.05, ymax = 2.75, fill = "#F2FAFD", colour = "grey35", linewidth = 0.35) +
  annotate("segment", x = 3.05, xend = 5.35, y = 2.25, yend = 2.25, colour = "#9FCBE2", linewidth = 0.45) +
  annotate("text", x = 4.20, y = 2.85, label = expression(mPEG[2*k]*"-PCL"[2*k]), size = 3.1) +
  annotate("curve", x = 3.35, xend = 4.10, y = 1.55, yend = 1.65, curvature = 0.35, colour = col_blue, linewidth = 0.55) +
  annotate("curve", x = 3.55, xend = 4.45, y = 1.90, yend = 1.80, curvature = -0.30, colour = col_blue, linewidth = 0.55) +
  annotate("curve", x = 4.15, xend = 5.05, y = 1.35, yend = 1.50, curvature = 0.30, colour = col_blue, linewidth = 0.55) +
  annotate("point", x = c(3.35, 4.10, 3.55, 4.45, 4.15, 5.05),
           y = c(1.55, 1.65, 1.90, 1.80, 1.35, 1.50), size = 2.0, colour = col_blue2) +
  annotate("point", x = c(3.95, 4.75), y = c(1.45, 2.0), size = 3.0, shape = 22, fill = "#73B66B", colour = "grey30", stroke = 0.25) +
  # syringe/pipette
  annotate("segment", x = 3.40, xend = 3.95, y = 3.75, yend = 2.65, colour = "grey35", linewidth = 0.65) +
  annotate("rect", xmin = 2.95, xmax = 3.55, ymin = 3.10, ymax = 3.70, angle = 0, fill = "#F7FBFD", colour = "grey35", linewidth = 0.35) +
  annotate("point", x = c(3.10, 3.25, 3.42), y = c(3.25, 3.45, 3.22),
           size = 1.8, shape = 21, fill = c("#078C82", "#D77A61", "#8D7AB8"), colour = "grey30", stroke = 0.25) +
  annotate("point", x = 4.10, y = 2.35, size = 2.2, shape = 21, fill = "#8FB3D9", colour = "#5F8EB4", stroke = 0.25) +
  annotate("text", x = 4.2, y = 0.55, label = "Nanoprecipitation", size = 2.8) +
  # arrow to nanoparticle
  annotate("segment", x = 5.75, xend = 6.75, y = 2.0, yend = 2.0,
           arrow = arrow(length = unit(0.16, "cm")), linewidth = 0.45, colour = col_dark) +
  # nanoparticle
  annotate("point", x = 8.0, y = 2.05, size = 44, shape = 21, fill = alpha(col_teal_light, 0.80), colour = alpha(col_teal, 0.70), stroke = 0.6) +
  annotate("point", x = 8.0, y = 2.05, size = 26, shape = 21, fill = alpha(col_cream, 0.85), colour = alpha(col_orange, 0.30), stroke = 0.4) +
  # corona chains
  annotate("segment", x = rep(8.0, 12), y = rep(2.05, 12),
           xend = 8.0 + 1.15 * cos(seq(0, 2*pi, length.out = 12)),
           yend = 2.05 + 1.15 * sin(seq(0, 2*pi, length.out = 12)),
           colour = col_blue, linewidth = 0.45) +
  annotate("point", x = 8.0 + 1.15 * cos(seq(0, 2*pi, length.out = 12)),
           y = 2.05 + 1.15 * sin(seq(0, 2*pi, length.out = 12)),
           size = 2.4, shape = 21, fill = "#A7CDE2", colour = col_blue, stroke = 0.35) +
  annotate("point", x = c(7.65, 7.95, 8.25, 8.10, 7.85), y = c(2.20, 2.45, 2.05, 1.72, 1.80),
           size = 3.4, shape = 21, stroke = 0.25,
           fill = c("#078C82", "#8D7AB8", "#D77A61", "#E08B3E", "#6EA9C8"), colour = "grey30") +
  annotate("text", x = 8.0, y = 0.55, label = "Multi-prodrug\nnanomedicine", size = 2.8, lineheight = 0.95)

# -----------------------------
# 5. Panel b: Single-prodrug score-size anchor
# -----------------------------
p_b <- ggplot(single_si_plot, aes(x = Score13, y = Single_size, fill = Score_group)) +
  geom_hline(yintercept = 30, linetype = 2, linewidth = 0.45, colour = "grey40") +
  geom_point(shape = 21, size = 3.0, colour = "black", stroke = 0.3) +
  geom_text(aes(label = gsub("-SI-C18", "", as.character(Candidate))),
            nudge_y = 5.2, size = 2.2, check_overlap = TRUE) +
  scale_fill_manual(values = pal_score, name = "13-point\nscore") +
  scale_x_continuous(breaks = c(6, 9, 13), limits = c(5.2, 13.8)) +
  scale_y_continuous(limits = c(0, 120), breaks = c(0, 30, 60, 90, 120)) +
  labs(
    title = "b",
    subtitle = "Single SI prodrug score-size anchor",
    x = "Segmented compatibility score",
    y = "Single-prodrug diameter (nm)"
  ) +
  base_theme +
  theme(legend.position = "right")

# -----------------------------
# 6. Panel c: DTX fraction heatmap
# -----------------------------
p_c <- ggplot(binary_long, aes(x = DTX_fraction, y = Candidate_label, fill = Size)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = Size_label), size = 2.35, colour = "black") +
  geom_hline(yintercept = 3.5, colour = "black", linewidth = 0.45) +
  scale_fill_gradientn(
    colours = c("#EAF4F8", "#9FCBE2", "#F4E2CE", "#D77A61", "#9E3D2B"),
    values = scales::rescale(c(16, 24, 30, 55, 106)),
    limits = c(16, 106.5),
    oob = scales::squish,
    name = "Size\n(nm)"
  ) +
  labs(
    title = "c",
    subtitle = "DTX-SI-C18 fraction-dependent binary co-assembly",
    x = "DTX-SI-C18 fraction",
    y = NULL
  ) +
  base_theme +
  theme(axis.text.y = element_text(size = 6.8), legend.position = "right")

# -----------------------------
# 7. Panel d: Outcome map
# -----------------------------
p_d <- ggplot(binary_long, aes(x = DTX_fraction, y = Candidate_label, fill = Outcome)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Outcome == "Sub-30 nm", "<30", ">30")), size = 2.35) +
  geom_hline(yintercept = 3.5, colour = "black", linewidth = 0.45) +
  scale_fill_manual(values = pal_outcome, name = "DLS outcome") +
  labs(
    title = "d",
    subtitle = "Sub-30 nm classification across binary combinations",
    x = "DTX-SI-C18 fraction",
    y = NULL
  ) +
  base_theme +
  theme(axis.text.y = element_text(size = 6.8), legend.position = "right")

# -----------------------------
# 8. Panel e: Size response curves
# -----------------------------
p_e <- ggplot(binary_long, aes(x = DTX_numeric, y = Size, group = Candidate, colour = Initial_group)) +
  geom_hline(yintercept = 30, linetype = 2, linewidth = 0.45, colour = "grey40") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8) +
  facet_wrap(~Initial_group, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = pal_group) +
  scale_x_continuous(breaks = c(0, 25, 50, 75), limits = c(-2, 77)) +
  labs(
    title = "e",
    subtitle = "Composition-dependent size response",
    x = "DTX-SI-C18 fraction (%)",
    y = "Hydrodynamic diameter (nm)",
    colour = NULL
  ) +
  base_theme +
  theme(legend.position = "none")

# -----------------------------
# 9. Panel f: Minimum DTX fraction required for sub-30 nm assembly
# -----------------------------
p_f <- ggplot(summary_df, aes(x = Candidate_label, y = rescue_fraction, fill = Initial_group)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = Rescue_label), vjust = -0.25, size = 2.25) +
  scale_fill_manual(values = pal_group, name = NULL) +
  scale_y_continuous(limits = c(0, 82), breaks = c(0, 25, 50, 75), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "f",
    subtitle = "Minimum DTX-SI-C18 fraction for sub-30 nm assembly",
    x = NULL,
    y = "DTX-SI-C18 fraction (%)"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), legend.position = "top")

# -----------------------------
# 10. Panel g: Rescue magnitude in oversized candidates
# -----------------------------
p_g <- ggplot(oversized_summary, aes(y = Candidate_label)) +
  geom_segment(aes(x = min_size, xend = baseline_size, yend = Candidate_label), linewidth = 0.7, colour = "grey45") +
  geom_point(aes(x = baseline_size), size = 2.6, colour = col_orange_dark) +
  geom_point(aes(x = min_size), size = 2.6, colour = col_teal) +
  geom_text(aes(x = baseline_size, label = sprintf("%.1f", baseline_size)),
            nudge_x = 4.0, size = 2.2, colour = col_orange_dark) +
  geom_text(aes(x = min_size, label = sprintf("%.1f", min_size)),
            nudge_x = -4.0, size = 2.2, colour = col_teal) +
  geom_vline(xintercept = 30, linetype = 2, linewidth = 0.45, colour = "grey40") +
  scale_x_continuous(limits = c(0, 118), breaks = c(0, 30, 60, 90, 120)) +
  labs(
    title = "g",
    subtitle = "Size reduction in initially oversized candidates",
    x = "Hydrodynamic diameter (nm)",
    y = NULL
  ) +
  base_theme +
  annotate("text", x = 82, y = 0.45, label = "0% DTX", colour = col_orange_dark, size = 2.3) +
  annotate("text", x = 18, y = 0.45, label = "minimum", colour = col_teal, size = 2.3)

# -----------------------------
# 11. Panel h: Multidrug co-assembly size summary
# -----------------------------
p_h <- ggplot(multi_df, aes(x = Size, y = Formulation_short, fill = factor(No_prodrugs))) +
  geom_vline(xintercept = 30, linetype = 2, linewidth = 0.45, colour = "grey40") +
  geom_col(width = 0.68, colour = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(xmin = Size - SD, xmax = Size + SD),
    orientation = "y",
    width = 0.16,
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_point(
    data = multi_replicates_plot,
    aes(x = Size, y = Formulation_short),
    inherit.aes = FALSE,
    position = position_jitter(width = 0, height = 0.055, seed = 123),
    shape = 21,
    size = 1.7,
    fill = "white",
    colour = "black",
    stroke = 0.3
  ) +
  geom_text(
    aes(x = Size / 2, label = sprintf("%.1f nm", Size)),
    size = 2.3,
    colour = "black"
  ) +
  scale_fill_manual(values = pal_multidrug, name = "No. of\nprodrugs") +
  scale_x_continuous(limits = c(0, 34), breaks = c(0, 10, 20, 30), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "h",
    subtitle = "Representative multi-prodrug co-assembly",
    x = "Hydrodynamic diameter (nm)",
    y = NULL
  ) +
  base_theme +
  theme(axis.text.y = element_text(size = 6.4), legend.position = "right")

# -----------------------------
# 12. Optional Panel i: Exploratory composition-weighted score map
# -----------------------------
p_i <- ggplot(score_mix, aes(x = DTX_fraction, y = Candidate_label, fill = Exploratory_weighted_score)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = Score_label), size = 2.2, colour = "black") +
  geom_hline(yintercept = 3.5, colour = "black", linewidth = 0.45) +
  scale_fill_gradientn(
    colours = c(col_orange, col_blue2, col_teal),
    values = scales::rescale(c(6, 10, 13)),
    limits = c(6, 13),
    oob = scales::squish,
    name = "Exploratory\nscore"
  ) +
  labs(
    title = "i",
    subtitle = "Exploratory weighted single-prodrug score (not a validated mixed score)",
    x = "DTX-SI-C18 fraction",
    y = NULL
  ) +
  base_theme +
  theme(axis.text.y = element_text(size = 6.8), legend.position = "right")

# -----------------------------
# 13. Optional Panel j: Conceptual formula for rigorous mixed score
# -----------------------------

p_j <- ggplot() +
  theme_void(base_size = 8) +
  labs(title = "j") +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0)) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 6), expand = FALSE) +
  annotate("text", x = 0.5, y = 5.35, hjust = 0, size = 3.0, fontface = "bold",
           label = "Mixed 13-point score: recommended calculation") +
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 4.05, ymax = 5.00, fill = "#F7FBFD", colour = "grey45", linewidth = 0.35) +
  annotate("text", x = 0.8, y = 4.62, hjust = 0, size = 2.55,
           label = "Carrier terms:  ΔG_mix,k = Σ x_i ΔG_i,k  (9 prodrug–carrier terms)") +
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 2.80, ymax = 3.80, fill = "#FFF7EF", colour = "grey45", linewidth = 0.35) +
  annotate("text", x = 0.8, y = 3.38, hjust = 0, size = 2.48,
           label = "Self-association terms:  ΔG_mix,k = Σ x_i²ΔG_ii,k + 2Σx_ix_jΔG_ij,k") +
  annotate("text", x = 0.8, y = 2.95, hjust = 0, size = 2.15, colour = "grey30",
           label = "Requires hetero pairwise MD for ΔG_ij,k; apply the same 13 single-prodrug windows.") +
  annotate("segment", x = 5.0, xend = 5.0, y = 2.45, yend = 1.80,
           arrow = arrow(length = unit(0.14, "cm")), linewidth = 0.45) +
  annotate("rect", xmin = 2.1, xmax = 7.9, ymin = 0.75, ymax = 1.75, fill = "#EEF7F5", colour = "grey45", linewidth = 0.35) +
  annotate("text", x = 5.0, y = 1.40, size = 2.60, fontface = "bold",
           label = "13 binary terms → mixed compatibility score") +
  annotate("text", x = 5.0, y = 1.02, size = 2.15, colour = "grey30",
           label = "Use as prediction only after hetero terms are computed")

# -----------------------------
# 14. Optional Panel k: MD descriptor template
# -----------------------------
md_long <- md_template %>%
  pivot_longer(cols = c(System_Rg_nm, Drug_cluster_Rg_nm, Drug_polymer_contact_nm2, Drug_drug_contact_index),
               names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = recode(Metric,
                    "System_Rg_nm"= "System Rg (nm)",
                    "Drug_cluster_Rg_nm"= "Drug-cluster Rg (nm)",
                    "Drug_polymer_contact_nm2"= "Drug-polymer contact (nm²)",
                    "Drug_drug_contact_index"= "Drug-drug contact index"),
    System = factor(System, levels = md_template$System)
  )

if (all(is.na(md_long$Value))) {
  p_k <- ggplot() +
    theme_void(base_size = 8) +
    labs(title = "k") +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0)) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 5), expand = FALSE) +
    annotate("rect", xmin = 0.4, xmax = 9.6, ymin = 0.8, ymax = 4.2,
             fill = "#F7F7F7", colour = "grey50", linewidth = 0.35) +
    annotate("text", x = 5, y = 3.2, size = 3.0, fontface = "bold",
             label = "MD descriptor panel placeholder") +
    annotate("text", x = 5, y = 2.45, size = 2.35,
             label = "Add real DTX/Abi multi-molecule MD values here") +
    annotate("text", x = 5, y = 1.85, size = 2.15, colour = "grey30",
             label = "Recommended metrics: system Rg, drug-cluster Rg,\ndrug-polymer contact area, and drug-drug association")
} else {
  p_k <- ggplot(md_long, aes(x = DTX_fraction, y = Value, colour = System, group = System)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.9) +
    facet_wrap(~Metric, scales = "free_y", ncol = 2) +
    labs(title = "k", subtitle = "DTX/Abi multi-molecule MD descriptors",
         x = "DTX-SI-C18 fraction (%)", y = NULL, colour = NULL) +
    base_theme +
    theme(legend.position = "bottom")
}

# -----------------------------
# 15. Export individual panels
# -----------------------------
panels <- list(
  coassembly_workflow = p_a,
  single_prodrug_score_size = p_b,
  DTX_fraction_size_heatmap = p_c,
  coassembly_outcome_map = p_d,
  size_response = p_e,
  rescue_fraction = p_f,
  rescue_magnitude = p_g,
  multidrug_coassembly_size = p_h,
  exploratory_score_map = p_i,
  mixed_score_formula = p_j,
  MD_template = p_k
)

for (nm in names(panels)) {
  save_panel(panels[[nm]], nm, width = 4.7, height = 3.45)
}

# -----------------------------
# 16. Assemble complete figures
# -----------------------------

figure5_main <- (p_a | p_b) /
  (p_c | p_d) /
  (p_e | p_f) /
  (p_g | p_h) +
  plot_layout(widths = c(1.10, 1.0), heights = c(0.98, 1.05, 1.02, 1.00))

figure5_wide <- (p_a | p_c | p_h) /
  (p_b | p_f | p_g) +
  plot_layout(heights = c(1.0, 1.0))

figure5_with_model <- (p_a | p_c | p_h) /
  (p_j | p_i | p_k) +
  plot_layout(heights = c(1.0, 1.0))

# Main manuscript layout.
ggsave(file.path(figure_dir, "coassembly_formulation.pdf"), figure5_main,
       width = 10.4, height = 15.4, useDingbats = FALSE)
ggsave(file.path(figure_dir, "coassembly_formulation.png"), figure5_main,
       width = 10.4, height = 15.4, dpi = 600)
ggsave(file.path(figure_dir, "coassembly_formulation.tiff"), figure5_main,
       width = 10.4, height = 15.4, dpi = 600, compression = "lzw")

# Compact wide layout.
ggsave(file.path(figure_dir, "coassembly_formulation_wide.pdf"), figure5_wide,
       width = 15.0, height = 7.6, useDingbats = FALSE)
ggsave(file.path(figure_dir, "coassembly_formulation_wide.png"), figure5_wide,
       width = 15.0, height = 7.6, dpi = 600)
ggsave(file.path(figure_dir, "coassembly_formulation_wide.tiff"), figure5_wide,
       width = 15.0, height = 7.6, dpi = 600, compression = "lzw")

# Mechanistic/model-support layout.
ggsave(file.path(figure_dir, "coassembly_model_support.pdf"), figure5_with_model,
       width = 15.0, height = 7.6, useDingbats = FALSE)
ggsave(file.path(figure_dir, "coassembly_model_support.png"), figure5_with_model,
       width = 15.0, height = 7.6, dpi = 600)
ggsave(file.path(figure_dir, "coassembly_model_support.tiff"), figure5_with_model,
       width = 15.0, height = 7.6, dpi = 600, compression = "lzw")

# -----------------------------
# 17. Export processed data tables
# -----------------------------
write_csv(binary_replicates, file.path(table_dir, "binary_DTX_fraction_size_replicates.csv"))
write_csv(binary_long, file.path(table_dir, "binary_DTX_fraction_size_long.csv"))
write_csv(summary_df, file.path(table_dir, "binary_rescue_summary.csv"))
write_csv(oversized_summary, file.path(table_dir, "oversized_rescue_summary.csv"))
write_csv(multi_replicates, file.path(table_dir, "multidrug_coassembly_size_replicates.csv"))
write_csv(multi_df, file.path(table_dir, "multidrug_coassembly_size.csv"))
write_csv(single_si_plot, file.path(table_dir, "single_SI_size_score_anchor.csv"))
write_csv(score_mix, file.path(table_dir, "exploratory_composition_weighted_score.csv"))
write_csv(md_template, file.path(table_dir, "MD_descriptor_template_fill_with_real_values.csv"))

message("Figure 5 outputs: ", normalizePath(figure_dir, winslash = "/", mustWork = FALSE))
message("Figure 5 tables: ", normalizePath(table_dir, winslash = "/", mustWork = FALSE))
message("Recommended main figure: coassembly_formulation.pdf")
message("Optional model-support figure: coassembly_model_support.pdf")
message("MD panel currently uses placeholders unless real MD values are added to md_template.")
