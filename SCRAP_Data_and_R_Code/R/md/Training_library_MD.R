# Training-library pairwise and multi-molecule MD analysis

# ============================================================
# ============================================================

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
  
  return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
}

script_dir <- get_script_dir()
project_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
input_root <- file.path(project_root, "data", "md", "training_library_MD")
pair_fig_dir <- file.path(project_root, "results", "figures", "md", "training_library_MD", "pairwise")
pair_table_dir <- file.path(project_root, "results", "tables", "md", "training_library_MD", "pairwise")
dir.create(pair_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pair_table_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
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

section_output_dir <- pair_fig_dir

message("Running section: training-library pairwise prodrug-polymer MD")
message("Outputs will be saved to: ", normalizePath(section_output_dir, winslash = "/", mustWork = FALSE))

# Pairwise prodrug-polymer MD analysis

# -----------------------------
# 0. Working directory
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
    if (rstudioapi::isAvailable()) {
      script_path <- rstudioapi::getActiveDocumentContext()$path
      if (!is.null(script_path) && nzchar(script_path)) {
        return(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
      }
    }
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
# -----------------------------
# 1. Package setup
# -----------------------------
required_pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "patchwork", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing R package(s): ", paste(missing_pkgs, collapse = ", "),
    ". Please install them first, for example: install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# -----------------------------
# 2. File paths and settings
# -----------------------------
input_file <- file.path(input_root, "MD_drug1polymer1.xlsx")
if (!file.exists(input_file)) {
  stop(
    "Input file was not found: ", input_file,
    "\nExpected location: data/md/training_library_MD/MD_drug1polymer1.xlsx"
  )
}

output_dir <- section_output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

production_start_ns <- 50

rmsd_smooth_window <- 31

sheet_map <- list(
  distance = "distance",
  rmsd     = "rmsd",
  rg       = "rg",
  hbonds   = "hbonds",
  area     = "area_drug+pol-com"
)

available_sheets <- readxl::excel_sheets(input_file)
missing_sheets <- unlist(sheet_map, use.names = FALSE)[!unlist(sheet_map, use.names = FALSE) %in% available_sheets]
if (length(missing_sheets) > 0) {
  stop(
    "Missing required sheet(s): ", paste(missing_sheets, collapse = ", "),
    "\nAvailable sheets: ", paste(available_sheets, collapse = ", ")
  )
}

# -----------------------------
# 3. Global plotting settings
# -----------------------------
compound_levels <- c("PTX", "P2-Me", "P2-C6", "P2-C10", "P2-C14", "P2-C18", "Ful", "Ful-C18")

compound_palette <- c(
  "PTX"    = "#C9B77D",
  "P2-Me"  = "#E69F00",
  "P2-C6"  = "#A6611A",
  "P2-C10" = "#D95F5F",
  "P2-C14" = "#6B8E23",
  "P2-C18" = "#008C95",
  "Ful"    = "#4A4A4A",
  "Ful-C18"= "#7B61A8"
)

base_family <- "sans"

nature_theme_clean <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(size = base_size + 1, face = "plain", hjust = 0.5,
                                margin = margin(b = 2)),
      axis.title = element_text(size = base_size, color = "black"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.ticks.length = unit(1.4, "mm"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.title = element_text(size = base_size - 0.2, color = "black"),
      legend.text = element_text(size = base_size - 1, color = "black"),
      legend.key.size = unit(3.0, "mm"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

nature_theme_box <- function(base_size = 7, base_family = "sans") {
  nature_theme_clean(base_size = base_size, base_family = base_family) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      axis.line = element_blank()
    )
}

# -----------------------------
# 4. Data import and processing
# -----------------------------
read_metric_long <- function(sheet_name) {
  raw_df <- readxl::read_excel(input_file, sheet = sheet_name)
  if (ncol(raw_df) < 2) stop("Sheet ", sheet_name, "does not contain enough columns.")
  
  names(raw_df)[1] <- "Time_ns"
  
  raw_df %>%
    tidyr::pivot_longer(
      cols = -Time_ns,
      names_to = "Mol",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Time_ns = suppressWarnings(as.numeric(Time_ns)),
      Value   = suppressWarnings(as.numeric(Value))
    ) %>%
    dplyr::filter(!is.na(Time_ns), !is.na(Value), Mol %in% compound_levels) %>%
    dplyr::mutate(Mol = factor(Mol, levels = compound_levels))
}

get_production_window <- function(metric_long, start_ns = production_start_ns) {
  prod_df <- metric_long %>% dplyr::filter(Time_ns >= start_ns)
  if (nrow(prod_df) == 0) {
    warning("No rows were found at Time_ns >= ", start_ns, "; full trajectory was used instead.")
    return(metric_long)
  }
  prod_df
}

summarize_metric <- function(sheet_name, metric_name, start_ns = production_start_ns) {
  read_metric_long(sheet_name) %>%
    get_production_window(start_ns = start_ns) %>%
    dplyr::group_by(Mol) %>%
    dplyr::summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      N_frames = sum(!is.na(Value)),
      Start_ns = min(Time_ns, na.rm = TRUE),
      End_ns = max(Time_ns, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Metric = metric_name,
      Mol = factor(Mol, levels = compound_levels)
    )
}

summary_distance <- summarize_metric(sheet_map$distance, "Distance")
summary_rg       <- summarize_metric(sheet_map$rg,       "Rg")
summary_hbonds   <- summarize_metric(sheet_map$hbonds,   "H-bonds")
summary_area     <- summarize_metric(sheet_map$area,     "Contact area")

summary_table <- dplyr::bind_rows(summary_distance, summary_rg, summary_hbonds, summary_area) %>%
  dplyr::select(Metric, Mol, Mean, SD, N_frames, Start_ns, End_ns)

write.csv(summary_table, file.path(pair_table_dir, "pairwise_summary_statistics.csv"), row.names = FALSE)

# -----------------------------
# 5. Helper functions
# -----------------------------
auto_ymax <- function(summary_df, floor_lim = NA_real_, mult = 1.18) {
  max_val <- max(summary_df$Mean + summary_df$SD, na.rm = TRUE)
  lim <- max_val * mult
  if (!is.na(floor_lim)) lim <- max(lim, floor_lim)
  max(pretty(c(0, lim), n = 4))
}

smooth_series <- function(x, k = 31) {
  if (length(x) < k || k <= 1) return(x)
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
}

z_score_safe <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

prepare_summary_for_plot <- function(summary_df, digits = 2) {
  summary_df %>%
    dplyr::mutate(
      ymin = Mean,
      ymax = Mean + SD,
      label = sprintf(paste0("%.", digits, "f"), Mean)
    )
}

plot_summary_bar <- function(summary_df, y_lab, floor_lim = NA_real_, y_breaks = waiver(),
                             title = NULL, digits = 2) {
  plot_df <- prepare_summary_for_plot(summary_df, digits = digits)
  y_lim <- auto_ymax(plot_df, floor_lim = floor_lim)
  label_offset <- y_lim * 0.03
  
  ggplot(plot_df, aes(x = Mol, y = Mean, fill = Mol)) +
    geom_col(width = 0.58, color = "black", linewidth = 0.30, alpha = 0.96) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.16, linewidth = 0.30, color = "black") +
    geom_text(
      aes(y = pmin(ymax + label_offset, y_lim * 0.97), label = label),
      size = 2.05, vjust = 0, color = "black"
    ) +
    scale_fill_manual(values = compound_palette, drop = FALSE) +
    scale_y_continuous(
      limits = c(0, y_lim),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(x = NULL, y = y_lab, title = title) +
    guides(fill = "none") +
    coord_cartesian(clip = "off") +
    nature_theme_clean(base_size = 7, base_family = base_family)
}

# -----------------------------
# 6. Static summary panels (bars only)
# -----------------------------
p_distance_bar <- plot_summary_bar(
  summary_distance, y_lab = "Distance (nm)", floor_lim = 4,
  y_breaks = seq(0, 4, 1), digits = 2
)

p_rg_bar <- plot_summary_bar(
  summary_rg, y_lab = "Rg (nm)", floor_lim = 3,
  y_breaks = seq(0, 3, 1), digits = 2
)

p_hbonds_bar <- plot_summary_bar(
  summary_hbonds, y_lab = "H-bonds", floor_lim = 0.8,
  y_breaks = seq(0, 0.8, 0.2), digits = 2
)

p_area_bar <- plot_summary_bar(
  summary_area, y_lab = expression("Contact area (nm"^2*")"), floor_lim = 18,
  y_breaks = seq(0, 20, 5), digits = 2
)

# -----------------------------
# 7. RMSD panel
# -----------------------------
rmsd_long <- read_metric_long(sheet_map$rmsd)

rmsd_smooth <- rmsd_long %>%
  dplyr::arrange(Mol, Time_ns) %>%
  dplyr::group_by(Mol) %>%
  dplyr::mutate(Value_smooth = smooth_series(Value, k = rmsd_smooth_window)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(Value_smooth))

p_rmsd <- ggplot(rmsd_smooth, aes(x = Time_ns, y = Value_smooth, color = Mol)) +
  geom_line(linewidth = 0.38, alpha = 0.98, lineend = "round") +
  scale_color_manual(values = compound_palette, drop = FALSE, name = "Molecule") +
  scale_x_continuous(
    limits = c(0, 200),
    breaks = seq(0, 200, 50),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Time (ns)", y = "Smoothed RMSD (nm)", title = NULL) +
  nature_theme_clean(base_size = 7, base_family = base_family) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    legend.key.height = unit(2.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )

# -----------------------------
# 8. Direction-aligned heatmap + carrier-integration overview
# -----------------------------
metric_direction <- c(
  "Distance"    = -1,
  "Rg"          = -1,
  "H-bonds"     =  1,
  "Contact area"=  1
)

heatmap_df <- summary_table %>%
  dplyr::mutate(
    Metric = factor(Metric, levels = c("Distance", "Rg", "H-bonds", "Contact area")),
    Mol = factor(Mol, levels = compound_levels),
    Direction = metric_direction[as.character(Metric)]
  ) %>%
  dplyr::group_by(Metric) %>%
  dplyr::mutate(
    Row_z = z_score_safe(Mean),
    Integration_score = Row_z * Direction
  ) %>%
  dplyr::ungroup()

p_metric_heatmap <- ggplot(heatmap_df, aes(x = Mol, y = Metric, fill = Integration_score)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", Mean)), size = 2.10, color = "black") +
  scale_fill_gradient2(
    low = "#4C78A8", mid = "#F7F7F7", high = "#C44E52",
    midpoint = 0,
    limits = c(-1.8, 1.8),
    oob = scales::squish,
    breaks = c(-1, 0, 1),
    labels = c("weaker", "0", "stronger"),
    name = "Relative\nintegration\nscore"
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  nature_theme_box(base_size = 7, base_family = base_family) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7),
    legend.position = "right",
    legend.title = element_text(size = 6.5, lineheight = 0.95),
    legend.text = element_text(size = 6),
    legend.key.height = unit(10, "mm"),
    legend.key.width = unit(3, "mm")
  )

integration_df <- summary_table %>%
  dplyr::select(Mol, Metric, Mean) %>%
  dplyr::mutate(
    Metric = dplyr::recode(
      Metric,
      "Contact area"= "Contact_area",
      "H-bonds"= "H_bonds"
    )
  ) %>%
  tidyr::pivot_wider(names_from = Metric, values_from = Mean) %>%
  dplyr::mutate(
    Mol = factor(Mol, levels = compound_levels),
    Mol_label = as.character(Mol)
  )

label_offsets <- data.frame(
  Mol_label = compound_levels,
  dx = c(-0.08,  0.10,  0.10, -0.22, -0.08,  0.08,  0.08,  0.12),
  dy = c(-0.45,  0.55, -0.35,  0.65,  0.55,  0.60,  0.30,  0.35),
  stringsAsFactors = FALSE
)

label_df <- integration_df %>%
  dplyr::left_join(label_offsets, by = "Mol_label") %>%
  dplyr::mutate(
    label_x = Distance + dx,
    label_y = Contact_area + dy
  )

x_range <- range(c(integration_df$Distance, label_df$label_x), na.rm = TRUE)
y_range <- range(c(integration_df$Contact_area, label_df$label_y), na.rm = TRUE)

p_integration_map <- ggplot(integration_df, aes(x = Distance, y = Contact_area)) +
  geom_segment(
    data = label_df,
    aes(x = Distance, y = Contact_area, xend = label_x, yend = label_y),
    inherit.aes = FALSE, color = "grey45", linewidth = 0.18
  ) +
  geom_point(aes(size = Rg, fill = Mol), shape = 21, color = "black", linewidth = 0.35, alpha = 0.96) +
  geom_text(
    data = label_df,
    aes(x = label_x, y = label_y, label = Mol_label),
    inherit.aes = FALSE, size = 2.15, color = "black"
  ) +
  scale_fill_manual(values = compound_palette, drop = FALSE, name = "Molecule") +
  scale_size_continuous(range = c(2.8, 6.8), name = "Rg (nm)", breaks = pretty(integration_df$Rg, n = 4)) +
  scale_x_continuous(limits = c(x_range[1] - 0.10, x_range[2] + 0.10),
                     breaks = pretty(integration_df$Distance, n = 5),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_continuous(limits = c(y_range[1] - 0.35, y_range[2] + 0.35),
                     breaks = pretty(integration_df$Contact_area, n = 5),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    x = "Distance (nm)",
    y = expression("Contact area (nm"^2*")"),
    title = NULL
  ) +
  nature_theme_clean(base_size = 7, base_family = base_family) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(1.0, "mm")
  )

fig_heatmap_overview <- (p_metric_heatmap | p_integration_map) +
  patchwork::plot_layout(widths = c(1.08, 1.25)) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold", color = "black"))

# -----------------------------
# 9. Figure assembly
# -----------------------------
fig_bars_2x2 <- (p_distance_bar | p_rg_bar) / (p_hbonds_bar | p_area_bar) +
  patchwork::plot_annotation(tag_levels = list(c("b", "c", "d", "e"))) &
  theme(plot.tag = element_text(size = 9, face = "bold", color = "black"))

fig_with_rmsd_bar_compact <- (p_rmsd + labs(tag = "a")) /
  fig_bars_2x2 +
  patchwork::plot_layout(heights = c(1.05, 1.35))

fig_with_rmsd_bar_compact_tall <- (p_rmsd + labs(tag = "a")) /
  ((p_distance_bar | p_rg_bar) / (p_hbonds_bar | p_area_bar) +
     patchwork::plot_annotation(tag_levels = list(c("b", "c", "d", "e"))) &
     theme(plot.tag = element_text(size = 9, face = "bold", color = "black"))) +
  patchwork::plot_layout(heights = c(1.0, 1.5))

# -----------------------------
# 10. Save outputs
# -----------------------------
save_figure <- function(plot_obj, file_stem, width_mm, height_mm, dpi = 600) {
  pdf_file  <- file.path(output_dir, paste0(file_stem, ".pdf"))
  png_file  <- file.path(output_dir, paste0(file_stem, ".png"))
  tiff_file <- file.path(output_dir, paste0(file_stem, ".tiff"))
  
  tryCatch(
    suppressMessages(
      ggsave(pdf_file, plot_obj, width = width_mm, height = height_mm,
             units = "mm", device = grDevices::cairo_pdf, bg = "white")
    ),
    error = function(e) suppressMessages(
      ggsave(pdf_file, plot_obj, width = width_mm, height = height_mm,
             units = "mm", device = "pdf", bg = "white")
    )
  )
  
  suppressMessages(
    ggsave(png_file, plot_obj, width = width_mm, height = height_mm,
           units = "mm", dpi = dpi, bg = "white")
  )
  suppressMessages(
    ggsave(tiff_file, plot_obj, width = width_mm, height = height_mm,
           units = "mm", dpi = dpi, compression = "lzw", bg = "white")
  )
}

save_figure(p_rmsd,         "pairwise_RMSD",  width_mm = 125, height_mm = 62, dpi = 600)
save_figure(p_distance_bar, "pairwise_distance",           width_mm = 65,  height_mm = 55, dpi = 600)
save_figure(p_rg_bar,       "pairwise_Rg",                 width_mm = 65,  height_mm = 55, dpi = 600)
save_figure(p_hbonds_bar,   "pairwise_Hbonds",             width_mm = 65,  height_mm = 55, dpi = 600)
save_figure(p_area_bar,         "pairwise_contact_area",        width_mm = 65,  height_mm = 55, dpi = 600)
save_figure(p_metric_heatmap,   "pairwise_metric_heatmap",          width_mm = 105, height_mm = 55, dpi = 600)
save_figure(p_integration_map,  "pairwise_carrier_integration_map", width_mm = 105, height_mm = 75, dpi = 600)

save_figure(fig_bars_2x2,                "pairwise_summary_bars", width_mm = 150, height_mm = 118, dpi = 600)
save_figure(fig_with_rmsd_bar_compact,   "pairwise_MD_support",   width_mm = 180, height_mm = 150, dpi = 600)
save_figure(fig_with_rmsd_bar_compact_tall, "pairwise_MD_support_tall", width_mm = 180, height_mm = 160, dpi = 600)
save_figure(fig_heatmap_overview, "pairwise_heatmap_overview", width_mm = 190, height_mm = 82, dpi = 600)

notes <- c(
  "Training-library pairwise MD output notes",
  paste0("Working directory: ", getwd()),
  paste0("Input file: ", input_file),
  paste0("Output directory: ", output_dir),
  paste0("Static summary metrics were calculated from Time_ns >= ", production_start_ns, "ns."),
  paste0("RMSD curves were smoothed using a centered rolling mean window of ", rmsd_smooth_window, "frames."),
  "1) Heatmap colour now uses one direction-aligned, row-normalized integration score.",
  "2) Distance and Rg were sign-reversed before colouring because lower values indicate closer/compact organization.",
  "3) H-bonds and contact area retain positive direction because higher values indicate stronger interaction in this overview.",
  "4) Raw mean values remain printed inside the tiles; only the fill colour is transformed.",
  "5) Upper-only SD bars and mean-value labels are retained for all bar charts.",
  "Recommended files:",
  "- pairwise_MD_support",
  "- pairwise_heatmap_overview"
)
writeLines(notes, file.path(pair_table_dir, "pairwise_notes.txt"))

# End of script

# ============================================================
# ============================================================

rm(list = ls())
graphics.off()

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
  
  return(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
}

script_dir <- get_script_dir()
project_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
input_root <- file.path(project_root, "data", "md", "training_library_MD")
multi_fig_dir <- file.path(project_root, "results", "figures", "md", "training_library_MD", "multi_molecule")
multi_table_dir <- file.path(project_root, "results", "tables", "md", "training_library_MD", "multi_molecule")
dir.create(multi_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(multi_table_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
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

section_output_dir <- multi_fig_dir

message("Running section: training-library multi-molecule prodrug-polymer MD")
message("Outputs will be saved to: ", normalizePath(section_output_dir, winslash = "/", mustWork = FALSE))

# ============================================================
# Multi-molecule prodrug-polymer MD analysis
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
install_missing <- TRUE

packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "stringr", "tibble",
  "ggplot2", "patchwork", "scales", "zoo", "grid", "readr",
  "viridisLite"
)

if (install_missing) {
  to_install <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
  }
}

invisible(lapply(packages, library, character.only = TRUE))

# -----------------------------
# 1. Paths and user-adjustable settings
# -----------------------------
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grepl(file_arg, cmd_args)])
  if (length(script_path) > 0 && nzchar(script_path[1])) {
    return(dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)))
  }
  
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)))
  }
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(active_path) && nzchar(active_path)) {
      return(dirname(normalizePath(active_path, winslash = "/", mustWork = FALSE)))
    }
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
input_xlsx <- file.path(input_root, "MD_drugspolymers.xlsx")
if (!file.exists(input_xlsx)) {
  stop(paste0(
    "Cannot find MD_drugspolymers.xlsx: ", input_xlsx
  ))
}

out_dir <- section_output_dir
fig_dir <- multi_fig_dir
panel_dir <- file.path(fig_dir, "individual_panels")
data_dir <- multi_table_dir
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

drug_order <- c("PTX", "P2-Me", "P2-C6", "P2-C10", "P2-C14", "P2-C18", "Ful", "Ful-C18")
taxane_order <- c("PTX", "P2-Me", "P2-C6", "P2-C10", "P2-C14", "P2-C18")
ful_order <- c("Ful", "Ful-C18")

assembly_tbl <- tibble::tibble(
  system = drug_order,
  assembly_class = factor(
    c("Precipitated", "Precipitated", "Sub-30 nm", "Sub-30 nm",
      "Sub-30 nm", "Sub-30 nm", "Sub-30 nm", "Oversized"),
    levels = c("Sub-30 nm", "Oversized", "Precipitated")
  )
)

drug_palette <- c(
  "PTX"     = "#111111",
  "P2-Me"  = "#7B61FF",
  "P2-C6"  = "#00A6D6",
  "P2-C10" = "#009E73",
  "P2-C14" = "#0072B2",
  "P2-C18" = "#332288",
  "Ful"    = "#CC79A7",
  "Ful-C18"= "#D55E00"
)

drug_linetype <- c(
  "PTX"     = "solid",
  "P2-Me"  = "solid",
  "P2-C6"  = "solid",
  "P2-C10" = "solid",
  "P2-C14" = "longdash",
  "P2-C18" = "solid",
  "Ful"    = "twodash",
  "Ful-C18"= "solid"
)

class_palette <- c(
  "Sub-30 nm"   = "#0072B2",
  "Oversized"   = "#D55E00",
  "Precipitated"= "#7F7F7F"
)

final150_start_ns <- 50
final50_start_ns <- 150
rdf_xmax <- 4.0
rmsd_smooth_window_ns <- 7.5
rdf_smooth_window_nm <- 0.035
rdf_plot_rmin <- 0.05

combined_width_in <- 15.2
combined_height_in <- 17.0
panel_dpi_png <- 600
panel_dpi_tiff <- 900

# -----------------------------
# 2. Plot theme and export helpers
# -----------------------------
theme_nature <- function(base_size = 8, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1.2, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.2, color = "grey25", hjust = 0),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1, color = "grey15"),
      axis.line = ggplot2::element_line(linewidth = 0.25, color = "grey20"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, color = "grey20"),
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.key.width = grid::unit(6, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.5),
      panel.grid.major = ggplot2::element_line(linewidth = 0.15, color = "grey91"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

theme_set(theme_nature())

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"

save_plot_all <- function(plot, prefix, width, height, dir = panel_dir) {
  pdf_file <- file.path(dir, paste0(prefix, ".pdf"))
  tiff_file <- file.path(dir, paste0(prefix, ".tiff"))
  png_file <- file.path(dir, paste0(prefix, ".png"))
  
  ggplot2::ggsave(
    filename = pdf_file, plot = plot, width = width, height = height,
    units = "in", device = pdf_device
  )
  ggplot2::ggsave(
    filename = tiff_file, plot = plot, width = width, height = height,
    units = "in", dpi = panel_dpi_tiff, compression = "lzw", bg = "white"
  )
  ggplot2::ggsave(
    filename = png_file, plot = plot, width = width, height = height,
    units = "in", dpi = panel_dpi_png, bg = "white"
  )
  
  invisible(c(pdf = pdf_file, tiff = tiff_file, png = png_file))
}

# -----------------------------
# 3. Data import helpers
# -----------------------------
clean_system_names <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("\\s+", "")
}

read_wide_numeric_sheet <- function(sheet_name) {
  dat <- readxl::read_excel(input_xlsx, sheet = sheet_name, .name_repair = "unique")
  first_col <- names(dat)[1]
  
  dat %>%
    dplyr::rename(x = dplyr::all_of(first_col)) %>%
    tidyr::pivot_longer(
      cols = -x,
      names_to = "system",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      system = clean_system_names(system),
      x = suppressWarnings(as.numeric(x)),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    dplyr::filter(!is.na(system), system != "", system %in% drug_order)
}

summarise_time_sheet <- function(sheet_name, metric, metric_label, module, unit,
                                 start_ns = 50, x_label = "Time (ns)") {
  raw_long <- read_wide_numeric_sheet(sheet_name) %>%
    dplyr::rename(time = x, raw_value = value) %>%
    dplyr::filter(!is.na(time), time >= start_ns)
  
  raw_long %>%
    dplyr::group_by(system) %>%
    dplyr::summarise(
      n = sum(!is.na(raw_value)),
      value = mean(raw_value, na.rm = TRUE),
      sd = if (dplyr::n() > 1 && sum(!is.na(raw_value)) > 1) stats::sd(raw_value, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      sd = dplyr::if_else(is.nan(sd), NA_real_, sd),
      metric = metric,
      metric_label = metric_label,
      module = module,
      unit = unit,
      source_sheet = sheet_name,
      window = paste0("raw Excel values with time >= ", start_ns, "ns"),
      x_label = x_label
    )
}

read_mmpbsa_sheet <- function() {
  dat <- readxl::read_excel(input_xlsx, sheet = "MMPBSA", .name_repair = "unique")
  
  dat %>%
    dplyr::slice(1) %>%
    dplyr::select(-1) %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = "system",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      system = clean_system_names(system),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    dplyr::filter(system %in% drug_order) %>%
    dplyr::mutate(
      sd = NA_real_,
      n = NA_integer_,
      metric = "mmpbsa_pol_drug",
      metric_label = "ΔG Pro–Pol",
      module = "Thermodynamic affinity",
      unit = "kcal/mol",
      source_sheet = "MMPBSA",
      window = "MM-PBSA value",
      x_label = NA_character_
    )
}

# -----------------------------
# 4. Summarize non-RDF metrics
# -----------------------------
metric_order_tbl <- tibble::tribble(
  ~metric,              ~metric_label,                   ~metric_label_heat,                       ~module,                    ~unit,
  "contact_area",       "Contact area",                  "Contact area\n(nm²)",                    "Polymer integration",      "nm²",
  "distance_drug_pol",  "COM distance",                  "COM distance\n(nm)",                    "Polymer integration",      "nm",
  "hbonds_drug_pol",    "H-bonds Pro–Pol",              "H-bonds\nPro–Pol",                    "Polymer integration",      "count",
  "hbonds_drug_pcl",    "H-bonds Pro–PCL",              "H-bonds\nPro–PCL",                    "Segment partitioning",     "count",
  "hbonds_drug_peg",    "H-bonds Pro–PEG",              "H-bonds\nPro–PEG",                    "Segment partitioning",     "count",
  "hbonds_drug_drug",   "H-bonds Pro–Pro",             "H-bonds\nPro–Pro",                   "Self-association",         "count",
  "mmpbsa_pol_drug",    "ΔG Pro–Pol",                   "ΔG Pro–Pol\n(kcal/mol)",             "Thermodynamic affinity",   "kcal/mol",
  "rmsd",               "RMSD",                          "RMSD\n(nm)",                          "End-state dynamics",       "nm",
  "rg_drugs",           "Rg drug cluster",               "Rg drug cluster\n(nm)",               "End-state dynamics",       "nm",
  "rg_system",          "Rg whole system",               "Rg whole system\n(nm)",               "End-state dynamics",       "nm"
)

metric_summary <- dplyr::bind_rows(
  summarise_time_sheet(
    sheet_name = "area_drug+pol-com",
    metric = "contact_area",
    metric_label = "Contact area",
    module = "Polymer integration",
    unit = "nm²",
    start_ns = final150_start_ns
  ),
  summarise_time_sheet(
    sheet_name = "drug-pol_distance",
    metric = "distance_drug_pol",
    metric_label = "COM distance",
    module = "Polymer integration",
    unit = "nm",
    start_ns = final150_start_ns
  ),
  summarise_time_sheet(
    sheet_name = "hbonds-drugpol",
    metric = "hbonds_drug_pol",
    metric_label = "H-bonds Pro–Pol",
    module = "Polymer integration",
    unit = "count",
    start_ns = final150_start_ns
  ),
  summarise_time_sheet(
    sheet_name = "hbonds-drugpcl",
    metric = "hbonds_drug_pcl",
    metric_label = "H-bonds Pro–PCL",
    module = "Segment partitioning",
    unit = "count",
    start_ns = final150_start_ns
  ),
  summarise_time_sheet(
    sheet_name = "hbonds-drugpeg",
    metric = "hbonds_drug_peg",
    metric_label = "H-bonds Pro–PEG",
    module = "Segment partitioning",
    unit = "count",
    start_ns = final150_start_ns
  ),
  summarise_time_sheet(
    sheet_name = "hbonds-drugdrug",
    metric = "hbonds_drug_drug",
    metric_label = "H-bonds Pro–Pro",
    module = "Self-association",
    unit = "count",
    start_ns = final150_start_ns
  ),
  read_mmpbsa_sheet(),
  summarise_time_sheet(
    sheet_name = "rmsd",
    metric = "rmsd",
    metric_label = "RMSD",
    module = "End-state dynamics",
    unit = "nm",
    start_ns = final50_start_ns
  ),
  summarise_time_sheet(
    sheet_name = "rg-drugs",
    metric = "rg_drugs",
    metric_label = "Rg drug cluster",
    module = "End-state dynamics",
    unit = "nm",
    start_ns = final50_start_ns
  ),
  summarise_time_sheet(
    sheet_name = "rg-system",
    metric = "rg_system",
    metric_label = "Rg whole system",
    module = "End-state dynamics",
    unit = "nm",
    start_ns = final50_start_ns
  )
) %>%
  dplyr::left_join(metric_order_tbl, by = c("metric", "metric_label", "module", "unit")) %>%
  dplyr::left_join(assembly_tbl, by = "system") %>%
  dplyr::mutate(
    system = factor(system, levels = drug_order),
    module = factor(
      module,
      levels = c("Polymer integration", "Segment partitioning",
                 "Self-association", "Thermodynamic affinity", "End-state dynamics")
    ),
    metric = factor(metric, levels = metric_order_tbl$metric),
    metric_label_heat = factor(metric_label_heat, levels = metric_order_tbl$metric_label_heat)
  )

readr::write_csv(
  metric_summary %>%
    dplyr::arrange(metric, system) %>%
    dplyr::select(system, assembly_class, module, metric, metric_label, value, sd, n, unit, source_sheet, window),
  file.path(data_dir, "multi_molecule_summary.csv")
)

format_metric_value <- function(metric, value) {
  dplyr::case_when(
    is.na(value) ~ "",
    metric == "contact_area"~ sprintf("%.0f", value),
    metric == "mmpbsa_pol_drug"~ sprintf("%.1f", value),
    metric %in% c("distance_drug_pol", "rmsd", "rg_drugs", "rg_system") ~ sprintf("%.2f", value),
    TRUE ~ sprintf("%.1f", value)
  )
}

make_one_sided_error <- function(value, sd) {
  has_sd <- !is.na(sd) & is.finite(sd) & sd > 0
  err_cap <- dplyr::case_when(
    !has_sd ~ NA_real_,
    value >= 0 ~ value + sd,
    TRUE ~ value - sd
  )
  
  tibble::tibble(
    has_sd = has_sd,
    err_start = dplyr::if_else(has_sd, value, NA_real_),
    err_cap = err_cap,
    err_low = pmin(value, err_cap, na.rm = TRUE),
    err_high = pmax(value, err_cap, na.rm = TRUE)
  ) %>%
    dplyr::mutate(
      err_low = dplyr::if_else(has_sd, err_low, NA_real_),
      err_high = dplyr::if_else(has_sd, err_high, NA_real_)
    )
}

heatmap_low_col <- "#F7F7F7"
heatmap_mid_col <- "#B9D8EA"
heatmap_high_col <- "#145C8A"

heatmap_df <- metric_summary %>%
  dplyr::mutate(
    metric_family = dplyr::case_when(
      stringr::str_detect(as.character(metric), "^hbonds_") ~ "H-bonds",
      as.character(metric) %in% c("rg_drugs", "rg_system") ~ "Rg",
      as.character(metric) == "contact_area"~ "Contact area",
      as.character(metric) == "distance_drug_pol"~ "COM distance",
      as.character(metric) == "mmpbsa_pol_drug"~ "ΔG magnitude",
      as.character(metric) == "rmsd"~ "RMSD",
      TRUE ~ as.character(metric)
    ),
    value_for_color = dplyr::if_else(as.character(metric) == "mmpbsa_pol_drug", -value, value)
  ) %>%
  dplyr::group_by(metric_family) %>%
  dplyr::mutate(
    metric_min = min(value_for_color, na.rm = TRUE),
    metric_max = max(value_for_color, na.rm = TRUE),
    value_norm = dplyr::if_else(
      is.finite(metric_min) & is.finite(metric_max) & abs(metric_max - metric_min) > .Machine$double.eps,
      (value_for_color - metric_min) / (metric_max - metric_min),
      0.5
    ),
    tile_label = format_metric_value(as.character(metric), value),
    text_col = dplyr::if_else(!is.na(value_norm) & value_norm > 0.72, "white", "black")
  ) %>%
  dplyr::ungroup()

readr::write_csv(
  heatmap_df %>%
    dplyr::arrange(metric, system) %>%
    dplyr::select(system, assembly_class, module, metric_family, metric, metric_label_heat, value, value_for_color, sd, n, unit,
                  metric_min, metric_max, value_norm, tile_label, source_sheet, window),
  file.path(data_dir, "multi_molecule_metric_heatmap_data.csv")
)

# -----------------------------
# 5. Read RDF data and calculate RDF summaries
# -----------------------------
rdf_sheets <- c(
  "Pro–Pol"= "rdf-drugpol",
  "Pro–PCL"= "rdf-drugpcl",
  "Pro–PEG"= "rdf-drugpeg",
  "Pro–Pro"= "rdf-drugs"
)

rdf_long <- purrr::imap_dfr(
  rdf_sheets,
  function(sheet_name, rdf_label) {
    read_wide_numeric_sheet(sheet_name) %>%
      dplyr::rename(r = x) %>%
      dplyr::mutate(
        rdf_metric = rdf_label,
        rdf_metric = factor(rdf_metric, levels = names(rdf_sheets)),
        system = factor(system, levels = drug_order)
      )
  }
) %>%
  dplyr::filter(!is.na(r), !is.na(value), r <= rdf_xmax)

rdf_peak_summary <- rdf_long %>%
  dplyr::group_by(system, rdf_metric) %>%
  dplyr::slice_max(order_by = value, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::rename(gmax = value, r_at_gmax = r)

rdf_auc_summary <- rdf_long %>%
  dplyr::arrange(system, rdf_metric, r) %>%
  dplyr::group_by(system, rdf_metric) %>%
  dplyr::summarise(
    auc_0_to_xmax = sum(diff(r) * zoo::rollmean(value, 2), na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(rdf_peak_summary, file.path(data_dir, "multi_molecule_RDF_peak_summary.csv"))
readr::write_csv(rdf_auc_summary, file.path(data_dir, "multi_molecule_RDF_AUC_summary.csv"))
readr::write_csv(rdf_long, file.path(data_dir, "multi_molecule_RDF_curves_raw.csv"))

rdf_dr <- median(diff(sort(unique(rdf_long$r))), na.rm = TRUE)
rdf_window_n <- max(5, round(rdf_smooth_window_nm / rdf_dr))
if (rdf_window_n %% 2 == 0) rdf_window_n <- rdf_window_n + 1

rdf_smooth <- rdf_long %>%
  dplyr::arrange(system, rdf_metric, r) %>%
  dplyr::group_by(system, rdf_metric) %>%
  dplyr::mutate(
    value_smooth = zoo::rollmean(value, k = rdf_window_n, fill = NA, align = "center")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(r >= rdf_plot_rmin)

readr::write_csv(rdf_smooth, file.path(data_dir, "multi_molecule_RDF_curves_smoothed.csv"))

# -----------------------------
# 6. Panel a: smoothed RMSD trajectories
# -----------------------------
rmsd_raw <- read_wide_numeric_sheet("rmsd") %>%
  dplyr::rename(time = x) %>%
  dplyr::mutate(system = factor(system, levels = drug_order)) %>%
  dplyr::arrange(system, time)

rmsd_dt <- median(diff(sort(unique(rmsd_raw$time))), na.rm = TRUE)
rmsd_window_n <- max(5, round(rmsd_smooth_window_ns / rmsd_dt))
if (rmsd_window_n %% 2 == 0) rmsd_window_n <- rmsd_window_n + 1

rmsd_smooth <- rmsd_raw %>%
  dplyr::group_by(system) %>%
  dplyr::arrange(time, .by_group = TRUE) %>%
  dplyr::mutate(
    value_smooth = zoo::rollmean(value, k = rmsd_window_n, fill = NA, align = "center")
  ) %>%
  dplyr::ungroup()

readr::write_csv(rmsd_smooth, file.path(data_dir, "multi_molecule_RMSD_smoothed_data.csv"))

p_rmsd <- ggplot2::ggplot(
  rmsd_smooth,
  ggplot2::aes(x = time, y = value_smooth, color = system, linetype = system)
) +
  ggplot2::geom_line(linewidth = 0.58, alpha = 0.98, na.rm = TRUE) +
  ggplot2::scale_color_manual(values = drug_palette, drop = FALSE) +
  ggplot2::scale_linetype_manual(values = drug_linetype, drop = FALSE) +
  ggplot2::coord_cartesian(xlim = c(0, 200)) +
  ggplot2::labs(
    title = "End-state dynamics from smoothed RMSD trajectories",
    subtitle = paste0("Rolling mean = ", rmsd_smooth_window_ns, "ns; raw trajectories are kept in the exported data file"),
    x = "Time (ns)",
    y = "RMSD (nm)",
    color = NULL,
    linetype = NULL
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(linewidth = 1.4, linetype = unname(drug_linetype[drug_order]))
    ),
    linetype = "none"
  ) +
  theme_nature(base_size = 8) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box.margin = ggplot2::margin(-3, 0, 0, 0),
    panel.grid.major = ggplot2::element_line(linewidth = 0.15, color = "grey91")
  )

# -----------------------------
# 7. Panel b: metric-wise heatmap with mean values in tiles
# -----------------------------
p_metric_heatmap <- ggplot2::ggplot(
  heatmap_df,
  ggplot2::aes(x = system, y = metric_label_heat, fill = value_norm)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.45, width = 0.96, height = 0.92) +
  ggplot2::geom_text(
    ggplot2::aes(label = tile_label, color = text_col),
    size = 2.35,
    fontface = "bold",
    na.rm = TRUE
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(module),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_fill_gradientn(
    colours = c(heatmap_low_col, heatmap_mid_col, heatmap_high_col),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("Low", "Mid", "High"),
    na.value = "grey94",
    name = "Relative within\neach metric family"
  ) +
  ggplot2::scale_color_manual(values = c("black"= "grey10", "white"= "white"), guide = "none") +
  ggplot2::labs(
    title = "Segmented non-RDF metric fingerprint",
    subtitle = "Color is scaled within each metric family; numbers indicate the original mean values",
    x = NULL,
    y = NULL
  ) +
  theme_nature(base_size = 8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = ggplot2::element_text(size = 6.8, color = "grey10"),
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, color = "grey20", size = 7.2),
    panel.spacing.y = grid::unit(2.2, "mm"),
    legend.position = "right",
    panel.grid = ggplot2::element_blank()
  )

# -----------------------------
# 8. Panel c: H-bonding summary as absolute values
# -----------------------------
hbond_base <- metric_summary %>%
  dplyr::filter(metric %in% c(
    "hbonds_drug_pol", "hbonds_drug_pcl", "hbonds_drug_peg", "hbonds_drug_drug"
  )) %>%
  dplyr::mutate(
    metric_label_bar = factor(
      metric_label,
      levels = c("H-bonds Pro–Pol", "H-bonds Pro–PCL", "H-bonds Pro–PEG", "H-bonds Pro–Pro"),
      labels = c("Pro–Pol", "Pro–PCL", "Pro–PEG", "Pro–Pro")
    )
  )

hbond_df <- dplyr::bind_cols(
  hbond_base,
  make_one_sided_error(hbond_base$value, hbond_base$sd)
) %>%
  dplyr::group_by(metric_label_bar) %>%
  dplyr::mutate(
    value_label = format_metric_value(as.character(metric), value),
    label_pad = pmax(0.06 * (max(c(value, err_high), na.rm = TRUE) - min(c(value, err_low), na.rm = TRUE)), 0.06),
    label_y = dplyr::if_else(has_sd, err_high + label_pad, value + label_pad)
  ) %>%
  dplyr::ungroup()

readr::write_csv(
  hbond_df %>%
    dplyr::arrange(metric, system) %>%
    dplyr::select(system, assembly_class, metric, metric_label_bar, value, sd, n, err_start, err_cap, err_low, err_high, unit, source_sheet, window),
  file.path(data_dir, "multi_molecule_Hbond_summary_data.csv")
)

p_hbond <- ggplot2::ggplot(
  hbond_df,
  ggplot2::aes(x = system, y = value, fill = assembly_class)
) +
  ggplot2::geom_col(width = 0.68, color = "grey25", linewidth = 0.18) +
  ggplot2::geom_linerange(
    ggplot2::aes(ymin = err_low, ymax = err_high),
    linewidth = 0.42,
    color = "grey10",
    na.rm = TRUE
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = err_cap, ymax = err_cap),
    width = 0.22,
    linewidth = 0.42,
    color = "grey10",
    na.rm = TRUE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = label_y, label = value_label),
    vjust = 0,
    size = 2.1,
    color = "grey10",
    na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~metric_label_bar, ncol = 2, scales = "free_y") +
  ggplot2::scale_fill_manual(values = class_palette, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
  ggplot2::labs(
    title = "H-bonding distinguishes polymer contacts from prodrug self-association",
    subtitle = "Mean ± s.d. calculated from raw Excel values; only the terminal one-sided s.d. bar is shown",
    x = NULL,
    y = "H-bonds",
    fill = "Assembly"
  ) +
  theme_nature(base_size = 8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.box.margin = ggplot2::margin(-3, 0, 0, 0),
    panel.grid.major.x = ggplot2::element_blank()
  )

# -----------------------------
# 9. Panel d: structural and thermodynamic summary
# -----------------------------
struct_metric_levels <- c(
  "Contact area", "COM distance", "ΔG Pro–Pol", "Rg drug cluster", "Rg whole system", "RMSD"
)

struct_base <- metric_summary %>%
  dplyr::filter(metric %in% c(
    "contact_area", "distance_drug_pol", "mmpbsa_pol_drug", "rg_drugs", "rg_system", "rmsd"
  )) %>%
  dplyr::mutate(
    metric_label_bar = factor(metric_label, levels = struct_metric_levels),
    value_label = format_metric_value(as.character(metric), value)
  )

struct_df <- dplyr::bind_cols(
  struct_base,
  make_one_sided_error(struct_base$value, struct_base$sd)
) %>%
  dplyr::group_by(metric_label_bar) %>%
  dplyr::mutate(
    label_pad = pmax(0.04 * (max(c(value, err_high), na.rm = TRUE) - min(c(value, err_low), na.rm = TRUE)), 0.06),
    label_y = dplyr::case_when(
      has_sd & value >= 0 ~ err_high + label_pad,
      has_sd & value < 0 ~ err_low - label_pad,
      value >= 0 ~ value + label_pad,
      TRUE ~ value - label_pad
    ),
    label_vjust = dplyr::if_else(value >= 0, 0, 1)
  ) %>%
  dplyr::ungroup()

readr::write_csv(
  struct_df %>%
    dplyr::arrange(metric, system) %>%
    dplyr::select(system, assembly_class, metric, metric_label_bar, value, sd, n, err_start, err_cap, err_low, err_high, unit, source_sheet, window),
  file.path(data_dir, "multi_molecule_structural_energy_summary_data.csv")
)

p_struct <- ggplot2::ggplot(
  struct_df,
  ggplot2::aes(x = system, y = value, fill = assembly_class)
) +
  ggplot2::geom_col(width = 0.68, color = "grey25", linewidth = 0.18) +
  ggplot2::geom_linerange(
    ggplot2::aes(ymin = err_low, ymax = err_high),
    linewidth = 0.42,
    color = "grey10",
    na.rm = TRUE
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = err_cap, ymax = err_cap),
    width = 0.22,
    linewidth = 0.42,
    color = "grey10",
    na.rm = TRUE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(y = label_y, label = value_label, vjust = label_vjust),
    size = 2.05,
    color = "grey10",
    na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~metric_label_bar, ncol = 3, scales = "free_y") +
  ggplot2::scale_fill_manual(values = class_palette, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.16, 0.18))) +
  ggplot2::labs(
    title = "Absolute structural and thermodynamic readouts across the eight systems",
    subtitle = "Mean ± s.d. calculated from raw Excel values where available; MM-PBSA is plotted without artificial s.d.",
    x = NULL,
    y = "Mean value",
    fill = "Assembly"
  ) +
  theme_nature(base_size = 8) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.box.margin = ggplot2::margin(-3, 0, 0, 0),
    panel.grid.major.x = ggplot2::element_blank()
  )

# -----------------------------
# 10. Panels e/f: RDF profiles separated as requested
# -----------------------------
rdf_theme <- theme_nature(base_size = 8) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box.margin = ggplot2::margin(-3, 0, 0, 0),
    panel.spacing.x = grid::unit(3, "mm"),
    panel.grid.major = ggplot2::element_line(linewidth = 0.12, color = "grey92")
  )

rdf_taxane_data <- rdf_smooth %>%
  dplyr::filter(system %in% taxane_order) %>%
  dplyr::mutate(system = factor(system, levels = taxane_order))

rdf_ful_data <- rdf_smooth %>%
  dplyr::filter(system %in% ful_order) %>%
  dplyr::mutate(system = factor(system, levels = ful_order))

readr::write_csv(rdf_taxane_data, file.path(data_dir, "multi_molecule_RDF_PTX_P2_data.csv"))
readr::write_csv(rdf_ful_data, file.path(data_dir, "multi_molecule_RDF_Ful_FulC18_data.csv"))

p_rdf_taxane <- ggplot2::ggplot(rdf_taxane_data, ggplot2::aes(x = r, y = value_smooth, color = system, linetype = system)) +
  ggplot2::geom_line(linewidth = 0.55, alpha = 0.98, na.rm = TRUE) +
  ggplot2::facet_wrap(~rdf_metric, nrow = 1, scales = "free_y") +
  ggplot2::scale_color_manual(values = drug_palette[taxane_order], drop = FALSE) +
  ggplot2::scale_linetype_manual(values = drug_linetype[taxane_order], drop = FALSE) +
  ggplot2::coord_cartesian(xlim = c(rdf_plot_rmin, rdf_xmax)) +
  ggplot2::labs(
    title = "Segment-resolved RDF profiles: PTX and P2 series",
    subtitle = paste0("Rolling mean = ", rdf_smooth_window_nm, "nm; r < ", rdf_plot_rmin, "nm is omitted to remove near-zero RDF spikes"),
    x = "r (nm)",
    y = "g(r)",
    color = NULL,
    linetype = NULL
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(linewidth = 1.3, linetype = unname(drug_linetype[taxane_order]))
    ),
    linetype = "none"
  ) +
  rdf_theme

p_rdf_ful <- ggplot2::ggplot(rdf_ful_data, ggplot2::aes(x = r, y = value_smooth, color = system, linetype = system)) +
  ggplot2::geom_line(linewidth = 0.65, alpha = 0.98, na.rm = TRUE) +
  ggplot2::facet_wrap(~rdf_metric, nrow = 1, scales = "free_y") +
  ggplot2::scale_color_manual(values = drug_palette[ful_order], drop = FALSE) +
  ggplot2::scale_linetype_manual(values = drug_linetype[ful_order], drop = FALSE) +
  ggplot2::coord_cartesian(xlim = c(rdf_plot_rmin, rdf_xmax)) +
  ggplot2::labs(
    title = "Segment-resolved RDF profiles: Ful versus Ful-C18",
    subtitle = paste0("Rolling mean = ", rdf_smooth_window_nm, "nm; Ful/Ful-C18 are plotted separately because Ful-C18 has a high drug–drug peak"),
    x = "r (nm)",
    y = "g(r)",
    color = NULL,
    linetype = NULL
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(linewidth = 1.4, linetype = unname(drug_linetype[ful_order]))
    ),
    linetype = "none"
  ) +
  rdf_theme

# -----------------------------
# 11. Export every individual subfigure
# -----------------------------
save_plot_all(p_rmsd,          "multi_molecule_RMSD",             width = 7.2,  height = 4.4)
save_plot_all(p_metric_heatmap,"multi_molecule_metric_heatmap",    width = 9.2,  height = 5.8)
save_plot_all(p_hbond,         "multi_molecule_Hbond_summary",            width = 7.2,  height = 5.0)
save_plot_all(p_struct,        "multi_molecule_structural_energy_summary", width = 10.2, height = 5.4)
save_plot_all(p_rdf_taxane,    "multi_molecule_RDF_PTX_P2",        width = 12.0, height = 3.8)
save_plot_all(p_rdf_ful,       "multi_molecule_RDF_Ful_FulC18",           width = 12.0, height = 3.8)

# -----------------------------
# 12. Assemble and export the combined figure
# -----------------------------
top_row <- p_rmsd | p_hbond
combined_plot <- (
  top_row /
    p_metric_heatmap /
    p_struct /
    p_rdf_taxane /
    p_rdf_ful
) +
  patchwork::plot_layout(
    heights = c(1.0, 1.12, 0.92, 0.70, 0.70)
  ) +
  patchwork::plot_annotation(
    tag_levels = "a",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(face = "bold", size = 10),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
  )

combined_pdf <- file.path(fig_dir, "multi_molecule_MD_support.pdf")
combined_tiff <- file.path(fig_dir, "multi_molecule_MD_support.tiff")
combined_png <- file.path(fig_dir, "multi_molecule_MD_support.png")

ggplot2::ggsave(
  filename = combined_pdf,
  plot = combined_plot,
  width = combined_width_in,
  height = combined_height_in,
  units = "in",
  device = pdf_device
)

ggplot2::ggsave(
  filename = combined_tiff,
  plot = combined_plot,
  width = combined_width_in,
  height = combined_height_in,
  units = "in",
  dpi = panel_dpi_tiff,
  compression = "lzw",
  bg = "white"
)

ggplot2::ggsave(
  filename = combined_png,
  plot = combined_plot,
  width = combined_width_in,
  height = combined_height_in,
  units = "in",
  dpi = panel_dpi_png,
  bg = "white"
)

# -----------------------------
# 13. Output log
# -----------------------------
message("Done.")
message("Code folder:   ", normalizePath(script_dir, mustWork = FALSE))
message("Output folder: ", normalizePath(out_dir, mustWork = FALSE))
message("Combined PDF:  ", normalizePath(combined_pdf, mustWork = FALSE))
message("Combined TIFF: ", normalizePath(combined_tiff, mustWork = FALSE))
message("Combined PNG:  ", normalizePath(combined_png, mustWork = FALSE))
message("Individual panels: ", normalizePath(panel_dir, mustWork = FALSE))
message("CSV data files:    ", normalizePath(data_dir, mustWork = FALSE))
