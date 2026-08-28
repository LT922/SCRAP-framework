# training-library DLS/PDI analysis
# Input: data/training/training_size_pdi.csv

rm(list = ls())
graphics.off()
set.seed(123)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# -----------------------------
# 1. Input and output
# -----------------------------
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grepl(file_arg, cmd_args)])
  
  if (length(script_path) > 0 && nzchar(script_path[1])) {
    return(dirname(normalizePath(script_path[1], winslash = "/", mustWork = TRUE)))
  }
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    context_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(context_path)) {
      return(dirname(normalizePath(context_path, winslash = "/", mustWork = TRUE)))
    }
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

input_file <- file.path(project_root, "data", "training", "training_size_pdi.csv")
figure_dir <- file.path(project_root, "results", "figures", "training", "training_library")
table_dir <- file.path(project_root, "results", "tables", "training", "training_library")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2. Read data and summarize
# -----------------------------
raw_df <- readr::read_csv(input_file, show_col_types = FALSE, na = c("", "NA", "N/A", "na", "NaN", "-"))

if (!"Mol"%in% names(raw_df)) {
  stop("The input file must contain a column named 'Mol'.")
}

size_cols <- grep("^size(\\.|$)", names(raw_df), value = TRUE, ignore.case = TRUE)
pdi_cols  <- grep("^pdi(\\.|$)",  names(raw_df), value = TRUE, ignore.case = TRUE)

if (length(size_cols) == 0 || length(pdi_cols) == 0) {
  stop("The input file must contain size replicate columns and pdi replicate columns.")
}


raw_df <- raw_df %>%
  mutate(
    across(all_of(size_cols), ~ readr::parse_number(as.character(.x))),
    across(all_of(pdi_cols),  ~ readr::parse_number(as.character(.x)))
  )

mol_order <- raw_df$Mol

summary_df <- raw_df %>%
  mutate(
    Mol = as.character(Mol),
    x_id = row_number(),
    size_mean = rowMeans(across(all_of(size_cols)), na.rm = TRUE),
    pdi_mean  = rowMeans(across(all_of(pdi_cols)),  na.rm = TRUE),
    size_sd = apply(select(., all_of(size_cols)), 1, sd, na.rm = TRUE),
    pdi_sd  = apply(select(., all_of(pdi_cols)),  1, sd, na.rm = TRUE),
    size_mean = if_else(is.nan(size_mean), NA_real_, size_mean),
    pdi_mean  = if_else(is.nan(pdi_mean),  NA_real_, pdi_mean),
    size_sd   = if_else(is.na(size_mean), NA_real_, size_sd),
    pdi_sd    = if_else(is.na(pdi_mean),  NA_real_, pdi_sd),
    Outcome = case_when(
      is.na(size_mean) ~ "Precipitated",
      size_mean <= 30 ~ "Carrier-scale",
      TRUE ~ "Oversized"
    ),
    Outcome = factor(Outcome, levels = c("Carrier-scale", "Oversized", "Precipitated")),
    Compatibility = if_else(Outcome == "Carrier-scale", "Compatible", "Incompatible")
  )

readr::write_csv(summary_df, file.path(table_dir, "training_size_pdi_summary.csv"))

# -----------------------------
# 3. Plot settings
# -----------------------------
carrier_col <- "#3E6F9F"
oversize_col <- "#C8735C"
precip_col <- "#D99582"
compatible_label_col <- "#82A3C1"
incompatible_label_col <- "#D8A090"

outcome_cols <- c(
  "Carrier-scale"= carrier_col,
  "Oversized"= oversize_col,
  "Precipitated"= precip_col
)

outcome_shapes <- c(
  "Carrier-scale"= 21,
  "Oversized"= 21,
  "Precipitated"= 24
)

# Plot missing sizes as precipitated triangles above the plotted size range.
max_size <- suppressWarnings(max(summary_df$size_mean, 30, na.rm = TRUE))
if (!is.finite(max_size)) max_size <- 30
y_top <- max(170, max_size * 1.25)
precip_y <- y_top * 0.86

plot_df <- summary_df %>%
  mutate(
    y_plot = if_else(Outcome == "Precipitated", precip_y, size_mean),
    y_start = if_else(Outcome == "Precipitated", 30, 0),
    size_label = if_else(Outcome == "Precipitated", "Precipitated", sprintf("%.1f", size_mean)),
    pdi_label = case_when(
      Outcome == "Precipitated"& is.na(pdi_mean) ~ "NA",
      is.na(pdi_mean) ~ "NA",
      TRUE ~ sprintf("%.2f", pdi_mean)
    )
  )

last_compatible <- suppressWarnings(max(plot_df$x_id[plot_df$Outcome == "Carrier-scale"], na.rm = TRUE))
if (!is.finite(last_compatible)) last_compatible <- 0
first_x <- min(plot_df$x_id)
last_x <- max(plot_df$x_id)
compatible_mid <- if (last_compatible > 0) mean(c(first_x, last_compatible)) else first_x
incompatible_mid <- if (last_compatible < last_x) mean(c(last_compatible + 1, last_x)) else last_x

# -----------------------------
# 4. Size lollipop plot
# -----------------------------
p_size <- ggplot(plot_df, aes(x = x_id)) +
  annotate("rect",
           xmin = 0.5, xmax = last_compatible + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "#EEF3F8", alpha = 0.75) +
  annotate("rect",
           xmin = last_compatible + 0.5, xmax = last_x + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "#FBF1ED", alpha = 0.45) +
  geom_hline(yintercept = 30, linetype = "dashed", linewidth = 0.45, colour = "grey35") +
  annotate("text", x = 0.75, y = 34, label = "30 nm threshold",
           hjust = 0, vjust = 0, size = 3.3, colour = "grey30") +
  geom_vline(xintercept = last_compatible + 0.5,
             linetype = "dotted", linewidth = 0.55, colour = "grey25") +
  geom_segment(aes(y = y_start, yend = y_plot, colour = Outcome),
               linewidth = 0.75, alpha = 0.55) +
  geom_errorbar(
    data = filter(plot_df, Outcome != "Precipitated"),
    aes(ymin = pmax(size_mean - size_sd, 0), ymax = size_mean + size_sd, colour = Outcome),
    width = 0.14, linewidth = 0.45, alpha = 0.80
  ) +
  geom_point(aes(y = y_plot, fill = Outcome, shape = Outcome),
             size = 3.0, colour = "black", stroke = 0.45) +
  annotate("text", x = compatible_mid, y = y_top * 0.97, label = "Compatible",
           colour = compatible_label_col, fontface = "bold", size = 5.2) +
  annotate("text", x = compatible_mid, y = y_top * 0.86, label = "Sub-30 nm",
           colour = "grey15", size = 4.2) +
  annotate("text", x = incompatible_mid, y = y_top * 0.97, label = "Incompatible",
           colour = incompatible_label_col, fontface = "bold", size = 5.2) +
  scale_fill_manual(values = outcome_cols, name = NULL) +
  scale_colour_manual(values = outcome_cols, name = NULL) +
  scale_shape_manual(values = outcome_shapes, name = NULL) +
  scale_x_continuous(
    limits = c(0.5, last_x + 0.5),
    breaks = plot_df$x_id,
    labels = plot_df$Mol,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, y_top),
    breaks = c(10, 30, 50, 100, 150),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(x = NULL, y = "Hydrodynamic diameter (nm)") +
  guides(
    colour = "none",
    fill = guide_legend(override.aes = list(shape = c(21, 21, 24), colour = "black", size = 3.5)),
    shape = "none"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.y = element_text(size = 11),
    plot.margin = margin(5.5, 8, 0, 8)
  )

# -----------------------------
# 5. PDI heatmap
# -----------------------------
pdi_max <- suppressWarnings(max(plot_df$pdi_mean, na.rm = TRUE))
if (!is.finite(pdi_max)) pdi_max <- 0.30

p_pdi <- ggplot(plot_df, aes(x = x_id, y = "PDI", fill = pdi_mean)) +
  geom_tile(colour = "white", linewidth = 0.4, height = 0.85) +
  geom_text(aes(label = pdi_label), size = 3.1, colour = "grey10") +
  scale_fill_gradientn(
    colours = c("#DCEAF5", "#F7F7F7", "#D37A61"),
    limits = c(0, max(0.30, pdi_max)),
    na.value = "#E6E6E6",
    name = "PDI"
  ) +
  scale_x_continuous(
    limits = c(0.5, last_x + 0.5),
    breaks = plot_df$x_id,
    labels = plot_df$Mol,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, colour = "grey10"),
    axis.text.y = element_text(size = 11, colour = "grey10"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(0, 8, 5.5, 8)
  )

# -----------------------------
# 6. Combine and export
# -----------------------------
final_plot <- p_size / p_pdi + plot_layout(heights = c(3.5, 0.9))

pdf_file <- file.path(figure_dir, "training_size_pdi.pdf")
png_file <- file.path(figure_dir, "training_size_pdi.png")

ggsave(pdf_file, final_plot, width = 12.5, height = 7.2, units = "in", device = cairo_pdf)
ggsave(png_file, final_plot, width = 12.5, height = 7.2, units = "in", dpi = 600, bg = "white")

message("Done. Figure outputs: ", normalizePath(figure_dir, winslash = "/", mustWork = FALSE))
message(" - ", pdf_file)
message(" - ", png_file)
message(" - ", file.path(table_dir, "training_size_pdi_summary.csv"))
