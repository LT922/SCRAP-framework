rm(list = ls())
graphics.off()
set.seed(123)


# ============================================================
# Molecular-dynamics analysis of DTX-SI-C18/Abi-SI-C18 co-assembly
# Input: data/md/DTX_Abi_coassembly_MD/<system>/ana/*.xvg
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr", "purrr",
  "patchwork", "scales", "tibble"
)
missing_packages <- packages[!packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}
invisible(lapply(packages, library, character.only = TRUE))

# -----------------------------
# 0b. Cross-platform graphics device and font settings
# -----------------------------
# Use the generic R font family rather than a named font such as Arial.
# This avoids Windows/RStudio PDF/TIFF errors such as:
#   grid.Call.graphics(...): font category error
plot_base_family <- "sans"

safe_pdf_device <- function(filename, width, height, ...) {
  grDevices::pdf(
    file = filename,
    width = width,
    height = height,
    family = plot_base_family,
    useDingbats = FALSE,
    ...
  )
}

# -----------------------------
# 1. Paths
# -----------------------------
get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && !is.null(ctx$path) && nzchar(ctx$path)) {
      return(dirname(normalizePath(ctx$path, winslash = "/", mustWork = FALSE)))
    }
  }
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grepl(file_arg, cmd_args)])
  if (length(script_path) > 0 && nzchar(script_path[1])) {
    return(dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

root_dir <- file.path(project_root, "data", "md", "DTX_Abi_coassembly_MD")
diagnostic_dir <- file.path(project_root, "results", "figures", "md", "DTX_Abi_coassembly_MD", "analysis")
out_dir <- file.path(project_root, "results", "tables", "md", "DTX_Abi_coassembly_MD", "analysis")
fig_dir <- diagnostic_dir
data_dir <- out_dir
converted_dir <- file.path(out_dir, "converted_xvg_csv")
summary_dir <- file.path(out_dir, "summary_tables")

for (d in c(diagnostic_dir, out_dir, converted_dir, summary_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(root_dir)) {
  stop("Simulation root folder not found: ", root_dir)
}
root_dir <- normalizePath(root_dir, winslash = "/", mustWork = TRUE)

systems <- tibble::tribble(
  ~system_dir,              ~system_label,                              ~system_type,
  "pclabisic18",            "Single Abi-SI-C18@PCL",                   "single Abi-SI-C18",
  "pcldtxsic18",            "Single DTX-SI-C18@PCL",                   "single DTX-SI-C18",
  "pcldtxsic18abisic18",    "Mixed DTX-SI-C18/Abi-SI-C18@PCL",         "mixed co-assembly"
) %>%
  mutate(
    ana_dir = file.path(root_dir, system_dir, "ana"),
    exists = file.exists(ana_dir)
  )

if (!all(systems$exists)) {
  warning("Some ana folders were not found:\n", paste(systems$ana_dir[!systems$exists], collapse = "\n"))
}

message("Root folder: ", root_dir)
message("Output folder: ", out_dir)

# -----------------------------
# 2. General plotting helpers
# -----------------------------
theme_nature <- function(base_size = 8, base_family = plot_base_family) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.30, colour = "black"),
      axis.ticks.length = unit(1.6, "mm"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size + 1),
      plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.size = unit(3.5, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      panel.spacing = unit(1.2, "lines")
    )
}

# Backward-compatible alias used by existing SASA plotting blocks.
nature_theme <- theme_nature


# -----------------------------
# 2b. Publication color palettes
# -----------------------------
# System-level palette: used when comparing the three simulation systems.
system_colors <- c(
  "Single Abi-SI-C18@PCL"= "#79B8AE",        # teal
  "Single DTX-SI-C18@PCL"= "#5B8FCF",        # blue
  "Mixed DTX-SI-C18/Abi-SI-C18@PCL"= "#D97968"# coral
)

# Mechanistic-role palette: used consistently across COM distance, contacts,
# H-bonds and SASA contact-area panels.
domain_colors <- c(
  "PCL core"= "#4C8EA0",
  "PEG/corona"= "#E0B678",
  "whole polymer"= "#B58AAE",
  "prodrug-prodrug"= "#D76F61",
  "other"= "#93A8B8"
)

# SASA panel uses longer class names; keep colors matched to domain_colors.
sasa_domain_colors <- c(
  "PCL-core contact area"= domain_colors[["PCL core"]],
  "PEG/corona contact area"= domain_colors[["PEG/corona"]],
  "DTX-Abi contact area"= domain_colors[["prodrug-prodrug"]],
  "whole-polymer contact area"= domain_colors[["whole polymer"]]
)

# Heatmap palette from the uploaded reference color strip: coral -> cream -> teal.
heatmap_colors <- c("#CF7467", "#F1D8C8", "#7DBDB7")
scale_fill_heatmap <- function(name = NULL) {
  ggplot2::scale_fill_gradientn(
    colours = heatmap_colors,
    values = scales::rescale(c(0, 0.5, 1)),
    name = name,
    na.value = "#F6F6F6"
  )
}

interaction_class <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "PCL") ~ "PCL core",
    stringr::str_detect(x, "PEG") ~ "PEG/corona",
    stringr::str_detect(x, "polymer") ~ "whole polymer",
    stringr::str_detect(x, "DTX-Abi|drug-drug|par-par|lig-lig|par-lig|lig-par") ~ "prodrug-prodrug",
    TRUE ~ "other"
  )
}

sasa_class <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "PCL") ~ "PCL-core contact area",
    stringr::str_detect(x, "PEG") ~ "PEG/corona contact area",
    stringr::str_detect(x, "DTX-Abi") ~ "DTX-Abi contact area",
    stringr::str_detect(x, "polymer") ~ "whole-polymer contact area",
    TRUE ~ "whole-polymer contact area"
  )
}

component_class <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "PCL") ~ "PCL core",
    stringr::str_detect(x, "PEG") ~ "PEG/corona",
    stringr::str_detect(x, "polymer|assembly") ~ "whole polymer",
    stringr::str_detect(x, "DTX|Abi|drug|prodrug|parent|ligand") ~ "prodrug-prodrug",
    TRUE ~ "other"
  )
}

metric_palette <- function(keys) {
  keys <- unique(as.character(keys))
  base <- c(
    system_colors,
    domain_colors,
    "single-system RMSD"= "#D97968",
    "polymer + total drug"= "#B58AAE",
    "total drug"= "#6EA8A1",
    "DTX-SI-C18"= "#5B8FCF",
    "Abi-SI-C18"= "#79B8AE",
    "single drug/prodrug"= "#D97968",
    "polymer"= "#B58AAE",
    "whole non-water assembly"= "#93A8B8",
    "PEG block"= "#E0B678",
    "PCL block"= "#4C8EA0",
    "DTX-PCL"= domain_colors[["PCL core"]],
    "Abi-PCL"= "#70AFC2",
    "total drug-PCL"= "#2E6F83",
    "DTX-PEG"= domain_colors[["PEG/corona"]],
    "Abi-PEG"= "#EDC891",
    "total drug-PEG"= "#C99A55",
    "DTX-Abi"= domain_colors[["prodrug-prodrug"]],
    "DTX-polymer"= domain_colors[["whole polymer"]],
    "Abi-polymer"= "#C79BC0",
    "total drug-polymer"= "#9D6F9A",
    "drug-polymer"= domain_colors[["whole polymer"]],
    "drug-PCL"= domain_colors[["PCL core"]],
    "drug-PEG"= domain_colors[["PEG/corona"]],
    "drug-drug"= domain_colors[["prodrug-prodrug"]],
    "polymer-polymer"= "#8FA6B3"
  )
  missing <- setdiff(keys, names(base))
  if (length(missing) > 0) {
    extra <- scales::hue_pal(l = 55, c = 80)(length(missing))
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[keys]
}

save_pdf <- function(plot, filename, width = 7.2, height = 4.8, dpi = 600) {
  # Export each figure as PDF plus a high-resolution raster copy.
  # The raster format is TIFF by default; if the local graphics device fails,
  # a PNG copy is written instead so every PDF has a paired raster file.
  pdf_path <- file.path(fig_dir, filename)
  stem <- tools::file_path_sans_ext(filename)
  tif_path <- file.path(fig_dir, paste0(stem, ".tif"))
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  
  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = safe_pdf_device,
    bg = "white",
    limitsize = FALSE
  )
  
  ok_tif <- tryCatch({
    ggplot2::ggsave(
      filename = tif_path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = "tiff",
      compression = "lzw",
      bg = "white",
      limitsize = FALSE
    )
    TRUE
  }, error = function(e) {
    warning("TIFF export failed for ", filename, ": ", conditionMessage(e), ". Writing PNG instead.")
    FALSE
  })
  
  if (!ok_tif) {
    ggplot2::ggsave(
      filename = png_path,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = "png",
      bg = "white",
      limitsize = FALSE
    )
  }
  invisible(pdf_path)
}

format_mean_label <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

add_bar_label_columns <- function(df, value_col = "mean", sd_col = "sd", digits = 2, pad_frac = 0.035) {
  if (nrow(df) == 0) return(df)
  value <- df[[value_col]]
  sdv <- if (sd_col %in% names(df)) df[[sd_col]] else rep(0, length(value))
  sdv[is.na(sdv) | !is.finite(sdv)] <- 0
  ymax <- suppressWarnings(max(value + sdv, na.rm = TRUE))
  pad <- ifelse(is.finite(ymax) && ymax > 0, ymax * pad_frac, 0.05)
  df %>%
    mutate(
      mean_label = format_mean_label(.data[[value_col]], digits = digits),
      label_y = .data[[value_col]] + sdv + pad
    )
}

bar_y_expand <- function() {
  ggplot2::expansion(mult = c(0.02, 0.22))
}

clean_label <- function(x) {
  x %>%
    stringr::str_replace_all("\\\\s", "_") %>%
    stringr::str_replace_all("\\\\N", "") %>%
    stringr::str_replace_all("[^A-Za-z0-9_+\\-()./ ]", "") %>%
    stringr::str_squish()
}

safe_filename <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

# -----------------------------
# 3. XVG parser
# -----------------------------
extract_quoted_value <- function(lines, pattern) {
  hit <- lines[stringr::str_detect(lines, pattern)]
  if (length(hit) == 0) return(NA_character_)
  val <- stringr::str_match(hit[1], '"(.*)"')[, 2]
  ifelse(is.na(val), NA_character_, val)
}

parse_xvg_metadata <- function(file) {
  lines <- readLines(file, warn = FALSE)
  title <- extract_quoted_value(lines, "^@\\s+title")
  x_label <- extract_quoted_value(lines, "^@\\s+xaxis\\s+label")
  y_label <- extract_quoted_value(lines, "^@\\s+yaxis\\s+label")
  
  legend_lines <- lines[stringr::str_detect(lines, "^@\\s+s[0-9]+\\s+legend")]
  legends <- character(0)
  if (length(legend_lines) > 0) {
    legend_tbl <- stringr::str_match(legend_lines, "^@\\s+s([0-9]+)\\s+legend\\s+\"(.*)\"")
    legends <- legend_tbl[, 3]
    legends <- legends[!is.na(legends)]
    legends <- clean_label(legends)
  }
  
  list(title = title, x_label = x_label, y_label = y_label, legends = legends)
}

read_xvg_long <- function(file, system_dir, system_label) {
  meta <- parse_xvg_metadata(file)
  lines <- readLines(file, warn = FALSE)
  numeric_lines <- lines[!stringr::str_detect(lines, "^\\s*[#@]") & nzchar(stringr::str_squish(lines))]
  if (length(numeric_lines) == 0) return(NULL)
  
  dat <- tryCatch(
    read.table(
      text = paste(numeric_lines, collapse = "\n"),
      header = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE,
      comment.char = ""
    ),
    error = function(e) NULL
  )
  if (is.null(dat) || ncol(dat) < 2) {
    warning("Could not parse numeric data: ", file)
    return(NULL)
  }
  
  # Keep only x and the first y column. Convert after subsetting and remove
  # non-numeric rows robustly. This prevents occasional malformed/filled rows
  # from entering plots as huge negative or non-physical values.
  dat <- dat[, 1:2, drop = FALSE]
  dat <- as.data.frame(lapply(dat, function(z) suppressWarnings(as.numeric(z))))
  dat <- dat[is.finite(dat[[1]]) & is.finite(dat[[2]]), , drop = FALSE]
  if (nrow(dat) == 0) {
    warning("No finite numeric x/y values after parsing: ", file)
    return(NULL)
  }
  n_y <- 1
  
  legends <- meta$legends
  if (length(legends) < n_y || is.na(legends[1]) || !nzchar(legends[1])) {
    legends <- "Y1"
  }
  legends <- legends[seq_len(n_y)]
  y_names <- make.names(legends, unique = TRUE)
  colnames(dat) <- c("x", y_names)
  
  x_label <- ifelse(is.na(meta$x_label), "", meta$x_label)
  y_label <- ifelse(is.na(meta$y_label), "", meta$y_label)
  file_stem <- stringr::str_remove(basename(file), "\\.xvg$")
  
  is_time <- stringr::str_detect(stringr::str_to_lower(x_label), "time") ||
    stringr::str_detect(stringr::str_to_lower(file_stem),
                        "rmsd|gyrate|dist|hbnum|sasa|area|temperature|pressure|density|potential|ntr|nclid")
  is_ps <- stringr::str_detect(stringr::str_to_lower(x_label), "ps") &&
    !stringr::str_detect(stringr::str_to_lower(x_label), "ns")
  is_ns <- stringr::str_detect(stringr::str_to_lower(x_label), "ns")
  
  x_converted <- dat$x
  x_unit <- "raw"
  if (is_time) {
    if (is_ps) {
      x_converted <- dat$x / 1000
      x_unit <- "ns"
    } else if (is_ns) {
      x_converted <- dat$x
      x_unit <- "ns"
    } else if (max(dat$x, na.rm = TRUE) > 1000) {
      x_converted <- dat$x / 1000
      x_unit <- "ns_inferred_from_ps"
    } else {
      x_converted <- dat$x
      x_unit <- "ns_assumed"
    }
  }
  
  dat_long <- dat %>%
    mutate(x_raw = x, x_value = x_converted) %>%
    select(-x) %>%
    tidyr::pivot_longer(
      cols = -c(x_raw, x_value),
      names_to = "series",
      values_to = "value"
    ) %>%
    mutate(
      series_index = match(series, y_names),
      series = stringr::str_replace_all(series, "\\.", ""),
      series = stringr::str_squish(series),
      file = basename(file),
      file_stem = file_stem,
      system_dir = system_dir,
      system = system_label,
      x_label = x_label,
      y_label = y_label,
      x_unit = x_unit,
      title = ifelse(is.na(meta$title), "", meta$title)
    ) %>%
    relocate(system_dir, system, file, file_stem, title, x_label, y_label, x_unit, series_index, series)
  
  dat_long
}

# -----------------------------
# 4. Convert all .xvg files to CSV
# -----------------------------
xvg_index <- systems %>%
  filter(exists) %>%
  mutate(xvg_files = purrr::map(ana_dir, ~list.files(.x, pattern = "\\.xvg$", full.names = TRUE))) %>%
  select(system_dir, system_label, system_type, ana_dir, xvg_files) %>%
  tidyr::unnest(xvg_files) %>%
  mutate(
    file = basename(xvg_files),
    file_stem = stringr::str_remove(file, "\\.xvg$"),
    converted_csv = file.path(converted_dir, system_dir, paste0(file_stem, ".csv"))
  )

if (nrow(xvg_index) == 0) stop("No .xvg files were found under the ana folders.")

all_xvg_list <- list()
metadata_rows <- list()

for (i in seq_len(nrow(xvg_index))) {
  row <- xvg_index[i, ]
  dir.create(dirname(row$converted_csv), recursive = TRUE, showWarnings = FALSE)
  dat_long <- read_xvg_long(row$xvg_files, row$system_dir, row$system_label)
  
  if (!is.null(dat_long)) {
    readr::write_csv(dat_long, row$converted_csv)
    all_xvg_list[[length(all_xvg_list) + 1]] <- dat_long
    
    metadata_rows[[length(metadata_rows) + 1]] <- dat_long %>%
      summarise(
        system_dir = first(system_dir),
        system = first(system),
        file = first(file),
        file_stem = first(file_stem),
        title = first(title),
        x_label = first(x_label),
        y_label = first(y_label),
        x_unit = first(x_unit),
        n_points = n(),
        n_series = n_distinct(series_index),
        .groups = "drop"
      )
  }
}

all_xvg <- dplyr::bind_rows(all_xvg_list)
xvg_metadata <- dplyr::bind_rows(metadata_rows)

readr::write_csv(xvg_index, file.path(summary_dir, "xvg_file_index.csv"))
readr::write_csv(xvg_metadata, file.path(summary_dir, "xvg_metadata.csv"))
readr::write_csv(all_xvg, file.path(summary_dir, "all_xvg_long.csv"))

xvg_value_diagnostics <- all_xvg %>%
  group_by(system_dir, system, file_stem, file, y_label) %>%
  summarise(
    n = sum(is.finite(value)),
    min_value = suppressWarnings(min(value, na.rm = TRUE)),
    max_value = suppressWarnings(max(value, na.rm = TRUE)),
    min_x = suppressWarnings(min(x_value, na.rm = TRUE)),
    max_x = suppressWarnings(max(x_value, na.rm = TRUE)),
    .groups = "drop"
  )
readr::write_csv(xvg_value_diagnostics, file.path(summary_dir, "xvg_value_diagnostics.csv"))

message("Converted XVG files: ", dplyr::n_distinct(paste(all_xvg$system_dir, all_xvg$file)))

# -----------------------------
# 5. Data selection helpers
# -----------------------------
get_xvg <- function(system_dir = NULL, file_stem = NULL, file_pattern = NULL, first_series = TRUE) {
  df <- all_xvg
  if (!is.null(system_dir)) df <- df %>% filter(.data$system_dir %in% system_dir)
  if (!is.null(file_stem)) df <- df %>% filter(.data$file_stem %in% file_stem)
  if (!is.null(file_pattern)) df <- df %>% filter(stringr::str_detect(.data$file_stem, file_pattern))
  if (nrow(df) == 0) return(tibble())
  if (first_series) df <- df %>% group_by(system_dir, system, file_stem) %>% filter(series_index == min(series_index, na.rm = TRUE)) %>% ungroup()
  df
}

summarise_final_window <- function(df, window_ns = 50) {
  if (nrow(df) == 0) return(tibble())
  df %>%
    filter(is.finite(x_value), is.finite(value)) %>%
    group_by(system_dir, system, file_stem, series, y_label) %>%
    mutate(
      max_x = max(x_value, na.rm = TRUE),
      min_keep = ifelse(max_x > window_ns, max_x - window_ns, quantile(x_value, 0.75, na.rm = TRUE))
    ) %>%
    filter(x_value >= min_keep) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      n = sum(!is.na(value)),
      final_window_start_ns = first(min_keep),
      final_window_end_ns = first(max_x),
      .groups = "drop"
    )
}

add_metric_label <- function(df, dictionary) {
  df %>%
    left_join(dictionary, by = "file_stem") %>%
    mutate(metric = ifelse(is.na(metric), file_stem, metric))
}

# factor() requires unique levels. Some file stems intentionally map to the
# same displayed metric, e.g. polymergyrate and gyrate_polymer both map to
# "polymer". Use this helper whenever dictionary-derived levels are used.
metric_levels <- function(dictionary) {
  unique(dictionary$metric)
}

# Use dictionary-only joins for plotting. This prevents unrelated XVG files or
# unmatched file stems from entering a plot as metric = NA. Structural metrics
# such as RMSD, Rg, distances, H-bonds, contacts and RDF are non-negative; any
# negative values are parsing artifacts or unrelated data and are removed here.
prepare_metric_df <- function(file_stems, dictionary, system_dir = NULL,
                              first_series = TRUE, nonnegative = TRUE,
                              max_value = NULL) {
  df <- get_xvg(system_dir = system_dir, file_stem = file_stems, first_series = first_series)
  if (nrow(df) == 0) return(tibble())
  
  dict_clean <- dictionary %>% distinct(file_stem, .keep_all = TRUE)
  
  df <- df %>%
    inner_join(dict_clean, by = "file_stem") %>%
    filter(is.finite(x_value), is.finite(value), !is.na(metric))
  
  if (nonnegative) {
    df <- df %>% filter(value >= 0)
  }
  if (!is.null(max_value)) {
    df <- df %>% filter(value <= max_value)
  }
  
  df %>% mutate(metric = factor(metric, levels = metric_levels(dictionary)))
}

# -----------------------------
# 6. Metric dictionaries
# -----------------------------
rmsd_dict <- tibble::tribble(
  ~file_stem,              ~metric,
  "rmsd",                  "single-system RMSD",
  "rmsd_polymer_drug",     "polymer + total drug",
  "rmsd_drug",             "total drug",
  "rmsd_dtxsic18",         "DTX-SI-C18",
  "rmsd_abisic18",         "Abi-SI-C18"
)

rg_dict <- tibble::tribble(
  ~file_stem,                ~metric,
  "allgyrate",               "whole non-water assembly",
  "druggyrate",              "single drug/prodrug",
  "polymergyrate",           "polymer",
  "gyrate_polymer_drug",     "polymer + total drug",
  "gyrate_polymer",          "polymer",
  "gyrate_total_drug",       "total drug",
  "gyrate_dtxsic18",         "DTX-SI-C18",
  "gyrate_abisic18",         "Abi-SI-C18",
  "gyrate_peg",              "PEG block",
  "gyrate_pcl",              "PCL block"
)

whole_distance_dict <- tibble::tribble(
  ~file_stem,                   ~metric,
  "drug-polymer-dist",          "total drug-polymer",
  "dtx-polymer-dist",           "DTX-polymer",
  "abi-polymer-dist",           "Abi-polymer",
  "dtx-abi-dist",               "DTX-Abi",
  "drug-pcl-dist",              "total drug-PCL",
  "dtx-pcl-dist",               "DTX-PCL",
  "abi-pcl-dist",               "Abi-PCL",
  "drug-peg-dist",              "total drug-PEG",
  "dtx-peg-dist",               "DTX-PEG",
  "abi-peg-dist",               "Abi-PEG"
)

hbond_dict <- tibble::tribble(
  ~file_stem,                     ~metric,
  "hbnum-dru-po",                 "drug-polymer",
  "hbnum-dru-pcl",                "drug-PCL",
  "hbnum-dru-peg",                "drug-PEG",
  "hbnum-drug-polymer",           "total drug-polymer",
  "hbnum-dtx-polymer",            "DTX-polymer",
  "hbnum-abi-polymer",            "Abi-polymer",
  "hbnum-dtx-abi",                "DTX-Abi",
  "hbnum-dtx-pcl",                "DTX-PCL",
  "hbnum-abi-pcl",                "Abi-PCL",
  "hbnum-dtx-peg",                "DTX-PEG",
  "hbnum-abi-peg",                "Abi-PEG",
  "hbnum-pol-pol",                "polymer-polymer",
  "hbnum-dru-dru",                "drug-drug"
)

contact_dict <- tibble::tribble(
  ~file_stem,                    ~metric,
  "dtx-abi-contacts",            "DTX-Abi",
  "dtx-polymer-contacts",        "DTX-polymer",
  "abi-polymer-contacts",        "Abi-polymer",
  "dtx-pcl-contacts",            "DTX-PCL",
  "abi-pcl-contacts",            "Abi-PCL",
  "dtx-peg-contacts",            "DTX-PEG",
  "abi-peg-contacts",            "Abi-PEG",
  "dtxlig-pcl-contacts",         "DTX lig-PCL",
  "dtxpar-pcl-contacts",         "DTX par-PCL",
  "abilig-pcl-contacts",         "Abi lig-PCL",
  "abipar-pcl-contacts",         "Abi par-PCL"
)

rdf_dict <- tibble::tribble(
  ~file_stem,          ~metric,
  "rdf-dtx-pcl",       "DTX-PCL",
  "rdf-abi-pcl",       "Abi-PCL",
  "rdf-dtx-peg",       "DTX-PEG",
  "rdf-abi-peg",       "Abi-PEG",
  "rdf-dtx-abi",       "DTX-Abi",
  "rdf-drug-pcl",      "total drug-PCL",
  "rdf-drug-peg",      "total drug-PEG"
)

segment_distance_dict <- tibble::tribble(
  ~file_stem,                 ~metric,             ~prodrug, ~segment_pair,
  "dtxpar-polymer-dist",      "DTX par-polymer",   "DTX",    "par-polymer",
  "dtxlig-polymer-dist",      "DTX lig-polymer",   "DTX",    "lig-polymer",
  "dtxpar-pcl-dist",          "DTX par-PCL",       "DTX",    "par-PCL",
  "dtxlig-pcl-dist",          "DTX lig-PCL",       "DTX",    "lig-PCL",
  "dtxpar-peg-dist",          "DTX par-PEG",       "DTX",    "par-PEG",
  "dtxlig-peg-dist",          "DTX lig-PEG",       "DTX",    "lig-PEG",
  "abipar-polymer-dist",      "Abi par-polymer",   "Abi",    "par-polymer",
  "abilig-polymer-dist",      "Abi lig-polymer",   "Abi",    "lig-polymer",
  "abipar-pcl-dist",          "Abi par-PCL",       "Abi",    "par-PCL",
  "abilig-pcl-dist",          "Abi lig-PCL",       "Abi",    "lig-PCL",
  "abipar-peg-dist",          "Abi par-PEG",       "Abi",    "par-PEG",
  "abilig-peg-dist",          "Abi lig-PEG",       "Abi",    "lig-PEG",
  "dtxpar-abipar-dist",       "DTX par-Abi par",   "DTX-Abi", "par-par",
  "dtxlig-abilig-dist",       "DTX lig-Abi lig",   "DTX-Abi", "lig-lig",
  "dtxpar-abilig-dist",       "DTX par-Abi lig",   "DTX-Abi", "par-lig",
  "dtxlig-abipar-dist",       "DTX lig-Abi par",   "DTX-Abi", "lig-par"
)
segment_hbond_dict <- segment_distance_dict %>% mutate(file_stem = paste0("hbnum-", stringr::str_remove(file_stem, "-dist$")))

# -----------------------------
# 7. RMSD
# -----------------------------
rmsd_df <- prepare_metric_df(rmsd_dict$file_stem, rmsd_dict, nonnegative = TRUE, max_value = 1000)

if (nrow(rmsd_df) > 0) {
  p_rmsd <- ggplot(rmsd_df, aes(x = x_value, y = value, colour = metric)) +
    geom_line(linewidth = 0.35, alpha = 0.95) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    labs(title = "Structural convergence of single-prodrug and mixed co-assembly systems", x = "Time (ns)", y = "RMSD (nm)") +
    scale_colour_manual(values = metric_palette(unique(as.character(rmsd_df$metric))), drop = FALSE) +
    theme_nature() + theme(legend.position = "right")
  save_pdf(p_rmsd, "MD_RMSD_time_series.pdf", width = 7.6, height = 6.0)
  readr::write_csv(rmsd_df, file.path(summary_dir, "MD_RMSD_time_series_data.csv"))
}

# -----------------------------
# 8. Radius of gyration
# -----------------------------
rg_df <- prepare_metric_df(rg_dict$file_stem, rg_dict, nonnegative = TRUE, max_value = 1000)

if (nrow(rg_df) > 0) {
  p_rg_time <- ggplot(rg_df, aes(x = x_value, y = value, colour = metric)) +
    geom_line(linewidth = 0.35, alpha = 0.95) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    labs(title = "Assembly compactness assessed by radius of gyration", x = "Time (ns)", y = "Rg (nm)") +
    scale_colour_manual(values = metric_palette(unique(as.character(rg_df$metric))), drop = FALSE) +
    theme_nature() + theme(legend.position = "right")
  save_pdf(p_rg_time, "MD_Rg_time_series.pdf", width = 7.6, height = 6.0)
  
  rg_summary <- summarise_final_window(rg_df) %>% add_metric_label(rg_dict) %>%
    mutate(metric = factor(metric, levels = metric_levels(rg_dict))) %>%
    add_bar_label_columns(digits = 2)
  p_rg_bar <- ggplot(rg_summary, aes(x = metric, y = mean, fill = system)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.72), width = 0.18, linewidth = 0.25) +
    geom_text(aes(y = label_y, label = mean_label), position = position_dodge(width = 0.72), hjust = 0, size = 2.4) +
    scale_fill_manual(values = system_colors, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Final-window Rg summary", x = NULL, y = "Rg (nm), mean +/- SD") +
    theme_nature() + theme(legend.position = "top")
  save_pdf(p_rg_bar, "MD_Rg_final_window.pdf", width = 7.4, height = 4.8)
  readr::write_csv(rg_df, file.path(summary_dir, "MD_Rg_time_series_data.csv"))
  readr::write_csv(rg_summary, file.path(summary_dir, "MD_Rg_final_window_summary.csv"))
}

# -----------------------------
# 9. COM distances
# -----------------------------
dist_df <- prepare_metric_df(whole_distance_dict$file_stem, whole_distance_dict, nonnegative = TRUE, max_value = 1000)

if (nrow(dist_df) > 0) {
  p_dist_time <- ggplot(dist_df, aes(x = x_value, y = value, colour = metric)) +
    geom_line(linewidth = 0.35, alpha = 0.95) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    labs(title = "COM distances define carrier integration and DTX-Abi proximity", x = "Time (ns)", y = "COM distance (nm)") +
    scale_colour_manual(values = metric_palette(unique(as.character(dist_df$metric))), drop = FALSE) +
    theme_nature() + theme(legend.position = "right")
  save_pdf(p_dist_time, "MD_COM_distance_time_series.pdf", width = 7.6, height = 6.0)
  
  dist_summary <- summarise_final_window(dist_df) %>% add_metric_label(whole_distance_dict) %>%
    mutate(metric = factor(metric, levels = metric_levels(whole_distance_dict))) %>%
    add_bar_label_columns(digits = 2)
  p_dist_bar <- ggplot(dist_summary, aes(x = metric, y = mean, fill = system)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.72), width = 0.18, linewidth = 0.25) +
    geom_text(aes(y = label_y, label = mean_label), position = position_dodge(width = 0.72), hjust = 0, size = 2.4) +
    scale_fill_manual(values = system_colors, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Final-window COM distance summary", x = NULL, y = "Distance (nm), mean +/- SD") +
    theme_nature() + theme(legend.position = "top")
  save_pdf(p_dist_bar, "MD_COM_distance_final_window.pdf", width = 7.5, height = 5.0)
  readr::write_csv(dist_df, file.path(summary_dir, "MD_COM_distance_time_series_data.csv"))
  readr::write_csv(dist_summary, file.path(summary_dir, "MD_COM_distance_final_window_summary.csv"))
}

# -----------------------------
# 10. Hydrogen bonds
# -----------------------------
hbond_df <- prepare_metric_df(hbond_dict$file_stem, hbond_dict, nonnegative = TRUE, max_value = 1e6)

if (nrow(hbond_df) > 0) {
  p_hb_time <- ggplot(hbond_df, aes(x = x_value, y = value, colour = metric)) +
    geom_line(linewidth = 0.35, alpha = 0.90) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    labs(title = "Segment-specific hydrogen-bond contribution during co-assembly", x = "Time (ns)", y = "Number of H-bonds") +
    scale_colour_manual(values = metric_palette(unique(as.character(hbond_df$metric))), drop = FALSE) +
    theme_nature() + theme(legend.position = "right")
  save_pdf(p_hb_time, "MD_Hbond_time_series.pdf", width = 7.6, height = 6.0)
  
  hbond_summary <- summarise_final_window(hbond_df) %>% add_metric_label(hbond_dict) %>%
    mutate(metric = factor(metric, levels = metric_levels(hbond_dict))) %>%
    add_bar_label_columns(digits = 2)
  p_hb_bar <- ggplot(hbond_summary, aes(x = metric, y = mean, fill = system)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), position = position_dodge(width = 0.72), width = 0.18, linewidth = 0.25) +
    geom_text(aes(y = label_y, label = mean_label), position = position_dodge(width = 0.72), hjust = 0, size = 2.4) +
    scale_fill_manual(values = system_colors, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Final-window H-bond summary", x = NULL, y = "H-bonds, mean +/- SD") +
    theme_nature() + theme(legend.position = "top")
  save_pdf(p_hb_bar, "MD_Hbond_final_window.pdf", width = 7.5, height = 5.2)
  readr::write_csv(hbond_df, file.path(summary_dir, "MD_Hbond_time_series_data.csv"))
  readr::write_csv(hbond_summary, file.path(summary_dir, "MD_Hbond_final_window_summary.csv"))
}

# -----------------------------
# 11. Contact numbers from gmx mindist -on
# -----------------------------
contact_df <- prepare_metric_df(contact_dict$file_stem, contact_dict, nonnegative = TRUE, max_value = 1e8)

if (nrow(contact_df) > 0) {
  contact_summary <- summarise_final_window(contact_df) %>% add_metric_label(contact_dict) %>%
    mutate(metric = factor(metric, levels = metric_levels(contact_dict))) %>%
    add_bar_label_columns(digits = 2) %>%
    mutate(class = interaction_class(as.character(metric)))
  p_contact <- ggplot(contact_summary, aes(x = metric, y = mean, fill = class)) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.18, linewidth = 0.25) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 2.4) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Near-contact numbers in the mixed DTX/Abi co-assembly", subtitle = "Contacts are based on the gmx mindist cutoff used during analysis, typically 0.6 nm.", x = NULL, y = "Contact number, mean +/- SD", fill = NULL) +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    theme_nature() + theme(legend.position = "none")
  save_pdf(p_contact, "MD_contact_number_final_window.pdf", width = 6.8, height = 4.8)
  readr::write_csv(contact_df, file.path(summary_dir, "MD_contact_number_time_series_data.csv"))
  readr::write_csv(contact_summary, file.path(summary_dir, "MD_contact_number_final_window_summary.csv"))
}

# -----------------------------
# 12. RDF curves and RDF peak summary
# -----------------------------
rdf_df <- prepare_metric_df(rdf_dict$file_stem, rdf_dict, nonnegative = TRUE, max_value = 1e8)

if (nrow(rdf_df) > 0) {
  p_rdf <- ggplot(rdf_df, aes(x = x_raw, y = value, colour = metric)) +
    geom_line(linewidth = 0.45, alpha = 0.95) +
    facet_wrap(~metric, scales = "free_y", ncol = 3) +
    labs(title = "RDF analysis resolves PCL anchoring, PEG exposure and DTX-Abi association", x = "r (nm)", y = "g(r)") +
    scale_colour_manual(values = metric_palette(unique(as.character(rdf_df$metric))), drop = FALSE) +
    theme_nature() + theme(legend.position = "none")
  save_pdf(p_rdf, "MD_RDF_mixed_core_shell_and_DTX_Abi.pdf", width = 7.5, height = 5.4)
  
  rdf_peaks <- rdf_df %>% group_by(system_dir, system, file_stem, metric) %>%
    summarise(peak_g = max(value, na.rm = TRUE), peak_r_nm = x_raw[which.max(value)][1], .groups = "drop")
  readr::write_csv(rdf_df, file.path(summary_dir, "MD_RDF_curve_data.csv"))
  readr::write_csv(rdf_peaks, file.path(summary_dir, "MD_RDF_peak_summary.csv"))
}

# -----------------------------
# 13. SASA and SASA-derived contact area
# -----------------------------
# SASA contact area is calculated directly from matched XVG files:
# The global XVG table is convenient for most readouts, but SASA contact area is
# sensitive to exact frame alignment. Previous versions could mix duplicated or
# unit-converted x values, producing many zeros or nearly identical mean values.
# Here SASA is read directly from exact .xvg files in each ana folder and aligned
# by frame order. For single-prodrug systems, the available area_* files are used:
#   contact area = (area_dru + area_pol - area_whole) / 2
# where area_whole is the SASA of the non-water drug+polymer assembly.
# For the mixed system, the sasa-* files generated for A, B and A+B are used.
# -----------------------------

xvg_path_from_stem <- function(system_dir, file_stem) {
  ana_dir <- systems$ana_dir[match(system_dir, systems$system_dir)]
  if (length(ana_dir) == 0 || is.na(ana_dir) || !dir.exists(ana_dir)) return(NA_character_)
  f <- file.path(ana_dir, paste0(file_stem, ".xvg"))
  if (file.exists(f)) normalizePath(f, winslash = "/", mustWork = TRUE) else NA_character_
}

read_xvg_first_xy_direct <- function(file) {
  if (is.na(file) || !file.exists(file)) return(tibble())
  meta <- parse_xvg_metadata(file)
  lines <- readLines(file, warn = FALSE)
  numeric_lines <- lines[!stringr::str_detect(lines, "^\\s*[#@]") & nzchar(stringr::str_squish(lines))]
  if (length(numeric_lines) == 0) return(tibble())
  
  dat <- tryCatch(
    read.table(
      text = paste(numeric_lines, collapse = "\n"),
      header = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE,
      comment.char = ""
    ),
    error = function(e) NULL
  )
  if (is.null(dat) || ncol(dat) < 2) return(tibble())
  
  dat <- dat[, 1:2, drop = FALSE]
  dat <- as.data.frame(lapply(dat, function(z) suppressWarnings(as.numeric(z))))
  dat <- dat[is.finite(dat[[1]]) & is.finite(dat[[2]]), , drop = FALSE]
  if (nrow(dat) == 0) return(tibble())
  colnames(dat) <- c("x_raw", "value")
  
  x_label <- ifelse(is.na(meta$x_label), "", meta$x_label)
  x_low <- stringr::str_to_lower(x_label)
  is_ps <- stringr::str_detect(x_low, "ps") && !stringr::str_detect(x_low, "ns")
  is_ns <- stringr::str_detect(x_low, "ns")
  
  # SASA files should be time series. Use the x-axis label when present; otherwise
  # keep raw values. This avoids falsely converting RDF-like x axes or other data.
  if (is_ps) {
    dat$x_value <- dat$x_raw / 1000
    x_unit <- "ns_from_ps"
  } else if (is_ns) {
    dat$x_value <- dat$x_raw
    x_unit <- "ns"
  } else {
    dat$x_value <- dat$x_raw
    x_unit <- "raw"
  }
  
  tibble(
    source_file = basename(file),
    file_stem = stringr::str_remove(basename(file), "\\.xvg$"),
    title = ifelse(is.na(meta$title), "", meta$title),
    x_label = x_label,
    y_label = ifelse(is.na(meta$y_label), "", meta$y_label),
    x_unit = x_unit,
    frame = seq_len(nrow(dat)),
    x_raw = dat$x_raw,
    x_value = dat$x_value,
    value = dat$value
  )
}

read_sasa_direct <- function(system_dir, file_stem, component_label = file_stem) {
  file <- xvg_path_from_stem(system_dir, file_stem)
  dat <- read_xvg_first_xy_direct(file)
  if (nrow(dat) == 0) return(tibble())
  dat %>%
    mutate(
      system_dir = system_dir,
      system = systems$system_label[match(system_dir, systems$system_dir)],
      component = component_label,
      file_stem_requested = file_stem
    ) %>%
    relocate(system_dir, system, component, file_stem_requested)
}

align_sasa_three_by_frame <- function(a, b, ab) {
  n_keep <- min(nrow(a), nrow(b), nrow(ab))
  if (!is.finite(n_keep) || n_keep <= 0) return(tibble())
  a <- a[seq_len(n_keep), , drop = FALSE]
  b <- b[seq_len(n_keep), , drop = FALSE]
  ab <- ab[seq_len(n_keep), , drop = FALSE]
  
  tibble(
    frame = seq_len(n_keep),
    x_value = a$x_value,
    x_raw = a$x_raw,
    value_a = a$value,
    value_b = b$value,
    value_ab = ab$value,
    file_a = first(a$file_stem),
    file_b = first(b$file_stem),
    file_ab = first(ab$file_stem),
    x_unit_a = first(a$x_unit),
    x_unit_b = first(b$x_unit),
    x_unit_ab = first(ab$x_unit),
    n_a = nrow(a),
    n_b = nrow(b),
    n_ab = nrow(ab)
  )
}

calculate_contact_area_direct <- function(system_dir, pair_label, file_a, file_b, file_ab) {
  a <- read_sasa_direct(system_dir, file_a, file_a)
  b <- read_sasa_direct(system_dir, file_b, file_b)
  ab <- read_sasa_direct(system_dir, file_ab, file_ab)
  if (nrow(a) == 0 || nrow(b) == 0 || nrow(ab) == 0) {
    message("Skipping SASA contact area for ", system_dir, "/ ", pair_label,
            "because one or more files are missing: ", file_a, ", ", file_b, ", ", file_ab)
    return(tibble())
  }
  
  aligned <- align_sasa_three_by_frame(a, b, ab)
  if (nrow(aligned) == 0) return(tibble())
  
  aligned %>%
    mutate(
      system_dir = system_dir,
      system = systems$system_label[match(system_dir, systems$system_dir)],
      pair = pair_label,
      contact_area_nm2_raw = (value_a + value_b - value_ab) / 2,
      # Do not hide large negative values; only remove tiny floating-point noise.
      contact_area_nm2 = ifelse(contact_area_nm2_raw < 0 & abs(contact_area_nm2_raw) < 1e-6, 0, contact_area_nm2_raw),
      negative_contact_area = contact_area_nm2 < 0
    ) %>%
    relocate(system_dir, system, pair, frame, x_value, contact_area_nm2, contact_area_nm2_raw)
}

# Direct raw SASA export. These curves help diagnose whether the contact area is
# meaningful before using it as a figure panel.
raw_sasa_defs <- tibble::tribble(
  ~system_dir,             ~component,                    ~file_stem,
  "pclabisic18",           "Abi prodrug",                 "area_dru",
  "pclabisic18",           "Abi parent",                  "area_par",
  "pclabisic18",           "Abi ligand",                  "area_lig",
  "pclabisic18",           "polymer",                     "area_pol",
  "pclabisic18",           "PEG block",                   "area_peg",
  "pclabisic18",           "drug+polymer assembly",       "area_whole",
  "pcldtxsic18",           "DTX prodrug",                 "area_dru",
  "pcldtxsic18",           "DTX parent",                  "area_par",
  "pcldtxsic18",           "DTX ligand",                  "area_lig",
  "pcldtxsic18",           "polymer",                     "area_pol",
  "pcldtxsic18",           "PEG block",                   "area_peg",
  "pcldtxsic18",           "drug+polymer assembly",       "area_whole",
  "pcldtxsic18abisic18",   "total prodrugs",              "sasa-drug",
  "pcldtxsic18abisic18",   "DTX-SI-C18",                  "sasa-dtx",
  "pcldtxsic18abisic18",   "Abi-SI-C18",                  "sasa-abi",
  "pcldtxsic18abisic18",   "polymer",                     "sasa-polymer",
  "pcldtxsic18abisic18",   "PCL block",                   "sasa-pcl",
  "pcldtxsic18abisic18",   "PEG block",                   "sasa-peg",
  "pcldtxsic18abisic18",   "total prodrugs+polymer",      "sasa-drug-polymer"
)

raw_sasa <- purrr::pmap_dfr(raw_sasa_defs, function(system_dir, component, file_stem) {
  read_sasa_direct(system_dir, file_stem, component)
}) %>%
  filter(is.finite(x_value), is.finite(value), value >= 0)

if (nrow(raw_sasa) > 0) {
  raw_sasa_diag <- raw_sasa %>%
    group_by(system_dir, system, component, file_stem_requested, source_file, x_label, y_label, x_unit) %>%
    summarise(
      n_points = n(),
      x_min = min(x_value, na.rm = TRUE),
      x_max = max(x_value, na.rm = TRUE),
      median_dx = suppressWarnings(median(diff(sort(unique(x_value))), na.rm = TRUE)),
      min_value = min(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE),
      mean_value = mean(value, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_csv(raw_sasa, file.path(summary_dir, "MD_raw_SASA_time_series_data.csv"))
  readr::write_csv(raw_sasa_diag, file.path(summary_dir, "MD_raw_SASA_file_axis_diagnostic.csv"))
  
  p_raw_sasa_time <- raw_sasa %>%
    filter(system_dir %in% c("pclabisic18", "pcldtxsic18", "pcldtxsic18abisic18")) %>%
    ggplot(aes(x = x_value, y = value, colour = component)) +
    geom_line(linewidth = 0.28, alpha = 0.85) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    labs(title = "Raw solvent-accessible surface area", x = "Time (ns)", y = "SASA (nm²)", colour = NULL) +
    scale_colour_manual(values = metric_palette(unique(as.character(raw_sasa$component))), drop = FALSE) +
    nature_theme() +
    theme(legend.position = "right")
  save_pdf(p_raw_sasa_time, "MD_raw_SASA_time_series.pdf", width = 8.8, height = 7.2)
  
  raw_sasa_summary <- raw_sasa %>%
    group_by(system_dir, system, component) %>%
    mutate(
      max_x = max(x_value, na.rm = TRUE),
      min_keep = ifelse(max_x > 50, max_x - 50, quantile(x_value, 0.75, na.rm = TRUE))
    ) %>%
    filter(x_value >= min_keep) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      n = sum(!is.na(value)),
      final_window_start_ns = first(min_keep),
      final_window_end_ns = first(max_x),
      .groups = "drop"
    ) %>%
    mutate(
      sd = ifelse(is.na(sd), 0, sd),
      label_y = mean + sd + 0.03 * max(mean + sd, na.rm = TRUE),
      mean_label = sprintf("%.1f", mean)
    )
  readr::write_csv(raw_sasa_summary, file.path(summary_dir, "MD_raw_SASA_final_window_summary.csv"))
  
  p_raw_sasa_bar <- raw_sasa_summary %>%
    ggplot(aes(x = reorder(component, mean), y = mean, fill = component_class(component))) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.20) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.22) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 2.4) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = bar_y_expand()) +
    labs(title = "Raw SASA final-window summary", x = NULL, y = "SASA (nm²), mean ± SD", fill = NULL) +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    nature_theme() +
    theme(legend.position = "none")
  save_pdf(p_raw_sasa_bar, "MD_raw_SASA_final_window.pdf", width = 8.4, height = 7.2)
}

# Contact-area definitions. Single systems use area_dru/area_pol/area_whole. The
# mixed system uses explicit A, B and A+B sasa-* files.
sasa_pairs <- tibble::tribble(
  ~system_dir,                ~pair,                         ~file_a,       ~file_b,          ~file_ab,
  "pclabisic18",              "Abi-polymer",                  "area_dru",   "area_pol",      "area_whole",
  "pcldtxsic18",              "DTX-polymer",                  "area_dru",   "area_pol",      "area_whole",
  "pcldtxsic18abisic18",      "total drug-polymer",           "sasa-drug",  "sasa-polymer",  "sasa-drug-polymer",
  "pcldtxsic18abisic18",      "DTX-polymer",                  "sasa-dtx",   "sasa-polymer",  "sasa-dtx-polymer",
  "pcldtxsic18abisic18",      "Abi-polymer",                  "sasa-abi",   "sasa-polymer",  "sasa-abi-polymer",
  "pcldtxsic18abisic18",      "DTX-Abi",                      "sasa-dtx",   "sasa-abi",      "sasa-dtx-abi",
  "pcldtxsic18abisic18",      "DTX-PCL",                      "sasa-dtx",   "sasa-pcl",      "sasa-dtx-pcl",
  "pcldtxsic18abisic18",      "Abi-PCL",                      "sasa-abi",   "sasa-pcl",      "sasa-abi-pcl",
  "pcldtxsic18abisic18",      "DTX-PEG",                      "sasa-dtx",   "sasa-peg",      "sasa-dtx-peg",
  "pcldtxsic18abisic18",      "Abi-PEG",                      "sasa-abi",   "sasa-peg",      "sasa-abi-peg"
)

sasa_contact <- purrr::pmap_dfr(sasa_pairs, ~calculate_contact_area_direct(..1, ..2, ..3, ..4, ..5))

if (nrow(sasa_contact) > 0) {
  readr::write_csv(sasa_contact, file.path(summary_dir, "MD_SASA_contact_area_time_series_data.csv"))
  
  sasa_contact_diag <- sasa_contact %>%
    group_by(system_dir, system, pair, file_a, file_b, file_ab, x_unit_a, x_unit_b, x_unit_ab) %>%
    summarise(
      n_points = n(),
      x_min = min(x_value, na.rm = TRUE),
      x_max = max(x_value, na.rm = TRUE),
      median_dx = suppressWarnings(median(diff(sort(unique(x_value))), na.rm = TRUE)),
      min_contact_area = min(contact_area_nm2, na.rm = TRUE),
      max_contact_area = max(contact_area_nm2, na.rm = TRUE),
      mean_contact_area = mean(contact_area_nm2, na.rm = TRUE),
      n_negative = sum(negative_contact_area, na.rm = TRUE),
      .groups = "drop"
    )
  readr::write_csv(sasa_contact_diag, file.path(summary_dir, "MD_SASA_contact_area_axis_diagnostic.csv"))
  
  # For statistical plotting, do not silently use non-physical negative areas.
  sasa_contact_plot_df <- sasa_contact %>%
    filter(is.finite(contact_area_nm2), contact_area_nm2 >= 0)
  
  p_sasa_time <- ggplot(sasa_contact_plot_df, aes(x = x_value, y = contact_area_nm2, colour = pair)) +
    geom_line(linewidth = 0.30, alpha = 0.85) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    labs(title = "SASA-derived interfacial contact area", x = "Time (ns)", y = "Contact area (nm²)", colour = NULL) +
    scale_colour_manual(values = metric_palette(unique(as.character(sasa_contact_plot_df$pair))), drop = FALSE) +
    nature_theme() +
    theme(legend.position = "right")
  save_pdf(p_sasa_time, "MD_SASA_contact_area_time_series.pdf", width = 8.8, height = 7.0)
  
  sasa_summary <- sasa_contact_plot_df %>%
    group_by(system_dir, system, pair) %>%
    mutate(
      max_x = max(x_value, na.rm = TRUE),
      min_keep = ifelse(max_x > 50, max_x - 50, quantile(x_value, 0.75, na.rm = TRUE))
    ) %>%
    filter(x_value >= min_keep) %>%
    summarise(
      mean = mean(contact_area_nm2, na.rm = TRUE),
      sd = sd(contact_area_nm2, na.rm = TRUE),
      n = sum(!is.na(contact_area_nm2)),
      final_window_start_ns = first(min_keep),
      final_window_end_ns = first(max_x),
      .groups = "drop"
    ) %>%
    mutate(
      sd = ifelse(is.na(sd), 0, sd),
      pair = factor(pair, levels = unique(sasa_pairs$pair)),
      label_y = mean + sd + 0.03 * max(mean + sd, na.rm = TRUE),
      mean_label = sprintf("%.1f", mean)
    )
  readr::write_csv(sasa_summary, file.path(summary_dir, "MD_SASA_contact_area_final_window_summary.csv"))
  
  p_sasa_bar <- sasa_summary %>%
    ggplot(aes(x = pair, y = mean, fill = sasa_class(as.character(pair)))) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.20) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.22) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 2.4) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = bar_y_expand()) +
    labs(title = "SASA-derived contact area", subtitle = "Calculated directly from exact XVG files; final-window mean ± SD.", x = NULL, y = "Contact area (nm²)", fill = NULL) +
    scale_fill_manual(values = sasa_domain_colors, drop = FALSE) +
    nature_theme() +
    theme(legend.position = "none")
  save_pdf(p_sasa_bar, "MD_SASA_contact_area_final_window.pdf", width = 8.4, height = 7.0)
}

# 14. Segment-specific COM-distance heatmap
# -----------------------------
segment_dist_df <- get_xvg(system_dir = "pcldtxsic18abisic18", file_stem = segment_distance_dict$file_stem, first_series = TRUE) %>%
  inner_join(segment_distance_dict, by = "file_stem") %>%
  filter(is.finite(value), value >= 0, value <= 1000)

if (nrow(segment_dist_df) > 0) {
  segment_dist_summary <- summarise_final_window(segment_dist_df) %>%
    select(system_dir, system, file_stem, mean, sd, n) %>%
    inner_join(segment_distance_dict, by = "file_stem") %>%
    mutate(segment_pair = factor(segment_pair, levels = unique(segment_distance_dict$segment_pair)))
  
  # Add whole-prodrug COM distances to the segment-resolved map.
  # These rows provide the missing "drug整体"reference, allowing direct comparison of
  # whole-prodrug localization with parent- and ligand-level localization.
  whole_drug_distance_dict <- tibble::tribble(
    ~file_stem,            ~metric,              ~prodrug,  ~segment_pair,
    "dtx-polymer-dist",    "DTX drug-polymer",   "DTX",     "drug-polymer",
    "dtx-pcl-dist",        "DTX drug-PCL",       "DTX",     "drug-PCL",
    "dtx-peg-dist",        "DTX drug-PEG",       "DTX",     "drug-PEG",
    "abi-polymer-dist",    "Abi drug-polymer",   "Abi",     "drug-polymer",
    "abi-pcl-dist",        "Abi drug-PCL",       "Abi",     "drug-PCL",
    "abi-peg-dist",        "Abi drug-PEG",       "Abi",     "drug-PEG",
    "dtx-abi-dist",        "DTX-Abi drug-drug",  "DTX-Abi", "drug-drug"
  )
  
  whole_drug_dist_df <- get_xvg(
    system_dir = "pcldtxsic18abisic18",
    file_stem = whole_drug_distance_dict$file_stem,
    first_series = TRUE
  ) %>%
    inner_join(whole_drug_distance_dict, by = "file_stem") %>%
    filter(is.finite(value), value >= 0, value <= 1000)
  
  if (nrow(whole_drug_dist_df) > 0) {
    whole_drug_dist_summary <- summarise_final_window(whole_drug_dist_df) %>%
      select(system_dir, system, file_stem, mean, sd, n) %>%
      inner_join(whole_drug_distance_dict, by = "file_stem")
    
    segment_dist_summary <- bind_rows(segment_dist_summary, whole_drug_dist_summary) %>%
      mutate(
        segment_pair = factor(
          segment_pair,
          levels = c(
            "drug-polymer", "par-polymer", "lig-polymer",
            "drug-PCL", "par-PCL", "lig-PCL",
            "drug-PEG", "par-PEG", "lig-PEG",
            "drug-drug", "par-par", "lig-lig", "par-lig", "lig-par"
          )
        )
      )
  }
  
  p_seg_dist <- ggplot(segment_dist_summary, aes(x = segment_pair, y = prodrug, fill = mean)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", mean)), size = 2.3) +
    scale_fill_heatmap("Distance\n(nm)") +
    labs(title = "Whole-prodrug and segment-specific COM distances in DTX/Abi mixed co-assembly", x = NULL, y = NULL) +
    theme_nature() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  save_pdf(p_seg_dist, "MD_whole_and_segment_COM_distance_heatmap.pdf", width = 9.4, height = 3.8)
  readr::write_csv(segment_dist_summary, file.path(summary_dir, "MD_whole_and_segment_COM_distance_summary.csv"))
}

# -----------------------------
# 15. Segment-specific H-bond heatmap
# -----------------------------
segment_hbond_df <- get_xvg(system_dir = "pcldtxsic18abisic18", file_stem = segment_hbond_dict$file_stem, first_series = TRUE) %>%
  inner_join(segment_hbond_dict, by = "file_stem") %>%
  filter(is.finite(value), value >= 0, value <= 1e6)

if (nrow(segment_hbond_df) > 0) {
  segment_hbond_summary <- summarise_final_window(segment_hbond_df) %>%
    select(system_dir, system, file_stem, mean, sd, n) %>%
    inner_join(segment_hbond_dict, by = "file_stem") %>%
    mutate(segment_pair = factor(segment_pair, levels = unique(segment_hbond_dict$segment_pair)))
  
  # Add whole-prodrug H-bond descriptors to keep the H-bond heatmap
  # consistent with the whole-drug + segment-resolved COM-distance heatmap.
  whole_drug_hbond_dict <- tibble::tribble(
    ~file_stem,             ~metric,              ~prodrug,  ~segment_pair,
    "hbnum-dtx-polymer",    "DTX drug-polymer",   "DTX",     "drug-polymer",
    "hbnum-dtx-pcl",        "DTX drug-PCL",       "DTX",     "drug-PCL",
    "hbnum-dtx-peg",        "DTX drug-PEG",       "DTX",     "drug-PEG",
    "hbnum-abi-polymer",    "Abi drug-polymer",   "Abi",     "drug-polymer",
    "hbnum-abi-pcl",        "Abi drug-PCL",       "Abi",     "drug-PCL",
    "hbnum-abi-peg",        "Abi drug-PEG",       "Abi",     "drug-PEG",
    "hbnum-dtx-abi",        "DTX-Abi drug-drug",  "DTX-Abi", "drug-drug"
  )
  
  whole_drug_hbond_df <- get_xvg(
    system_dir = "pcldtxsic18abisic18",
    file_stem = whole_drug_hbond_dict$file_stem,
    first_series = TRUE
  ) %>%
    inner_join(whole_drug_hbond_dict, by = "file_stem") %>%
    filter(is.finite(value), value >= 0, value <= 1e6)
  
  if (nrow(whole_drug_hbond_df) > 0) {
    whole_drug_hbond_summary <- summarise_final_window(whole_drug_hbond_df) %>%
      select(system_dir, system, file_stem, mean, sd, n) %>%
      inner_join(whole_drug_hbond_dict, by = "file_stem")
    
    segment_hbond_summary <- bind_rows(segment_hbond_summary, whole_drug_hbond_summary) %>%
      mutate(
        segment_pair = factor(
          segment_pair,
          levels = c(
            "drug-polymer", "par-polymer", "lig-polymer",
            "drug-PCL", "par-PCL", "lig-PCL",
            "drug-PEG", "par-PEG", "lig-PEG",
            "drug-drug", "par-par", "lig-lig", "par-lig", "lig-par"
          )
        )
      )
  }
  
  p_seg_hb <- ggplot(segment_hbond_summary, aes(x = segment_pair, y = prodrug, fill = mean)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", mean)), size = 2.3) +
    scale_fill_heatmap("H-bonds") +
    labs(title = "Whole-prodrug and segment-specific hydrogen bonding in DTX/Abi mixed co-assembly", x = NULL, y = NULL) +
    theme_nature() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  save_pdf(p_seg_hb, "MD_whole_and_segment_Hbond_heatmap.pdf", width = 9.4, height = 3.8)
  readr::write_csv(segment_hbond_summary, file.path(summary_dir, "MD_whole_and_segment_Hbond_summary.csv"))
}

# -----------------------------
# 16. Cluster-size output, if available
# -----------------------------
cluster_dict <- tibble::tribble(
  ~file_stem,     ~metric,
  "size",         "single-system clusters",
  "size_drug",    "total drug clusters",
  "size_dtx",     "DTX clusters",
  "size_abi",     "Abi clusters"
)
cluster_df <- prepare_metric_df(cluster_dict$file_stem, cluster_dict, nonnegative = TRUE, max_value = 1e8)

if (nrow(cluster_df) > 0) {
  p_cluster <- ggplot(cluster_df, aes(x = x_raw, y = value, colour = metric)) +
    geom_line(linewidth = 0.35, alpha = 0.95) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    labs(title = "Cluster-size output from gmx cluster", x = "Raw x value from size.xvg", y = "Cluster size / population") +
    scale_colour_manual(values = metric_palette(unique(as.character(cluster_df$metric))), drop = FALSE) +
    theme_nature() + theme(legend.position = "right")
  save_pdf(p_cluster, "MD_cluster_size_distribution.pdf", width = 7.6, height = 5.6)
  readr::write_csv(cluster_df, file.path(summary_dir, "MD_cluster_size_distribution_data.csv"))
}

# -----------------------------
# 17. MD quality control
# -----------------------------
qc_dict <- tibble::tribble(
  ~file_stem,       ~metric,
  "potential",      "Potential energy",
  "temperature",    "Temperature",
  "pressure",       "Pressure",
  "density",        "Density"
)
qc_df <- prepare_metric_df(qc_dict$file_stem, qc_dict, nonnegative = FALSE, max_value = NULL)

if (nrow(qc_df) > 0) {
  p_qc <- ggplot(qc_df, aes(x = x_value, y = value, colour = system)) +
    geom_line(linewidth = 0.30, alpha = 0.90) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    labs(title = "MD quality-control readouts", x = "Time (ns)", y = NULL) +
    scale_colour_manual(values = system_colors, drop = FALSE) +
    theme_nature() + theme(legend.position = "top")
  save_pdf(p_qc, "MD_quality_control.pdf", width = 7.2, height = 5.2)
  readr::write_csv(qc_df, file.path(summary_dir, "MD_quality_control_data.csv"))
}

# -----------------------------
# 18. Integrated wide multi-panel Figure
# -----------------------------
# The individual PDFs above are useful for checking each metric. This section
# builds one wide integrated Figure using the most manuscript-relevant panels:
#   A, structural convergence; B, assembly compactness; C, core/corona distances;
#   D, RDF; E, H-bonds; F, near-contact numbers; G, SASA contact area;
#   H-I, segment-resolved distance and H-bond maps.
# The figure is intentionally wide to fit a Figure-style layout.

compact_theme <- function(base_size = 7) {
  theme_nature(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 0.8, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.8),
      legend.text = element_text(size = base_size - 1.0),
      legend.key.size = unit(2.6, "mm"),
      strip.text = element_text(face = "bold", size = base_size - 0.2),
      panel.spacing = unit(0.8, "lines"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

panel_letter <- function(p, lab) {
  p +
    labs(tag = lab) +
    theme(
      plot.tag = element_text(face = "bold", size = 11),
      plot.tag.position = c(0.01, 0.99)
    )
}

integrated_panels <- list()

# A. Structural convergence: use only a few interpretable curves to avoid a crowded panel.
if (exists("rmsd_df") && nrow(rmsd_df) > 0) {
  rmsd_keep <- c("single-system RMSD", "polymer + total drug", "total drug", "DTX-SI-C18", "Abi-SI-C18")
  pA <- rmsd_df %>%
    filter(metric %in% rmsd_keep) %>%
    mutate(metric = factor(as.character(metric), levels = rmsd_keep)) %>%
    ggplot(aes(x = x_value, y = value, colour = metric)) +
    geom_line(linewidth = 0.25, alpha = 0.95) +
    facet_wrap(~system, scales = "free_y", ncol = 1) +
    labs(title = "Structural convergence", x = "Time (ns)", y = "RMSD (nm)") +
    compact_theme() +
    theme(legend.position = "bottom")
  integrated_panels[["A"]] <- panel_letter(pA, "a")
}

# B. Final-window Rg summary.
if (exists("rg_summary") && nrow(rg_summary) > 0) {
  rg_keep <- c("single drug/prodrug", "polymer + total drug", "total drug", "DTX-SI-C18", "Abi-SI-C18", "PEG block", "PCL block")
  pB <- rg_summary %>%
    filter(as.character(metric) %in% rg_keep) %>%
    mutate(metric = factor(as.character(metric), levels = rev(rg_keep))) %>%
    ggplot(aes(x = metric, y = mean, fill = system)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62, colour = "black", linewidth = 0.18) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
                  position = position_dodge(width = 0.72), width = 0.16, linewidth = 0.18) +
    geom_text(aes(y = label_y, label = mean_label), position = position_dodge(width = 0.72), hjust = 0, size = 1.8) +
    scale_fill_manual(values = system_colors, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Assembly compactness", x = NULL, y = "Rg (nm)") +
    compact_theme() +
    theme(legend.position = "bottom")
  integrated_panels[["B"]] <- panel_letter(pB, "b")
}

# C. Mixed-system COM distances: focus on PCL anchoring, PEG exposure and DTX-Abi proximity.
if (exists("dist_summary") && nrow(dist_summary) > 0) {
  dist_keep <- c("DTX-PCL", "Abi-PCL", "DTX-PEG", "Abi-PEG", "DTX-Abi", "DTX-polymer", "Abi-polymer")
  pC <- dist_summary %>%
    filter(system_dir == "pcldtxsic18abisic18", as.character(metric) %in% dist_keep) %>%
    mutate(metric = factor(as.character(metric), levels = rev(dist_keep)), class = interaction_class(as.character(metric))) %>%
    ggplot(aes(x = metric, y = mean, fill = class)) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.18) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.18) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 1.8) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Core/corona localization", x = NULL, y = "COM distance (nm)", fill = NULL) +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    compact_theme() +
    theme(legend.position = "none")
  integrated_panels[["C"]] <- panel_letter(pC, "c")
}

# D. RDF: compact facet panel for local enrichment patterns.
if (exists("rdf_df") && nrow(rdf_df) > 0) {
  rdf_keep <- c("DTX-PCL", "Abi-PCL", "DTX-PEG", "Abi-PEG", "DTX-Abi")
  pD <- rdf_df %>%
    filter(system_dir == "pcldtxsic18abisic18", as.character(metric) %in% rdf_keep) %>%
    mutate(metric = factor(as.character(metric), levels = rdf_keep)) %>%
    ggplot(aes(x = x_raw, y = value)) +
    geom_line(linewidth = 0.35, colour = "#4C8EA0", alpha = 0.95) +
    facet_wrap(~metric, scales = "free_y", ncol = 3) +
    labs(title = "RDF-resolved local interactions", x = "r (nm)", y = "g(r)") +
    compact_theme() +
    theme(legend.position = "none")
  integrated_panels[["D"]] <- panel_letter(pD, "d")
}

# E. Mixed-system H-bond summary.
if (exists("hbond_summary") && nrow(hbond_summary) > 0) {
  hb_keep <- c("DTX-polymer", "Abi-polymer", "DTX-Abi", "DTX-PCL", "Abi-PCL", "DTX-PEG", "Abi-PEG")
  pE <- hbond_summary %>%
    filter(system_dir == "pcldtxsic18abisic18", as.character(metric) %in% hb_keep) %>%
    mutate(metric = factor(as.character(metric), levels = rev(hb_keep)), class = interaction_class(as.character(metric))) %>%
    ggplot(aes(x = metric, y = mean, fill = class)) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.18) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.18) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 1.8) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Segment-specific hydrogen-bond contribution", x = NULL, y = "H-bonds", fill = NULL) +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    compact_theme() +
    theme(legend.position = "none")
  integrated_panels[["E"]] <- panel_letter(pE, "e")
}

# F. Mixed-system contact numbers.
if (exists("contact_summary") && nrow(contact_summary) > 0) {
  contact_keep <- c("DTX-PCL", "Abi-PCL", "DTX-PEG", "Abi-PEG", "DTX-Abi", "DTX lig-PCL", "Abi lig-PCL")
  pF <- contact_summary %>%
    filter(as.character(metric) %in% contact_keep) %>%
    mutate(metric = factor(as.character(metric), levels = rev(contact_keep)), class = interaction_class(as.character(metric))) %>%
    ggplot(aes(x = metric, y = mean, fill = class)) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.18) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.18) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 1.8) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Near-contact enrichment", x = NULL, y = "Contacts", fill = NULL) +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    compact_theme() +
    theme(legend.position = "none")
  integrated_panels[["F"]] <- panel_letter(pF, "f")
}

# G. SASA-derived contact area.
if (exists("sasa_summary") && nrow(sasa_summary) > 0) {
  sasa_keep <- c("DTX-PCL", "Abi-PCL", "DTX-PEG", "Abi-PEG", "DTX-Abi", "DTX-polymer", "Abi-polymer")
  pG <- sasa_summary %>%
    filter(as.character(pair) %in% sasa_keep) %>%
    mutate(pair = factor(as.character(pair), levels = rev(sasa_keep))) %>%
    ggplot(aes(x = pair, y = mean, fill = sasa_class(as.character(pair)))) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.18) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.18) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 1.8) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Interfacial contact area", x = NULL, y = "Area (nm²)", fill = NULL) +
    scale_fill_manual(values = sasa_domain_colors, drop = FALSE) +
    compact_theme() +
    theme(legend.position = "none")
  integrated_panels[["G"]] <- panel_letter(pG, "g")
}

# H. Segment-specific distance heatmap.
if (exists("segment_dist_summary") && nrow(segment_dist_summary) > 0) {
  pH <- segment_dist_summary %>%
    filter(segment_pair %in% c("par-polymer", "lig-polymer", "par-PCL", "lig-PCL", "par-PEG", "lig-PEG", "par-par", "lig-lig", "par-lig", "lig-par")) %>%
    ggplot(aes(x = segment_pair, y = prodrug, fill = mean)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", mean)), size = 1.8) +
    scale_fill_heatmap("nm") +
    labs(title = "Whole-prodrug and segment-resolved distance", x = NULL, y = NULL) +
    compact_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  integrated_panels[["H"]] <- panel_letter(pH, "h")
}

# I. Segment-specific H-bond heatmap.
if (exists("segment_hbond_summary") && nrow(segment_hbond_summary) > 0) {
  pI <- segment_hbond_summary %>%
    filter(segment_pair %in% c("par-polymer", "lig-polymer", "par-PCL", "lig-PCL", "par-PEG", "lig-PEG", "par-par", "lig-lig", "par-lig", "lig-par")) %>%
    ggplot(aes(x = segment_pair, y = prodrug, fill = mean)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", mean)), size = 1.8) +
    scale_fill_heatmap("H-bonds") +
    labs(title = "Whole-prodrug and segment-resolved H-bonds", x = NULL, y = NULL) +
    compact_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  integrated_panels[["I"]] <- panel_letter(pI, "i")
}

# Export one wide integrated figure. The design keeps the figure horizontal:
# row 1 emphasizes global assembly metrics; row 2 emphasizes local and segmental mechanisms.
if (length(integrated_panels) >= 2) {
  # Preserve biological narrative order even if some panels are missing.
  panel_order <- intersect(c("A", "B", "C", "D", "E", "F", "G", "H", "I"), names(integrated_panels))
  integrated_panels <- integrated_panels[panel_order]
  
  integrated_wide <- patchwork::wrap_plots(integrated_panels, ncol = 4, guides = "collect") +
    patchwork::plot_annotation(
      title = "Compatibility-focused molecular dynamics analysis of DTX-SI-C18/Abi-SI-C18 co-assembly",
      subtitle = "Final-window summaries emphasize PCL anchoring, PEG/corona exposure, controlled DTX-Abi association and segment-resolved interaction balance.",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, hjust = 0),
        plot.margin = margin(5, 5, 5, 5)
      )
    )
  
  save_pdf(integrated_wide, "MD_coassembly_integrated_wide.pdf", width = 18.5, height = 10.5)
  
  # Also export a slightly taller version that is easier to inspect on screen.
  save_pdf(integrated_wide, "MD_coassembly_integrated_wide_large.pdf", width = 20, height = 12)
}


# -----------------------------
# 19. Rev7 clear integrated Figure with explicit system and component labels
# -----------------------------
# Rationale:
# - The legacy overview overlays several systems and uses labels such as
#   "single drug/prodrug", which mean different molecules in different folders.
# - Rev7 therefore builds a clearer, manuscript-oriented wide Figure. Each panel
#   either facets by system or focuses only on the mixed co-assembly system.
# - Contact counts are reported as x10^3 atom contacts to avoid unreadable large
#   numbers. SASA-derived contact area is retained as a diagnostic output, but it
#   is not forced into the main integrated Figure if the summary appears nearly
#   degenerate or dominated by very large SD values.

system_order_clear <- c(
  "Single Abi-SI-C18@PCL",
  "Single DTX-SI-C18@PCL",
  "Mixed DTX-SI-C18/Abi-SI-C18@PCL"
)

system_colors_clear <- system_colors

system_short_label <- function(x) {
  dplyr::case_when(
    x == "Single Abi-SI-C18@PCL"~ "Abi single",
    x == "Single DTX-SI-C18@PCL"~ "DTX single",
    x == "Mixed DTX-SI-C18/Abi-SI-C18@PCL"~ "DTX+Abi mixed",
    TRUE ~ x
  )
}

metric_clear_label <- function(system_dir, metric) {
  dplyr::case_when(
    system_dir == "pclabisic18"& metric == "single drug/prodrug"~ "Abi-SI-C18 prodrug",
    system_dir == "pcldtxsic18"& metric == "single drug/prodrug"~ "DTX-SI-C18 prodrug",
    system_dir == "pcldtxsic18abisic18"& metric == "total drug"~ "DTX+Abi total prodrugs",
    system_dir == "pcldtxsic18abisic18"& metric == "polymer + total drug"~ "polymer + DTX+Abi",
    metric == "whole non-water assembly"~ "whole non-water assembly",
    metric == "polymer"~ "polymer",
    metric == "DTX-SI-C18"~ "DTX-SI-C18 prodrug",
    metric == "Abi-SI-C18"~ "Abi-SI-C18 prodrug",
    metric == "PEG block"~ "PEG block",
    metric == "PCL block"~ "PCL block",
    TRUE ~ as.character(metric)
  )
}

mixed_metric_class <- function(metric) {
  dplyr::case_when(
    stringr::str_detect(metric, "PCL") ~ "PCL core",
    stringr::str_detect(metric, "PEG") ~ "PEG/corona",
    stringr::str_detect(metric, "polymer") ~ "whole polymer",
    stringr::str_detect(metric, "DTX-Abi") ~ "prodrug-prodrug",
    TRUE ~ "other"
  )
}

bar_label_columns_custom <- function(df, value_col = "mean", sd_col = "sd", digits = 2, pad_frac = 0.04) {
  if (nrow(df) == 0) return(df)
  value <- df[[value_col]]
  sdv <- if (sd_col %in% names(df)) df[[sd_col]] else rep(0, length(value))
  sdv[is.na(sdv) | !is.finite(sdv)] <- 0
  ymax <- suppressWarnings(max(value + sdv, na.rm = TRUE))
  pad <- ifelse(is.finite(ymax) && ymax > 0, ymax * pad_frac, 0.05)
  df %>%
    mutate(
      mean_label = format_mean_label(.data[[value_col]], digits = digits),
      label_y = .data[[value_col]] + sdv + pad
    )
}

# A small text panel is added so that readers immediately know what the three
# systems mean. This avoids repeated ambiguity in legends.
p_system_map <- ggplot() +
  annotate("text", x = 0, y = 1.00, hjust = 0, vjust = 1, size = 3.0, fontface = "bold",
           label = "System definitions") +
  annotate("text", x = 0, y = 0.78, hjust = 0, vjust = 1, size = 2.45,
           label = "Abi single: 1 Abi-SI-C18 + 1 mPEG-PCL\nDTX single: 1 DTX-SI-C18 + 1 mPEG-PCL\nDTX+Abi mixed: 10 DTX-SI-C18 + 10 Abi-SI-C18 + 30 mPEG-PCL") +
  annotate("text", x = 0, y = 0.33, hjust = 0, vjust = 1, size = 2.25,
           label = "Panels b-c compare final-window means.\nPanels d-h focus on the mixed co-assembly and separate PCL-core anchoring, PEG/corona exposure and DTX-Abi association.") +
  xlim(0, 1) + ylim(0, 1) +
  theme_void(base_family = plot_base_family) +
  theme(plot.margin = margin(5, 5, 5, 5))
p_system_map <- panel_letter(p_system_map, "a")

clear_panels <- list(A = p_system_map)

# B. Rg summary with system-specific component labels. Each facet is one system,
# so "single drug/prodrug"is no longer ambiguous.
if (exists("rg_summary") && nrow(rg_summary) > 0) {
  rg_keep_clear <- c(
    "Abi-SI-C18 prodrug", "DTX-SI-C18 prodrug", "DTX+Abi total prodrugs",
    "polymer + DTX+Abi", "polymer", "PEG block", "PCL block", "whole non-water assembly"
  )
  rg_clear <- rg_summary %>%
    mutate(
      system = factor(system, levels = system_order_clear),
      system_short = system_short_label(as.character(system)),
      component = metric_clear_label(system_dir, as.character(metric))
    ) %>%
    filter(component %in% rg_keep_clear) %>%
    mutate(component = factor(component, levels = rev(rg_keep_clear))) %>%
    bar_label_columns_custom(digits = 2)
  
  p_rg_clear <- ggplot(rg_clear, aes(x = component, y = mean, fill = system)) +
    geom_col(width = 0.60, colour = "black", linewidth = 0.18) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.18) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 1.85) +
    scale_fill_manual(values = system_colors_clear, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    facet_wrap(~system_short, nrow = 1, scales = "free_y") +
    labs(title = "Assembly compactness", x = NULL, y = "Rg (nm)") +
    compact_theme() +
    theme(legend.position = "none")
  
  save_pdf(p_rg_clear, "MD_Rg_final_window_clear.pdf", width = 10.5, height = 3.2)
  readr::write_csv(rg_clear, file.path(summary_dir, "MD_Rg_final_window_clear_summary.csv"))
  clear_panels[["B"]] <- panel_letter(p_rg_clear, "b")
}

# C. Mixed COM distances. This panel uses only the mixed system and labels every
# interaction explicitly.
if (exists("dist_summary") && nrow(dist_summary) > 0) {
  dist_keep_clear <- c("total drug-PCL", "DTX-PCL", "Abi-PCL", "total drug-PEG", "DTX-PEG", "Abi-PEG", "DTX-Abi", "DTX-polymer", "Abi-polymer")
  dist_clear <- dist_summary %>%
    filter(system_dir == "pcldtxsic18abisic18", as.character(metric) %in% dist_keep_clear) %>%
    mutate(
      interaction = as.character(metric),
      class = mixed_metric_class(interaction),
      interaction = factor(interaction, levels = rev(dist_keep_clear))
    ) %>%
    bar_label_columns_custom(digits = 2)
  
  p_dist_clear <- ggplot(dist_clear, aes(x = interaction, y = mean, fill = class)) +
    geom_col(width = 0.60, colour = "black", linewidth = 0.18) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.18) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 1.85) +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Mixed-system localization", x = NULL, y = "COM distance (nm)") +
    compact_theme() +
    theme(legend.position = "bottom")
  
  save_pdf(p_dist_clear, "MD_mixed_COM_distance.pdf", width = 6.4, height = 4.0)
  readr::write_csv(dist_clear, file.path(summary_dir, "MD_mixed_COM_distance_summary.csv"))
  clear_panels[["C"]] <- panel_letter(p_dist_clear, "c")
}

# D. Mixed RDF curves. Faceting by interaction removes the need for a line legend.
if (exists("rdf_df") && nrow(rdf_df) > 0) {
  rdf_keep_clear <- c("DTX-PCL", "Abi-PCL", "DTX-PEG", "Abi-PEG", "DTX-Abi", "total drug-PCL", "total drug-PEG")
  rdf_clear <- rdf_df %>%
    filter(system_dir == "pcldtxsic18abisic18", as.character(metric) %in% rdf_keep_clear) %>%
    mutate(interaction = factor(as.character(metric), levels = rdf_keep_clear))
  
  p_rdf_clear <- ggplot(rdf_clear, aes(x = x_raw, y = value)) +
    geom_line(linewidth = 0.35, colour = "#4C8EA0", alpha = 0.95) +
    facet_wrap(~interaction, scales = "free_y", ncol = 4) +
    labs(title = "Local enrichment by RDF", x = "r (nm)", y = "g(r)") +
    compact_theme() +
    theme(legend.position = "none")
  
  save_pdf(p_rdf_clear, "MD_mixed_RDF.pdf", width = 9.5, height = 4.2)
  clear_panels[["D"]] <- panel_letter(p_rdf_clear, "d")
}

# E. Mixed H-bond summary.
if (exists("hbond_summary") && nrow(hbond_summary) > 0) {
  hb_keep_clear <- c("DTX-polymer", "Abi-polymer", "DTX-Abi", "DTX-PCL", "Abi-PCL", "DTX-PEG", "Abi-PEG", "total drug-polymer")
  hb_clear <- hbond_summary %>%
    filter(system_dir == "pcldtxsic18abisic18", as.character(metric) %in% hb_keep_clear) %>%
    mutate(
      interaction = as.character(metric),
      class = mixed_metric_class(interaction),
      interaction = factor(interaction, levels = rev(hb_keep_clear))
    ) %>%
    bar_label_columns_custom(digits = 2)
  
  p_hb_clear <- ggplot(hb_clear, aes(x = interaction, y = mean, fill = class)) +
    geom_col(width = 0.60, colour = "black", linewidth = 0.18) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.18) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 1.85) +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Segment-specific hydrogen-bond contribution", x = NULL, y = "H-bonds", fill = NULL) +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    compact_theme() +
    theme(legend.position = "none")
  
  save_pdf(p_hb_clear, "MD_mixed_Hbond.pdf", width = 6.2, height = 3.8)
  readr::write_csv(hb_clear, file.path(summary_dir, "MD_mixed_Hbond_summary.csv"))
  clear_panels[["E"]] <- panel_letter(p_hb_clear, "e")
}

# F. Mixed contact numbers, scaled to x10^3 atom contacts.
if (exists("contact_summary") && nrow(contact_summary) > 0) {
  contact_keep_clear <- c("DTX-polymer", "Abi-polymer", "DTX-PCL", "Abi-PCL", "DTX-PEG", "Abi-PEG", "DTX-Abi", "DTX lig-PCL", "Abi lig-PCL", "DTX par-PCL", "Abi par-PCL")
  contact_clear <- contact_summary %>%
    filter(as.character(metric) %in% contact_keep_clear) %>%
    mutate(
      interaction = as.character(metric),
      class = mixed_metric_class(interaction),
      mean_raw = mean,
      sd_raw = sd,
      mean = mean / 1000,
      sd = sd / 1000,
      interaction = factor(interaction, levels = rev(contact_keep_clear))
    ) %>%
    bar_label_columns_custom(digits = 1)
  
  p_contact_clear <- ggplot(contact_clear, aes(x = interaction, y = mean, fill = class)) +
    geom_col(width = 0.60, colour = "black", linewidth = 0.18) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.18) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 1.85) +
    scale_fill_manual(values = domain_colors, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(title = "Near-contact enrichment", subtitle = "Atom contacts within the gmx mindist cutoff, shown as x10^3.", x = NULL, y = "Contacts (x10³)") +
    compact_theme() +
    theme(legend.position = "none")
  
  save_pdf(p_contact_clear, "MD_mixed_contact_number.pdf", width = 6.2, height = 4.2)
  readr::write_csv(contact_clear, file.path(summary_dir, "MD_mixed_contact_number_summary.csv"))
  clear_panels[["F"]] <- panel_letter(p_contact_clear, "f")
}

# G. SASA-derived contact area, mixed-system focused panel.
if (exists("sasa_summary") && nrow(sasa_summary) > 0) {
  sasa_keep_clear <- c("DTX-PCL", "Abi-PCL", "DTX-PEG", "Abi-PEG", "DTX-Abi", "DTX-polymer", "Abi-polymer")
  sasa_clear <- sasa_summary %>%
    filter(system == "Mixed DTX-SI-C18/Abi-SI-C18@PCL", as.character(pair) %in% sasa_keep_clear) %>%
    mutate(
      class = dplyr::case_when(
        stringr::str_detect(as.character(pair), "PCL") ~ "PCL-core contact area",
        stringr::str_detect(as.character(pair), "PEG") ~ "PEG/corona contact area",
        stringr::str_detect(as.character(pair), "DTX-Abi") ~ "DTX-Abi contact area",
        TRUE ~ "whole-polymer contact area"
      ),
      pair = factor(as.character(pair), levels = rev(sasa_keep_clear))
    )
  if (nrow(sasa_clear) > 0) {
    p_sasa_contact_clear <- ggplot(sasa_clear, aes(x = pair, y = mean, fill = class)) +
      geom_col(width = 0.62, colour = "black", linewidth = 0.18) +
      geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.18) +
      geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 1.8) +
      scale_fill_manual(values = sasa_domain_colors, drop = FALSE) +
      scale_y_continuous(expand = bar_y_expand()) +
      coord_flip(clip = "off") +
      labs(title = "SASA-derived contact area", subtitle = "Calculated directly from matched SASA XVG files.", x = NULL, y = "Area (nm²)", fill = NULL) +
      compact_theme() +
      theme(legend.position = "right")
    save_pdf(p_sasa_contact_clear, "MD_mixed_SASA_contact_area.pdf", width = 6.8, height = 4.2)
    readr::write_csv(sasa_clear, file.path(summary_dir, "MD_mixed_SASA_contact_area.csv"))
    clear_panels[["G"]] <- panel_letter(p_sasa_contact_clear, "g")
  }
}

# H. Segment-specific COM distance heatmap.
if (exists("segment_dist_summary") && nrow(segment_dist_summary) > 0) {
  seg_keep_clear <- c(
    "drug-polymer", "par-polymer", "lig-polymer",
    "drug-PCL", "par-PCL", "lig-PCL",
    "drug-PEG", "par-PEG", "lig-PEG",
    "drug-drug", "par-par", "lig-lig", "par-lig", "lig-par"
  )
  segment_dist_clear <- segment_dist_summary %>%
    filter(segment_pair %in% seg_keep_clear) %>%
    mutate(
      prodrug = dplyr::recode(prodrug, "DTX"= "DTX-SI-C18", "Abi"= "Abi-SI-C18", "DTX-Abi"= "DTX-Abi cross-pair"),
      prodrug = factor(prodrug, levels = c("Abi-SI-C18", "DTX-SI-C18", "DTX-Abi cross-pair")),
      segment_pair = factor(segment_pair, levels = seg_keep_clear)
    )
  
  p_seg_dist_clear <- ggplot(segment_dist_clear, aes(x = segment_pair, y = prodrug, fill = mean)) +
    geom_tile(colour = "white", linewidth = 0.30) +
    geom_text(aes(label = sprintf("%.2f", mean)), size = 2.0) +
    scale_fill_heatmap("Distance\n(nm)") +
    labs(title = "Whole-prodrug and segment-resolved localization", x = NULL, y = NULL) +
    compact_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
  save_pdf(p_seg_dist_clear, "MD_whole_and_segment_distance_heatmap_clear.pdf", width = 9.2, height = 3.4)
  readr::write_csv(segment_dist_clear, file.path(summary_dir, "MD_whole_and_segment_distance_summary_clear.csv"))
  clear_panels[["H"]] <- panel_letter(p_seg_dist_clear, "h")
}

# I. Segment-specific H-bond heatmap.
if (exists("segment_hbond_summary") && nrow(segment_hbond_summary) > 0) {
  seg_keep_clear <- c(
    "drug-polymer", "par-polymer", "lig-polymer",
    "drug-PCL", "par-PCL", "lig-PCL",
    "drug-PEG", "par-PEG", "lig-PEG",
    "drug-drug", "par-par", "lig-lig", "par-lig", "lig-par"
  )
  segment_hbond_clear <- segment_hbond_summary %>%
    filter(segment_pair %in% seg_keep_clear) %>%
    mutate(
      prodrug = dplyr::recode(prodrug, "DTX"= "DTX-SI-C18", "Abi"= "Abi-SI-C18", "DTX-Abi"= "DTX-Abi cross-pair"),
      prodrug = factor(prodrug, levels = c("Abi-SI-C18", "DTX-SI-C18", "DTX-Abi cross-pair")),
      segment_pair = factor(segment_pair, levels = seg_keep_clear)
    )
  
  p_seg_hbond_clear <- ggplot(segment_hbond_clear, aes(x = segment_pair, y = prodrug, fill = mean)) +
    geom_tile(colour = "white", linewidth = 0.30) +
    geom_text(aes(label = sprintf("%.2f", mean)), size = 2.0) +
    scale_fill_heatmap("H-bonds") +
    labs(title = "Whole-prodrug and segment-resolved H-bonds", x = NULL, y = NULL) +
    compact_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
  save_pdf(p_seg_hbond_clear, "MD_whole_and_segment_Hbond_heatmap_clear.pdf", width = 9.2, height = 3.4)
  readr::write_csv(segment_hbond_clear, file.path(summary_dir, "MD_whole_and_segment_Hbond_summary_clear.csv"))
  clear_panels[["I"]] <- panel_letter(p_seg_hbond_clear, "i")
}

# Optional SASA diagnostic. It is exported separately because degenerate means or
# large SD can make this readout less interpretable than Rg/distance/RDF/H-bond/contact data.
if (exists("sasa_summary") && nrow(sasa_summary) > 0) {
  sasa_diag <- sasa_summary %>%
    mutate(sd_to_mean = ifelse(mean > 0, sd / mean, NA_real_)) %>%
    summarise(
      n_pairs = n(),
      n_unique_mean = n_distinct(round(mean, 3)),
      max_sd_to_mean = suppressWarnings(max(sd_to_mean, na.rm = TRUE)),
      .groups = "drop"
    )
  readr::write_csv(sasa_diag, file.path(summary_dir, "MD_SASA_contact_area_diagnostic.csv"))
}

# Export the clearer integrated Figure. The layout is wide and separates global
# comparison from mixed-system mechanism.
if (length(clear_panels) >= 3) {
  order_clear <- intersect(c("A", "B", "C", "D", "E", "F", "G", "H", "I"), names(clear_panels))
  clear_panels <- clear_panels[order_clear]
  
  clear_integrated <- patchwork::wrap_plots(clear_panels, ncol = 4, guides = "collect") +
    patchwork::plot_annotation(
      title = "Molecular dynamics analysis of DTX-SI-C18/Abi-SI-C18 compatibility-guided co-assembly",
      subtitle = "Explicit system labels distinguish single-prodrug controls from the mixed co-assembly; mixed-system panels resolve PCL-core anchoring, PEG/corona exposure, prodrug-prodrug association and segment-level interaction balance.",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 9, hjust = 0),
        plot.margin = margin(5, 5, 5, 5)
      )
    )
  
  save_pdf(clear_integrated, "MD_coassembly_integrated_clear.pdf", width = 19.5, height = 10.6)
  save_pdf(clear_integrated, "MD_coassembly_integrated_clear_large.pdf", width = 21, height = 12.0)
}

# A compact interpretation table is exported for figure-legend drafting.
interpretation_rows <- tibble::tibble()
if (exists("rg_clear") && nrow(rg_clear) > 0) {
  interpretation_rows <- bind_rows(
    interpretation_rows,
    rg_clear %>% transmute(panel = "Rg", system, readout = component, mean, sd, unit = "nm")
  )
}
if (exists("dist_clear") && nrow(dist_clear) > 0) {
  interpretation_rows <- bind_rows(
    interpretation_rows,
    dist_clear %>% transmute(panel = "COM distance", system = "Mixed DTX-SI-C18/Abi-SI-C18@PCL", readout = as.character(interaction), mean, sd, unit = "nm")
  )
}
if (exists("hb_clear") && nrow(hb_clear) > 0) {
  interpretation_rows <- bind_rows(
    interpretation_rows,
    hb_clear %>% transmute(panel = "H-bonds", system = "Mixed DTX-SI-C18/Abi-SI-C18@PCL", readout = as.character(interaction), mean, sd, unit = "number")
  )
}
if (exists("contact_clear") && nrow(contact_clear) > 0) {
  interpretation_rows <- bind_rows(
    interpretation_rows,
    contact_clear %>% transmute(panel = "Contacts", system = "Mixed DTX-SI-C18/Abi-SI-C18@PCL", readout = as.character(interaction), mean = mean_raw, sd = sd_raw, unit = "atom contacts")
  )
}
if (exists("sasa_clear") && nrow(sasa_clear) > 0) {
  interpretation_rows <- bind_rows(
    interpretation_rows,
    sasa_clear %>% transmute(panel = "SASA contact area", system = "Mixed DTX-SI-C18/Abi-SI-C18@PCL", readout = as.character(pair), mean, sd, unit = "nm^2")
  )
}
if (nrow(interpretation_rows) > 0) {
  readr::write_csv(interpretation_rows, file.path(summary_dir, "MD_coassembly_interpretation_table.csv"))
}

# -----------------------------
# 20. Legacy combined overview, retained for quick checking
# -----------------------------
combined_panels <- list()
if (exists("p_rg_bar")) combined_panels[[length(combined_panels) + 1]] <- p_rg_bar + ggtitle("Rg")
if (exists("p_dist_bar")) combined_panels[[length(combined_panels) + 1]] <- p_dist_bar + ggtitle("COM distance")
if (exists("p_hb_bar")) combined_panels[[length(combined_panels) + 1]] <- p_hb_bar + ggtitle("H-bonds")
if (exists("p_contact")) combined_panels[[length(combined_panels) + 1]] <- p_contact + ggtitle("Contacts")
if (exists("p_sasa_bar")) combined_panels[[length(combined_panels) + 1]] <- p_sasa_bar + ggtitle("SASA contact area")
if (exists("p_seg_dist")) combined_panels[[length(combined_panels) + 1]] <- p_seg_dist + ggtitle("Segment distance")

if (length(combined_panels) >= 2) {
  combined_plot <- patchwork::wrap_plots(combined_panels, ncol = 2) +
    patchwork::plot_annotation(title = "Drugspolymers MD overview: compatibility-guided co-assembly metrics")
  save_pdf(combined_plot, "MD_legacy_overview.pdf", width = 12, height = 10)
}

# -----------------------------
# 21. Run log
# -----------------------------
run_log <- tibble::tibble(
  item = c("script_dir", "root_dir", "out_dir", "n_xvg_files", "n_converted_series_rows", "integrated_panels"),
  value = c(script_dir, root_dir, out_dir, nrow(xvg_index), nrow(all_xvg), paste(c(names(integrated_panels), paste0("clear:", names(clear_panels))), collapse = ","))
)
readr::write_csv(run_log, file.path(out_dir, "run_log.csv"))

message("Done.")
message("Diagnostic figures saved to: ", fig_dir)
message("Integrated wide Figure: ", file.path(fig_dir, "MD_coassembly_integrated_wide.pdf"))
message("Converted CSV and summary tables saved to: ", out_dir)


# ============================================================
# 22. Organized DTX/Abi co-assembly MD figure sets
# ============================================================
# This section rewrites the organized figure output after the full XVG analysis.
# It creates:
#   1) main-text Figure 5 MD module: concise, direct mechanistic evidence;
#   2) Extended Data module: secondary mechanistic support;
#   3) SI module: trajectory diagnostics and auxiliary transparency.
# Each panel and composite figure is exported as PDF plus TIFF, with PNG fallback.

MolecularDynamic_root <- file.path(project_root, "results", "figures")
MolecularDynamic_main_dir <- file.path(project_root, "results", "figures", "md", "DTX_Abi_coassembly_MD", "summary")
MolecularDynamic_ext_dir  <- file.path(project_root, "results", "figures", "md", "DTX_Abi_coassembly_MD", "mechanistic_support")
MolecularDynamic_si_dir   <- file.path(project_root, "results", "figures", "md", "DTX_Abi_coassembly_MD", "trajectory_support")
for (d in c(MolecularDynamic_main_dir, MolecularDynamic_ext_dir, MolecularDynamic_si_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

save_MolecularDynamic_plot <- function(plot, folder, stem, width, height, dpi = 600) {
  pdf_file <- file.path(folder, paste0(stem, ".pdf"))
  tif_file <- file.path(folder, paste0(stem, ".tif"))
  png_file <- file.path(folder, paste0(stem, ".png"))
  
  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = safe_pdf_device,
    bg = "white",
    limitsize = FALSE
  )
  
  ok_tif <- tryCatch({
    ggplot2::ggsave(
      filename = tif_file,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = "tiff",
      compression = "lzw",
      bg = "white",
      limitsize = FALSE
    )
    TRUE
  }, error = function(e) {
    warning("TIFF export failed for ", stem, ": ", conditionMessage(e), ". Writing PNG instead.")
    FALSE
  })
  
  if (!ok_tif) {
    ggplot2::ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = "png",
      bg = "white",
      limitsize = FALSE
    )
  }
  invisible(pdf_file)
}

MolecularDynamic_label_bar <- function(df, value_col = "mean", sd_col = "sd", digits = 1, pad_frac = 0.045) {
  if (nrow(df) == 0) return(df)
  value <- df[[value_col]]
  sdv <- if (sd_col %in% names(df)) df[[sd_col]] else rep(0, length(value))
  sdv[is.na(sdv) | !is.finite(sdv)] <- 0
  ymax <- suppressWarnings(max(value + sdv, na.rm = TRUE))
  pad <- ifelse(is.finite(ymax) && ymax > 0, ymax * pad_frac, 0.05)
  df %>%
    mutate(
      mean_label = ifelse(is.na(.data[[value_col]]), "", sprintf(paste0("%.", digits, "f"), .data[[value_col]])),
      label_y = .data[[value_col]] + sdv + pad
    )
}

MolecularDynamic_main_colors <- c(
  "single prodrug domain"= "#79B8AE",
  "mixed prodrug domain"= "#D97968",
  "single carrier+prodrug assembly"= "#8FB9C9",
  "mixed carrier+prodrug assembly"= "#B58AAE",
  "single-system interface"= "#79B8AE",
  "mixed total prodrug-polymer interface"= "#B58AAE",
  "mixed individual prodrug-polymer interface"= "#C39BC0",
  "PCL-core contact area"= "#4C8EA0",
  "PEG/corona contact area"= "#E0B678",
  "DTX-Abi contact area"= "#D76F61"
)

# -----------------------------
# 22.1 Summary MD module
# -----------------------------
main_panels <- list()
main_panel_meta <- tibble::tibble()

# Main panel A: include both single and mixed prodrug-domain Rg and the whole carrier+prodrug assembly Rg.
# This directly addresses the user's concern that the main compactness panel must include
# single-system carrier+drug assembly comparisons.
if (exists("rg_clear") && nrow(rg_clear) > 0) {
  rg_main <- rg_clear %>%
    mutate(
      system_short = as.character(system_short),
      component_chr = as.character(component),
      display = dplyr::case_when(
        system_short == "Abi single"& component_chr == "Abi-SI-C18 prodrug"~ "Abi prodrug\n(single)",
        system_short == "DTX single"& component_chr == "DTX-SI-C18 prodrug"~ "DTX prodrug\n(single)",
        system_short == "DTX+Abi mixed"& component_chr == "Abi-SI-C18 prodrug"~ "Abi prodrug\n(+DTX mixed)",
        system_short == "DTX+Abi mixed"& component_chr == "DTX-SI-C18 prodrug"~ "DTX prodrug\n(+Abi mixed)",
        system_short == "DTX+Abi mixed"& component_chr == "DTX+Abi total prodrugs"~ "DTX+Abi total\nprodrugs",
        system_short == "Abi single"& component_chr == "whole non-water assembly"~ "carrier+Abi\n(single)",
        system_short == "DTX single"& component_chr == "whole non-water assembly"~ "carrier+DTX\n(single)",
        system_short == "DTX+Abi mixed"& component_chr == "polymer + DTX+Abi"~ "carrier+DTX+Abi\n(mixed)",
        TRUE ~ NA_character_
      ),
      class = dplyr::case_when(
        stringr::str_detect(display, "carrier\\+") & stringr::str_detect(display, "mixed") ~ "mixed carrier+prodrug assembly",
        stringr::str_detect(display, "carrier\\+") ~ "single carrier+prodrug assembly",
        stringr::str_detect(display, "mixed") | stringr::str_detect(display, "DTX\\+Abi total") ~ "mixed prodrug domain",
        TRUE ~ "single prodrug domain"
      )
    ) %>%
    filter(!is.na(display)) %>%
    mutate(display = factor(display, levels = rev(c(
      "Abi prodrug\n(single)",
      "Abi prodrug\n(+DTX mixed)",
      "DTX prodrug\n(single)",
      "DTX prodrug\n(+Abi mixed)",
      "DTX+Abi total\nprodrugs",
      "carrier+Abi\n(single)",
      "carrier+DTX\n(single)",
      "carrier+DTX+Abi\n(mixed)"
    )))) %>%
    MolecularDynamic_label_bar(digits = 2, pad_frac = 0.04)
  
  p_main_rg <- ggplot(rg_main, aes(x = display, y = mean, fill = class)) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.22) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.15, linewidth = 0.22) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 2.45) +
    scale_fill_manual(values = MolecularDynamic_main_colors, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(
      title = "Assembly compactness",
      subtitle = "Final-window Rg of prodrug domains and carrier+prodrug assemblies",
      x = NULL,
      y = "Rg (nm)",
      fill = NULL
    ) +
    compact_theme(base_size = 7) +
    theme(legend.position = "bottom", legend.box = "vertical")
  
  main_panels[["A"]] <- panel_letter(p_main_rg, "a")
  save_MolecularDynamic_plot(main_panels[["A"]], MolecularDynamic_main_dir, "MD_assembly_compactness", 5.4, 4.5)
}

# Main panel B: include single-system drug-polymer interfaces plus mixed total and pair-specific interfaces.
# This panel is designed to show: (i) single-system carrier interface context;
# (ii) mixed total prodrug-polymer interface; (iii) Abi remains carrier/PCL-associated;
# (iv) DTX-Abi direct interface is limited.
if (exists("sasa_summary") && nrow(sasa_summary) > 0) {
  sasa_main <- sasa_summary %>%
    mutate(pair_chr = as.character(pair), system_chr = as.character(system)) %>%
    filter(
      (system_chr == "Single Abi-SI-C18@PCL"& pair_chr == "Abi-polymer") |
        (system_chr == "Single DTX-SI-C18@PCL"& pair_chr == "DTX-polymer") |
        (system_chr == "Mixed DTX-SI-C18/Abi-SI-C18@PCL"& pair_chr %in% c(
          "total drug-polymer", "Abi-polymer", "DTX-polymer", "Abi-PCL", "DTX-PCL", "Abi-PEG", "DTX-PEG", "DTX-Abi"
        ))
    ) %>%
    mutate(
      display = dplyr::case_when(
        system_chr == "Single Abi-SI-C18@PCL"& pair_chr == "Abi-polymer"~ "Abi-polymer\n(single)",
        system_chr == "Single DTX-SI-C18@PCL"& pair_chr == "DTX-polymer"~ "DTX-polymer\n(single)",
        pair_chr == "total drug-polymer"~ "DTX+Abi-polymer\n(mixed total)",
        pair_chr == "Abi-polymer"~ "Abi-polymer\n(mixed)",
        pair_chr == "DTX-polymer"~ "DTX-polymer\n(mixed)",
        pair_chr == "Abi-PCL"~ "Abi-PCL\n(mixed)",
        pair_chr == "DTX-PCL"~ "DTX-PCL\n(mixed)",
        pair_chr == "Abi-PEG"~ "Abi-PEG\n(mixed)",
        pair_chr == "DTX-PEG"~ "DTX-PEG\n(mixed)",
        pair_chr == "DTX-Abi"~ "DTX-Abi\n(mixed)",
        TRUE ~ pair_chr
      ),
      class = dplyr::case_when(
        stringr::str_detect(display, "single") ~ "single-system interface",
        pair_chr == "total drug-polymer"~ "mixed total prodrug-polymer interface",
        pair_chr %in% c("Abi-polymer", "DTX-polymer") ~ "mixed individual prodrug-polymer interface",
        stringr::str_detect(pair_chr, "PCL") ~ "PCL-core contact area",
        stringr::str_detect(pair_chr, "PEG") ~ "PEG/corona contact area",
        pair_chr == "DTX-Abi"~ "DTX-Abi contact area",
        TRUE ~ "mixed individual prodrug-polymer interface"
      ),
      display = factor(display, levels = rev(c(
        "Abi-polymer\n(single)",
        "DTX-polymer\n(single)",
        "DTX+Abi-polymer\n(mixed total)",
        "Abi-polymer\n(mixed)",
        "DTX-polymer\n(mixed)",
        "Abi-PCL\n(mixed)",
        "DTX-PCL\n(mixed)",
        "Abi-PEG\n(mixed)",
        "DTX-PEG\n(mixed)",
        "DTX-Abi\n(mixed)"
      )))
    ) %>%
    MolecularDynamic_label_bar(digits = 1, pad_frac = 0.035)
  
  p_main_sasa <- ggplot(sasa_main, aes(x = display, y = mean, fill = class)) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.22) +
    geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.15, linewidth = 0.22) +
    geom_text(aes(y = label_y, label = mean_label), hjust = 0, size = 2.35) +
    scale_fill_manual(values = MolecularDynamic_main_colors, drop = FALSE) +
    scale_y_continuous(expand = bar_y_expand()) +
    coord_flip(clip = "off") +
    labs(
      title = "Prodrug-carrier interfacial contact",
      subtitle = "SASA-derived contact area; single controls and mixed interfaces are shown together",
      x = NULL,
      y = "Contact area (nm²)",
      fill = NULL
    ) +
    compact_theme(base_size = 7) +
    theme(legend.position = "bottom", legend.box = "vertical")
  
  main_panels[["B"]] <- panel_letter(p_main_sasa, "b")
  save_MolecularDynamic_plot(main_panels[["B"]], MolecularDynamic_main_dir, "MD_interfacial_contact_area", 6.0, 5.0)
}

# Main panel C: whole-prodrug and segment-resolved localization heatmap.
if (exists("segment_dist_clear") && nrow(segment_dist_clear) > 0) {
  seg_order <- c(
    "drug-polymer", "par-polymer", "lig-polymer",
    "drug-PCL", "par-PCL", "lig-PCL",
    "drug-PEG", "par-PEG", "lig-PEG",
    "drug-drug", "par-par", "lig-lig", "par-lig", "lig-par"
  )
  prodrug_order <- c("DTX-Abi cross-pair", "DTX-SI-C18", "Abi-SI-C18")
  segment_main <- segment_dist_clear %>%
    mutate(
      segment_pair = factor(as.character(segment_pair), levels = seg_order),
      prodrug = dplyr::recode(as.character(prodrug), "DTX-Abi"= "DTX-Abi cross-pair", "DTX"= "DTX-SI-C18", "Abi"= "Abi-SI-C18"),
      prodrug = factor(prodrug, levels = prodrug_order)
    )
  
  p_main_segdist <- ggplot(segment_main, aes(x = segment_pair, y = prodrug, fill = mean)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", mean)), size = 2.25) +
    scale_fill_heatmap("COM distance\n(nm)") +
    labs(
      title = "Whole-prodrug and segment-resolved localization",
      subtitle = "Whole-drug terms precede parent/ligand segment terms for each carrier domain",
      x = NULL,
      y = NULL
    ) +
    compact_theme(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
  main_panels[["C"]] <- panel_letter(p_main_segdist, "c")
  save_MolecularDynamic_plot(main_panels[["C"]], MolecularDynamic_main_dir, "MD_segment_localization", 8.6, 3.8)
}

if (length(main_panels) > 0) {
  main_order <- intersect(c("A", "B", "C"), names(main_panels))
  main_widths <- c("A"= 1.05, "B"= 1.15, "C"= 1.75)[main_order]
  main_figure_MolecularDynamic <- patchwork::wrap_plots(main_panels[main_order], nrow = 1, widths = main_widths, guides = "collect") +
    patchwork::plot_annotation(
      title = "DTX-SI-C18 supports carrier-integrated incorporation of Abi-SI-C18 in mixed co-assembly",
      subtitle = "Final-window MD readouts compare single-prodrug controls with the mixed system and resolve prodrug-carrier interfaces and segment-level localization.",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(size = 8.5, hjust = 0)
      )
    )
  save_MolecularDynamic_plot(main_figure_MolecularDynamic, MolecularDynamic_main_dir, "MD_coassembly", 20.0, 5.7)
}

# -----------------------------
# 22.2 Mechanistic-support panels
# -----------------------------
ext_panels <- list()

if (exists("p_dist_clear")) {
  ext_panels[["A"]] <- panel_letter(p_dist_clear + labs(title = "Core/corona localization"), "a")
  save_MolecularDynamic_plot(ext_panels[["A"]], MolecularDynamic_ext_dir, "MD_COM_localization", 5.2, 3.8)
}
if (exists("p_rdf_clear")) {
  ext_panels[["B"]] <- panel_letter(p_rdf_clear + labs(title = "RDF-resolved local enrichment"), "b")
  save_MolecularDynamic_plot(ext_panels[["B"]], MolecularDynamic_ext_dir, "MD_RDF_local_enrichment", 7.4, 3.6)
}
if (exists("p_contact_clear")) {
  ext_panels[["C"]] <- panel_letter(p_contact_clear + labs(title = "Near-contact enrichment"), "c")
  save_MolecularDynamic_plot(ext_panels[["C"]], MolecularDynamic_ext_dir, "MD_near_contact_enrichment", 5.4, 4.0)
}
if (exists("p_hb_clear")) {
  ext_panels[["D"]] <- panel_letter(p_hb_clear + labs(title = "Auxiliary H-bond contribution"), "d")
  save_MolecularDynamic_plot(ext_panels[["D"]], MolecularDynamic_ext_dir, "MD_Hbond_contribution", 5.2, 3.6)
}
if (exists("p_seg_hbond_clear")) {
  ext_panels[["E"]] <- panel_letter(p_seg_hbond_clear + labs(title = "Whole-prodrug and segment-resolved H-bonds"), "e")
  save_MolecularDynamic_plot(ext_panels[["E"]], MolecularDynamic_ext_dir, "MD_segment_Hbonds", 6.2, 3.4)
}
if (exists("p_sasa_bar")) {
  ext_panels[["F"]] <- panel_letter(p_sasa_bar + labs(title = "Single- and mixed-system interfacial area"), "f")
  save_MolecularDynamic_plot(ext_panels[["F"]], MolecularDynamic_ext_dir, "MD_SASA_single_mixed_comparison", 6.0, 5.1)
}

if (length(ext_panels) > 0) {
  ext_order <- intersect(c("A", "B", "C", "D", "E", "F"), names(ext_panels))
  extended_figure_MolecularDynamic <- patchwork::wrap_plots(ext_panels[ext_order], ncol = 3, guides = "collect") +
    patchwork::plot_annotation(
      title = "Extended molecular dynamics readouts for DTX-SI-C18/Abi-SI-C18 mixed co-assembly",
      subtitle = "Secondary metrics resolve core/corona localization, local enrichment, atom-contact enrichment and segment-specific hydrogen-bond contributions.",
      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0), plot.subtitle = element_text(size = 8.5, hjust = 0))
    )
  save_MolecularDynamic_plot(extended_figure_MolecularDynamic, MolecularDynamic_ext_dir, "MD_extended_support", 16.8, 9.6)
}

# -----------------------------
# 22.3 Trajectory-support panels
# -----------------------------
si_panels <- list()
if (exists("p_rmsd")) {
  si_panels[["A"]] <- panel_letter(p_rmsd + labs(title = "RMSD time-series"), "a")
  save_MolecularDynamic_plot(si_panels[["A"]], MolecularDynamic_si_dir, "MD_RMSD_time_series", 6.2, 4.8)
}
if (exists("p_rg_time")) {
  si_panels[["B"]] <- panel_letter(p_rg_time + labs(title = "Rg time-series"), "b")
  save_MolecularDynamic_plot(si_panels[["B"]], MolecularDynamic_si_dir, "MD_Rg_time_series", 6.2, 4.8)
}
if (exists("p_dist_time")) {
  si_panels[["C"]] <- panel_letter(p_dist_time + labs(title = "COM-distance time-series"), "c")
  save_MolecularDynamic_plot(si_panels[["C"]], MolecularDynamic_si_dir, "MD_COM_distance_time_series", 6.2, 4.8)
}
if (exists("p_hb_time")) {
  si_panels[["D"]] <- panel_letter(p_hb_time + labs(title = "H-bond time-series"), "d")
  save_MolecularDynamic_plot(si_panels[["D"]], MolecularDynamic_si_dir, "MD_Hbond_time_series", 6.2, 4.8)
}
if (exists("p_raw_sasa_time")) {
  si_panels[["E"]] <- panel_letter(p_raw_sasa_time + labs(title = "Raw SASA time-series"), "e")
  save_MolecularDynamic_plot(si_panels[["E"]], MolecularDynamic_si_dir, "MD_raw_SASA_time_series", 6.8, 5.0)
}
if (exists("p_raw_sasa_bar")) {
  si_panels[["F"]] <- panel_letter(p_raw_sasa_bar + labs(title = "Raw SASA final-window summary"), "f")
  save_MolecularDynamic_plot(si_panels[["F"]], MolecularDynamic_si_dir, "MD_raw_SASA_final_window", 6.8, 5.2)
}
if (exists("p_sasa_time")) {
  si_panels[["G"]] <- panel_letter(p_sasa_time + labs(title = "SASA-contact-area time-series"), "g")
  save_MolecularDynamic_plot(si_panels[["G"]], MolecularDynamic_si_dir, "MD_SASA_contact_area_time_series", 6.8, 5.0)
}
if (exists("p_cluster")) {
  si_panels[["H"]] <- panel_letter(p_cluster + labs(title = "Cluster-size diagnostic"), "h")
  save_MolecularDynamic_plot(si_panels[["H"]], MolecularDynamic_si_dir, "MD_cluster_size_diagnostic", 6.2, 4.8)
}
if (exists("p_qc")) {
  si_panels[["I"]] <- panel_letter(p_qc + labs(title = "MD quality-control readouts"), "i")
  save_MolecularDynamic_plot(si_panels[["I"]], MolecularDynamic_si_dir, "MD_quality_control", 6.2, 4.8)
}

if (length(si_panels) > 0) {
  si_order <- intersect(LETTERS[1:9], names(si_panels))
  si_figure_MolecularDynamic <- patchwork::wrap_plots(si_panels[si_order], ncol = 3, guides = "collect") +
    patchwork::plot_annotation(
      title = "Supplementary molecular dynamics diagnostics for DTX-SI-C18/Abi-SI-C18 co-assembly",
      subtitle = "Auxiliary time-series, raw SASA, cluster and QC readouts support trajectory assessment and data transparency.",
      theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0), plot.subtitle = element_text(size = 8.5, hjust = 0))
    )
  save_MolecularDynamic_plot(si_figure_MolecularDynamic, MolecularDynamic_si_dir, "MD_diagnostics", 17.0, 13.0)
}

# -----------------------------
# 22.4 figure.log: English legends and Chinese translations
# -----------------------------
figure_log_lines <- c(
  "figure.log for drugspolymers260606MolecularDynamic.r",
  "============================================================",
  "",
  "SUMMARY MD MODULE",
  "------------------------------------------------------------",
  "Large figure file: MD_coassembly.pdf / .tif or .png",
  "English legend:",
  "Molecular dynamics analysis of DTX-SI-C18-assisted Abi-SI-C18 incorporation in mixed co-assembly. (a) Final-window radius of gyration (Rg) compares single-prodrug domains, the corresponding single carrier+prodrug assemblies, and the mixed DTX-SI-C18/Abi-SI-C18 assembly. The mixed prodrug domain and carrier+prodrug assembly remain compact, indicating that adding DTX-SI-C18 does not induce pronounced expansion of Abi-SI-C18-containing assemblies. (b) SASA-derived interfacial contact area compares single-system drug-polymer interfaces with the total and pair-specific interfaces in the mixed assembly. In the mixed system, Abi-SI-C18 and DTX-SI-C18 maintain larger polymer/PCL interfaces than the direct DTX-Abi interface, supporting carrier-integrated incorporation rather than dominant DTX-Abi self-association. (c) Whole-prodrug and segment-resolved COM-distance heatmap. Whole-drug rows provide overall localization references, while parent and ligand/modifier rows resolve segment-level contributions. DTX-SI-C18 is more carrier-proximal than Abi-SI-C18, particularly through its ligand/modifier segment, supporting its role as the stronger anchoring component in the mixed co-assembly.",
  "中文说明:",
  "DTX-SI-C18 促进 Abi-SI-C18 纳入混合共组装体系的分子动力学分析。(a) 终末窗口回转半径比较单一前药区域、相应的单一 carrier+prodrug 组装体，以及 DTX-SI-C18/Abi-SI-C18 混合组装体。混合前药区域和 carrier+prodrug 组装体仍保持紧凑，说明加入 DTX-SI-C18 后没有导致含 Abi-SI-C18 组装体明显膨胀。(b) 基于 SASA 的界面接触面积，将单一体系 drug-polymer 界面与混合体系的总界面和分子对界面进行比较。在混合体系中，Abi-SI-C18 和 DTX-SI-C18 与 polymer/PCL 的界面大于直接 DTX-Abi 界面，支持载体整合型纳入，而不是 DTX-Abi 自聚集主导。(c) 整体前药和分段 COM 距离热图。drug 整体行提供总体定位参照，parent 和 ligand/modifier 行解析片段水平贡献。DTX-SI-C18 比 Abi-SI-C18 更靠近载体，尤其是 ligand/modifier 片段，支持其在混合共组装中作为更强锚定组分。",
  "",
  "Main panel a: MD_assembly_compactness.pdf / .tif or .png",
  "English: Final-window Rg of single and mixed prodrug domains and carrier+prodrug assemblies. This panel provides the structural compactness context for DTX-SI-C18-assisted Abi-SI-C18 incorporation.",
  "中文: 单一和混合前药区域以及 carrier+prodrug 组装体的终末窗口 Rg。该图为 DTX-SI-C18 辅助 Abi-SI-C18 纳入提供结构紧密性的背景。",
  "",
  "Main panel b: MD_interfacial_contact_area.pdf / .tif or .png",
  "English: SASA-derived interfacial contact area comparing single drug-polymer interfaces with mixed total and pair-specific interfaces. Larger prodrug-polymer/PCL interfaces than DTX-Abi interface support carrier-integrated co-assembly.",
  "中文: 基于 SASA 计算的界面接触面积，比较单一 drug-polymer 界面以及混合体系中的总界面和分子对界面。prodrug-polymer/PCL 界面大于 DTX-Abi 界面，支持载体整合型共组装。",
  "",
  "Main panel c: MD_segment_localization.pdf / .tif or .png",
  "English: Whole-prodrug and segment-resolved localization by COM distance. Whole-drug terms show overall positioning, while parent and ligand/modifier terms identify segment-level anchoring differences between DTX-SI-C18 and Abi-SI-C18.",
  "中文: 基于 COM 距离的整体前药和分段定位分析。drug 整体项显示总体定位，parent 和 ligand/modifier 项显示 DTX-SI-C18 与 Abi-SI-C18 的片段锚定差异。",
  "",
  "MECHANISTIC SUPPORT MODULE",
  "------------------------------------------------------------",
  "Large figure file: MD_extended_support.pdf / .tif or .png",
  "English legend:",
  "Extended MD readouts supporting the DTX-SI-C18/Abi-SI-C18 mixed co-assembly mechanism. COM distances, RDF curves, atom-contact counts, H-bond counts and segment-resolved H-bonds provide secondary evidence for core/corona localization, local enrichment, carrier-associated contacts and auxiliary polar stabilization. These panels support the main-text conclusion but are not required for the most concise Figure 5 narrative.",
  "中文说明:",
  "支持 DTX-SI-C18/Abi-SI-C18 混合共组装机制的扩展 MD 读数。COM 距离、RDF 曲线、原子接触数、氢键数和分段氢键提供关于 core/corona 定位、局部富集、载体相关接触和辅助极性稳定的次级证据。这些图支持正文结论，但不是 Figure 5 最简洁叙事所必需。",
  "",
  "TRAJECTORY SUPPORT MODULE",
  "------------------------------------------------------------",
  "Large figure file: MD_diagnostics.pdf / .tif or .png",
  "English legend:",
  "Supplementary MD diagnostics and auxiliary readouts. RMSD, Rg, COM-distance and H-bond time-series support trajectory assessment. Raw SASA curves and final-window raw SASA values verify correct SASA parsing before contact-area calculation. SASA-contact-area time-series, cluster-size output and MD quality-control readouts are provided for transparency and internal validation.",
  "中文说明:",
  "补充 MD 诊断和辅助读数。RMSD、Rg、COM 距离和氢键时间序列用于评估轨迹。Raw SASA 曲线和终末窗口 raw SASA 用于验证 SASA 读取正确性，再进行 contact area 计算。SASA contact area 时间序列、cluster size 输出和 MD 质量控制读数用于数据透明性和内部验证。"
)

writeLines(figure_log_lines, con = file.path(summary_dir, "MD_figure_legend.log"), useBytes = TRUE)
# Figure legend log is written once to the summary-table directory.

MolecularDynamic_manifest <- tibble::tibble(
  figure_tier = c(
    rep("summary", 4),
    rep("mechanistic_support", 7),
    rep("trajectory_support", 10)
  ),
  figure_name = c(
    "MD_coassembly",
    "MD_assembly_compactness",
    "MD_interfacial_contact_area",
    "MD_segment_localization",
    "MD_extended_support",
    "MD_COM_localization",
    "MD_RDF_local_enrichment",
    "MD_near_contact_enrichment",
    "MD_Hbond_contribution",
    "MD_segment_Hbonds",
    "MD_SASA_single_mixed_comparison",
    "MD_diagnostics",
    "MD_RMSD_time_series",
    "MD_Rg_time_series",
    "MD_COM_distance_time_series",
    "MD_Hbond_time_series",
    "MD_raw_SASA_time_series",
    "MD_raw_SASA_final_window",
    "MD_SASA_contact_area_time_series",
    "MD_cluster_size_diagnostic",
    "MD_quality_control"
  ),
  output_folder = c(
    rep(MolecularDynamic_main_dir, 4),
    rep(MolecularDynamic_ext_dir, 7),
    rep(MolecularDynamic_si_dir, 10)
  )
)
readr::write_csv(MolecularDynamic_manifest, file.path(summary_dir, "MD_figure_manifest.csv"))

message("DTX/Abi co-assembly MD summary outputs: ", MolecularDynamic_main_dir)
message("DTX/Abi co-assembly MD legend log: ", file.path(summary_dir, "MD_figure_legend.log"))
message("All figures generated through save_pdf() or save_MolecularDynamic_plot() were exported as PDF plus TIFF, with PNG fallback if TIFF export failed.")

# ============================================================
# End of MolecularDynamic organized figure section
# ============================================================
