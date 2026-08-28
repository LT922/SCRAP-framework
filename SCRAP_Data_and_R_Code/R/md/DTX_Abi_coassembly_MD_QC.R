# DTX-SI-C18/Abi-SI-C18 co-assembly MD quality-control and extended analysis

rm(list = ls())
graphics.off()
set.seed(123)

# -----------------------------
# 0. Packages
# -----------------------------
packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr", "purrr",
  "tibble", "patchwork", "scales"
)
missing_packages <- packages[!packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE, repos = "https://cloud.r-project.org")
}
invisible(lapply(packages, library, character.only = TRUE))

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
rename <- dplyr::rename

# -----------------------------
# 1. User-adjustable settings
# -----------------------------
MAKE_INDIVIDUAL_XVG_PLOTS <- TRUE
FINAL_WINDOW_NS <- 50
OUTPUT_DPI <- 600
BASE_FAMILY <- "sans"

SMOOTH_RMSD_DISPLAY <- TRUE
RMSD_SMOOTH_WINDOW_NS <- 3

systems <- tibble::tribble(
  ~system_dir,             ~system_label,                                      ~short_label,      ~composition,
  "pcldtxsic18",           "DTX-SI-C18 single-prodrug system",                "DTX single",      "30 mPEG2k-PCL2k + 20 DTX-SI-C18",
  "pclabisic18",           "Abi-SI-C18 single-prodrug system",                "Abi single",      "30 mPEG2k-PCL2k + 20 Abi-SI-C18",
  "pcldtxsic18abisic18",   "DTX-SI-C18/Abi-SI-C18 mixed-prodrug system",      "DTX+Abi mixed",   "30 mPEG2k-PCL2k + 10 DTX-SI-C18 + 10 Abi-SI-C18"
)

system_levels <- c("DTX single", "Abi single", "DTX+Abi mixed")

system_colors <- c(
  "DTX single"= "#4F7DBA",
  "Abi single"= "#70AFA3",
  "DTX+Abi mixed"= "#D56B5F"
)

interaction_colors <- c(
  "PCL anchoring"= "#4C8EA0",
  "PEG contact"= "#E0B678",
  "polymer integration"= "#B58AAE",
  "DTX-Abi association"= "#D76F61",
  "self/cluster"= "#8FA6B3",
  "compactness"= "#76A9C2",
  "other"= "#9A9A9A"
)

heatmap_cols <- c("#D97968", "#F2E4D8", "#78BDB8")

# -----------------------------
# 2. Paths
# -----------------------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

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
project_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

root_dir <- file.path(project_root, "data", "md", "DTX_Abi_coassembly_MD")
fig_dir <- file.path(project_root, "results", "figures", "md", "DTX_Abi_coassembly_MD", "QC_support")
diagnostic_dir <- file.path(fig_dir, "trajectory_diagnostics")
fig_individual_dir <- file.path(diagnostic_dir, "all_individual_xvg")
table_dir <- file.path(project_root, "results", "tables", "md", "DTX_Abi_coassembly_MD", "QC_support")
converted_dir <- file.path(table_dir, "converted_xvg_csv")
summary_dir <- file.path(table_dir, "summary_tables")

for (d in c(fig_dir, diagnostic_dir, fig_individual_dir, table_dir, converted_dir, summary_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(root_dir)) {
  stop("Simulation root folder not found: ", root_dir, call. = FALSE)
}
root_dir <- normalizePath(root_dir, winslash = "/", mustWork = TRUE)

locate_system_folder <- function(root, system_dir) {
  candidates <- c(
    file.path(root, system_dir),
    file.path(root, system_dir, "ana"),
    file.path(root, "pcl260606", system_dir),
    file.path(root, "pcl260606", system_dir, "ana")
  )
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  has_xvg <- purrr::map_lgl(candidates, function(p) {
    dir.exists(p) && length(list.files(p, pattern = "\\.xvg$", full.names = TRUE, recursive = FALSE)) > 0
  })
  if (any(has_xvg)) return(candidates[which(has_xvg)[1]])
  candidates[dir.exists(candidates)][1] %||% NA_character_
}

systems <- systems %>%
  dplyr::mutate(
    ana_dir = purrr::map_chr(system_dir, ~locate_system_folder(root_dir, .x)),
    exists = !is.na(ana_dir) & dir.exists(ana_dir)
  )

if (!all(systems$exists)) {
  warning("Some system folders were not found or contained no XVG files:\n",
          paste(systems$system_dir[!systems$exists], collapse = ", "))
}

message("Root folder: ", root_dir)
message("Figure folder: ", fig_dir)
print(systems)

# -----------------------------
# 3. Plot helpers
# -----------------------------
clean_file_stem <- function(x) stringr::str_remove(basename(x), "\\.xvg$")

safe_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9_+\\-./]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

metric_group <- function(file_stem) {
  dplyr::case_when(
    stringr::str_detect(file_stem, "^temperature$|^pressure$|^density$|^potential$") ~ "QC",
    stringr::str_detect(file_stem, "^rmsd") & !stringr::str_detect(file_stem, "rmsd-dist") ~ "RMSD",
    stringr::str_detect(file_stem, "gyrate") ~ "Rg",
    stringr::str_detect(file_stem, "^rdf-") ~ "RDF",
    stringr::str_detect(file_stem, "^hbnum-") ~ "H-bond",
    stringr::str_detect(file_stem, "contacts") ~ "Contact number",
    stringr::str_detect(file_stem, "mindist") ~ "Minimum distance",
    stringr::str_detect(file_stem, "cogdist") ~ "COG distance",
    stringr::str_detect(file_stem, "-dist$|drug-polymer-dist|drug-pcl-dist|drug-peg-dist") ~ "COM distance",
    stringr::str_detect(file_stem, "^sasa-|^area_") ~ "SASA",
    stringr::str_detect(file_stem, "^size|^ntr|^nclid|cluster") ~ "Cluster",
    TRUE ~ "Other"
  )
}

metric_label <- function(x) {
  y <- x
  y <- stringr::str_replace_all(y, "^hbnum-", "H-bond: ")
  y <- stringr::str_replace_all(y, "^rdf-", "RDF: ")
  y <- stringr::str_replace_all(y, "^sasa-", "SASA: ")
  y <- stringr::str_replace_all(y, "^area_", "SASA area: ")
  y <- stringr::str_replace_all(y, "-cogdist", "COG distance")
  y <- stringr::str_replace_all(y, "-mindist", "minimum distance")
  y <- stringr::str_replace_all(y, "-contacts", "contacts")
  y <- stringr::str_replace_all(y, "-dist", "distance")
  y <- stringr::str_replace_all(y, "druggyrate", "prodrug Rg")
  y <- stringr::str_replace_all(y, "allgyrate", "whole assembly Rg")
  y <- stringr::str_replace_all(y, "polymergyrate", "polymer Rg")
  y <- stringr::str_replace_all(y, "gyrate_", "Rg: ")
  y <- stringr::str_replace_all(y, "dtxsic18", "DTX-SI-C18")
  y <- stringr::str_replace_all(y, "abisic18", "Abi-SI-C18")
  y <- stringr::str_replace_all(y, "dtx", "DTX")
  y <- stringr::str_replace_all(y, "abi", "Abi")
  y <- stringr::str_replace_all(y, "dru|drug", "Pro")
  y <- stringr::str_replace_all(y, "par", "Par")
  y <- stringr::str_replace_all(y, "lig", "Mod")
  y <- stringr::str_replace_all(y, "pcl", "PCL")
  y <- stringr::str_replace_all(y, "peg", "PEG")
  y <- stringr::str_replace_all(y, "polymer", "Pol")
  y <- stringr::str_replace_all(y, "_", "")
  y <- stringr::str_replace_all(y, "-", "–")
  stringr::str_squish(y)
}

interaction_class <- function(label) {
  x <- stringr::str_to_lower(label)
  dplyr::case_when(
    stringr::str_detect(x, "pcl") ~ "PCL anchoring",
    stringr::str_detect(x, "peg") ~ "PEG contact",
    stringr::str_detect(x, "polymer|pol") ~ "polymer integration",
    stringr::str_detect(x, "dtx.*abi|abi.*dtx") ~ "DTX-Abi association",
    stringr::str_detect(x, "cluster|self|drug-drug") ~ "self/cluster",
    stringr::str_detect(x, "rg|compact") ~ "compactness",
    TRUE ~ "other"
  )
}

theme_pub <- function(base_size = 8) {
  ggplot2::theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.30, colour = "black"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size),
      axis.title = ggplot2::element_text(colour = "black", size = base_size + 1),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1.5, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key.size = grid::unit(4, "mm"),
      panel.spacing = grid::unit(1.0, "lines")
    )
}

save_plot <- function(p, stem, width = 7.0, height = 4.5, folder = fig_dir, dpi = OUTPUT_DPI) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(folder, paste0(stem, ".pdf"))
  png_file <- file.path(folder, paste0(stem, ".png"))
  ggplot2::ggsave(pdf_file, p, width = width, height = height, units = "in",
                  device = grDevices::pdf, family = BASE_FAMILY, useDingbats = FALSE,
                  bg = "white", limitsize = FALSE)
  ggplot2::ggsave(png_file, p, width = width, height = height, units = "in",
                  dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(c(pdf_file, png_file))
}

mean_label_df <- function(df, value_col = "mean", sd_col = "sd", digits = 2, pad_frac = 0.04) {
  if (nrow(df) == 0) return(df)
  y <- df[[value_col]]
  s <- if (sd_col %in% names(df)) df[[sd_col]] else rep(0, length(y))
  s[!is.finite(s) | is.na(s)] <- 0
  ymax <- suppressWarnings(max(y + s, na.rm = TRUE))
  pad <- ifelse(is.finite(ymax) && ymax > 0, ymax * pad_frac, 0.05)
  df %>%
    dplyr::mutate(
      mean_label = sprintf(paste0("%.", digits, "f"), .data[[value_col]]),
      label_y = .data[[value_col]] + s + pad
    )
}

bar_expand <- function() ggplot2::expansion(mult = c(0.02, 0.22))

# -----------------------------
# 4. XVG parser
# -----------------------------
extract_xvg_meta <- function(lines, pattern) {
  hit <- lines[stringr::str_detect(lines, pattern)]
  if (length(hit) == 0) return(NA_character_)
  val <- stringr::str_match(hit[1], '"(.*)"')[, 2]
  ifelse(is.na(val), NA_character_, val)
}

parse_xvg_metadata <- function(file) {
  lines <- readLines(file, warn = FALSE)
  legend_lines <- lines[stringr::str_detect(lines, "^@\\s+s[0-9]+\\s+legend")]
  legends <- character(0)
  if (length(legend_lines) > 0) {
    tmp <- stringr::str_match(legend_lines, '^@\\s+s([0-9]+)\\s+legend\\s+"(.*)"')
    legends <- tmp[, 3]
    legends <- legends[!is.na(legends)]
  }
  list(
    title = extract_xvg_meta(lines, "^@\\s+title"),
    x_label = extract_xvg_meta(lines, "^@\\s+xaxis\\s+label"),
    y_label = extract_xvg_meta(lines, "^@\\s+yaxis\\s+label"),
    legends = legends
  )
}

read_xvg_long <- function(file, system_dir, system_label, short_label, composition) {
  meta <- parse_xvg_metadata(file)
  lines <- readLines(file, warn = FALSE)
  numeric_lines <- lines[!stringr::str_detect(lines, "^\\s*[#@]") & nzchar(stringr::str_squish(lines))]
  if (length(numeric_lines) == 0) return(NULL)
  dat <- tryCatch(
    utils::read.table(text = paste(numeric_lines, collapse = "\n"), header = FALSE,
                      stringsAsFactors = FALSE, fill = TRUE, comment.char = ""),
    error = function(e) NULL
  )
  if (is.null(dat) || ncol(dat) < 2) return(NULL)
  dat <- as.data.frame(lapply(dat, function(z) suppressWarnings(as.numeric(z))))
  dat <- dat[is.finite(dat[[1]]), , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)
  n_y <- ncol(dat) - 1
  y_names <- meta$legends
  if (length(y_names) < n_y) y_names <- c(y_names, paste0("Y", seq_len(n_y - length(y_names))))
  if (length(y_names) == 0) y_names <- paste0("Y", seq_len(n_y))
  y_names <- y_names[seq_len(n_y)]
  y_names <- make.names(y_names, unique = TRUE)
  colnames(dat) <- c("x_raw", y_names)
  file_stem <- clean_file_stem(file)
  xlab <- ifelse(is.na(meta$x_label), "", meta$x_label)
  ylab <- ifelse(is.na(meta$y_label), "", meta$y_label)
  x_low <- stringr::str_to_lower(xlab)
  stem_low <- stringr::str_to_lower(file_stem)
  is_time_file <- stringr::str_detect(stem_low, "rmsd|gyrate|dist|cogdist|mindist|contacts|hbnum|sasa|area_|temperature|pressure|density|potential|ntr|nclid|size")
  is_ps <- stringr::str_detect(x_low, "ps") && !stringr::str_detect(x_low, "ns")
  is_ns <- stringr::str_detect(x_low, "ns")
  x_value <- dat$x_raw
  x_unit <- "raw"
  if (is_time_file && !stringr::str_detect(stem_low, "^rdf-")) {
    if (is_ps) {
      x_value <- dat$x_raw / 1000
      x_unit <- "ns_from_ps"
    } else if (is_ns) {
      x_value <- dat$x_raw
      x_unit <- "ns"
    } else if (max(dat$x_raw, na.rm = TRUE) > 1000) {
      x_value <- dat$x_raw / 1000
      x_unit <- "ns_inferred_from_ps"
    } else {
      x_unit <- "raw_or_ns"
    }
  }
  dat %>%
    dplyr::mutate(x_value = x_value) %>%
    tidyr::pivot_longer(cols = -c(x_raw, x_value), names_to = "series", values_to = "value") %>%
    dplyr::mutate(
      value = suppressWarnings(as.numeric(value)),
      file = basename(file),
      file_stem = file_stem,
      metric_group = metric_group(file_stem),
      metric_label = metric_label(file_stem),
      system_dir = system_dir,
      system = system_label,
      system_short = short_label,
      composition = composition,
      x_label = xlab,
      y_label = ylab,
      x_unit = x_unit,
      title = ifelse(is.na(meta$title), "", meta$title)
    ) %>%
    dplyr::filter(is.finite(x_raw), is.finite(x_value), is.finite(value)) %>%
    dplyr::relocate(system_dir, system, system_short, composition, file, file_stem,
                    metric_group, metric_label, series, x_raw, x_value, value)
}

# -----------------------------
# 5. Convert all XVG files
# -----------------------------
xvg_index <- systems %>%
  dplyr::filter(exists) %>%
  dplyr::mutate(xvg_files = purrr::map(ana_dir, ~list.files(.x, pattern = "\\.xvg$", full.names = TRUE, recursive = FALSE))) %>%
  tidyr::unnest(xvg_files) %>%
  dplyr::mutate(
    file = basename(xvg_files),
    file_stem = clean_file_stem(file),
    metric_group = metric_group(file_stem),
    metric_label = metric_label(file_stem)
  )

if (nrow(xvg_index) == 0) stop("No XVG files were found in the selected system folders.")

all_list <- list()
for (i in seq_len(nrow(xvg_index))) {
  row <- xvg_index[i, ]
  dat <- read_xvg_long(row$xvg_files, row$system_dir, row$system_label, row$short_label, row$composition)
  if (!is.null(dat) && nrow(dat) > 0) {
    sys_conv_dir <- file.path(converted_dir, row$system_dir)
    dir.create(sys_conv_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(dat, file.path(sys_conv_dir, paste0(row$file_stem, ".csv")))
    all_list[[length(all_list) + 1]] <- dat
  }
}

all_xvg <- dplyr::bind_rows(all_list)
readr::write_csv(xvg_index, file.path(summary_dir, "xvg_file_index.csv"))
readr::write_csv(all_xvg, file.path(summary_dir, "all_xvg_long.csv"))

metadata <- all_xvg %>%
  dplyr::group_by(system_dir, system, system_short, file, file_stem, metric_group, metric_label, x_label, y_label, x_unit) %>%
  dplyr::summarise(
    n_points = dplyr::n(),
    x_min = min(x_value, na.rm = TRUE),
    x_max = max(x_value, na.rm = TRUE),
    value_min = min(value, na.rm = TRUE),
    value_max = max(value, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(metadata, file.path(summary_dir, "xvg_metadata_and_value_ranges.csv"))

message("Parsed XVG files: ", dplyr::n_distinct(paste(all_xvg$system_dir, all_xvg$file_stem)))

# -----------------------------
# 6. Summary helpers
# -----------------------------
get_xvg_exact <- function(dict, systems_keep = NULL, value_min = NULL, value_max = NULL) {
  if (nrow(dict) == 0) return(tibble::tibble())
  df <- all_xvg %>% dplyr::inner_join(dict, by = c("system_dir", "file_stem"))
  if (!is.null(systems_keep)) df <- df %>% dplyr::filter(system_dir %in% systems_keep)
  if (!is.null(value_min)) df <- df %>% dplyr::filter(value >= value_min)
  if (!is.null(value_max)) df <- df %>% dplyr::filter(value <= value_max)
  df
}

get_metric <- function(groups = NULL, pattern = NULL, systems_keep = NULL) {
  df <- all_xvg
  if (!is.null(groups)) df <- df %>% dplyr::filter(metric_group %in% groups)
  if (!is.null(pattern)) df <- df %>% dplyr::filter(stringr::str_detect(file_stem, pattern))
  if (!is.null(systems_keep)) df <- df %>% dplyr::filter(system_dir %in% systems_keep)
  df
}

final_window_summary <- function(df, window_ns = FINAL_WINDOW_NS, nonnegative = TRUE, group_cols = NULL) {
  if (nrow(df) == 0) return(tibble::tibble())
  out <- df %>% dplyr::filter(is.finite(x_value), is.finite(value))
  if (nonnegative) out <- out %>% dplyr::filter(value >= 0)
  default_cols <- c("system_dir", "system", "system_short", "composition", "file_stem", "metric_group", "metric_label", "series")
  group_cols <- group_cols %||% intersect(default_cols, names(out))
  out %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(
      max_x = max(x_value, na.rm = TRUE),
      min_keep = ifelse(max_x > window_ns, max_x - window_ns, stats::quantile(x_value, 0.75, na.rm = TRUE))
    ) %>%
    dplyr::filter(x_value >= min_keep) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE),
      median = stats::median(value, na.rm = TRUE),
      n = sum(is.finite(value)),
      final_window_start = first(min_keep),
      final_window_end = first(max_x),
      .groups = "drop"
    ) %>%
    dplyr::mutate(sd = ifelse(is.na(sd), 0, sd))
}

smooth_time_series <- function(df, window_ns = 3) {
  if (nrow(df) == 0) return(df)
  df %>%
    dplyr::arrange(system_dir, file_stem, series, x_value) %>%
    dplyr::group_by(system_dir, system_short, file_stem, metric_label, series) %>%
    dplyr::mutate(
      median_dx = suppressWarnings(stats::median(diff(sort(unique(x_value))), na.rm = TRUE)),
      smooth_n = dplyr::case_when(
        is.finite(median_dx) & median_dx > 0 ~ max(3L, as.integer(round(window_ns / median_dx))),
        TRUE ~ 11L
      ),
      smooth_n = ifelse(smooth_n %% 2 == 0, smooth_n + 1L, smooth_n),
      value_smooth = as.numeric(stats::filter(value, rep(1 / dplyr::first(smooth_n), dplyr::first(smooth_n)), sides = 2)),
      value_smooth = ifelse(is.finite(value_smooth), value_smooth, value)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-median_dx, -smooth_n)
}

plot_time_series_faceted <- function(df, title, y_lab = NULL, file_stem, width = 8.5, height = 6.0, free_y = TRUE,
                                     smooth = FALSE, smooth_window_ns = 3,
                                     show_raw_when_smooth = TRUE) {
  if (nrow(df) == 0) return(NULL)
  y_lab <- y_lab %||% unique(df$y_label)[1] %||% "Value"
  group_aes <- ggplot2::aes(x = x_value, colour = metric_label, group = interaction(file_stem, series))
  if (smooth) {
    df_plot <- smooth_time_series(df, window_ns = smooth_window_ns)
    p <- ggplot2::ggplot(df_plot, group_aes)
    if (show_raw_when_smooth) {
      p <- p + ggplot2::geom_line(ggplot2::aes(y = value), linewidth = 0.18, alpha = 0.18)
    }
    p <- p + ggplot2::geom_line(ggplot2::aes(y = value_smooth), linewidth = 0.55, alpha = 0.95)
  } else {
    df_plot <- df
    p <- ggplot2::ggplot(df_plot, group_aes) +
      ggplot2::geom_line(ggplot2::aes(y = value), linewidth = 0.32, alpha = 0.9)
  }
  p <- p +
    ggplot2::facet_wrap(~system_short, scales = ifelse(free_y, "free_y", "fixed"), ncol = 1) +
    ggplot2::labs(title = title, x = "Time (ns)", y = y_lab) +
    ggplot2::theme(legend.position = "right") +
    theme_pub(8)
  save_plot(p, file_stem, width, height)
  p
}

plot_curated_bar <- function(summary_df, title, y_lab, file_stem, width = 8.0, height = 5.0,
                             digits = 2, facet = TRUE, fill_by = "system_short") {
  if (nrow(summary_df) == 0) return(NULL)
  df <- summary_df %>%
    mean_label_df(digits = digits) %>%
    dplyr::mutate(
      display = factor(display, levels = rev(unique(display))),
      system_short = factor(system_short, levels = system_levels),
      class = ifelse(is.na(class), interaction_class(display), class)
    )
  aes_fill <- if (fill_by == "class") ggplot2::aes(fill = class) else ggplot2::aes(fill = system_short)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = display, y = mean)) +
    ggplot2::geom_col(aes_fill, width = 0.62, colour = "black", linewidth = 0.22) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(mean - sd, 0), ymax = mean + sd), width = 0.16, linewidth = 0.22) +
    ggplot2::geom_text(ggplot2::aes(y = label_y, label = mean_label), hjust = 0, size = 2.25) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(expand = bar_expand()) +
    ggplot2::labs(title = title, x = NULL, y = y_lab)
  if (fill_by == "class") {
    p <- p + ggplot2::scale_fill_manual(values = interaction_colors, drop = FALSE)
  } else {
    p <- p + ggplot2::scale_fill_manual(values = system_colors, drop = FALSE)
  }
  if (facet) p <- p + ggplot2::facet_wrap(~system_short, scales = "free_y", ncol = 1)
  p <- p + theme_pub(8) + ggplot2::theme(legend.position = "top")
  save_plot(p, file_stem, width, height)
  p
}

# -----------------------------
# 7. Individual plot for each XVG file
# -----------------------------
if (MAKE_INDIVIDUAL_XVG_PLOTS) {
  individual_keys <- all_xvg %>%
    dplyr::distinct(system_dir, system_short, file_stem, metric_group, metric_label, x_label, y_label, x_unit)
  for (i in seq_len(nrow(individual_keys))) {
    key <- individual_keys[i, ]
    df_i <- all_xvg %>% dplyr::filter(system_dir == key$system_dir, file_stem == key$file_stem)
    if (nrow(df_i) == 0) next
    xlab <- ifelse(stringr::str_detect(key$x_unit, "ns"), "Time (ns)", ifelse(nzchar(key$x_label), key$x_label, "x"))
    ylab <- ifelse(nzchar(key$y_label), key$y_label, "Value")
    p_i <- ggplot2::ggplot(df_i, ggplot2::aes(x = x_value, y = value, colour = series)) +
      ggplot2::geom_line(linewidth = 0.35, alpha = 0.9) +
      ggplot2::labs(title = paste0(key$system_short, ": ", key$metric_label), subtitle = key$file_stem, x = xlab, y = ylab) +
      theme_pub(8) + ggplot2::theme(legend.position = ifelse(dplyr::n_distinct(df_i$series) > 1, "right", "none"))
    save_plot(p_i, paste0(safe_name(key$system_dir), "__", safe_name(key$file_stem)), 5.4, 3.8, folder = fig_individual_dir)
  }
}

# -----------------------------
# 8. Full diagnostic grouped figures and curated summary figures
# -----------------------------
plot_registry <- list()

qc_df <- get_metric(groups = "QC")
if (nrow(qc_df) > 0) {
  p_qc <- ggplot2::ggplot(qc_df, ggplot2::aes(x = x_value, y = value, colour = system_short)) +
    ggplot2::geom_line(linewidth = 0.30, alpha = 0.9) +
    ggplot2::facet_wrap(~metric_label, scales = "free_y", ncol = 2) +
    ggplot2::scale_colour_manual(values = system_colors, drop = FALSE) +
    ggplot2::labs(title = "MD quality-control readouts", x = "Time (ns)", y = NULL) +
    theme_pub(8) + ggplot2::theme(legend.position = "top")
  save_plot(p_qc, "MD_quality_control", 7.4, 5.2)
  plot_registry[["QC"]] <- p_qc
  readr::write_csv(qc_df, file.path(summary_dir, "MD_quality_control_data.csv"))
}

rmsd_df <- get_metric(groups = "RMSD") %>% dplyr::filter(value >= 0, value < 1000)
if (nrow(rmsd_df) > 0) {
  p_rmsd <- plot_time_series_faceted(
    rmsd_df,
    "RMSD time-series",
    "RMSD (nm)",
    "RMSD_time_series_smoothed",
    8.0,
    5.8,
    smooth = SMOOTH_RMSD_DISPLAY,
    smooth_window_ns = RMSD_SMOOTH_WINDOW_NS,
    show_raw_when_smooth = TRUE
  )
  p_rmsd_raw <- plot_time_series_faceted(
    rmsd_df,
    "RMSD time-series, raw traces",
    "RMSD (nm)",
    "RMSD_time_series_raw",
    8.0,
    5.8,
    smooth = FALSE
  )
  rmsd_summary_all <- final_window_summary(rmsd_df)
  readr::write_csv(rmsd_summary_all, file.path(summary_dir, "RMSD_final_window_summary_all.csv"))
}

rg_dict <- tibble::tribble(
  ~system_dir,             ~file_stem,                ~display,                         ~class,
  "pcldtxsic18",           "druggyrate",              "DTX prodrug Rg",                 "compactness",
  "pcldtxsic18",           "polymergyrate",           "polymer Rg",                     "polymer integration",
  "pcldtxsic18",           "allgyrate",               "whole assembly Rg",              "polymer integration",
  "pclabisic18",           "druggyrate",              "Abi prodrug Rg",                 "compactness",
  "pclabisic18",           "polymergyrate",           "polymer Rg",                     "polymer integration",
  "pclabisic18",           "allgyrate",               "whole assembly Rg",              "polymer integration",
  "pcldtxsic18abisic18",   "gyrate_dtxsic18",         "DTX prodrug Rg",                 "compactness",
  "pcldtxsic18abisic18",   "gyrate_abisic18",         "Abi prodrug Rg",                 "compactness",
  "pcldtxsic18abisic18",   "gyrate_total_drug",       "total prodrug Rg",               "compactness",
  "pcldtxsic18abisic18",   "gyrate_polymer",          "polymer Rg",                     "polymer integration",
  "pcldtxsic18abisic18",   "gyrate_polymer_drug",     "whole assembly Rg",              "polymer integration"
)
rg_df <- get_xvg_exact(rg_dict, value_min = 0, value_max = 1000)
if (nrow(rg_df) > 0) {
  p_rg_time <- plot_time_series_faceted(rg_df, "Curated Rg time-series", "Rg (nm)", "Rg_time_series_curated", 8.0, 5.6)
  rg_summary <- final_window_summary(rg_df, group_cols = c("system_dir", "system", "system_short", "composition", "file_stem", "display", "class")) %>%
    dplyr::mutate(display = factor(display, levels = unique(rg_dict$display)))
  readr::write_csv(rg_summary, file.path(summary_dir, "Rg_final_window_summary_curated.csv"))
  p_rg_bar <- plot_curated_bar(rg_summary, "Assembly compactness from final-window Rg", "Rg (nm), mean ± SD", "Rg_final_window_summary", 7.2, 5.4, digits = 2, facet = TRUE, fill_by = "class")
  plot_registry[["Rg"]] <- p_rg_bar
}

com_dict <- tibble::tribble(
  ~system_dir,             ~file_stem,              ~display,                    ~class,
  "pcldtxsic18",           "drug-polymer-dist",     "DTX-Pol distance",          "polymer integration",
  "pcldtxsic18",           "drug-pcl-dist",         "DTX-PCL distance",          "PCL anchoring",
  "pcldtxsic18",           "drug-peg-dist",         "DTX-PEG distance",          "PEG contact",
  "pclabisic18",           "drug-polymer-dist",     "Abi-Pol distance",          "polymer integration",
  "pclabisic18",           "drug-pcl-dist",         "Abi-PCL distance",          "PCL anchoring",
  "pclabisic18",           "drug-peg-dist",         "Abi-PEG distance",          "PEG contact",
  "pcldtxsic18abisic18",   "drug-polymer-dist",     "total Pro-Pol distance",    "polymer integration",
  "pcldtxsic18abisic18",   "dtx-polymer-dist",      "DTX-Pol distance",          "polymer integration",
  "pcldtxsic18abisic18",   "abi-polymer-dist",      "Abi-Pol distance",          "polymer integration",
  "pcldtxsic18abisic18",   "dtx-pcl-dist",          "DTX-PCL distance",          "PCL anchoring",
  "pcldtxsic18abisic18",   "abi-pcl-dist",          "Abi-PCL distance",          "PCL anchoring",
  "pcldtxsic18abisic18",   "dtx-peg-dist",          "DTX-PEG distance",          "PEG contact",
  "pcldtxsic18abisic18",   "abi-peg-dist",          "Abi-PEG distance",          "PEG contact",
  "pcldtxsic18abisic18",   "dtx-abi-dist",          "DTX-Abi distance",          "DTX-Abi association"
)
com_df_all <- get_metric(groups = "COM distance") %>% dplyr::filter(value >= 0, value < 1000)
if (nrow(com_df_all) > 0) {
  com_summary_all <- final_window_summary(com_df_all)
  readr::write_csv(com_summary_all, file.path(summary_dir, "COM_distance_final_window_summary_all.csv"))
}
com_df <- get_xvg_exact(com_dict, value_min = 0, value_max = 1000)
if (nrow(com_df) > 0) {
  p_com_time <- plot_time_series_faceted(com_df, "Key COM distance time-series", "COM distance (nm)", "COM_distance_time_series_key", 8.2, 5.8)
  com_summary <- final_window_summary(com_df, group_cols = c("system_dir", "system", "system_short", "composition", "file_stem", "display", "class"))
  readr::write_csv(com_summary, file.path(summary_dir, "COM_distance_final_window_summary_key.csv"))
  p_com_bar <- plot_curated_bar(com_summary, "Final-window key COM distances", "Distance (nm), mean ± SD", "COM_distance_final_window_summary", 7.6, 6.2, digits = 2, facet = TRUE, fill_by = "class")
  plot_registry[["COM"]] <- p_com_bar
}

cog_df <- get_metric(groups = "COG distance") %>% dplyr::filter(value >= 0, value < 1000)
if (nrow(cog_df) > 0) {
  p_cog <- plot_time_series_faceted(cog_df, "COG distance time-series", "COG distance (nm)", "COG_distance_time_series", 8.3, 6.2)
  cog_summary <- final_window_summary(cog_df)
  readr::write_csv(cog_summary, file.path(summary_dir, "COG_distance_final_window_summary_all.csv"))
}

mindist_df <- get_metric(groups = "Minimum distance") %>% dplyr::filter(value >= 0, value < 1000)
if (nrow(mindist_df) > 0) {
  p_mindist <- plot_time_series_faceted(mindist_df, "Minimum-distance time-series", "Minimum distance (nm)", "minimum_distance_time_series", 8.2, 6.0)
  mindist_summary <- final_window_summary(mindist_df)
  readr::write_csv(mindist_summary, file.path(summary_dir, "minimum_distance_final_window_summary_all.csv"))
}

contact_dict <- tibble::tribble(
  ~system_dir,             ~file_stem,                 ~display,             ~class,
  "pcldtxsic18abisic18",   "dtx-polymer-contacts",     "DTX-Pol contacts",   "polymer integration",
  "pcldtxsic18abisic18",   "abi-polymer-contacts",     "Abi-Pol contacts",   "polymer integration",
  "pcldtxsic18abisic18",   "dtx-pcl-contacts",         "DTX-PCL contacts",   "PCL anchoring",
  "pcldtxsic18abisic18",   "abi-pcl-contacts",         "Abi-PCL contacts",   "PCL anchoring",
  "pcldtxsic18abisic18",   "dtx-peg-contacts",         "DTX-PEG contacts",   "PEG contact",
  "pcldtxsic18abisic18",   "abi-peg-contacts",         "Abi-PEG contacts",   "PEG contact",
  "pcldtxsic18abisic18",   "dtx-abi-contacts",         "DTX-Abi contacts",   "DTX-Abi association",
  "pcldtxsic18abisic18",   "dtxlig-pcl-contacts",      "DTX Mod-PCL contacts", "PCL anchoring",
  "pcldtxsic18abisic18",   "abilig-pcl-contacts",      "Abi Mod-PCL contacts", "PCL anchoring"
)
contact_df_all <- get_metric(groups = "Contact number") %>% dplyr::filter(value >= 0, value < 1e9)
if (nrow(contact_df_all) > 0) {
  contact_summary_all <- final_window_summary(contact_df_all)
  readr::write_csv(contact_summary_all, file.path(summary_dir, "contact_number_final_window_summary_all.csv"))
}
contact_df <- get_xvg_exact(contact_dict, value_min = 0, value_max = 1e9)
if (nrow(contact_df) > 0) {
  contact_summary <- final_window_summary(contact_df, group_cols = c("system_dir", "system", "system_short", "composition", "file_stem", "display", "class")) %>%
    dplyr::mutate(mean_k = mean / 1000, sd_k = sd / 1000) %>%
    dplyr::rename(mean_raw = mean, sd_raw = sd) %>%
    dplyr::mutate(mean = mean_k, sd = sd_k)
  readr::write_csv(contact_summary, file.path(summary_dir, "contact_number_final_window_summary_key.csv"))
  p_contact_bar <- plot_curated_bar(contact_summary, "Final-window key near-contact numbers", "Contact number (×10³), mean ± SD", "contact_number_final_window_summary", 6.8, 4.6, digits = 1, facet = FALSE, fill_by = "class")
  plot_registry[["Contact"]] <- p_contact_bar
}

hbond_dict <- tibble::tribble(
  ~system_dir,             ~file_stem,              ~display,                  ~class,
  "pcldtxsic18",           "hbnum-dru-po",          "DTX-Pol H-bonds",         "polymer integration",
  "pcldtxsic18",           "hbnum-dru-pcl",         "DTX-PCL H-bonds",         "PCL anchoring",
  "pcldtxsic18",           "hbnum-dru-peg",         "DTX-PEG H-bonds",         "PEG contact",
  "pclabisic18",           "hbnum-dru-po",          "Abi-Pol H-bonds",         "polymer integration",
  "pclabisic18",           "hbnum-dru-pcl",         "Abi-PCL H-bonds",         "PCL anchoring",
  "pclabisic18",           "hbnum-dru-peg",         "Abi-PEG H-bonds",         "PEG contact",
  "pcldtxsic18abisic18",   "hbnum-drug-polymer",    "total Pro-Pol H-bonds",   "polymer integration",
  "pcldtxsic18abisic18",   "hbnum-dtx-polymer",     "DTX-Pol H-bonds",         "polymer integration",
  "pcldtxsic18abisic18",   "hbnum-abi-polymer",     "Abi-Pol H-bonds",         "polymer integration",
  "pcldtxsic18abisic18",   "hbnum-dtx-pcl",         "DTX-PCL H-bonds",         "PCL anchoring",
  "pcldtxsic18abisic18",   "hbnum-abi-pcl",         "Abi-PCL H-bonds",         "PCL anchoring",
  "pcldtxsic18abisic18",   "hbnum-dtx-peg",         "DTX-PEG H-bonds",         "PEG contact",
  "pcldtxsic18abisic18",   "hbnum-abi-peg",         "Abi-PEG H-bonds",         "PEG contact",
  "pcldtxsic18abisic18",   "hbnum-dtx-abi",         "DTX-Abi H-bonds",         "DTX-Abi association"
)
hbond_df_all <- get_metric(groups = "H-bond") %>% dplyr::filter(value >= 0, value < 1e6)
if (nrow(hbond_df_all) > 0) {
  hbond_summary_all <- final_window_summary(hbond_df_all)
  readr::write_csv(hbond_summary_all, file.path(summary_dir, "Hbond_final_window_summary_all.csv"))
}
hbond_df <- get_xvg_exact(hbond_dict, value_min = 0, value_max = 1e6)
if (nrow(hbond_df) > 0) {
  p_hb_time <- plot_time_series_faceted(hbond_df, "Key hydrogen-bond time-series", "Number of H-bonds", "Hbond_time_series_key", 8.2, 5.8)
  hbond_summary <- final_window_summary(hbond_df, group_cols = c("system_dir", "system", "system_short", "composition", "file_stem", "display", "class"))
  readr::write_csv(hbond_summary, file.path(summary_dir, "Hbond_final_window_summary_key.csv"))
  p_hb_bar <- plot_curated_bar(hbond_summary, "Final-window key hydrogen bonds", "H-bonds, mean ± SD", "Hbond_final_window_summary", 7.6, 6.0, digits = 2, facet = TRUE, fill_by = "class")
  plot_registry[["Hbond"]] <- p_hb_bar
}

rdf_dict <- tibble::tribble(
  ~system_dir,             ~file_stem,        ~display,             ~class,
  "pcldtxsic18abisic18",   "rdf-dtx-pcl",    "DTX-PCL",           "PCL anchoring",
  "pcldtxsic18abisic18",   "rdf-abi-pcl",    "Abi-PCL",           "PCL anchoring",
  "pcldtxsic18abisic18",   "rdf-dtx-peg",    "DTX-PEG",           "PEG contact",
  "pcldtxsic18abisic18",   "rdf-abi-peg",    "Abi-PEG",           "PEG contact",
  "pcldtxsic18abisic18",   "rdf-dtx-abi",    "DTX-Abi",           "DTX-Abi association",
  "pcldtxsic18abisic18",   "rdf-drug-pcl",   "total Pro-PCL",     "PCL anchoring",
  "pcldtxsic18abisic18",   "rdf-drug-peg",   "total Pro-PEG",     "PEG contact"
)
rdf_df_all <- get_metric(groups = "RDF") %>% dplyr::filter(value >= 0, value < 1e9)
if (nrow(rdf_df_all) > 0) {
  readr::write_csv(rdf_df_all, file.path(summary_dir, "RDF_curve_data_all.csv"))
}
rdf_df <- get_xvg_exact(rdf_dict, value_min = 0, value_max = 1e9)
if (nrow(rdf_df) > 0) {
  rdf_df <- rdf_df %>% dplyr::mutate(display = factor(display, levels = rdf_dict$display))
  p_rdf <- ggplot2::ggplot(rdf_df, ggplot2::aes(x = x_raw, y = value, colour = class)) +
    ggplot2::geom_line(linewidth = 0.42, alpha = 0.95) +
    ggplot2::facet_wrap(~display, scales = "fixed", ncol = 4) +
    ggplot2::coord_cartesian(ylim = c(0, 9), expand = FALSE) +
    ggplot2::scale_y_continuous(breaks = c(0, 3, 6, 9)) +
    ggplot2::scale_colour_manual(values = interaction_colors, drop = FALSE) +
    ggplot2::labs(title = "Radial distribution functions", x = "r (nm)", y = "g(r)") +
    theme_pub(8) + ggplot2::theme(legend.position = "none")
  save_plot(p_rdf, "RDF_curves_fixed_y_0_9", 8.8, 4.8)
  save_plot(p_rdf, "RDF_curves", 8.8, 4.8)
  rdf_peaks <- rdf_df %>%
    dplyr::group_by(system_dir, system_short, file_stem, display, class) %>%
    dplyr::summarise(peak_g = max(value, na.rm = TRUE), peak_r_nm = x_raw[which.max(value)][1], .groups = "drop")
  readr::write_csv(rdf_df, file.path(summary_dir, "RDF_curve_data_key.csv"))
  readr::write_csv(rdf_peaks, file.path(summary_dir, "RDF_peak_summary_key.csv"))
  plot_registry[["RDF"]] <- p_rdf
}

sasa_df <- get_metric(groups = "SASA") %>% dplyr::filter(value >= 0, value < 1e9)
if (nrow(sasa_df) > 0) {
  p_sasa_time <- plot_time_series_faceted(sasa_df, "Raw SASA time-series", "SASA / area (nm²)", "raw_SASA_time_series", 8.2, 6.2)
  sasa_raw_summary <- final_window_summary(sasa_df)
  readr::write_csv(sasa_raw_summary, file.path(summary_dir, "raw_SASA_final_window_summary_all.csv"))
}

cluster_df <- get_metric(groups = "Cluster") %>% dplyr::filter(value >= 0, value < 1e9)
if (nrow(cluster_df) > 0) {
  p_cluster <- ggplot2::ggplot(cluster_df, ggplot2::aes(x = x_raw, y = value, colour = metric_label)) +
    ggplot2::geom_line(linewidth = 0.35, alpha = 0.9) +
    ggplot2::facet_wrap(~system_short, scales = "free_y", ncol = 1) +
    ggplot2::labs(title = "Cluster-size and clustering diagnostics", x = "Raw x value", y = "Cluster output value") +
    theme_pub(8) + ggplot2::theme(legend.position = "right")
  save_plot(p_cluster, "cluster_size_and_clustering_diagnostics", 8.0, 5.8)
  readr::write_csv(cluster_df, file.path(summary_dir, "cluster_size_and_clustering_data.csv"))
}

# -----------------------------
# 9. SASA-derived contact area, using explicit matched files
# -----------------------------
read_component_sasa <- function(system_id, stem) {
  all_xvg %>%
    dplyr::filter(system_dir == system_id, file_stem == stem) %>%
    dplyr::arrange(x_value) %>%
    dplyr::group_by(file_stem) %>%
    dplyr::mutate(frame = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(frame, x_value, value)
}

calc_contact_area <- function(system_id, pair, a_stem, b_stem, ab_stem) {
  a <- read_component_sasa(system_id, a_stem)
  b <- read_component_sasa(system_id, b_stem)
  ab <- read_component_sasa(system_id, ab_stem)
  if (nrow(a) == 0 || nrow(b) == 0 || nrow(ab) == 0) return(tibble::tibble())
  n <- min(nrow(a), nrow(b), nrow(ab))
  a <- a[seq_len(n), ]; b <- b[seq_len(n), ]; ab <- ab[seq_len(n), ]
  sys_row <- systems %>% dplyr::filter(system_dir == system_id) %>% dplyr::slice(1)
  tibble::tibble(
    system_dir = system_id,
    system = sys_row$system_label,
    system_short = sys_row$short_label,
    composition = sys_row$composition,
    pair = pair,
    x_value = a$x_value,
    contact_area_nm2_raw = (a$value + b$value - ab$value) / 2,
    file_a = a_stem,
    file_b = b_stem,
    file_ab = ab_stem
  ) %>%
    dplyr::mutate(contact_area_nm2 = ifelse(contact_area_nm2_raw < 0 & abs(contact_area_nm2_raw) < 1e-6, 0, contact_area_nm2_raw)) %>%
    dplyr::filter(is.finite(contact_area_nm2))
}

sasa_pairs <- tibble::tribble(
  ~system_dir,             ~pair,                  ~a,             ~b,              ~ab,                  ~class,
  "pcldtxsic18",           "DTX-Pol",              "area_dru",    "area_pol",     "area_whole",        "polymer integration",
  "pclabisic18",           "Abi-Pol",              "area_dru",    "area_pol",     "area_whole",        "polymer integration",
  "pcldtxsic18abisic18",   "total Pro-Pol",        "sasa-drug",   "sasa-polymer", "sasa-drug-polymer", "polymer integration",
  "pcldtxsic18abisic18",   "DTX-Pol",              "sasa-dtx",    "sasa-polymer", "sasa-dtx-polymer",  "polymer integration",
  "pcldtxsic18abisic18",   "Abi-Pol",              "sasa-abi",    "sasa-polymer", "sasa-abi-polymer",  "polymer integration",
  "pcldtxsic18abisic18",   "DTX-Abi",              "sasa-dtx",    "sasa-abi",     "sasa-dtx-abi",      "DTX-Abi association",
  "pcldtxsic18abisic18",   "DTX-PCL",              "sasa-dtx",    "sasa-pcl",     "sasa-dtx-pcl",      "PCL anchoring",
  "pcldtxsic18abisic18",   "Abi-PCL",              "sasa-abi",    "sasa-pcl",     "sasa-abi-pcl",      "PCL anchoring",
  "pcldtxsic18abisic18",   "DTX-PEG",              "sasa-dtx",    "sasa-peg",     "sasa-dtx-peg",      "PEG contact",
  "pcldtxsic18abisic18",   "Abi-PEG",              "sasa-abi",    "sasa-peg",     "sasa-abi-peg",      "PEG contact"
)

sasa_contact <- purrr::pmap_dfr(sasa_pairs[, c("system_dir", "pair", "a", "b", "ab")], ~calc_contact_area(..1, ..2, ..3, ..4, ..5))
if (nrow(sasa_contact) > 0) {
  sasa_contact <- sasa_contact %>% dplyr::left_join(sasa_pairs %>% dplyr::select(system_dir, pair, class), by = c("system_dir", "pair"))
  readr::write_csv(sasa_contact, file.path(summary_dir, "SASA_contact_area_time_series_data.csv"))
  sasa_contact_plot <- sasa_contact %>% dplyr::filter(contact_area_nm2 >= 0)
  p_sasa_contact_time <- ggplot2::ggplot(sasa_contact_plot, ggplot2::aes(x = x_value, y = contact_area_nm2, colour = pair)) +
    ggplot2::geom_line(linewidth = 0.32, alpha = 0.9) +
    ggplot2::facet_wrap(~system_short, scales = "free_y", ncol = 1) +
    ggplot2::labs(title = "SASA-derived interfacial contact area", x = "Time (ns)", y = "Contact area (nm²)") +
    theme_pub(8) + ggplot2::theme(legend.position = "right")
  save_plot(p_sasa_contact_time, "SASA_contact_area_time_series", 8.0, 5.8)
  
  sasa_contact_summary <- sasa_contact_plot %>%
    dplyr::group_by(system_dir, system_short, pair, class) %>%
    dplyr::mutate(max_x = max(x_value, na.rm = TRUE), min_keep = ifelse(max_x > FINAL_WINDOW_NS, max_x - FINAL_WINDOW_NS, stats::quantile(x_value, 0.75, na.rm = TRUE))) %>%
    dplyr::filter(x_value >= min_keep) %>%
    dplyr::summarise(mean = mean(contact_area_nm2, na.rm = TRUE), sd = sd(contact_area_nm2, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(sd = ifelse(is.na(sd), 0, sd), display = pair)
  readr::write_csv(sasa_contact_summary, file.path(summary_dir, "SASA_contact_area_final_window_summary.csv"))
  p_sasa_contact_bar <- plot_curated_bar(sasa_contact_summary, "Final-window SASA-derived contact area", "Contact area (nm²), mean ± SD", "SASA_contact_area_final_window_summary", 7.6, 5.8, digits = 1, facet = TRUE, fill_by = "class")
  plot_registry[["SASAcontact"]] <- p_sasa_contact_bar
}

# -----------------------------
# 10. Segment-resolved heatmaps using explicit dictionaries
# -----------------------------
segment_order <- c(
  "Pro-Pol", "Par-Pol", "Mod-Pol",
  "Pro-PCL", "Par-PCL", "Mod-PCL",
  "Pro-PEG", "Par-PEG", "Mod-PEG",
  "DTX-Abi", "Par-Par", "Mod-Mod", "Par-Mod", "Mod-Par"
)

segment_dist_dict <- tibble::tribble(
  ~system_dir,             ~file_stem,                 ~row,                 ~segment_pair,
  "pcldtxsic18abisic18",   "drug-polymer-dist",        "total prodrugs",     "Pro-Pol",
  "pcldtxsic18abisic18",   "dtx-polymer-dist",         "DTX-SI-C18",         "Pro-Pol",
  "pcldtxsic18abisic18",   "abi-polymer-dist",         "Abi-SI-C18",         "Pro-Pol",
  "pcldtxsic18abisic18",   "dtxpar-polymer-dist",      "DTX-SI-C18",         "Par-Pol",
  "pcldtxsic18abisic18",   "dtxlig-polymer-dist",      "DTX-SI-C18",         "Mod-Pol",
  "pcldtxsic18abisic18",   "abipar-polymer-dist",      "Abi-SI-C18",         "Par-Pol",
  "pcldtxsic18abisic18",   "abilig-polymer-dist",      "Abi-SI-C18",         "Mod-Pol",
  "pcldtxsic18abisic18",   "drug-pcl-dist",            "total prodrugs",     "Pro-PCL",
  "pcldtxsic18abisic18",   "dtx-pcl-dist",             "DTX-SI-C18",         "Pro-PCL",
  "pcldtxsic18abisic18",   "abi-pcl-dist",             "Abi-SI-C18",         "Pro-PCL",
  "pcldtxsic18abisic18",   "dtxpar-pcl-dist",          "DTX-SI-C18",         "Par-PCL",
  "pcldtxsic18abisic18",   "dtxlig-pcl-dist",          "DTX-SI-C18",         "Mod-PCL",
  "pcldtxsic18abisic18",   "abipar-pcl-dist",          "Abi-SI-C18",         "Par-PCL",
  "pcldtxsic18abisic18",   "abilig-pcl-dist",          "Abi-SI-C18",         "Mod-PCL",
  "pcldtxsic18abisic18",   "drug-peg-dist",            "total prodrugs",     "Pro-PEG",
  "pcldtxsic18abisic18",   "dtx-peg-dist",             "DTX-SI-C18",         "Pro-PEG",
  "pcldtxsic18abisic18",   "abi-peg-dist",             "Abi-SI-C18",         "Pro-PEG",
  "pcldtxsic18abisic18",   "dtxpar-peg-dist",          "DTX-SI-C18",         "Par-PEG",
  "pcldtxsic18abisic18",   "dtxlig-peg-dist",          "DTX-SI-C18",         "Mod-PEG",
  "pcldtxsic18abisic18",   "abipar-peg-dist",          "Abi-SI-C18",         "Par-PEG",
  "pcldtxsic18abisic18",   "abilig-peg-dist",          "Abi-SI-C18",         "Mod-PEG",
  "pcldtxsic18abisic18",   "dtx-abi-dist",             "DTX-Abi cross-pair", "DTX-Abi",
  "pcldtxsic18abisic18",   "dtxpar-abipar-dist",       "DTX-Abi cross-pair", "Par-Par",
  "pcldtxsic18abisic18",   "dtxlig-abilig-dist",       "DTX-Abi cross-pair", "Mod-Mod",
  "pcldtxsic18abisic18",   "dtxpar-abilig-dist",       "DTX-Abi cross-pair", "Par-Mod",
  "pcldtxsic18abisic18",   "dtxlig-abipar-dist",       "DTX-Abi cross-pair", "Mod-Par"
)

seg_dist_df <- get_xvg_exact(segment_dist_dict, value_min = 0, value_max = 1000)
if (nrow(seg_dist_df) > 0) {
  seg_dist_heat <- final_window_summary(seg_dist_df, group_cols = c("system_dir", "file_stem", "row", "segment_pair")) %>%
    dplyr::mutate(
      sd = ifelse(is.na(sd), 0, sd),
      segment_pair = factor(segment_pair, levels = segment_order),
      row = factor(row, levels = c("total prodrugs", "DTX-SI-C18", "Abi-SI-C18", "DTX-Abi cross-pair")),
      label_mean_sd = sprintf("%.2f\n± %.2f", mean, sd),
      label_mean_sd_inline = sprintf("%.2f ± %.2f", mean, sd)
    )
  p_seg_dist <- ggplot2::ggplot(seg_dist_heat, ggplot2::aes(x = segment_pair, y = row, fill = mean)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = label_mean_sd), size = 1.85, lineheight = 0.82) +
    ggplot2::scale_fill_gradientn(colours = heatmap_cols, name = "Distance\n(nm)", na.value = "#F5F5F5") +
    ggplot2::labs(title = "Whole-prodrug and segment-resolved localization in the mixed DTX/Abi system", x = NULL, y = NULL) +
    theme_pub(8) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "right")
  save_plot(p_seg_dist, "segment_resolved_distance_heatmap_mixed_system", 10.6, 4.2)
  plot_registry[["SegmentDist"]] <- p_seg_dist
  readr::write_csv(seg_dist_heat, file.path(summary_dir, "segment_resolved_distance_heatmap_summary.csv"))
}

segment_hb_dict <- segment_dist_dict %>%
  dplyr::mutate(file_stem = paste0("hbnum-", stringr::str_remove(file_stem, "-dist$"))) %>%
  dplyr::filter(!file_stem %in% c("hbnum-drug-pcl", "hbnum-drug-peg"))

segment_hb_dict <- segment_hb_dict %>%
  dplyr::mutate(file_stem = dplyr::case_when(
    file_stem == "hbnum-drug-polymer"~ "hbnum-drug-polymer",
    TRUE ~ file_stem
  ))

seg_hb_df <- get_xvg_exact(segment_hb_dict, value_min = 0, value_max = 1e6)
if (nrow(seg_hb_df) > 0) {
  seg_hb_heat <- final_window_summary(seg_hb_df, group_cols = c("system_dir", "file_stem", "row", "segment_pair")) %>%
    dplyr::mutate(
      sd = ifelse(is.na(sd), 0, sd),
      segment_pair = factor(segment_pair, levels = segment_order),
      row = factor(row, levels = c("total prodrugs", "DTX-SI-C18", "Abi-SI-C18", "DTX-Abi cross-pair")),
      label_mean_sd = sprintf("%.2f\n± %.2f", mean, sd),
      label_mean_sd_inline = sprintf("%.2f ± %.2f", mean, sd)
    )
  p_seg_hb <- ggplot2::ggplot(seg_hb_heat, ggplot2::aes(x = segment_pair, y = row, fill = mean)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = label_mean_sd), size = 1.85, lineheight = 0.82) +
    ggplot2::scale_fill_gradientn(colours = heatmap_cols, name = "H-bonds", na.value = "#F5F5F5") +
    ggplot2::labs(title = "Whole-prodrug and segment-resolved hydrogen bonds in the mixed DTX/Abi system", x = NULL, y = NULL) +
    theme_pub(8) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "right")
  save_plot(p_seg_hb, "segment_resolved_Hbond_heatmap_mixed_system", 10.6, 4.2)
  plot_registry[["SegmentHB"]] <- p_seg_hb
  readr::write_csv(seg_hb_heat, file.path(summary_dir, "segment_resolved_Hbond_heatmap_summary.csv"))
}

# -----------------------------
# 11. Clean integrated SI overview
# -----------------------------
panel_letter <- function(p, lab) {
  p + ggplot2::labs(tag = lab) +
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 12), plot.tag.position = c(0.01, 0.99))
}

p_system_text <- ggplot2::ggplot() +
  ggplot2::annotate("text", x = 0, y = 1.0, hjust = 0, vjust = 1, size = 3.0, fontface = "bold", label = "MD system definitions") +
  ggplot2::annotate("text", x = 0, y = 0.78, hjust = 0, vjust = 1, size = 2.55,
                    label = paste0(
                      "DTX single: 30 mPEG2k-PCL2k + 20 DTX-SI-C18\n",
                      "Abi single: 30 mPEG2k-PCL2k + 20 Abi-SI-C18\n",
                      "DTX+Abi mixed: 30 mPEG2k-PCL2k + 10 DTX-SI-C18 + 10 Abi-SI-C18"
                    )) +
  ggplot2::annotate("text", x = 0, y = 0.28, hjust = 0, vjust = 1, size = 2.25,
                    label = paste0(
                      "Final-window summaries use the last ", FINAL_WINDOW_NS,
                      "ns when available; otherwise the final quartile is used.\n",
                      "Detailed all-file plots and converted CSV files are provided separately."
                    )) +
  ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
  ggplot2::theme_void(base_family = BASE_FAMILY)
p_system_text <- panel_letter(p_system_text, "a")

selected_panels <- list()
if (exists("p_rg_bar")) selected_panels[["b"]] <- panel_letter(p_rg_bar + ggplot2::theme(legend.position = "none"), "b")
if (exists("p_sasa_contact_bar")) selected_panels[["c"]] <- panel_letter(p_sasa_contact_bar + ggplot2::theme(legend.position = "none"), "c")
if (exists("p_com_bar")) selected_panels[["d"]] <- panel_letter(p_com_bar + ggplot2::theme(legend.position = "none"), "d")
if (exists("p_rdf")) selected_panels[["e"]] <- panel_letter(p_rdf, "e")
if (exists("p_hb_bar")) selected_panels[["f"]] <- panel_letter(p_hb_bar + ggplot2::theme(legend.position = "none"), "f")
if (exists("p_contact_bar")) selected_panels[["g"]] <- panel_letter(p_contact_bar + ggplot2::theme(legend.position = "none"), "g")
if (exists("p_seg_dist")) selected_panels[["h"]] <- panel_letter(p_seg_dist, "h")
if (exists("p_seg_hb")) selected_panels[["i"]] <- panel_letter(p_seg_hb, "i")

if (length(selected_panels) > 0) {
  overview <- patchwork::wrap_plots(c(list(p_system_text), selected_panels), ncol = 3, guides = "collect") +
    patchwork::plot_annotation(
      title = "Supplementary MD support for DTX-SI-C18/Abi-SI-C18 mixed co-assembly",
      subtitle = "Curated readouts compare single-prodrug controls with the mixed formulation and resolve compactness, prodrug-carrier contact, PCL/PEG localization, DTX-Abi association and segment-level interactions.",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
                             plot.subtitle = ggplot2::element_text(size = 9, hjust = 0))
    )
  save_plot(overview, "MD_support_overview", 18.0, 11.0)
}

qc_panels <- list()
if (exists("p_qc")) qc_panels[["QC"]] <- p_qc
if (exists("p_rmsd")) qc_panels[["RMSD"]] <- p_rmsd
if (exists("p_rg_time")) qc_panels[["Rg"]] <- p_rg_time
if (exists("p_com_time")) qc_panels[["COM"]] <- p_com_time
if (exists("p_hb_time")) qc_panels[["Hbond"]] <- p_hb_time
if (exists("p_sasa_time")) qc_panels[["SASA"]] <- p_sasa_time
if (length(qc_panels) > 1) {
  diag_overview <- patchwork::wrap_plots(qc_panels, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = "MD trajectory diagnostics and auxiliary time-series",
      subtitle = "These panels support trajectory stability and transparent inspection of converted XVG outputs.",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
                             plot.subtitle = ggplot2::element_text(size = 9, hjust = 0))
    )
  save_plot(diag_overview, "MD_diagnostics_time_series", 15.5, 13.0, folder = diagnostic_dir)
}

# -----------------------------
# 12. Figure log and run log
# -----------------------------
figure_log <- c(
  "DTX-SI-C18/Abi-SI-C18 co-assembly MD output guide",
  "============================================================",
  paste0("Input root: ", root_dir),
  paste0("Figure folder: ", fig_dir),
  "",
  "System definitions:",
  "DTX single: 30 mPEG2k-PCL2k + 20 DTX-SI-C18.",
  "Abi single: 30 mPEG2k-PCL2k + 20 Abi-SI-C18.",
  "DTX+Abi mixed: 30 mPEG2k-PCL2k + 10 DTX-SI-C18 + 10 Abi-SI-C18.",
  "",
  "1. Rg, COM-distance and H-bond panels are curated using explicit file dictionaries to avoid unreadable all-file bar charts.",
  "2. Full all-file summaries are still exported as CSV tables, and all individual XVG plots remain available in figures/all_individual_xvg.",
  "3. Segment-resolved distance and H-bond heatmaps are generated from explicit mixed-system file mappings rather than inferred labels.",
  "4. The integrated SI overview uses key mechanistic readouts for concise interpretation.",
  "",
  "Recommended supplementary MD narrative:",
  "Use 15_integrated_SI_MD_support_overview as the main supplementary MD support figure. Use 16_integrated_MD_diagnostics_time_series for trajectory/QC transparency if space permits.",
  "",
  "Suggested legend:",
  "Supplementary MD analysis. Molecular dynamics support for DTX-SI-C18-assisted incorporation of Abi-SI-C18 in mixed prodrug-polymer co-assembly. The simulations compared single-prodrug systems containing 30 mPEG2k-PCL2k and 20 DTX-SI-C18 or Abi-SI-C18 molecules with a mixed system containing 30 mPEG2k-PCL2k, 10 DTX-SI-C18 and 10 Abi-SI-C18 molecules. XVG outputs from GROMACS were converted and analysed to quantify assembly compactness, prodrug-carrier distances, hydrogen bonding, contact numbers, RDF profiles and SASA-derived interfacial contact areas. Final-window values were calculated from the last 50 ns when available. Curated summary panels focus on key prodrug-polymer, PCL, PEG and DTX-Abi interactions, while full data tables and all individual XVG plots are provided for transparency."
)
writeLines(figure_log, file.path(summary_dir, "MD_support_figure_guide.txt"), useBytes = TRUE)

run_log <- tibble::tibble(
  item = c("script_dir", "root_dir", "fig_dir", "n_xvg_files_indexed", "n_parsed_rows", "make_individual_xvg_plots"),
  value = c(script_dir, root_dir, fig_dir, nrow(xvg_index), nrow(all_xvg), as.character(MAKE_INDIVIDUAL_XVG_PLOTS))
)
readr::write_csv(run_log, file.path(summary_dir, "MD_support_run_log.csv"))

message("Done.")
message("Figures saved to: ", fig_dir)
message("Converted CSV files saved to: ", converted_dir)
message("Summary tables saved to: ", summary_dir)
message("Figure guide saved to: ", file.path(summary_dir, "MD_support_figure_guide.txt"))

# ============================================================
# ============================================================
