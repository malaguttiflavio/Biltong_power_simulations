#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(RColorBrewer)
})

# Regenerate PNG/PDF figures from *_long.csv files in output/cattle/tabs
# This script reads each long CSV, filters ITT series, maps series names to
# friendly labels, and writes PNG and PDF files into output/cattle/figs.

project_root <- normalizePath(".")
 tabs_dir <- file.path(project_root, "output", "cattle", "tabs")
 figs_dir <- file.path(project_root, "output", "cattle", "figs")
 if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

 long_files <- list.files(tabs_dir, pattern = "_long\\.csv$", full.names = TRUE)
 if (length(long_files) == 0) {
   message("No *_long.csv files found in: ", tabs_dir)
   quit(status = 0)
 }

 get_palette <- function(n, pal_name = "Dark2") {
   fallback_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
   if (!requireNamespace("RColorBrewer", quietly = TRUE)) return(fallback_colors[seq_len(min(n, length(fallback_colors)))])
   max_in_set <- suppressWarnings(RColorBrewer::brewer.pal.info[pal_name, "maxcolors"])
   if (is.na(max_in_set)) return(fallback_colors[seq_len(min(n, length(fallback_colors)))])
   k <- min(max_in_set, max(3, n))
   RColorBrewer::brewer.pal(k, pal_name)[seq_len(n)]
 }

for (f in long_files) {
  message("Processing ", basename(f))
  long_res <- readr::read_csv(f, show_col_types = FALSE)
  if (!"series" %in% names(long_res) || !"power" %in% names(long_res)) {
    message("  skipping (not a long power CSV): ", basename(f)); next
  }

  # Use filename stem to build output names and try to read companion meta CSV
  stem <- sub("_long\\.csv$", "", basename(f))
  # e.g., cattle_power_left_night_left_night_long.csv -> cattle_power_left_night_left_night
  png_path <- file.path(figs_dir, paste0(stem, ".png"))
  pdf_path <- file.path(figs_dir, paste0(stem, ".pdf"))

  # Map series labels to user-friendly names
  plot_res <- long_res %>%
    filter(grepl("^power_itt", series)) %>%
    mutate(series = case_when(
      series == "power_itt_T1_vs_C" ~ "ITT: T1 vs Control",
      series == "power_itt_T2_vs_C" ~ "ITT: T2 vs Control",
      series == "power_itt_T1_vs_T2" ~ "ITT: T1 vs T2",
      TRUE ~ series
    ))

  if (nrow(plot_res) == 0) {
    message("  no ITT series found; skipping plot for ", stem); next
  }

  series_levels <- unique(plot_res$series)
  cols <- setNames(get_palette(length(series_levels), "Dark2"), series_levels)

  # Attempt to read meta CSV for caption; fall back to a short caption
  meta_path <- file.path(tabs_dir, paste0(stem, "_meta.csv"))
  if (file.exists(meta_path)) {
    meta_df <- readr::read_csv(meta_path, show_col_types = FALSE)
    # Build a compact caption from key/value rows (first few rows)
    foot_lines <- meta_df %>%
      mutate(line = glue::glue("{parameter} = {value}")) %>%
      pull(line)
    footnote_text <- paste(foot_lines, collapse = "\n")
  } else {
    footnote_text <- glue::glue("Source: {basename(f)}\nGenerated: {Sys.time()}")
  }

  plt <- ggplot(plot_res, aes(x = sweep_value, y = power, color = series)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = cols) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
    theme_minimal(base_size = 12) +
    labs(
      title = glue::glue("Power vs Sweep: {stem}"),
      y = "Power",
      color = "Legend",
      caption = footnote_text
    ) +
    theme(
      plot.caption = element_text(hjust = 0, size = 8),
      plot.caption.position = "plot"
    )

  # x label: try to infer from column names
  xlab <- if ("sweep_value" %in% names(plot_res)) "Sweep value" else "x"
  plt <- plt + labs(x = xlab)

  tryCatch({
    ggsave(png_path, plt, width = 8, height = 5, dpi = 150)
    ggsave(pdf_path, plt, width = 8, height = 5)
    message("  wrote: ", png_path, ", ", pdf_path)
  }, error = function(e) message("  failed to write plot for ", stem, ": ", e$message))
}

message("Done.")
