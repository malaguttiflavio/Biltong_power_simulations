# Example script using the updated power simulation function
# This exports all results to Google Drive

# Load the main simulation function
# Note: Make sure Biltong_power_simulations.r is in your working directory or adjust path
# source('Biltong_power_simulations.r')

# For now, let's use a simplified version with Google Drive export
suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

# Simple power simulation example with Google Drive export
run_power_analysis_example <- function() {
  
  # Output directory (Google Drive)
  output_dir <- "/Users/flavioamalagutti/Library/CloudStorage/GoogleDrive-fmalagutti@ucsb.edu/Shared drives/emlab/projects/current-projects/arnhold-biltong/biltong_rct_design/Power_simulations"
  
  # Create directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Example power simulation data (replace with actual simulation results)
  effect_sizes <- seq(1.0, 1.5, by = 0.1)
  
  # Simulate power analysis results
  results <- tibble(
    sweep_param = "effect_meeting",
    sweep_value = effect_sizes,
    experiment_type = "community_level",
    participation_dist = "negbin",
    power_ITT_participation = pmin(0.05 + 0.7 * (effect_sizes - 1), 0.95),
    power_TOT_participation = pmin(0.03 + 0.65 * (effect_sizes - 1), 0.90),
    power_ITT_VCU = pmin(0.08 + 0.75 * (effect_sizes - 1), 0.98),
    power_TOT_VCU = pmin(0.06 + 0.70 * (effect_sizes - 1), 0.92)
  )
  
  # Create figure
  plt <- results |>
    pivot_longer(cols = starts_with("power_"), names_to = "metric", values_to = "power") |>
    mutate(series = case_when(
      metric == "power_ITT_participation" ~ "ITT - Participation",
      metric == "power_TOT_participation" ~ "TOT - Participation",
      metric == "power_ITT_VCU" ~ "ITT - VCU", 
      metric == "power_TOT_VCU" ~ "TOT - VCU",
      TRUE ~ metric
    )) |>
    ggplot(aes(x = sweep_value, y = power, color = series)) +
    geom_line(linewidth = 1) + 
    geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = 0.9, linetype = "dotted", alpha = 0.7) +
    labs(
      x = "Effect Size (Meeting Rate Ratio)",
      y = "Statistical Power",
      title = "Biltong RCT Power Analysis: Effect Size vs Statistical Power",
      subtitle = "Community-level randomization | Negative binomial participation",
      color = "Estimand",
      caption = "Dashed line: 80% power threshold; Dotted line: 90% power threshold"
    ) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  # Save files to Google Drive
  outfile_stem <- "biltong_power_analysis"
  sweep_param <- "effect_meeting"
  
  csv_path <- file.path(output_dir, glue("{outfile_stem}_{sweep_param}.csv"))
  png_path <- file.path(output_dir, glue("{outfile_stem}_{sweep_param}.png"))
  pdf_path <- file.path(output_dir, glue("{outfile_stem}_{sweep_param}.pdf"))
  
  write_csv(results, csv_path)
  ggsave(png_path, plt, width = 12, height = 8, dpi = 300)
  ggsave(pdf_path, plt, width = 12, height = 8)
  
  cat("Power analysis completed!\n")
  cat("Files saved to Google Drive:\n")
  cat("CSV:", csv_path, "\n")
  cat("PNG:", png_path, "\n") 
  cat("PDF:", pdf_path, "\n")
  
  return(list(results = results, csv = csv_path, png = png_path, pdf = pdf_path))
}

# Run the example
cat("Running Biltong RCT Power Analysis...\n")
output <- run_power_analysis_example()
cat("Analysis complete!\n")
print(output$results)