# ---------------- Dynamic Targeting Power - Master ---------------- #

suppressPackageStartupMessages({
  library(tidyverse)
  library(furrr)
  library(progressr)
  library(glue)
})
handlers(global = TRUE)
plan(multisession, workers = 8)

# Source engine relative to this file
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    ca <- commandArgs(trailingOnly = FALSE)
    idx <- grep("^--file=", ca)
    if (length(idx) > 0) {
      return(dirname(normalizePath(sub("^--file=", "", ca[idx][1]))))
    }
    getwd()
  }
)
source(file.path(script_dir, "dynamic_targeting_power_engine.r"))
project_root <- normalizePath(file.path(script_dir, ".."))
output_dir_base <- file.path(project_root, "output", "dynamic")

# ---------- Configure scenario ---------- #
config <- list(
  sims = 300,
  alpha = 0.05,
  n_communities = 120,
  avg_ind_obs_per_comm = 80,
  arms = c("control","T1","T2"),
  # Period 1 (unconstrained)
  alloc_ratios_p1 = c(control = 1/3, T1 = 1/3, T2 = 1/3),
  take_up_p1 = c(control = 0, T1 = 1, T2 = 1),
  # Eligibility targeting (top share by Y1 within community)
  target_top_share = 0.2,
  # Period 2 (targeted)
  alloc_ratios_p2_eligible = c(control = 0.1, T1 = 0.45, T2 = 0.45),
  alloc_ratios_p2_ineligible = c(control = 1.0, T1 = 0.0, T2 = 0.0),
  take_up_p2 = c(control = 0, T1 = 1, T2 = 1),
  # Outcome DGP
  y_base_mean = 10,
  y_sd = 5,
  ability_sd = 3,
  icc_y = 0.05,
  effect_p1 = c(T1 = 0, T2 = 0),
  effect_p2 = c(T1 = 2, T2 = 3)
)

# ---------- Run sweep over period-2 effect magnitude ---------- #
out <- simulate_dynamic_power(
  config,
  sweep_param = "effect_p2",
  sweep_values = seq(0, 6, by = 0.5),
  outfile_stem = "dynamic_power",
  parallel_layer = "inner",
  seed = 42,
  output_dir_base = output_dir_base
)

cat("Dynamic targeting sweep artifacts:\n",
    out$csv, "\n",
    out$png, "\n",
    out$pdf, "\n")

plan(sequential)
