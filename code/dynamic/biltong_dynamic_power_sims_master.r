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
"# Source shared utils first, then the engine"
source(file.path(script_dir, "..", "common", "biltong_power_utils.r"))
source(file.path(script_dir, "biltong_dynamic_power_sims_engine.r"))
project_root <- normalizePath(file.path(script_dir, "..", ".."))
output_dir_base <- file.path(project_root, "output", "dynamic")

# ---------- Configure scenario ---------- #
config <- list(
  sims = 300,
  alpha = 0.05,
  n_communities = 120,
  avg_ind_obs_per_comm = 80,
  arms = c("control", "T1", "T2"),

  # Quinn-style: two periods, targeted period-2 assignment to top performers
  T_periods = 2,
  fixed_control = FALSE,           # no anchor by default for Quinn
  # Community-level randomization (Option A)
  randomize_at = "community",     # assign one arm per community each period
  use_eligibility = FALSE,         # eligibility targeting not used in community-level assignment
  adapt_rule = "softmax",         # multi-armed bandit rule for period-2 community choice
  temperature = 1.0,
  epsilon = 0.1,
  min_arm_prob = 0.05,
  min_control_prob = 0.1,
  use_ipw = FALSE,

  # Period 1 allocation
  alloc_ratios_p1 = c(control = 1/3, T1 = 1/3, T2 = 1/3),

  # Period 2 bandit adaptation occurs at community level; eligibility ratios are not used when randomize_at="community"
  target_top_share = 0.2,
  alloc_ratios_eligible = c(control = 0.10, T1 = 0.45, T2 = 0.45),
  alloc_ratios_ineligible = c(control = 1.00, T1 = 0.00, T2 = 0.00),

  # Take-up by period
  take_up_by_period = list(
    t1 = c(control = 0, T1 = 1, T2 = 1),
    t2 = c(control = 0, T1 = 1, T2 = 1)
  ),

  # Outcome DGP
  y_base_mean = 10,
  y_sd = 5,
  ability_sd = 3,
  icc_y = 0.05,
  # Community characteristic + interaction with treatment
  community_char_sd = 1.0,
  y_char_coef = 0.0,                     # baseline effect of the community characteristic on Y
  trt_char_interaction = c(T1 = 0.0, T2 = 0.0),  # interaction strength by arm (set >0 or <0 as needed)

  # Period effects (t1 baseline zero; t2 swept below)
  effects_by_period = list(
    t1 = c(T1 = 0, T2 = 0),
    t2 = c(T1 = 2, T2 = 3)
  )
)

# ---------- Run sweep over last-period effect magnitude ---------- #
out <- simulate_dynamic_power(
  config,
  sweep_param = "effect_last_period",
  sweep_values = seq(0, 6, by = 0.5),
  outfile_stem = "dynamic_power",
  parallel_layer = "inner",
  seed = 42,
  output_dir_base = output_dir_base
)

cat("Dynamic (Quinn two-period) sweep artifacts:\n",
    out$csv, "\n",
    out$png, "\n",
    out$pdf, "\n")

plan(sequential)
