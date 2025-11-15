# ------------------------------ GRASS RCT Power Simulation ------------------------------ #
# Master script to run power simulations using sourced engine and utilities

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(AER)
  library(glue)
  library(parallel)
  library(furrr)
  library(progressr)
})
handlers(global = TRUE)
plan(multisession, workers = 14)
# plan(sequential) for no parallelization

# Source utility and engine files (robust to being sourced from any working directory)
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
source(file.path(script_dir, "..", "common", "biltong_power_utils.r"))
source(file.path(script_dir, "biltong_static_power_sims_engine.r"))

# Resolve project root relative to this script and set output directories (static artifacts)
project_root <- normalizePath(file.path(script_dir, "..", ".."))
tabs_dir <- file.path(project_root, "output", "static", "tabs")
figs_dir <- file.path(project_root, "output", "static", "figs")
if (!dir.exists(tabs_dir)) dir.create(tabs_dir, recursive = TRUE)
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# ------------------------------ Specify Scenario ------------------------------ #
n_communities            <- 160
avg_ind_obs_per_comm     <- 100
sd_indiv_per_comm        <- 10
soc_outcome_baseline_pct <- 0.35 # x% of the individual in the community are counted 
assoc_area               <- 2000
base_fire_avg_pct        <- 0.2 
base_fire_sd_pct         <- 0.1
arm_mode                 <- "Mult"  # change to "Mult" to use the multi-arm parameters below

config <- list(

  # General parameters
  sims                 = 300,                         # Number of Monte Carlo simulation repetitions
  alpha                = 0.05,                        # Significance level for power calculations (maybe should be 0.025 for two-sided test)
  hc_type              = "HC1",                       # Type of heteroskedasticity-consistent SE for inference # nolint 

  experiment_type      = "community_level",            # Unit of randomization: "community_level" or "individual_within_community"
  n_communities        = n_communities,                # Number of clusters (communities or associations)
  avg_ind_obs_per_comm = avg_ind_obs_per_comm,         # Average number of individuals per cluster
  sd_indiv_per_comm    = sd_indiv_per_comm,            # Std dev for individuals per cluster # nolint
  
  # Stratification variables and treatment-effect heterogeneity (additive ATE by stratum)
  year_in_program    = sample(1:3, n_communities, TRUE), # Stratifier: years since joining program
  ngo_id             = sample(1:6, n_communities, TRUE), # Stratifier: NGO identifier
  tribe_id           = sample(1:5, n_communities, TRUE), # Stratifier: tribal/community group ID

  # If *_ate is NULL, the engine draws per-level effects from N(0, *_ate_var).
  year_in_program_ate     = NULL,  # optionally c(`1`=0.0,`2`=0.02,`3`=-0.01)
  ngo_id_ate              = NULL,  # optionally a named vector per NGO id
  tribe_id_ate            = NULL,  # optionally a named vector per tribe id
  year_in_program_ate_var = 0.50,  # variance for draws on log-scale (soc) and level-scale (env) for N(0, v)
  ngo_id_ate_var          = 0.50,
  tribe_id_ate_var        = 0.50,

  # Econometric model parameters
  cluster_se         = TRUE,  # whether to use cluster-robust SEs at the community level, otherwise: heteroskedasticity-robust 
  cluster_fe_yn      = FALSE, # if TRUE, the model is Y ~ arm + factor(time) + factor(community_id); otherwise Y ~ arm + factor(time).
                              # In a pure community_level RCT where treatment assignment is constant within each community and does not change over time, including community FE in the ITT model makes the arm indicator collinear with the fixed effects

  # Single-arm parameters
  arms_single                 = c("control","T1"),
  alloc_ratios_single         = c(control = 0.50, T1 = 0.50),
  take_up_single              = c(control = 0.00, T1 = 1.00),
  soc_outcome_ate_pct_single  = c(T1 = 1.05),
  env_outcome_ate_single      = c(T1 = 10),

  # Multi-arm parameters
  arms_multi                  = c("control","T1","T2"),
  alloc_ratios_multi          = c(control = 0.34, T1 = 0.33, T2 = 0.33),
  take_up_multi               = c(control = 0.00, T1 = 1.00, T2 = 1.00),
  soc_outcome_ate_pct_multi   = c(T1 = 1.05, T2 = 1.45),
  env_outcome_ate_multi       = c(T1 = 10,  T2 = 10),        # Per-arm additive effects on environmental outcome (omit 'control')
  
  # Socioeconomic outcome parameters: log and counts
  soc_outcome_T         = 3,                            # Number of time points for socio-economic outcome
  soc_outcome_T_months  = c(1,6,12),                    # Timeline (in months) for socio-economic observations
  soc_outcome_dist      = "negbin",                     # Distribution for socio-economic ("poisson","negbin","none" for deterministic)
  soc_outcome_base_mean = avg_ind_obs_per_comm * soc_outcome_baseline_pct, # Base expected outcome mean per cluster
  soc_outcome_theta     = 6,                            # Dispersion parameter for negative binomial distribution, doesnt matter for others
  soc_outcome_ICC       = 0.05,                         # Intra-cluster correlation for socio-economic outcome
  soc_outcome_AR1_rho   = 0.20,                         # Autocorrelation across time for socio-economic outcome: e_t = \rho e_{t-1} + \varepsilon_t; rho = 1 means random walk, rho = 0 means no correlation only variance\white noise
  soc_outcome_AR1_var   = 0.5,                          # NEW: Variance of AR(1) innovation for socio-economic outcome (on log scale). If omitted, defaults to 1 in engine.

  # Environmental outcome parameters: level changes
  env_outcome_T         = 6,                            # Number of time points for environmental outcome
  env_outcome_T_months  = c(1, 2, 3, 4, 5, 6),          # Timeline (in months) for environmental observations
  env_outcome_base_mean = assoc_area*base_fire_avg_pct, # Base mean environmental outcome per cluster
  env_outcome_base_sd   = assoc_area*base_fire_sd_pct,  # Standard deviation of environmental outcome
  env_outcome_ICC       = 0.03,                         # Intra-cluster correlation for environmental outcome
  env_outcome_AR1_rho   = 0.50,                         # Autocorrelation across time for environmental outcome: e_t = \rho e_{t-1} + \varepsilon_t; rho = 1 means random walk, rho = 0 means no correlation only variance\white noise. For observational data, it is likely high
  env_outcome_AR1_var   = 1.00                          # NEW: Variance of AR(1) innovation for environmental outcome. If omitted, engine defaults to env_outcome_base_sd^2.

  )

# Map the chosen arm set into the standard fields expected by the engine
if (toupper(arm_mode) %in% c("SINGLE","S")) {
  config$arms           <- config$arms_single
  config$alloc_ratios   <- config$alloc_ratios_single
  config$take_up        <- config$take_up_single
  # Engine now reads soc_outcome_ate_pct_* directly; no legacy mapping needed
  # Engine now reads env_outcome_ate_* directly; no legacy mapping needed
} else if (toupper(arm_mode) %in% c("MULT","MULTI","M")) {
  config$arms           <- config$arms_multi
  config$alloc_ratios   <- config$alloc_ratios_multi
  config$take_up        <- config$take_up_multi
  # Engine now reads soc_outcome_ate_pct_* directly; no legacy mapping needed
  # Engine now reads env_outcome_ate_* directly; no legacy mapping needed
} else {
  stop("arm_mode must be 'Single' or 'Mult'")
}

# Enable stratified randomization by default across the three community-level characteristics
config$stratify_by    <- c("year_in_program","ngo_id","tribe_id")
config$stratify_exact <- TRUE

  
# ------------------------------ Run Simulations (separate sweeps) ------------------------------ #
# Toggle to also run the original single-scenario sweeps (kept for reference)
run_baseline_sweeps <- FALSE
# Note on sweep_param keys (single vs multi arm):
# - We consistently use the multi-arm keys ("soc_outcome_ate_pct_multi", "env_outcome_ate_multi").
# - The engine auto-resolves to the appropriate container based on your config: if only a single
#   treatment arm is present, it maps the sweep to the single-arm container; with multiple arms,
#   it uses the multi-arm container.
# - This lets the same master script work unchanged for both single- and multi-arm scenarios.
# - With sweep_each_arm = TRUE, the engine runs independent sweeps per treatment arm.

if (isTRUE(run_baseline_sweeps)) {
  # 1) Socio-economic: sweep over soc_outcome_ate_pct_multi
  out_soc <- simulate_power(
    config,
    sweep_param    = "soc_outcome_ate_pct_multi",
    sweep_arm      = NULL,    # If NULL, sweeps all arms
    sweep_each_arm = TRUE,
    sweep_values   = seq(1.05, 1.4, by = 0.05),
    parallel_layer = "inner",
    seed = 123
  )

  cat("Socio-economic sweep artifacts:\n",
      out_soc$csv, "\n",
      out_soc$long_csv, "\n",
      out_soc$soc_outcome_png, "\n",
      out_soc$soc_outcome_pdf, "\n")

  # 2) Environmental: sweep over env_outcome_ate_multi
  out_env <- simulate_power(
    config,
    # Use multi-arm key universally; engine auto-resolves to single- or multi- container as needed
    sweep_param    = "env_outcome_ate_multi",
    sweep_arm      = NULL,
    sweep_each_arm = TRUE,
    sweep_values   = seq(5, 30, by = 5),
    parallel_layer = "inner",
    seed = 456
  )

  cat("Environmental sweep artifacts:\n",
      out_env$csv, "\n",
      out_env$long_csv, "\n",
      out_env$env_outcome_png, "\n",
      out_env$env_outcome_pdf, "\n")
}

# ------------------------------ Sensitivity Grid (ICC × AR1 × theta) ------------------------------ #
# We vary: ICC (applied to both outcomes), AR1 rho (applied to both outcomes), and theta
# (dispersion for Negative Binomial) which applies only to socio-economic outcome.
# Assumptions:
# - We set soc_outcome_dist = "negbin" so theta is relevant; adjust values below as needed.
# - The same ICC and AR1 values are used for both outcomes for each scenario.

power_target <- 0.80   # Target power for MDE computation

# Define small grids (edit as needed)
icc_vals   <- c(0.02, 0.05, 0.10)
ar1_vals   <- c(0.20, 0.50, 0.80)
theta_vals <- c(4, 6, 10)

# Ensure NEGBIN is used for socio-economic outcome when running the grid
config$soc_outcome_dist <- "negbin"

grid <- expand.grid(ICC = icc_vals, AR1 = ar1_vals, theta = theta_vals, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

compute_mde <- function(long_df, target_power = 0.80) {
  # long_df expected columns: outcome, estimator, arm, sweep_value, power
  long_df %>%
    dplyr::filter(!is.na(arm)) %>%
    dplyr::group_by(outcome, estimator, arm) %>%
    dplyr::summarise(
      MDE = {
        v <- sweep_value[power >= target_power]
        if (length(v) == 0) NA_real_ else min(v, na.rm = TRUE)
      },
      .groups = "drop"
    )
}

mde_accum <- list()

for (i in seq_len(nrow(grid))) {
  scn <- grid[i, ]

  # Scenario-tag for filenames (safe for paths): replace '.' with 'p'
  fmt_num <- function(x) gsub("\\.", "p", sprintf("%.3f", x))
  tag <- glue::glue("ICC{fmt_num(scn$ICC)}_AR1{fmt_num(scn$AR1)}_theta{fmt_num(scn$theta)}")

  config_scn <- config
  # Apply scenario parameters
  config_scn$soc_outcome_ICC     <- scn$ICC
  config_scn$env_outcome_ICC     <- scn$ICC
  config_scn$soc_outcome_AR1_rho <- scn$AR1
  config_scn$env_outcome_AR1_rho <- scn$AR1
  config_scn$soc_outcome_theta   <- scn$theta
  config_scn$soc_outcome_dist    <- "negbin"

  # Socio-economic sweep (multiplicative ATEs)
  out_soc <- simulate_power(
    config_scn,
    sweep_param    = "soc_outcome_ate_pct_multi",
    sweep_arm      = NULL,
    sweep_each_arm = TRUE,
    sweep_values   = seq(1.05, 1.40, by = 0.05),
    parallel_layer = "inner",
    seed = 123,
    outfile_stem   = glue::glue("biltong_power_{tag}")
  )

  # Environmental sweep (additive ATEs)
  out_env <- simulate_power(
    config_scn,
    sweep_param    = "env_outcome_ate_multi",
    sweep_arm      = NULL,
    sweep_each_arm = TRUE,
    sweep_values   = seq(5, 30, by = 5),
    parallel_layer = "inner",
    seed = 456,
    outfile_stem   = glue::glue("biltong_power_{tag}")
  )

  # Read long tables for MDE computation (structure is stable and convenient)
  soc_long <- tryCatch(readr::read_csv(out_soc$long_csv, show_col_types = FALSE), error = function(e) NULL)
  env_long <- tryCatch(readr::read_csv(out_env$long_csv, show_col_types = FALSE), error = function(e) NULL)

  mde_soc <- if (!is.null(soc_long)) compute_mde(soc_long, target_power = power_target) else dplyr::tibble()
  mde_env <- if (!is.null(env_long)) compute_mde(env_long, target_power = power_target) else dplyr::tibble()

  mde_scn <- dplyr::bind_rows(mde_soc, mde_env) %>%
    dplyr::mutate(
      scenario_id = tag,
      ICC = scn$ICC,
      AR1 = scn$AR1,
      theta = scn$theta
    ) %>%
    dplyr::relocate(scenario_id, ICC, AR1, theta)

  mde_accum[[length(mde_accum) + 1]] <- mde_scn
}

if (length(mde_accum) > 0) {
  mde_all <- dplyr::bind_rows(mde_accum)
  mde_outfile <- file.path(tabs_dir, "biltong_power_sensitivity_mde.csv")
  readr::write_csv(mde_all, mde_outfile)
  cat("Saved aggregated MDE table to:\n", mde_outfile, "\n")
}

# ------------------------------ Cleanup ------------------------------ #
plan(sequential)