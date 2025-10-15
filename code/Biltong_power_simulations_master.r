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
source(file.path(script_dir, "biltong_power_simulations_utils.r"))
source(file.path(script_dir, "biltong_power_simulations_engine.r"))

# ------------------------------ Specify Scenario ------------------------------ #
n_communities            <- 130
avg_ind_obs_per_comm     <- 100
sd_indiv_per_comm        <- 10
soc_outcome_baseline_pct <- 0.35
sims                     <- 300
alpha                    <- 0.05
assoc_area               <- 2000
base_fire_avg_pct        <- 0.2 
base_fire_sd_pct         <- 0.1
arm_mode                 <- "Single"  # change to "Mult" to use the multi-arm parameters below

config <- list(

  # General parameters
  sims                 = sims,                         # Number of Monte Carlo simulation repetitions
  alpha                = alpha,                        # Significance level for power calculations
  hc_type              = "HC1",                      # Type of heteroskedasticity-consistent SE for inference # nolint 

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
  year_in_program_ate_var = 0.00,  # variance for draws on log-scale (soc) and level-scale (env)
  ngo_id_ate_var          = 0.00,
  tribe_id_ate_var        = 0.00

  # Econometric model parameters
  cluster_se         = TRUE, # whether to use cluster-robust SEs at the community level, otherwise: heteroskedasticity-robust 
  cluster_fe_yn      = FALSE, #if TRUE, the model is Y ~ arm + factor(time) + factor(community_id); otherwise Y ~ arm + factor(time).
                              #In a pure community_level RCT where treatment assignment is constant within each community and does not change over time, including community FE in the ITT model makes the arm indicator collinear with the fixed effects

  # Single-arm parameters
  arms_single                   = c("control","T1"),
  alloc_ratios_single           = c(control = 0.50, T1 = 0.50),
  take_up_single                = c(control = 0.00, T1 = 1.00),
  ate_pct_single                = c(T1 = 1.05),
  env_outcome_delta_mean_single = c(T1 = 15),

  # Multi-arm parameters
  arms_multi                   = c("control","T1","T2"),
  alloc_ratios_multi           = c(control = 0.34, T1 = 0.33, T2 = 0.33),
  take_up_multi                = c(control = 0.05, T1 = 0.50, T2 = 0.85),
  ate_pct_multi                = c(T1 = 1.20, T2 = 1.45),
  env_outcome_delta_mean_multi = c(T1 = 15,  T2 = 10),        # Per-arm additive effects on environmental outcome (omit 'control')
  
  # Socioeconomic outcome parameters: log and counts
  soc_outcome_T         = 3,                            # Number of time points for socio-economic outcome
  soc_outcome_T_months  = c(1,6,12),                    # Timeline (in months) for socio-economic observations
  soc_outcome_dist      = "none",                       # Distribution for socio-economic ("poisson","negbin","none" for deterministic)
  soc_outcome_base_mean = avg_ind_obs_per_comm * soc_outcome_baseline_pct, # Base expected outcome mean per cluster
  soc_outcome_theta     = 6,                            # Dispersion parameter for negative binomial distribution, doesnt matter for others
  soc_outcome_ICC       = 0.00,                         # Intra-cluster correlation for socio-economic outcome
  soc_outcome_AR1_rho   = 0.00,                         # Autocorrelation across time for socio-economic outcome: e_t = \rho e_{t-1} + \varepsilon_t; rho = 1 means random walk, rho = 0 means no correlation only variance\white noise
  soc_outcome_AR1_var   = 0.00,                         # NEW: Variance of AR(1) innovation for socio-economic outcome (on log scale). If omitted, defaults to 1 in engine.

  # Environmental outcome parameters: level changes
  env_outcome_T         = 6,                            # Number of time points for environmental outcome
  env_outcome_T_months  = c(1, 2, 3, 4, 5, 6),          # Timeline (in months) for environmental observations
  env_outcome_base_mean = assoc_area*base_fire_avg_pct, # Base mean environmental outcome per cluster
  env_outcome_base_sd   = assoc_area*base_fire_sd_pct,  # Standard deviation of environmental outcome
  env_outcome_ICC       = 0.00,                         # Intra-cluster correlation for environmental outcome
  env_outcome_AR1_rho   = 0.00,                         # Autocorrelation across time for environmental outcome: e_t = \rho e_{t-1} + \varepsilon_t; rho = 1 means random walk, rho = 0 means no correlation only variance\white noise
  env_outcome_AR1_var   = 0.00,                         # NEW: Variance of AR(1) innovation for environmental outcome. If omitted, engine defaults to env_outcome_base_sd^2.

)

# Map the chosen arm set into the standard fields expected by the engine
if (toupper(arm_mode) %in% c("SINGLE","S")) {
  config$arms           <- config$arms_single
  config$alloc_ratios   <- config$alloc_ratios_single
  config$take_up        <- config$take_up_single
  config$ate_pct        <- config$ate_pct_single
  config$env_outcome_delta_mean <- config$env_outcome_delta_mean_single
} else if (toupper(arm_mode) %in% c("MULT","MULTI","M")) {
  config$arms           <- config$arms_multi
  config$alloc_ratios   <- config$alloc_ratios_multi
  config$take_up        <- config$take_up_multi
  config$ate_pct        <- config$ate_pct_multi
  config$env_outcome_delta_mean <- config$env_outcome_delta_mean_multi
} else {
  stop("arm_mode must be 'Single' or 'Mult'")
}

# Enable stratified randomization by default across the three community-level characteristics
config$stratify_by    <- c("year_in_program","ngo_id","tribe_id")
config$stratify_exact <- TRUE

# ------------------------------ Run Simulation ------------------------------ #
out <- simulate_power(
  config,
  sweep_param    = "ate_pct",
  sweep_arm      = NULL,    # If NULL, sweeps all arms
  sweep_each_arm = TRUE,
  sweep_values   = seq(1.05, 1.4, by = 0.05),
  parallel_layer = "inner",
  seed = 123
  )

  # Parallelization layer options: "inner" (default), "outer", "both", or "none".
  # "inner" is generally fastest and avoids oversubscription.

print(out$results)
print(out$csv)
print(out$long_csv)
print(out$soc_outcome_png)
print(out$soc_outcome_pdf)
print(out$env_outcome_png)
print(out$env_outcome_pdf)

# ------------------------------ Cleanup ------------------------------ #
plan(sequential)