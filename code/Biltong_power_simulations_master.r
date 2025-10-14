# ------------------------------ GRASS RCT Power Simulation ------------------------------ #

# Header
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
n_communities        <- 130
avg_indiv_per_comm   <- 100
sd_indiv_per_comm    <- 10
baseline_outcome_pct <- 0.35
sims                 <- 300
alpha                <- 0.05
hc_type              <- "HC1"  # Type of heteroskedasticity-consistent SE for inference # nolint 

config <- list(

  experiment_type    = "community_level",            # Unit of randomization: "community_level" or "individual_within_community"
  n_communities      = n_communities,                # Number of clusters (communities or associations)
  avg_indiv_per_comm = avg_indiv_per_comm,           # Average number of individuals per cluster
  sd_indiv_per_comm  = sd_indiv_per_comm,            # Std dev for individuals per cluster # nolint
  
  arms               = c("control","T1"),
  alloc_ratios       = c(control = 0.50, T1 = 0.50),
  take_up            = c(control = 0.00, T1 = 1.00),
  ate_pct            = c(T1 = 1.05),
  VCU_delta_mean     = c(T1 = 15),

  # arms               = c("control","T1","T2"),
  # alloc_ratios       = c(control = 0.34, T1 = 0.33, T2 = 0.33),
  # take_up            = c(control = 0.05, T1 = 0.50, T2 = 0.85),
  # ate_pct            = c(T1 = 1.20, T2 = 1.45),
  # VCU_delta_mean     = c(T1 = 15,  T2 = 10),        # Per-arm additive effects on VCU (omit 'control')
  
  cluster_se         = TRUE, # whether to use cluster-robust SEs at the community level, otherwise: heteroskedasticity-robust 
  use_cluster_fe     = FALSE, #if TRUE, the model is Y ~ arm + factor(time) + factor(community_id); otherwise Y ~ arm + factor(time).
                              #In a pure community_level RCT where treatment assignment is constant within each community and does not change over time, including community FE in the ITT model makes the arm indicator collinear with the fixed effects

  year_in_program    = sample(1:3, n_communities, TRUE), # Stratifier: years since joining program
  ngo_id             = sample(1:7, n_communities, TRUE), # Stratifier: NGO identifier
  tribe_id           = sample(1:5, n_communities, TRUE), # Stratifier: tribal/community group ID

  T_outcome          = 3,                            # Number of time points for (formerly meeting) outcome
  months_outcome     = c(1,6,12),                    # Timeline (in months) for outcome observations
  outcome_dist       = "negbin",                     # Distribution for outcome ("poisson" or "negbin")
  outcome_base_mean  = avg_indiv_per_comm * baseline_outcome_pct, # Base expected outcome mean per cluster
  theta_outcome      = 6,                            # Dispersion parameter for negative binomial distribution
  ICC_outcome        = 0.05,                         # Intra-cluster correlation for outcome
  AR1_outcome        = 0.3,                          # Autocorrelation across time for outcome

  T_VCU              = 3,                            # Number of time points for VCU outcome
  months_VCU         = c(1,6,12),                    # Timeline (in months) for VCU observations
  VCU_base_mean      = 220,                          # Base mean VCUs per cluster
  VCU_base_sd        = 70,                           # Standard deviation of VCU outcome
  ICC_VCU            = 0.1,                          # Intra-cluster correlation for VCU outcome
  AR1_VCU            = 0.4,                          # Autocorrelation across time for VCU outcome

  sims               = sims,                         # Number of Monte Carlo simulation repetitions
  alpha              = alpha,                        # Significance level for power calculations
  hc_type            = hc_type                       # Type of heteroskedasticity-consistent SE for inference # nolint

)

# ------------------------------ Run Simulation ------------------------------ #
out <- simulate_power(
  config,
  sweep_param    = "n_communities",
  sweep_arm      = NULL,    # If NULL, sweeps all arms
  sweep_each_arm = TRUE,
  sweep_values   = seq(100, 300, by = 20),
  # Parallelization layer options: "inner" (default), "outer", "both", or "none".
  # "inner" is generally fastest and avoids oversubscription.
  parallel_layer = "inner",
  seed = 123
  )

print(out$results)
print(out$csv)
print(out$long_csv)
print(out$participation_png)
print(out$participation_pdf)
print(out$VCU_png)
print(out$VCU_pdf)

# ------------------------------ Cleanup ------------------------------ #
plan(sequential)