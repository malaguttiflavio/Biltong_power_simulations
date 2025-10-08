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

# Source utility and engine files
source("Biltong_power_simulations_utils.r")
source("Biltong_power_simulations_engine.r")

# ------------------------------ Specify Scenario ------------------------------ #
n_communities       <- 130
avg_indiv_per_comm  <- 100
pct_participation   <- 0.35

config <- list(
  experiment_type    = "community_level",            # Unit of randomization: "community_level" or "individual_within_community"
  n_communities      = n_communities,                # Number of clusters (communities or associations)
  avg_indiv_per_comm = avg_indiv_per_comm,           # Average number of individuals per cluster
  sd_indiv_per_comm  = 10,                           # Std dev for individuals per cluster # nolint
  alloc_ratio        = 0.5,                          # Treatment allocation ratio

  year_in_program    = sample(1:3, n_communities, TRUE), # Stratifier: years since joining program
  ngo_id             = sample(1:7, n_communities, TRUE), # Stratifier: NGO identifier
  tribe_id           = sample(1:5, n_communities, TRUE), # Stratifier: tribal/community group ID

  T_meeting          = 3,                            # Number of time points for meeting attendance outcome
  times_meeting      = c(1,6,12),                    # Timeline (in months) for meeting attendance observations
  participation_dist = "negbin",                     # Distribution for meeting attendance ("poisson" or "negbin")
  meeting_base_mean  = avg_indiv_per_comm * pct_participation, # Base expected meeting attendance per cluster
  theta_meeting      = 6,                            # Dispersion parameter for negative binomial distribution
  ICC_meeting        = 0.05,                         # Intra-cluster correlation for meeting attendance
  AR1_meeting        = 0.3,                          # Autocorrelation across time for meeting attendance

  T_VCU              = 3,                            # Number of time points for VCU outcome
  times_VCU          = c(1,6,12),                    # Timeline (in months) for VCU observations
  VCU_base_mean      = 220,                          # Base mean VCUs per cluster
  VCU_base_sd        = 70,                           # Standard deviation of VCU outcome
  ICC_VCU            = 0.1,                          # Intra-cluster correlation for VCU outcome
  AR1_VCU            = 0.4,                          # Autocorrelation across time for VCU outcome

  meeting_rate_ratio = 1.15,                         # Treatment effect multiplier on meeting attendance rate
  VCU_delta_mean     = 25,                           # Absolute treatment effect on VCU outcome

  take_up_T          = 0.9,                          # Take-up rate in treatment group
  take_up_C          = 0.05,                         # Take-up rate in control group (noncompliance)

  sims               = 300,                          # Number of Monte Carlo simulation repetitions
  alpha              = 0.05,                         # Significance level for power calculations
  hc_type            = "HC1"                         # Type of heteroskedasticity-consistent SE for inference # nolint
)

# ------------------------------ Run Simulation ------------------------------ #
out <- simulate_power(config,
                      sweep_param = "effect_meeting",
                      sweep_values = seq(1.0, 1.3, by = 0.05),
                      outfile_stem = "grass_power",
                      seed = 123)

print(out$results)
print(out$png)
print(out$pdf)

# ------------------------------ Cleanup ------------------------------ #
plan(sequential)