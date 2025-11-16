# ------------------------------ Cattle Tracking Power Simulation - Master ------------------------------ #

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(furrr)
  library(progressr)
})

handlers(global = TRUE)
# Ensure a terminal-friendly progress bar handler is active
handlers(handler_txtprogressbar)
# Parallelism: adjust workers to your machine. The engine uses parallel at the inner loop.
plan(multisession, workers = 12)

# Robust sourcing of common and cattle utilities/engine
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    ca <- commandArgs(trailingOnly = FALSE)
    idx <- grep("^--file=", ca)
    if (length(idx) > 0) return(dirname(normalizePath(sub("^--file=", "", ca[idx][1]))))
    getwd()
  }
)
source(file.path(script_dir, "..", "common", "biltong_power_utils.r"))
source(file.path(script_dir, "cattle_tracking_power_utils.r"))
source(file.path(script_dir, "cattle_tracking_power_engine.r"))

# Resolve project root (engine ensures output dirs under output/cattle)
project_root <- normalizePath(file.path(script_dir, "..", ".."))

# ------------------------------ Configuration ------------------------------ #

# Assignment & compliance
# - arms: include "control" first; names must match throughout.
# - alloc: randomization probabilities across arms; should sum to 1.
# - takeup: probability that an association assigned to arm actually takes up treatment (D=1). With 1.0, ITT ≈ TOT.
arms  <- c("control","T1","T2")
alloc <- c(control = 0.34, T1 = 0.33, T2 = 0.33)
takeup <- c(control = 0.0, T1 = 1.0, T2 = 1.0)
n_associations = 130              # number of associations (clusters); main driver of power

# THE FOLLOWING CONTROL NOISE IN THE DGP
# We feed it sigmas, that is standard deviations (SD).
# in N(0, sigma), 68% of the draws lie within ± sigma from the mean (0). 
sigma_across_assoc     = 1 
sigma_within_assoc     = 0.750
sigma_within_cow       = 0.5

ar1_cow_event_persistence = 0.8
ar1_cow_event_sigma = 0.12

ar1_cow_month_persistence = 0.4
ar1_cow_month_sigma = 0.12

# ------------------------------ Run Simulation ------------------------------ #
config <- list(

  # Monte Carlo and inference settings
  sims = 50,           # number of Monte Carlo repetitions (default); reduce (e.g. 6) for quick smoke tests
  alpha = 0.05,         # test size for two-sided tests
  hc_type = "HC1",      # HC type when cluster_se = FALSE
  cluster_se = TRUE,    # TRUE: cluster-robust SEs at association level (recommended; assignment is at association)
  estimate_tot = FALSE, # TRUE: run TOT (IV) estimations in addition to ITT; FALSE: ITT only (faster)

  # Design sizes
  n_associations = n_associations,   # number of associations (clusters); main driver of power
  cows_per_association = 300,        # population size per association (sampling frame)
  cows_tagged_per_association = 10,  # number tagged/observed per association (measurement intensity)
  months_T = 12,                     # number of months observed per association/cow
  events_per_month_E = 30,           # number of tracking events tracked per month : event = data point (used when analysis_mode = "event")

  # Randomization
  arms = arms,             # must include "control"; treatment arms like "T1","T2"
  alloc_ratios = alloc,    # assignment probabilities per arm
  take_up = takeup,        # per-arm compliance probabilities; control typically 0

  # Association-level stratifiers
  ngo_id          = sample(1:6, n_associations, TRUE), # used for stratified randomization and TE heterogeneity (if enabled); also included as controls in regressions
  tribe_id        = sample(1:5, n_associations, TRUE), # idem
  year_in_program = sample(1:3, n_associations, TRUE), # idem

  # VARIANCE ACROSS ASSOCIATIONS: NGO, TRIBE, YEAR
  # Stratifier baseline variance (link scale) - creates baseline heterogeneity across strata
  strat_sigma_ngo   = sigma_across_assoc,  # variance of baseline shifts by NGO stratum
  strat_sigma_tribe = sigma_across_assoc,  # variance of baseline shifts by tribe stratum
  strat_sigma_yr    = sigma_across_assoc,  # variance of baseline shifts by year-in-program

  # VARIANCE WITHIN ASSOCIATIONS
  # Variance components (ICC components on link scale)
  icc_sigma_assoc = sigma_within_assoc,  # between-association SD (link scale)
  icc_sigma_cow   = sigma_within_cow,    # between-cow SD (link scale)
  
  # AR(1) temporal components
  ar1_cow_event_rho   = ar1_cow_event_persistence,  # event-level AR(1) autocorrelation (events within month)
  ar1_cow_event_sigma = ar1_cow_event_sigma,        # event-level AR innovation (new shock) SD (link scale)
  ar1_cow_month_rho   = ar1_cow_month_persistence,  # month-level AR(1) autocorrelation (across months)
  ar1_cow_month_sigma = ar1_cow_month_sigma,        # month-level AR innovation (new shock) SD (link scale)

  # Missingness
  missingness_p_cow_month = 0.0,  # probability a cow-month observation is missing (MCAR)

  # Baseline intercepts (link scale)
  mu_baseline = list(
    distance     = log(3.0),         # baseline expected mean (km) at event level (log-link)
    resting      = log(60.0),        # baseline minutes (log-link)
    pasture      = log(120.0),       # baseline minutes (log-link)
    home         = log(180.0),       # baseline minutes (log-link)
    left_morning = qlogis(0.95),     # baseline probability on logit scale
    left_night   = qlogis(0.10)      # baseline probability on logit scale
    ),

  # Treatment effects (multiplicative): per-arm multiplier for continuous; probability-scale multiplier for binary
  effect_mult = c(T1 = 1.025, T2 = 1.035), # e.g., distance: 1.10 means +10% multiplicative effect; binary: percent-change on baseline probability

  # Expected direction for signed power calculation: "positive", "negative", or "auto"
  # - "positive": treatment increases outcome (effect_mult > 1.0 expected to give positive coefficients)
  # - "negative": treatment decreases outcome (effect_mult > 1.0 expected to give negative coefficients)
  # - "auto": determine from effect_mult values (> 1.0 = positive, < 1.0 = negative)
  expected_direction = "positive",   # change to "negative" for outcomes where treatment should reduce the outcome

  # Analysis configuration
  analysis_mode        = "month",    # "month": aggregate events to month-level; "event": analyze each event
  month_aggregate_mode = "sum"       # for continuous outcomes under month mode: aggregate by mean (or sum)
  
  )

# ------------------------------ Run a sweep ------------------------------ #

# Choose the outcome and sweep values depending on type
# Example: distance multiplier sweep (multiplicative effect on mean distance)
# Set the outcome in config (engine reads from config$outcome_selected)
# choose one of: distance, resting, pasture, home (continuous, log-normal); left_morning, left_night (binary, logistic)
outcome_selected = "left_morning"  
config$outcome_selected <- outcome_selected
config$print_sim_numbers <- FALSE # set TRUE for literal 1..sims lines (sequential)
config$sim_progress <- TRUE       # show per-simulation progress within each sweep step (parallel-friendly)

out <- simulate_cattle_power(
  config,
  sweep_param = "cows_tagged_per_association",            
  sweep_values = c(1, 10, 25, 50, 100), 
  parallel_layer = "outer",
  outfile_stem = "cattle_power_cows_tagged",
  output_root = project_root,
  seed = 45
  )

message(glue::glue(
  "Artifacts saved to:\n{out$csv}\n{out$long_csv}\n{out$png}\n{out$pdf}\n"
))

# ------------------------------ Cleanup ------------------------------ #
plan(sequential)
