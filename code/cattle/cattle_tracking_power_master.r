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

config <- list(

  # Monte Carlo and inference settings
  sims = 2,            # number of Monte Carlo repetitions (default); reduce (e.g. 6) for quick smoke tests
  alpha = 0.05,        # test size for two-sided tests
  hc_type = "HC1",     # HC type when cluster_se = FALSE
  cluster_se = TRUE,   # TRUE: cluster-robust SEs at association level (recommended; assignment is at association)

  # Design sizes
  n_associations = 160,              # number of associations (clusters); main driver of power
  cows_per_association = 300,        # population size per association (sampling frame)
  cows_tagged_per_association = 100, # number tagged/observed per association (measurement intensity)
  months_T = 12,                     # number of months observed per association/cow
  events_per_month_E = 30,           # number of tracking events tracked per month : event = data point (used when analysis_mode = "event")

  # Randomization
  arms = arms,             # must include "control"; treatment arms like "T1","T2"
  alloc_ratios = alloc,    # assignment probabilities per arm
  take_up = takeup,        # per-arm compliance probabilities; control typically 0

  # Association-level stratifiers
  ngo_id          = sample(1:6, 160, TRUE), # used for stratified randomization and TE heterogeneity (if enabled)
  tribe_id        = sample(1:5, 160, TRUE), # idem
  year_in_program = sample(1:3, 160, TRUE),

  # Variance components
  sigma_assoc = 0.25,  # between-association SD on the link scale (random intercept per association)
  sigma_cow   = 0.25,  # between-cow SD on the link scale (random intercept per cow)
  rho_ar1     = 0.5,   # AR(1) correlation across months within cow/association
  sigma_ar    = 0.5,   # AR(1) innovation SD on the link scale

  # Heterogeneity variances (link scale)
  het_var_ngo   = 0.5,  # variance of treatment-effect heterogeneity by NGO stratum
  het_var_tribe = 0.5,  # variance of TE heterogeneity by tribe stratum
  het_var_year  = 0.5,  # variance of TE heterogeneity by year-in-program

  # Missingness
  missingness_p_cow_month = 0.0,  # probability a cow-month observation is missing (MCAR)

  # Baseline intercepts (link scale)
  mu_baseline = list(
    distance     = log(3.0),         # baseline expected mean (km) at event level (log-link)
    resting      = log(60.0),        # baseline minutes (log-link)
    pasture      = log(120.0),       # baseline minutes (log-link)
    home         = log(180.0),       # baseline minutes (log-link)
    left_morning = qlogis(0.90),     # baseline probability on logit scale
    left_night   = qlogis(0.10)      # baseline probability on logit scale
  ),

  # Treatment effects (multiplicative): per-arm multiplier for continuous; probability-scale multiplier for binary
  effect_mult = c(T1 = 1.05, T2 = 1.20), # e.g., distance: 1.10 means +10% multiplicative effect; binary: percent-change on baseline probability

  # Analysis configuration
  analysis_mode        = "month",    # "month": aggregate events to month-level; "event": analyze each event
  month_aggregate_mode = "sum"       # for continuous outcomes under month mode: aggregate by mean (or sum)
)

# ------------------------------ Run a sweep ------------------------------ #

# Choose the outcome and sweep values depending on type
# Example: distance multiplier sweep (multiplicative effect on mean distance)
# Set the outcome in config (engine reads from config$outcome_selected)
# choose one of: distance, resting, pasture, home (continuous, log-normal); left_morning, left_night (binary, logistic)
outcome_selected = "left_night"  
config$outcome_selected <- outcome_selected
config$print_sim_numbers <- FALSE # set TRUE for literal 1..sims lines (sequential)
config$sim_progress <- TRUE       # show per-simulation progress within each sweep step (parallel-friendly)

out <- simulate_cattle_power(
  config,
  sweep_param = "cows_tagged_per_association",
  sweep_values = c(1, seq(from=5, to=155, by=25)),  # number of cows tagged per association (integers; ≤ cows_per_association)
  sweep_each_arm = FALSE,                   # ignored for design-size sweep; set FALSE for clarity
  parallel_layer = "outer",
  outfile_stem = "cattle_power_left_night",
  output_root = project_root,
  seed = 42
  )

message(glue::glue(
  "Artifacts saved to:\n{out$csv}\n{out$long_csv}\n{out$png}\n{out$pdf}\n"
))

# ------------------------------ Cleanup ------------------------------ #
plan(sequential)
