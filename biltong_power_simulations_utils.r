# ---------------------------------------------------------------------
# BILTONG POWER SIMULATIONS - UTILITY FUNCTIONS
#
# This file defines helper functions used throughout the power 
# simulation engine. These functions are modular and reusable, 
# providing functionality for:
#  - Generating AR(1) time series errors
#  - Simulating count data (Poisson or Negative Binomial)
#  - Stratified treatment assignment across clusters
#  - Cluster-robust statistical testing (HC1)
#
# The term "utilities" reflects their general-purpose nature—they
# support the simulation logic but are not tied to any specific 
# RCT design or outcome structure. 
#
# ---------------------------------------------------------------------
# ------------------ Utility Functions ------------------ #

library(data.table)

# Robust SE coeftest wrapper for lm and ivreg
robust_test <- function(fit, hc_type = "HC1") {
  ct <- coeftest(fit, vcov = sandwich::vcovHC(fit, type = hc_type))
  tibble(term = rownames(ct), estimate = ct[,1], std_error = ct[,2],
         statistic = ct[,3], p_value = ct[,4])
}

# Create balanced assignment within strata according to allocation ratio
stratified_assign <- function(df, treat_share) {
  setDT(df)
  df[, Z := 0L]
  df[, Z := {
    n <- .N
    if (n == 0) return(integer(0))
    n_t <- round(treat_share * n)
    n_t <- pmin(n_t, n)
    if (n_t > 0) {
      idx <- sample.int(n, n_t)
      z <- integer(n); z[idx] <- 1L; z
    } else {
      integer(n)
    }
  }, by = .(year_in_program, ngo_id, tribe_id)]
  setDF(df)
  df
}

# Generate AR(1) errors for a unit across provided times
ar1_series <- function(times_vec, rho, sigma = 1) {
  Tn <- length(times_vec)
  if (Tn == 0) return(numeric(0))
  sigma <- pmax(sigma, 1e-8)
  if (!is.finite(sigma)) sigma <- 1
  rho <- pmax(pmin(rho, 0.999), -0.999)
  e <- numeric(Tn)
  e[1] <- rnorm(1, 0, sigma / sqrt(1 - rho^2 + 1e-9))
  if (Tn >= 2) {
    for (t in 2:Tn) e[t] <- rho * e[t-1] + rnorm(1, 0, sigma)
  }
  e
}

# Draw counts from Poisson or NegBin given mean mu
r_count <- function(n, mu, family = c("poisson","negbin"), theta = 1) {
  family <- match.arg(family)
  if(length(n) == 1 && length(mu) > 1) n <- length(mu)
  if(length(mu) == 1 && n > 1) mu <- rep(mu, n)
  if (family == "poisson") return(rpois(n, lambda = pmax(mu, 1e-8)))
  prob <- theta / (theta + pmax(mu, 1e-8))
  rnbinom(n, size = theta, prob = prob)
}