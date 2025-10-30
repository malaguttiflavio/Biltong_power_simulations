# ---------------------------------------------------------------------
# BILTONG POWER SIMULATIONS - COMMON UTILITY FUNCTIONS
# ---------------------------------------------------------------------

# Robust SE coeftest wrapper for lm and ivreg
# - If `cluster` is provided, uses sandwich::vcovCL(fit, cluster=...)
# - Otherwise, uses sandwich::vcovHC(fit, type = hc_type)
robust_test <- function(fit, cluster = NULL, hc_type = "HC1") {
  vc <- tryCatch({
    if (is.null(cluster)) {
      sandwich::vcovHC(fit, type = hc_type)
    } else {
      sandwich::vcovCL(fit, cluster = cluster)
    }
  }, error = function(e) {
    # Fallback defensively to HC1
    sandwich::vcovHC(fit, type = hc_type)
  })
  ct <- lmtest::coeftest(fit, vcov. = vc)
  tibble::tibble(term = rownames(ct),
                 estimate = ct[,1], std_error = ct[,2],
                 statistic = ct[,3], p_value = ct[,4])
}

# Balanced assignment within strata by treat_share (binary)
stratified_assign <- function(df, treat_share) {
  DT <- data.table::as.data.table(df)
  DT[, Z := 0L]
  DT[, Z := {
    n <- .N
    if (n <= 0) return(integer(0))
    n_t <- round(treat_share * n)
    n_t <- pmin(n_t, n)
    if (n_t > 0) {
      idx <- sample.int(n, n_t)
      z <- integer(n); z[idx] <- 1L; z
    } else {
      integer(n)
    }
  }, by = .(year_in_program, ngo_id, tribe_id)]
  as.data.frame(DT)
}

# Generate AR(1) errors for a unit across provided times
ar1_series <- function(times_vec, rho, sigma = 1) {
  Tn <- length(times_vec)
  if (Tn == 0) return(numeric(0))
  sigma <- pmax(sigma, 1e-8)
  if (!is.finite(sigma)) sigma <- 1
  rho <- pmax(pmin(rho, 0.999), -0.999)
  e <- numeric(Tn)
  e[1] <- stats::rnorm(1, 0, sigma / sqrt(1 - rho^2 + 1e-9))
  if (Tn >= 2) {
    for (t in 2:Tn) e[t] <- rho * e[t-1] + stats::rnorm(1, 0, sigma)
  }
  e
}

# Draw counts from Poisson or NegBin given mean mu
r_count <- function(n, mu, family = c("poisson","negbin"), theta = 1) {
  # 'none' handled upstream; if passed here just return mu
  if (length(family) == 1 && tolower(family) == "none") return(mu)
  family <- match.arg(family)
  if (length(n) == 1 && length(mu) > 1) n <- length(mu)
  if (length(mu) == 1 && n > 1) mu <- rep(mu, n)
  if (family == "poisson") return(stats::rpois(n, lambda = pmax(mu, 1e-8)))
  prob <- theta / (theta + pmax(mu, 1e-8))
  stats::rnbinom(n, size = theta, prob = prob)
}
