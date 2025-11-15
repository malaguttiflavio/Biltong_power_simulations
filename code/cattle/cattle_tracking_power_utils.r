# ---------------------------------------------------------------------
# CATTLE TRACKING POWER SIMULATIONS - UTILITY FUNCTIONS
# ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Use common robust test utilities from the static module
# (Assumes this script is sourced from the cattle master/engine which will also source the common utils.)

# Safe AR(1) with arbitrary times (integers) and guardrails
ar1_series_int <- function(Tn, rho, sigma = 1) {
  if (Tn <= 0) return(numeric(0))
  rho <- pmax(pmin(rho, 0.999), -0.999)
  sigma <- pmax(sigma, 1e-8)
  e <- numeric(Tn)
  e[1] <- rnorm(1, 0, sigma / sqrt(1 - rho^2 + 1e-9))
  if (Tn >= 2) for (t in 2:Tn) e[t] <- rho * e[t-1] + rnorm(1, 0, sigma)
  e
}

# Assign arms to associations, optionally stratified by association-level factors
assign_arms_assoc <- function(assoc_df, arms, alloc_ratios, strat_vars = NULL) {
  stopifnot(length(arms) >= 2, arms[1] == "control")
  if (is.null(names(alloc_ratios))) stop("alloc_ratios must be named by arms")
  alloc_ratios <- alloc_ratios[arms]
  alloc_ratios <- alloc_ratios / sum(alloc_ratios)

  DT <- as.data.table(assoc_df)
  if (!is.null(strat_vars) && all(strat_vars %in% names(DT))) {
    byvars <- strat_vars
  } else {
    byvars <- NULL
  }

  if (is.null(byvars)) {
    n <- nrow(DT)
    k <- length(arms)
    # target counts per arm
    target <- round(alloc_ratios * n)
    # adjust to sum to n
    while (sum(target) != n) {
      diff <- n - sum(target)
      j <- sample(seq_len(k), 1)
      target[j] <- target[j] + sign(diff)
    }
    DT[, arm := sample(rep(arms, times = target))]
  } else {
    DT[, arm := NA_character_]
    DT[, arm := {
      n <- .N
      k <- length(arms)
      target <- round(alloc_ratios * n)
      while (sum(target) != n) {
        diff <- n - sum(target)
        j <- sample(seq_len(k), 1)
        target[j] <- target[j] + sign(diff)
      }
      sample(rep(arms, times = target))
    }, by = byvars]
  }
  as_tibble(DT)
}

# Build event panel for tagged cows
build_event_panel <- function(assoc_tbl, cows_tagged_per_assoc, months_T, events_per_month_E) {
  # Expand cows per association
  cow_frame <- assoc_tbl |>
    rowwise() |>
    mutate(cow_id = list(seq_len(cows_tagged_per_assoc))) |>
    unnest(cow_id) |>
    ungroup()

  # Months and events
  panel <- cow_frame |>
    crossing(month_index = seq_len(months_T), event_index = seq_len(events_per_month_E)) |>
    arrange(association_id, cow_id, month_index, event_index) |>
    mutate(event_time_index = (month_index - 1L) * events_per_month_E + event_index)

  panel
}

# Apply cow-month missingness (removes all events for a cow-month)
apply_cow_month_missingness <- function(panel, p_miss) {
  if (is.null(p_miss) || p_miss <= 0) return(panel)
  set.seed(NULL)
  miss_df <- panel |>
    distinct(association_id, cow_id, month_index) |>
    mutate(miss = rbinom(n(), 1L, prob = p_miss) == 1L)
  panel |>
    left_join(miss_df, by = c("association_id","cow_id","month_index")) |>
    filter(!miss) |>
    select(-miss)
}

# Logistic inverse with clamping
inv_logit_clamp <- function(x, eps = 1e-6) {
  p <- 1 / (1 + exp(-x))
  pmin(pmax(p, eps), 1 - eps)
}
