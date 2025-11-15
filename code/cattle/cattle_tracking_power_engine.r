# ------------------------------ Cattle Tracking Power Simulator ----------------------------- #

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(AER)
  library(glue)
  library(furrr)
  library(progressr)
  library(RColorBrewer)
})

# Expect robust_test and helpers from common utils
# Source order: master -> common utils -> cattle utils -> engine

# Helper: robust vcov (cluster or HC)
get_vcov <- function(fit, cluster = NULL, hc_type = "HC1") {
  if (is.null(cluster)) {
    sandwich::vcovHC(fit, type = hc_type)
  } else {
    sandwich::vcovCL(fit, cluster = cluster)
  }
}

# Wald test for linear contrast R beta = r (typically r = 0)
wald_contrast <- function(fit, R, r = 0, cluster = NULL, hc_type = "HC1") {
  V <- tryCatch(
    get_vcov(fit, cluster = cluster, hc_type = hc_type),
    error = function(e) {
      message("Warning: cluster vcov failed (", conditionMessage(e), "); falling back to HC robust vcov.")
      get_vcov(fit, cluster = NULL, hc_type = hc_type)
    }
  )
  b <- coef(fit)
  # Ensure V has row/col names aligned with coefficients to allow safe subsetting
  vn <- names(b)
  if (is.null(rownames(V)) || is.null(colnames(V))) {
    if (!is.null(vn) && length(vn) == nrow(V)) {
      dimnames(V) <- list(vn, vn)
    }
  }
  # Align contrast to vcov/coef names; drop any terms not present
  if (is.null(vn)) {
    # Fallback: dimensions must match; otherwise return NA
    if (length(R) != length(b) || nrow(V) != length(b)) return(list(stat = NA_real_, p = NA_real_))
  } else {
    # Keep only entries present in both R and V/b
    if (!is.null(names(R))) {
      keep <- intersect(vn, names(R))
      if (length(keep) == 0) return(list(stat = NA_real_, p = NA_real_))
      vn <- keep
      b  <- b[keep]
      V  <- V[keep, keep, drop = FALSE]
      R  <- R[keep]
    } else {
      # Unnamed R; require matching length
      if (length(R) != length(b) || nrow(V) != length(b)) return(list(stat = NA_real_, p = NA_real_))
    }
  }
  # Guard: singular or empty
  if (length(R) == 0 || any(!is.finite(R)) || any(!is.finite(b)) || any(!is.finite(V))) {
    return(list(stat = NA_real_, p = NA_real_))
  }
  num <- sum(R * b) - r
  den <- tryCatch(as.numeric(t(matrix(R, ncol = 1)) %*% V %*% matrix(R, ncol = 1)), error = function(e) NA_real_)
  if (!is.finite(den) || den <= 0) return(list(stat = NA_real_, p = NA_real_))
  stat <- (num^2) / den
  p <- 1 - pchisq(stat, df = 1)
  list(stat = stat, p = p)
}

# Palette helper
get_palette <- function(n, pal_name = "Dark2") {
  fallback_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) return(fallback_colors[seq_len(min(n, length(fallback_colors)))])
  max_in_set <- suppressWarnings(RColorBrewer::brewer.pal.info[pal_name, "maxcolors"])
  if (is.na(max_in_set)) return(fallback_colors[seq_len(min(n, length(fallback_colors)))])
  k <- min(max_in_set, max(3, n))
  RColorBrewer::brewer.pal(k, pal_name)[seq_len(n)]
}

# Main function
simulate_cattle_power <- function(
  config,
  sweep_param = c("effect_mult", "cows_tagged_per_association"), # effect_mult: multiplier for continuous; probability-scale multiplier for binary; cows_tagged_per_association: sample size sweep
  sweep_values = NULL,
  sweep_each_arm = TRUE,
  parallel_layer = c("inner","outer","none"),
  outfile_stem = "cattle_power",
  seed = 1,
  output_root = NULL
) {
  set.seed(seed)
  sweep_param <- match.arg(sweep_param)
  parallel_layer <- match.arg(parallel_layer)
  parallel_inner <- parallel_layer == "inner"
  parallel_outer <- parallel_layer == "outer"
  print_sim_numbers <- isTRUE(config$print_sim_numbers %||% FALSE)
  sim_progress <- isTRUE(config$sim_progress %||% FALSE)
  # Printing sim numbers and inner parallelism don't mix well; prefer sequential when printing
  if (print_sim_numbers) parallel_inner <- FALSE

  # If the caller asked for progress, set a terminal-friendly progress handler so
  # `progressr` / `furrr` show a usable progress bar in non-interactive shells.
  if (sim_progress) {
    tryCatch({
      progressr::handlers(progressr::handler_txtprogressbar())
    }, error = function(e) {
      # Fallback quietly if the handler isn't available in this environment
      message("Note: could not set txtprogressbar handler for progressr: ", conditionMessage(e))
    })
  }

  # Unpack config with defaults
  arms <- config$arms; stopifnot(!is.null(arms), length(arms) >= 2, arms[1] == "control")
  alloc_ratios <- config$alloc_ratios; stopifnot(!is.null(alloc_ratios))
  alloc_ratios <- alloc_ratios[arms] / sum(alloc_ratios[arms])
  take_up <- config$take_up; if (is.null(take_up)) take_up <- setNames(c(0, rep(1, length(arms) - 1)), arms)

  n_assoc <- config$n_associations
  cows_per_assoc <- config$cows_per_association
  cows_tagged <- config$cows_tagged_per_association
  months_T <- config$months_T
  events_E <- config$events_per_month_E

  # Stratifiers (association-level)
  ngo <- config$ngo_id; tribe <- config$tribe_id; year <- config$year_in_program
  if (is.null(ngo))  ngo  <- sample(1:6, n_assoc, TRUE)
  if (is.null(tribe)) tribe <- sample(1:5, n_assoc, TRUE)
  if (is.null(year)) year <- sample(1:3, n_assoc, TRUE)

  # Outcome choice and link
  outcome_selected <- config$outcome_selected # e.g., "distance","resting","pasture","home","left_morning","left_night"
  stopifnot(!is.null(outcome_selected))
  is_binary <- outcome_selected %in% c("left_morning","left_night")

  # Effects containers (multiplicative)
  eff_mult <- config$effect_mult # named by treatment arms (exclude control)
  if (is.null(eff_mult)) eff_mult <- setNames(rep(1, length(arms) - 1), arms[-1])

  # Heterogeneity variances (link-scale)
  het_var <- list(
    ngo = config$het_var_ngo %||% 0,
    tribe = config$het_var_tribe %||% 0,
    year = config$het_var_year %||% 0
  )

  # Variance components
  sigma_assoc <- config$sigma_assoc %||% 0.2
  sigma_cow   <- config$sigma_cow   %||% 0.2
  rho_ar1     <- config$rho_ar1     %||% 0.5
  sigma_ar    <- config$sigma_ar    %||% 0.5
  # Optional month-level AR(1) (across months per cow). Defaults to 0 (disabled).
  rho_month   <- config$rho_month   %||% 0
  sigma_month <- config$sigma_month %||% 0

  # Missingness
  p_miss <- config$missingness_p_cow_month %||% 0

  # Analysis settings
  analysis_mode <- match.arg(config$analysis_mode %||% "month", c("month","event"))
  agg_mode_cont <- match.arg(config$month_aggregate_mode %||% "mean", c("mean","sum"))
  alpha <- config$alpha %||% 0.05
  sims  <- config$sims  %||% 300
  hc_type <- config$hc_type %||% "HC1"
  cluster_se <- isTRUE(config$cluster_se %||% TRUE)

  # Output dirs (relative to repo structure using static conventions)
  script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
  project_root <- normalizePath(file.path(script_dir, "..", ".."))
  # Allow caller (master) to override where outputs are written. When running
  # from other contexts R may compute a different script_dir; passing
  # `output_root` forces outputs into the provided project root.
  if (!is.null(output_root)) {
    project_root <- normalizePath(output_root)
  }
  tabs_dir <- file.path(project_root, "output", "cattle", "tabs"); if (!dir.exists(tabs_dir)) dir.create(tabs_dir, recursive = TRUE)
  figs_dir <- file.path(project_root, "output", "cattle", "figs"); if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

  # Association frame and arm assignment
  assoc_df <- tibble(
    association_id = seq_len(n_assoc),
    ngo_id = ngo,
    tribe_id = tribe,
    year_in_program = year
  )
  assoc_df <- assign_arms_assoc(assoc_df, arms = arms, alloc_ratios = alloc_ratios, strat_vars = c("ngo_id","tribe_id","year_in_program"))

  # One Monte Carlo run
  one_run <- function(cfg_local, capture_data = FALSE) {
    # Allow sweeps to override design sizes from cfg_local
  cows_tagged_local <- cfg_local$cows_tagged_per_association %||% cows_tagged
  # Do not exceed the available cow population per association
  cows_tagged_local <- as.integer(pmin(cows_tagged_local, cows_per_assoc))
    months_T_local    <- cfg_local$months_T %||% months_T
    events_E_local    <- cfg_local$events_per_month_E %||% events_E
    # Tag cows per association
  tagged_counts <- rep(cows_tagged_local, n_assoc)
    cow_map <- assoc_df |>
      rowwise() |>
      mutate(cow_id = list(seq_len(tagged_counts[association_id]))) |>
      unnest(cow_id) |>
      ungroup()

    # Build full event panel for tagged cows
    # Build panel inherits columns from assoc_df (including arm and stratifiers) via assoc_tbl
    # Avoid re-joining 'arm' to prevent duplicate columns (arm.x/arm.y). Only join extra stratifiers if needed.
    panel <- build_event_panel(assoc_tbl = assoc_df, cows_tagged_per_assoc = cows_tagged_local, months_T = months_T_local, events_per_month_E = events_E_local) |>
      left_join(assoc_df |> select(association_id, ngo_id, tribe_id, year_in_program), by = "association_id")

    # Cow random intercepts
    cow_re <- cow_map |>
      mutate(u_cow = rnorm(n(), 0, sigma_cow)) |>
      select(association_id, cow_id, u_cow)

    # Association random intercepts
    assoc_re <- assoc_df |>
      transmute(association_id, u_assoc = rnorm(n(), 0, sigma_assoc))

    # Heterogeneity by strata on link scale (per association, then applied to treated D=1)
    tau_strata <- assoc_df |>
      transmute(
        association_id,
        tau_strata = rnorm(n(), 0, sqrt(het_var$ngo)) + rnorm(n(), 0, sqrt(het_var$tribe)) + rnorm(n(), 0, sqrt(het_var$year))
      )

    # Merge random components
    panel <- panel |>
      left_join(cow_re,  by = c("association_id","cow_id")) |>
      left_join(assoc_re, by = "association_id") |>
      left_join(tau_strata, by = "association_id")

    # Take-up at cow-month (default equals 1 for treated arms, 0 for control)
    # Construct D at event level via join from cow-month mapping
    cow_month <- panel |>
      distinct(association_id, cow_id, month_index, arm) |>
      mutate(D = dplyr::case_when(
        arm == "control" ~ 0,
        arm == arms[2]    ~ take_up[[arms[2]]],
        length(arms) >= 3 & arm == arms[3] ~ take_up[[arms[3]]],
        TRUE ~ 0
      ))

    # Apply cow-month missingness
    cow_month <- cow_month |>
      mutate(drop = rbinom(n(), 1L, p_miss) == 1L) |>
      filter(!drop) |>
      select(-drop)

    panel <- panel |>
      inner_join(cow_month, by = c("association_id","cow_id","month_index","arm"))

    # AR(1) per cow across global event time
    panel <- panel |>
      arrange(association_id, cow_id, event_time_index) |>
      group_by(association_id, cow_id) |>
      # Event-level AR(1) across events for this cow (length = months_T * events_E)
      mutate(e_ar = ar1_series_int(n(), rho_ar1, sigma = sigma_ar)) |>
      # Month-level AR(1): generate an AR(1) series across months (length = months_T)
      # and index it by month_index so the same month-level effect applies to all
      # events within that month.
      mutate(month_ar = {
        mseries <- ar1_series_int(months_T, rho_month, sigma = sigma_month)
        # Ensure month_index indexes into mseries (month_index in 1:months_T)
        mseries[month_index]
      }) |>
      ungroup()

    # Baseline mean on link scale
    mu0 <- cfg_local$mu_baseline[[outcome_selected]] %||% 0

    # Link-scale components and linear predictor
    if (is_binary) {
      # Odds ratio effects on logit link
      or_map <- setNames(rep(1, length(arms)), arms)
      for (a in setdiff(arms, "control")) or_map[[a]] <- eff_mult[[a]] %||% 1

      # Build components on link scale (without treatment)
      comp_mu0   <- rep(mu0, nrow(panel))
      comp_assoc <- panel$u_assoc
      comp_cow   <- panel$u_cow
      comp_ar_event <- panel$e_ar
      comp_ar_month <- panel$month_ar
      comp_ar <- comp_ar_event + comp_ar_month

      # Linear predictor without treatment effect
      linpred_no_trt <- comp_mu0 + comp_assoc + comp_cow + comp_ar
      baseline_p <- inv_logit_clamp(linpred_no_trt)

      # Treatment effects are applied multiplicatively on the probability scale.
      # User-specified eff_mult gives the multiplier (e.g., 1.05 = +5%).
      # Heterogeneity (tau_strata) is applied as an additive shift on the probability scale
      # and only affects treated observations (D == 1). We clip final probabilities to [0,1].
      trt_mult_map <- setNames(rep(1, length(arms)), arms)
      for (a in setdiff(arms, "control")) trt_mult_map[[a]] <- eff_mult[[a]] %||% 1

      # Apply multiplier then additive heterogeneity (tau_strata stored per association)
      comp_mult_effect <- trt_mult_map[panel$arm]
      # Ensure comp_mult_effect is numeric vector aligned with panel rows
      comp_mult_effect <- as.numeric(comp_mult_effect)
      # Additive heterogeneity per association (assumed to be on probability scale)
      comp_tau_add <- panel$tau_strata

      # Compute observed probability: baseline for control, baseline * multiplier + additive tau for treated
      comp_prob <- ifelse(panel$D == 1, baseline_p * comp_mult_effect + comp_tau_add, baseline_p)
      # Clip probabilities to [0, 1]
      comp_prob <- pmin(pmax(comp_prob, 0), 1)

      y_event <- rbinom(nrow(panel), 1L, comp_prob)

      df <- panel |>
        mutate(
          y = y_event,
          comp_mu0 = comp_mu0,
          comp_assoc = comp_assoc,
          comp_cow = comp_cow,
          comp_ar = comp_ar,
          comp_tau = comp_tau_add * D, # show tau only when applied
          comp_trt = comp_mult_effect * D,
          comp_linpred = linpred_no_trt,
          comp_prob = comp_prob,
          baseline_prob = baseline_p
        )
    } else {
      # Log-normal multiplicative effects (log-link)
      mult_map <- setNames(rep(1, length(arms)), arms)
      for (a in setdiff(arms, "control")) mult_map[[a]] <- eff_mult[[a]] %||% 1

      comp_mu0   <- rep(mu0, nrow(panel))
      comp_assoc <- panel$u_assoc
      comp_cow   <- panel$u_cow
  comp_ar_event    <- panel$e_ar
  comp_ar_month    <- panel$month_ar
  comp_ar          <- comp_ar_event + comp_ar_month
      comp_tau   <- panel$tau_strata * panel$D
      comp_trt   <- log(mult_map[panel$arm]) * panel$D

      linpred <- comp_mu0 + comp_assoc + comp_cow + comp_ar + comp_tau + comp_trt
      comp_expected_mean <- exp(linpred)
      y_event <- comp_expected_mean # no extra iid noise beyond AR(1); adjust if needed

      df <- panel |>
        mutate(
          y = y_event,
          comp_mu0 = comp_mu0,
          comp_assoc = comp_assoc,
          comp_cow = comp_cow,
          comp_ar = comp_ar,
          comp_tau = comp_tau,
          comp_trt = comp_trt,
          comp_linpred = linpred,
          comp_expected_mean = comp_expected_mean
        )
    }

    # Aggregate to month if needed and prepare analysis frame
    if (analysis_mode == "month") {
      if (is_binary) {
        dfm <- df |>
          group_by(association_id, cow_id, month_index, arm) |>
          summarise(
            y = mean(y),
            D = mean(D),
            comp_mu0 = mean(comp_mu0),
            comp_assoc = mean(comp_assoc),
            comp_cow = mean(comp_cow),
            comp_ar = mean(comp_ar),
            comp_tau = mean(comp_tau),
            comp_trt = mean(comp_trt),
            comp_linpred = mean(comp_linpred),
            comp_prob = mean(comp_prob),
            .groups = "drop"
          )
      } else if (agg_mode_cont == "mean") {
        dfm <- df |>
          group_by(association_id, cow_id, month_index, arm) |>
          summarise(
            y = mean(y),
            D = mean(D),
            comp_mu0 = mean(comp_mu0),
            comp_assoc = mean(comp_assoc),
            comp_cow = mean(comp_cow),
            comp_ar = mean(comp_ar),
            comp_tau = mean(comp_tau),
            comp_trt = mean(comp_trt),
            comp_linpred = mean(comp_linpred),
            comp_expected_mean = mean(comp_expected_mean),
            .groups = "drop"
          )
      } else { # sum
        dfm <- df |>
          group_by(association_id, cow_id, month_index, arm) |>
          summarise(
            y = sum(y),
            D = mean(D),
            comp_mu0 = mean(comp_mu0),
            comp_assoc = mean(comp_assoc),
            comp_cow = mean(comp_cow),
            comp_ar = mean(comp_ar),
            comp_tau = mean(comp_tau),
            comp_trt = mean(comp_trt),
            comp_linpred = mean(comp_linpred),
            comp_expected_mean = mean(comp_expected_mean),
            .groups = "drop"
          )
      }
      df_use <- dfm |>
        mutate(month_f = factor(month_index), assoc_f = factor(association_id))
    } else {
      # Event-level analysis; add factors directly
      df_use <- df |>
        mutate(month_f = factor(month_index), assoc_f = factor(association_id))
    }

    # Arm indicators (control vs T1 vs T2)
    df_use <- df_use |>
      mutate(
        Z_T1 = as.integer(arm == arms[2]),
        Z_T2 = if (length(arms) >= 3) as.integer(arm == arms[3]) else 0L,
        any_treat = as.integer(arm != "control")
      )

    # Construct multi-arm treatment receipt indicators based on actual D; fallback if missing
    if (!"D" %in% names(df_use)) df_use <- df_use |> mutate(D = any_treat)
    df_use <- df_use |>
      mutate(D1 = D * Z_T1, D2 = D * Z_T2, Z1 = Z_T1, Z2 = Z_T2)

    # ITT models
    if (is_binary && analysis_mode == "event") {
      # LPM at event level
      fit_itt <- lm(y ~ Z_T1 + Z_T2 + month_f, data = df_use)
    } else if (is_binary && analysis_mode == "month") {
      fit_itt <- lm(y ~ Z_T1 + Z_T2 + month_f, data = df_use)
    } else {
      # Continuous: log transform for stability
      y_tr <- log(pmax(df_use$y, 1e-8))
      fit_itt <- lm(y_tr ~ Z_T1 + Z_T2 + month_f, data = df_use)
    }

  cluster <- df_use$association_id
  vc <- get_vcov(fit_itt, cluster = if (cluster_se) cluster else NULL, hc_type = hc_type)
    ct <- lmtest::coeftest(fit_itt, vcov. = vc)

    # Extract p-values for T1 vs C, T2 vs C, and T1 vs T2 (ITT)
    pv_T1 <- tryCatch(ct["Z_T1", "Pr(>|t|)"], error = function(e) NA_real_)
    pv_T2 <- tryCatch(ct["Z_T2", "Pr(>|t|)"], error = function(e) NA_real_)
    # T1 - T2 contrast
    R <- rep(0, length(coef(fit_itt))); names(R) <- names(coef(fit_itt))
    if ("Z_T1" %in% names(R)) R["Z_T1"] <- 1
    if ("Z_T2" %in% names(R)) R["Z_T2"] <- -1
    wt <- wald_contrast(fit_itt, R = R, r = 0, cluster = cluster)
    pv_T1T2 <- wt$p

    # TOT models using actual treatment receipt

    # T1 vs Control sample
    df_T1C <- df_use |>
      filter(arm %in% c("control", arms[2])) |>
      mutate(Z = as.integer(arm == arms[2]))
    if (nrow(df_T1C) > 0) {
      if (is_binary) {
        fit_tot_T1 <- AER::ivreg(y ~ D + month_f | Z + month_f, data = df_T1C)
      } else {
        y_tr_T1 <- log(pmax(df_T1C$y, 1e-8))
        fit_tot_T1 <- AER::ivreg(y_tr_T1 ~ D + month_f | Z + month_f, data = df_T1C)
      }
  vc1 <- get_vcov(fit_tot_T1, cluster = if (cluster_se) df_T1C$association_id else NULL, hc_type = hc_type)
      ct1 <- lmtest::coeftest(fit_tot_T1, vcov. = vc1)
      pv_tot_T1 <- tryCatch(ct1["D", "Pr(>|t|)"], error = function(e) NA_real_)
    } else pv_tot_T1 <- NA_real_

    # T2 vs Control sample (if exists)
    if (length(arms) >= 3) {
      df_T2C <- df_use |>
        filter(arm %in% c("control", arms[3])) |>
        mutate(Z = as.integer(arm == arms[3]))
      if (nrow(df_T2C) > 0) {
        if (is_binary) {
          fit_tot_T2 <- AER::ivreg(y ~ D + month_f | Z + month_f, data = df_T2C)
        } else {
          y_tr_T2 <- log(pmax(df_T2C$y, 1e-8))
          fit_tot_T2 <- AER::ivreg(y_tr_T2 ~ D + month_f | Z + month_f, data = df_T2C)
        }
  vc2 <- get_vcov(fit_tot_T2, cluster = if (cluster_se) df_T2C$association_id else NULL, hc_type = hc_type)
        ct2 <- lmtest::coeftest(fit_tot_T2, vcov. = vc2)
        pv_tot_T2 <- tryCatch(ct2["D", "Pr(>|t|)"], error = function(e) NA_real_)
      } else pv_tot_T2 <- NA_real_
    } else pv_tot_T2 <- NA_real_

    # Joint multi-arm TOT across full sample: y ~ D1 + D2 + FE | Z1 + Z2 + FE
    if (is_binary) {
      fit_tot_joint <- AER::ivreg(y ~ D1 + D2 + month_f | Z1 + Z2 + month_f, data = df_use)
    } else {
      y_tr_joint <- log(pmax(df_use$y, 1e-8))
      fit_tot_joint <- AER::ivreg(y_tr_joint ~ D1 + D2 + month_f | Z1 + Z2 + month_f, data = df_use)
    }
  vcj <- get_vcov(fit_tot_joint, cluster = if (cluster_se) df_use$association_id else NULL, hc_type = hc_type)
    ctj <- lmtest::coeftest(fit_tot_joint, vcov. = vcj)
    pv_totJ_T1 <- tryCatch(ctj["D1", "Pr(>|t|)"], error = function(e) NA_real_)
    pv_totJ_T2 <- tryCatch(ctj["D2", "Pr(>|t|)"], error = function(e) NA_real_)

    pvals <- tibble(
      p_itt_T1_vs_C = pv_T1,
      p_itt_T2_vs_C = pv_T2,
      p_itt_T1_vs_T2 = pv_T1T2,
      p_tot_T1_vs_C = pv_tot_T1,
      p_tot_T2_vs_C = pv_tot_T2,
      p_tot_joint_T1_vs_C = pv_totJ_T1,
      p_tot_joint_T2_vs_C = pv_totJ_T2
    )
    if (isTRUE(capture_data)) {
      return(list(pvals = pvals, df = df_use, outcome = outcome_selected))
    } else {
      return(pvals)
    }
  }

  stopifnot(!is.null(sweep_values), length(sweep_values) >= 1)
  trt_arms <- setdiff(arms, "control")

  # Progress over sweep values (and arms if applicable)
  results_list <- list()
  total_steps <- if (sweep_param == "effect_mult" && sweep_each_arm) length(trt_arms) * length(sweep_values) else length(sweep_values)

  # Function to run Monte Carlo sims for a given sweep value
  run_for_val <- function(val, which_arm = NULL) {
    cfg_local <- config
    if (sweep_param == "effect_mult") {
      # Set per-arm multiplicative effect on the selected outcome (probability-scale multiplier for binary)
      eff <- eff_mult
      if (is.null(which_arm)) {
        for (a in trt_arms) eff[[a]] <- val
      } else {
        eff[[which_arm]] <- val
      }
      cfg_local$effect_mult <- eff
    } else if (sweep_param == "cows_tagged_per_association") {
      # Set number of tagged cows per association (rounded to integer >= 1)
      cfg_local$cows_tagged_per_association <- as.integer(pmax(1, round(val)))
    }

    # Monte Carlo
    if (parallel_inner) {
      # Use furrr's built-in progress to avoid cross-process progressor binding issues
      sim_res <- furrr::future_map(
        1:sims,
        ~ one_run(cfg_local, capture_data = FALSE),
        .options = furrr::furrr_options(seed = TRUE),
        .progress = sim_progress
      )
    } else if (print_sim_numbers) {
      sim_res <- vector("list", sims)
      for (i in seq_len(sims)) {
        # Print a concise simulation counter so users can see progress even when
        # not using a full progress bar. Use message() so output appears in logs
        # and is flush-friendly across environments.
        message(sprintf("Sim %d/%d (sweep=%s)", i, sims, as.character(val)))
        sim_res[[i]] <- one_run(cfg_local, capture_data = FALSE)
      }
    } else {
      sim_res <- purrr::map(1:sims, ~ one_run(cfg_local, capture_data = FALSE))
    }

    dfp <- bind_rows(sim_res)
    tibble(
      sweep_value = val,
      power_itt_T1_vs_C = mean(dfp$p_itt_T1_vs_C < alpha, na.rm = TRUE),
      power_itt_T2_vs_C = mean(dfp$p_itt_T2_vs_C < alpha, na.rm = TRUE),
      power_itt_T1_vs_T2 = mean(dfp$p_itt_T1_vs_T2 < alpha, na.rm = TRUE),
      power_tot_T1_vs_C = mean(dfp$p_tot_T1_vs_C < alpha, na.rm = TRUE),
      power_tot_T2_vs_C = mean(dfp$p_tot_T2_vs_C < alpha, na.rm = TRUE),
      power_tot_joint_T1_vs_C = mean(dfp$p_tot_joint_T1_vs_C < alpha, na.rm = TRUE),
      power_tot_joint_T2_vs_C = mean(dfp$p_tot_joint_T2_vs_C < alpha, na.rm = TRUE)
    )
  }

  progressr::with_progress({
    p <- progressr::progressor(along = total_steps)

    if (parallel_outer) {
      # Outer loop parallelization: distribute sweep values across workers
      if (sweep_param == "effect_mult" && sweep_each_arm) {
        # Create task list: combination of arms and values
        tasks <- expand.grid(arm = trt_arms, val = sweep_values, stringsAsFactors = FALSE)
        results_list <- furrr::future_map2(
          tasks$val, tasks$arm,
          function(v, a) {
            p(message = glue::glue("sweep {sweep_param}={format(v, digits=4)} arm={a}"))
            run_for_val(v, which_arm = a) |> mutate(swept_arm = a)
          },
          .options = furrr::furrr_options(seed = TRUE),
          .progress = FALSE  # Manual progress via p()
        )
      } else {
        # Single sweep over values
        results_list <- furrr::future_map(
          sweep_values,
          function(v) {
            p(message = glue::glue("sweep {sweep_param}={format(v, digits=4)}"))
            run_for_val(v, which_arm = NULL)
          },
          .options = furrr::furrr_options(seed = TRUE),
          .progress = FALSE  # Manual progress via p()
        )
      }
    } else {
      # Sequential outer loop (inner may still be parallel)
      if (sweep_param == "effect_mult" && sweep_each_arm) {
        for (a in trt_arms) {
          for (v in sweep_values) {
            p(message = glue::glue("sweep {sweep_param}={format(v, digits=4)} arm={a}"))
            results_list[[length(results_list) + 1]] <- run_for_val(v, which_arm = a) |>
              mutate(swept_arm = a)
          }
        }
      } else {
        for (v in sweep_values) {
          p(message = glue::glue("sweep {sweep_param}={format(v, digits=4)}"))
          results_list[[length(results_list) + 1]] <- run_for_val(v, which_arm = NULL)
        }
      }
    }
  })

  results <- bind_rows(results_list)

  # Save outputs
  csv_path <- file.path(tabs_dir, glue::glue("{outfile_stem}_{outcome_selected}.csv"))
  readr::write_csv(results, csv_path)

  # Long format for downstream (ensure swept_arm exists for consistent ordering)
  if (!("swept_arm" %in% names(results))) results <- results |> mutate(swept_arm = NA_character_)
  long_res <- results |>
    pivot_longer(cols = starts_with("power_"), names_to = "series", values_to = "power") |>
    arrange(sweep_value, swept_arm)
  long_csv_path <- file.path(tabs_dir, glue::glue("{outfile_stem}_{outcome_selected}_long.csv"))
  readr::write_csv(long_res, long_csv_path)

  # Persist a small metadata table describing the run (CSV is canonical).
  # Build a small metadata table and write it alongside the exported CSV as a companion
  # key/value CSV; if available, also write an Excel workbook with two sheets (long + meta).
  meta <- list(
    outfile_stem = outfile_stem,
    outcome = outcome_selected,
    seed = seed,
    sims = sims,
    n_associations = n_assoc,
    cows_per_association = cows_per_assoc,
    cows_tagged_per_association = cows_tagged,
    months_T = months_T,
    events_per_month_E = events_E,
    analysis_mode = analysis_mode,
    month_aggregate_mode = agg_mode_cont,
    cluster_se = cluster_se,
    arms = paste(arms, collapse = ","),
    alloc_ratios = paste(names(alloc_ratios), "=", alloc_ratios, collapse = ", "),
    take_up = paste(names(take_up), "=", take_up, collapse = ", "),
    sigma_assoc = sigma_assoc,
    sigma_cow = sigma_cow,
    rho_ar1 = rho_ar1,
    sigma_ar = sigma_ar,
    rho_month = rho_month,
    sigma_month = sigma_month,
    het_var_ngo = het_var$ngo,
    het_var_tribe = het_var$tribe,
    het_var_year = het_var$year,
    missingness_p_cow_month = p_miss,
    effect_multipliers = paste(names(eff_mult), "=", eff_mult, collapse = ", "),
    # Additional provenance and simulation parameters used during estimation
    alpha = alpha,
    hc_type = hc_type,
    parallel_layer = parallel_layer,
    sweep_param = sweep_param,
    sweep_values = paste(sweep_values, collapse = ","),
    sweep_each_arm = sweep_each_arm,
    output_root = ifelse(is.null(output_root), "", output_root),
    print_sim_numbers = print_sim_numbers,
    sim_progress = sim_progress,
    agg_mode_cont = agg_mode_cont,
    # Summarize stratifier assignments (unique levels or samples)
    ngo_levels = paste(unique(ngo), collapse = ","),
    tribe_levels = paste(unique(tribe), collapse = ","),
    year_levels = paste(unique(year), collapse = ",")
  )
  meta_vec <- unlist(meta)
  meta_df <- tibble::tibble(parameter = names(meta_vec), value = as.character(meta_vec))
  meta_path <- file.path(tabs_dir, glue::glue("{outfile_stem}_{outcome_selected}_meta.csv"))
  # Robustly attempt to write meta CSV; fall back to write.table if readr fails
  tryCatch(
    {
      readr::write_csv(meta_df, meta_path)
      message("Wrote meta CSV: ", meta_path)
    },
    error = function(e) {
      message("Warning: failed to write meta CSV with readr: ", conditionMessage(e), " — attempting fallback write.table")
      tryCatch(
        write.table(meta_df, file = meta_path, sep = ",", row.names = FALSE, col.names = TRUE),
        error = function(e2) message("Fallback write.table also failed: ", conditionMessage(e2))
      )
    }
  )
  # Also write a copy of the meta CSV to the repository-relative output path (cwd-based)
  repo_meta_path <- file.path(getwd(), "output", "cattle", "tabs", glue::glue("{outfile_stem}_{outcome_selected}_meta.csv"))
  tryCatch({
    readr::write_csv(meta_df, repo_meta_path)
    message("Wrote repo-local meta CSV: ", repo_meta_path)
  }, error = function(e) {
    message("Warning: failed to write repo-local meta CSV: ", conditionMessage(e))
  })

  # If writexl is available, write a small workbook with two sheets (long + meta) for Excel users
  xlsx_path <- file.path(tabs_dir, glue::glue("{outfile_stem}_{outcome_selected}.xlsx"))
  if (requireNamespace("writexl", quietly = TRUE)) {
    tryCatch(
      writexl::write_xlsx(list(long = long_res, meta = meta_df), path = xlsx_path),
      error = function(e) NULL
    )
  }

  # Read the tidy CSV back from disk and use that as the plotting data source so the CSV is canonical
  long_res <- readr::read_csv(long_csv_path, show_col_types = FALSE)

  # Plot (only ITT series in the figure)
  plot_res <- long_res |>
    filter(grepl("^power_itt", series)) |>
    mutate(series = case_when(
      series == "power_itt_T1_vs_C" ~ "ITT: T1 vs Control",
      series == "power_itt_T2_vs_C" ~ "ITT: T2 vs Control",
      series == "power_itt_T1_vs_T2" ~ "ITT: T1 vs T2",
      TRUE ~ series
    ))

  series_levels <- unique(plot_res$series)
  cols <- setNames(get_palette(length(series_levels), "Dark2"), series_levels)
  # Build a left-justified, multi-line footnote with main simulation parameters (grouped, verbose labels)
  # Grouping helps readers quickly find design, stratifiers, variance, missingness and effect settings.
  foot_lines <- c(
    # Run-level
    glue::glue("Run: sims={sims}    |    Outcome={outcome_selected}    |    Sweep={sweep_param}"),

    # Design
    glue::glue("Design: associations={n_assoc}, cows_per_association={cows_per_assoc}, cows_tagged={cows_tagged}, months_T={months_T}, events_per_month_E={events_E}"),

    # Stratifiers
    glue::glue("Stratifiers: ngo_levels={length(unique(ngo))}, tribe_levels={length(unique(tribe))}, year_levels={length(unique(year))}"),

    # Variance components & temporal correlation
    glue::glue("Variance: sigma_assoc={format(sigma_assoc, digits=3)}, sigma_cow={format(sigma_cow, digits=3)}, rho_ar1={format(rho_ar1, digits=3)}, sigma_ar={format(sigma_ar, digits=3)}"),

    # Heterogeneity of treatment effects
    glue::glue("Heterogeneity (TE var): ngo={format(het_var$ngo, digits=3)}, tribe={format(het_var$tribe, digits=3)}, year={format(het_var$year, digits=3)}"),

    # Missingness and inference
    glue::glue("Missingness: p_cow_month={format(p_miss, digits=3)}    |    Cluster SE={cluster_se}    |    HC type={hc_type}"),

    # Treatment effects
    glue::glue("Effect multipliers: {paste0(names(eff_mult), '=', eff_mult, collapse = ', ')}"),

    # Parallelism / reproducibility
    glue::glue("Parallel layer: {parallel_layer}    |    seed={seed}")
  )
  footnote_text <- paste(foot_lines, collapse = "\n")

  plt <- ggplot(plot_res, aes(x = sweep_value, y = power, color = series)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = cols) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
    theme_minimal(base_size = 12) +
    labs(
      title = glue::glue("Power vs Effect for {outcome_selected}"),
      y = "Power",
      color = "Legend",
      caption = footnote_text
    ) +
    theme(
      plot.caption = element_text(hjust = 0, size = 8),
      plot.caption.position = "plot"
    )
  # Dynamic X label based on sweep_param
  xlab <- if (sweep_param == "effect_mult") "Effect (multiplier; probability-scale for binary)" else if (sweep_param == "cows_tagged_per_association") "Cows tagged per association" else "Sweep value"
  plt <- plt + labs(x = xlab)
  png_path <- file.path(figs_dir, glue::glue("{outfile_stem}_{outcome_selected}.png"))
  pdf_path <- file.path(figs_dir, glue::glue("{outfile_stem}_{outcome_selected}.pdf"))
  ggsave(png_path, plt, width = 8, height = 5, dpi = 150)
  ggsave(pdf_path, plt, width = 8, height = 5)
  last_val <- sweep_values[length(sweep_values)]
  last_arm <- if (sweep_param == "effect_mult" && isTRUE(sweep_each_arm)) tail(trt_arms, 1) else NULL
  cfg_last <- config
  if (sweep_param == "effect_mult") {
    eff <- eff_mult
    if (is.null(last_arm)) {
      for (a in trt_arms) eff[[a]] <- last_val
    } else {
      eff[[last_arm]] <- last_val
    }
    cfg_last$effect_mult <- eff
  } else if (sweep_param == "cows_tagged_per_association") {
    cfg_last$cows_tagged_per_association <- as.integer(pmax(1, round(last_val)))
  }
  last_run <- one_run(cfg_last, capture_data = TRUE)
  if (is.list(last_run) && !is.null(last_run$df)) {
    last_data_csv <- file.path(tabs_dir, glue::glue("{outfile_stem}_{outcome_selected}_last_sim_data.csv"))
    readr::write_csv(last_run$df, last_data_csv)
  }

  message("Exported results to:\n", csv_path, "\n", long_csv_path, "\n", png_path, "\n", pdf_path)

  list(results = results, csv = csv_path, long_csv = long_csv_path, png = png_path, pdf = pdf_path)
}
