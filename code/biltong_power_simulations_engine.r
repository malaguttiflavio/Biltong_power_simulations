# ------------------------------ Core Simulator ----------------------------- #

# Parallelized using {furrr} for multicore processing
# Uses 14 cores based on MacBook Pro M4 Pro
# Supports multi-arm designs: arms[1] is control; per-arm effects via soc_outcome_ate_pct_* and env_outcome_ate_*
simulate_power <- function(
  config,
  sweep_param = c(
    # generic
  "n_communities","avg_ind_obs_per_comm","alloc_ratio",
    # canonical soc/env names
    "soc_outcome_T","soc_outcome_ICC","soc_outcome_AR1_rho","soc_outcome_AR1_var",
    "env_outcome_T","env_outcome_ICC","env_outcome_AR1_rho","env_outcome_AR1_var",
    # convenience aliases from master script
    "soc_outcome_ate_pct_single","env_outcome_ate_single"
  ),
  sweep_arm = NULL,
  sweep_values = NULL,
  sweep_each_arm = FALSE,  # when TRUE and sweep_param is arm-level (e.g. soc_outcome_ate_pct_*), run an independent sweep for each treatment arm
  parallel_layer = c("inner","outer","both","none"), # choose which layer to parallelize
  outfile_stem = "biltong_power",
  seed = 1
) {
  set.seed(seed)
  parallel_layer <- match.arg(parallel_layer)
  parallel_inner <- parallel_layer %in% c("inner","both")
  parallel_outer <- parallel_layer %in% c("outer","both")
  if (identical(parallel_layer, "both")) {
    message("[note] parallel_layer='both' enables nested parallelism; ensure your future plan avoids oversubscription.")
  }
  # Helper for palettes with graceful fallback
  get_palette <- function(n, pal_name = "Dark2", fallback_colors = NULL) {
    if (is.null(fallback_colors)) fallback_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      return(fallback_colors[seq_len(min(n, length(fallback_colors)))] )
    }
    max_in_set <- suppressWarnings(RColorBrewer::brewer.pal.info[pal_name, "maxcolors"])
    if (is.na(max_in_set)) return(fallback_colors[seq_len(min(n, length(fallback_colors)))])
    k <- min(max_in_set, max(3, n))
    cols <- RColorBrewer::brewer.pal(k, pal_name)[seq_len(n)]
    cols
  }
  sweep_param <- match.arg(sweep_param)
  # ---------------------- Basic multi-arm validation & normalization ---------------------- #
  normalize_config <- function(cfg) {
    # Ensure arms include control as first element
    stopifnot(!is.null(cfg$arms), length(cfg$arms) >= 2, cfg$arms[1] == "control")
    treatment_arms <- setdiff(cfg$arms, "control")

    # Helper: normalize allocation ratios for 2 or 3 arms (or more if extended later)
    process_alloc <- function(cfg) {
      arms <- cfg$arms
      trt_arms <- setdiff(arms, "control")
      if (is.null(cfg$alloc_ratios)) {
        cfg$alloc_ratios <- setNames(rep(1/length(arms), length(arms)), arms)
        return(cfg)
      }
      ar <- cfg$alloc_ratios
      # If unnamed, try to assign
      if (is.null(names(ar))) {
        if (length(ar) == length(arms)) {
          names(ar) <- arms
        } else if (length(ar) == length(trt_arms)) {
          names(ar) <- trt_arms
          rem <- 1 - sum(ar)
          if (rem <= 0) stop("Cannot infer control allocation; remainder is non-positive")
          ar <- c(control = rem, ar)
        } else {
          stop("Unnamed alloc_ratios length (", length(ar), ") does not match arms (", length(arms), ") or treatment arms (", length(trt_arms), ")")
        }
      }
      # If control missing, attempt to infer
      if (!"control" %in% names(ar)) {
        missing_trt <- setdiff(trt_arms, names(ar))
        if (length(missing_trt) > 0) stop("alloc_ratios missing treatment arms: ", paste(missing_trt, collapse=","))
        rem <- 1 - sum(ar)
        if (rem <= 0) stop("Cannot infer control allocation; non-positive remainder after summing provided treatment allocations")
        ar <- c(control = rem, ar)
      }
      # Reorder & trim to specified arms (ignore any extraneous names)
      ar <- ar[arms]
      if (any(is.na(ar))) stop("alloc_ratios produced NA after aligning to arms")
      if (any(!is.finite(ar))) stop("alloc_ratios contain non-finite values")
      total <- sum(ar)
      # If total not close to 1 treat as weights (e.g., percents or counts) and normalize
      if (abs(total - 1) > 1e-6) {
        ar <- ar / total
        total <- 1
      }
      # Guard against extremely small or zero allocations
      if (any(ar <= 0)) stop("All alloc_ratios must be strictly positive after normalization")
      # Final precise normalization against drift
      ar <- ar / sum(ar)
      cfg$alloc_ratios <- ar
      cfg
    }
    cfg <- process_alloc(cfg)

    # Default stratification settings: if not provided, and canonical stratifier columns exist,
    # use all three to stratify randomization at the community level (applies to both experiment types).
    if (is.null(cfg$stratify_by)) {
      if (all(c("year_in_program","ngo_id","tribe_id") %in% names(cfg))) {
        cfg$stratify_by <- c("year_in_program","ngo_id","tribe_id")
      }
    }
    if (is.null(cfg$stratify_exact)) cfg$stratify_exact <- TRUE

    # alloc_ratios: must be named and cover all arms, sum to ~1
    # (Already normalized above)

    # take_up: must include all arms; provide defaults if absent
    if (is.null(cfg$take_up)) {
      cfg$take_up <- setNames(rep(1, length(cfg$arms)), cfg$arms)
      cfg$take_up["control"] <- 0
    } else {
      if (is.null(names(cfg$take_up))) stop("take_up must be a *named* numeric vector")
      missing_tu <- setdiff(cfg$arms, names(cfg$take_up))
      if (length(missing_tu) > 0) {
        # Assume full take-up for missing treatment arms, 0 for control if missing
        add_vals <- ifelse(missing_tu == "control", 0, 1)
        cfg$take_up <- c(cfg$take_up, setNames(add_vals, missing_tu))
      }
    }

    # Use only soc_* and env_* names; legacy aliasing removed
    # Ensure environmental ATE present under current names; default to zeros per treatment arm
    if (is.null(cfg$env_outcome_ate_single) && is.null(cfg$env_outcome_ate_multi)) {
      cfg$env_outcome_ate_single <- setNames(rep(0, length(treatment_arms)), treatment_arms)
    } else {
      # Normalize any provided container to be named and exclude control
      if (!is.null(cfg$env_outcome_ate_single)) {
        if (is.null(names(cfg$env_outcome_ate_single)) && length(cfg$env_outcome_ate_single) == 1 && length(treatment_arms) == 1) {
          cfg$env_outcome_ate_single <- setNames(cfg$env_outcome_ate_single, treatment_arms)
        }
        missing_vdm <- setdiff(treatment_arms, names(cfg$env_outcome_ate_single))
        if (length(missing_vdm) > 0) cfg$env_outcome_ate_single <- c(cfg$env_outcome_ate_single, setNames(rep(0, length(missing_vdm)), missing_vdm))
        cfg$env_outcome_ate_single <- cfg$env_outcome_ate_single[setdiff(names(cfg$env_outcome_ate_single), "control")]
      }
      if (!is.null(cfg$env_outcome_ate_multi)) {
        if (is.null(names(cfg$env_outcome_ate_multi)) && length(cfg$env_outcome_ate_multi) == 1 && length(treatment_arms) == 1) {
          cfg$env_outcome_ate_multi <- setNames(cfg$env_outcome_ate_multi, treatment_arms)
        }
        missing_vdm2 <- setdiff(treatment_arms, names(cfg$env_outcome_ate_multi))
        if (length(missing_vdm2) > 0) cfg$env_outcome_ate_multi <- c(cfg$env_outcome_ate_multi, setNames(rep(0, length(missing_vdm2)), missing_vdm2))
        cfg$env_outcome_ate_multi <- cfg$env_outcome_ate_multi[setdiff(names(cfg$env_outcome_ate_multi), "control")]
      }
    }

    # Legacy/VCU keys no longer supported and not referenced

    # Defaults for strata treatment-effect heterogeneity (variances default to 0 = no extra heterogeneity)
    if (is.null(cfg$year_in_program_ate))      cfg$year_in_program_ate <- NULL  # vector or named values per level; optional
    if (is.null(cfg$ngo_id_ate))               cfg$ngo_id_ate <- NULL
    if (is.null(cfg$tribe_id_ate))             cfg$tribe_id_ate <- NULL
    if (is.null(cfg$year_in_program_ate_var))  cfg$year_in_program_ate_var <- 0
    if (is.null(cfg$ngo_id_ate_var))           cfg$ngo_id_ate_var <- 0
    if (is.null(cfg$tribe_id_ate_var))         cfg$tribe_id_ate_var <- 0

  # Cluster FE flag: use 'cluster_fe_yn' exclusively

  # Optional diagnostics toggle (prints/writes per-stratum arm counts and TE draws for last run)
  if (is.null(cfg$print_strata_diag)) cfg$print_strata_diag <- FALSE

    cfg
  }
  config <- normalize_config(config)
  # Default: include stratifier controls in regressions unless disabled
  if (is.null(config$use_strata_controls)) config$use_strata_controls <- TRUE
  # Require ICCs to be explicitly defined
  if (is.null(config$soc_outcome_ICC)) stop("soc_outcome_ICC must be defined in config (separate from AR1 parameters)")
  if (is.null(config$env_outcome_ICC)) stop("env_outcome_ICC must be defined in config (separate from AR1 parameters)")
  # Validate ICC ranges
  if (!is.numeric(config$soc_outcome_ICC) || any(config$soc_outcome_ICC < 0 | config$soc_outcome_ICC >= 1)) {
    stop("soc_outcome_ICC must be in [0, 1)")
  }
  if (!is.numeric(config$env_outcome_ICC) || any(config$env_outcome_ICC < 0 | config$env_outcome_ICC >= 1)) {
    stop("env_outcome_ICC must be in [0, 1)")
  }
  # Defaults for new inference options
  if (is.null(config$cluster_se))     config$cluster_se     <- TRUE
  if (is.null(config$cluster_fe_yn)) config$cluster_fe_yn <- FALSE  # If TRUE keeps factor(community_id) in model

  if (!is.null(sweep_arm)) {
    if (is.null(config$arms)) stop("config$arms must be set for multi-arm sweep")
    if (!(sweep_arm %in% setdiff(config$arms, config$arms[1])))
      stop("sweep_arm must be one of the non-control arms")
  }
  if (is.null(sweep_values) || length(sweep_values) < 1) stop("Provide sweep_values")

  # Strict check for stratifiers
  stopifnot(
    !is.null(config$year_in_program),
    !is.null(config$ngo_id),
    !is.null(config$tribe_id),
    length(config$year_in_program) == config$n_communities,
    length(config$ngo_id)          == config$n_communities,
    length(config$tribe_id)        == config$n_communities
  )

  # Build sample frames
  make_sample <- function(cfg_local) {
    # Community roster
    # Harmonize stratifier vectors with n_communities in case a sweep changed n_communities
    # If lengths mismatch: replicate singletons or sample with replacement to match size
    nC <- cfg_local$n_communities
    yip <- cfg_local$year_in_program
    if (length(yip) != nC) {
      yip <- if (length(yip) == 1) rep(yip, nC) else sample(yip, nC, replace = TRUE)
    }
    ngo <- cfg_local$ngo_id
    if (length(ngo) != nC) {
      ngo <- if (length(ngo) == 1) rep(ngo, nC) else sample(ngo, nC, replace = TRUE)
    }
    tri <- cfg_local$tribe_id
    if (length(tri) != nC) {
      tri <- if (length(tri) == 1) rep(tri, nC) else sample(tri, nC, replace = TRUE)
    }

    comm <- tibble(
      community_id    = 1:nC,
      year_in_program = yip,
      ngo_id          = ngo,
      tribe_id        = tri
    )
    # Individuals per community
    if (cfg_local$sd_indiv_per_comm > 0) {
      n_i <- pmax(1L, round(rnorm(cfg_local$n_communities, cfg_local$avg_ind_obs_per_comm, cfg_local$sd_indiv_per_comm)))
    } else {
      n_i <- rep(cfg_local$avg_ind_obs_per_comm, cfg_local$n_communities)
    }
    comm$N_indiv <- n_i

    if (cfg_local$experiment_type == "community_level") {
      # Participation panel: community-time
  df_soc <- comm |>
        rowwise() |>
  mutate(times = list(cfg_local$soc_outcome_T_months)) |>
        unnest(times) |>
        mutate(time = times, .keep = "unused") |>
        ungroup()
  # Environmental panel: community-time
  df_env <- comm |>
        rowwise() |>
  mutate(times = list(cfg_local$env_outcome_T_months)) |>
        unnest(times) |>
        mutate(time = times, .keep = "unused") |>
        ungroup()
    } else {
      # individual_within_community
  df_soc <- comm |>
        rowwise() |>
  mutate(indiv = list(1:N_indiv), times = list(cfg_local$soc_outcome_T_months)) |>
        unnest(c(indiv, times)) |>
        rename(time = times) |>
        ungroup()
  # Environmental remains at community-time
  df_env <- comm |>
        rowwise() |>
  mutate(times = list(cfg_local$env_outcome_T_months)) |>
        unnest(times) |>
        mutate(time = times, .keep = "unused") |>
        ungroup()
    }

    list(comm = comm, df_soc = df_soc, df_env = df_env)
  }

  # Helpers to retrieve per-arm effects using config variable names directly
  get_soc_ate_pct <- function(cfgx) {
    # Prefer scenario-specific names; fallback to canonical if present
    if (!is.null(cfgx$soc_outcome_ate_pct_single)) return(cfgx$soc_outcome_ate_pct_single)
    if (!is.null(cfgx$soc_outcome_ate_pct_multi))  return(cfgx$soc_outcome_ate_pct_multi)
    NULL
  }
  get_env_add <- function(cfgx) {
    if (!is.null(cfgx$env_outcome_ate_single)) return(cfgx$env_outcome_ate_single)
    if (!is.null(cfgx$env_outcome_ate_multi))  return(cfgx$env_outcome_ate_multi)
    # Fallback default: zeros per treatment arm if neither provided
    trt_arms <- setdiff(cfgx$arms %||% character(), "control")
    if (length(trt_arms) > 0) return(setNames(rep(0, length(trt_arms)), trt_arms))
    NULL
  }

  # Helper: assign arms within strata for balanced randomization
  assign_arms_stratified <- function(comm_df, arms, alloc_ratios, strat_vars = NULL, exact = TRUE) {
    # Ensure required columns exist
    if (!is.null(strat_vars) && length(strat_vars) > 0) {
      missing <- setdiff(strat_vars, names(comm_df))
      if (length(missing) > 0) {
        warning("Stratification variables missing in comm_df: ", paste(missing, collapse=","), "; falling back to unstratified assignment.")
        strat_vars <- NULL
      }
    }
    # Guarantee names and order
    probs <- as.numeric(alloc_ratios[arms])
    if (any(is.na(probs))) stop("alloc_ratios must include all arms for assignment")

    # Split by strata (or single group)
    if (!is.null(strat_vars) && length(strat_vars) > 0) {
      groups <- split(comm_df, f = interaction(comm_df[, strat_vars], drop = TRUE))
    } else {
      groups <- list(all = comm_df)
    }

    out_list <- lapply(groups, function(gdf) {
      n <- nrow(gdf)
      if (n == 0) return(gdf)
      if (exact) {
        # Target counts per arm within stratum using largest remainders method
        target <- probs * n
        base_ct <- floor(target)
        remainder <- n - sum(base_ct)
        if (remainder > 0) {
          frac <- target - base_ct
          # Avoid all-zero frac
          if (sum(frac) <= 1e-12) frac <- rep(1, length(frac))
          add_idx <- sample.int(length(arms), size = remainder, replace = FALSE, prob = frac)
          add_tab <- tabulate(add_idx, nbins = length(arms))
          ct <- base_ct + add_tab
        } else {
          ct <- base_ct
        }
        # Create vector of assignments and shuffle
        arm_vec <- rep(arms, times = ct)
        # In rare cases due to rounding, length may be < n; pad by sampling arms by probs
        if (length(arm_vec) < n) {
          arm_vec <- c(arm_vec, sample(arms, size = n - length(arm_vec), replace = TRUE, prob = probs))
        } else if (length(arm_vec) > n) {
          arm_vec <- arm_vec[seq_len(n)]
        }
        arm_vec <- sample(arm_vec, size = n, replace = FALSE)
      } else {
        arm_vec <- sample(arms, size = n, replace = TRUE, prob = probs)
      }
      gdf$arm <- factor(arm_vec, levels = arms)
      gdf
    })
    assigned <- dplyr::bind_rows(out_list)
    assigned[, c("community_id","arm")]
  }

  # One Monte Carlo run returning p-values for ITT and TOT for both outcomes, per treatment arm
  one_run <- function(config_run, capture_data = FALSE) {
    smp <- make_sample(config_run)

    # ----------------------------- Assignment ----------------------------- #
    # Prefer multi-arm assignment using arms + alloc_ratios; fallback to binary alloc_ratio
    if (!is.null(config_run$alloc_ratios) && !is.null(config_run$arms)) {
      stopifnot(abs(sum(as.numeric(config_run$alloc_ratios[config_run$arms])) - 1) < 1e-8)
      arm_comm <- assign_arms_stratified(
        smp$comm,
        arms = config_run$arms,
        alloc_ratios = config_run$alloc_ratios,
        strat_vars = config_run$stratify_by,
        exact = isTRUE(config_run$stratify_exact)
      )
    } else {
      prob <- if (!is.null(config_run$alloc_ratio)) config_run$alloc_ratio else 0.5
      # Use the same stratified helper with binary fallback arms
      bin_arms <- c("control","T1")
      bin_ar <- c(control = 1 - prob, T1 = prob)
      arm_comm <- assign_arms_stratified(
        smp$comm,
        arms = bin_arms,
        alloc_ratios = bin_ar,
        strat_vars = config_run$stratify_by,
        exact = isTRUE(config_run$stratify_exact)
      )
      if (is.null(config_run$arms)) config_run$arms <- bin_arms
      if (is.null(config_run$alloc_ratios)) config_run$alloc_ratios <- bin_ar
    }

    # ---------------------- Strata treatment-effect heterogeneity ---------------------- #
    # Build per-community additive treatment effects based on stratification variables.
    # For soc_outcome (log scale), these are added to log(mu) when D==1; for env_outcome (level), added to mean when D==1.
    gen_level_effects <- function(levels_vec, vec_param, var_param) {
      lvls <- sort(unique(levels_vec))
      # If user provided explicit vector, prefer it; attempt name-based match first
      if (!is.null(vec_param)) {
        if (!is.null(names(vec_param))) {
          # Match by names; coerce to character for robust matching
          nm <- names(vec_param)
          map <- setNames(rep(0, length(lvls)), as.character(lvls))
          intersecting <- intersect(as.character(lvls), nm)
          map[intersecting] <- vec_param[intersecting]
          return(map)
        } else if (length(vec_param) == length(lvls)) {
          return(setNames(as.numeric(vec_param), as.character(lvls)))
        }
      }
      # Otherwise, draw from N(0, var)
      v <- if (is.null(var_param)) 0 else as.numeric(var_param)
      sdv <- if (v <= 0) 0 else sqrt(v)
      setNames(stats::rnorm(length(lvls), mean = 0, sd = sdv), as.character(lvls))
    }
    yr_map <- gen_level_effects(smp$comm$year_in_program, config_run$year_in_program_ate, config_run$year_in_program_ate_var)
    ngo_map <- gen_level_effects(smp$comm$ngo_id,          config_run$ngo_id_ate,          config_run$ngo_id_ate_var)
    tri_map <- gen_level_effects(smp$comm$tribe_id,        config_run$tribe_id_ate,        config_run$tribe_id_ate_var)
    tau_comm_df <- tibble(
      community_id = smp$comm$community_id,
      tau_strata = as.numeric(yr_map[as.character(smp$comm$year_in_program)]) +
                   as.numeric(ngo_map[as.character(smp$comm$ngo_id)]) +
                   as.numeric(tri_map[as.character(smp$comm$tribe_id)])
    )
    # Optional diagnostics: per-stratum arm counts and realized TE draws
    diag_strata_counts <- NULL
    diag_te_draws <- NULL
    if (isTRUE(config_run$print_strata_diag)) {
      if (!is.null(config_run$stratify_by) && length(config_run$stratify_by) > 0) {
        tmp_comm <- smp$comm |> left_join(arm_comm, by = "community_id")
        grp_vars <- c(config_run$stratify_by, "arm")
        diag_strata_counts <- tmp_comm |>
          group_by(across(all_of(grp_vars))) |>
          summarize(n = n(), .groups = "drop") |>
          arrange(across(all_of(config_run$stratify_by)), arm)
      }
      diag_te_draws <- bind_rows(
        tibble(variable = "year_in_program", level = names(yr_map), effect = as.numeric(yr_map)),
        tibble(variable = "ngo_id",          level = names(ngo_map), effect = as.numeric(ngo_map)),
        tibble(variable = "tribe_id",        level = names(tri_map), effect = as.numeric(tri_map))
      )
    }

    if (config_run$experiment_type == "community_level") {
      df_soc <- smp$df_soc |> left_join(arm_comm, by = "community_id")
      df_env <- smp$df_env |> left_join(arm_comm, by = "community_id")
    } else {
      # individual_within_community: inherit the community arm
      df_soc <- smp$df_soc |> left_join(arm_comm, by = "community_id")
      df_env <- smp$df_env |> left_join(arm_comm, by = "community_id")
    }
    # Attach strata TE to both panels
    df_soc <- df_soc |> left_join(tau_comm_df, by = "community_id")
    df_env <- df_env |> left_join(tau_comm_df, by = "community_id")

    # ----------------------------- Compliance ----------------------------- #
    tu_vec <- config_run$take_up
    if (is.null(tu_vec)) {
      tu_vec <- setNames(rep(1, length(config_run$arms)), config_run$arms)
      if ("control" %in% names(tu_vec)) tu_vec["control"] <- 0
    }
    if (!all(config_run$arms %in% names(tu_vec))) stop("take_up must be named for all arms")

    if (config_run$experiment_type == "community_level") {
      Dc <- tibble(community_id = smp$comm$community_id) |>
        left_join(arm_comm, by = "community_id") |>
        mutate(D = rbinom(n(), 1, tu_vec[as.character(arm)])) |>
        select(community_id, D)
      df_soc <- df_soc |> left_join(Dc, by = "community_id")
      df_env <- df_env |> left_join(Dc, by = "community_id")
    } else {
      df_soc <- df_soc |> mutate(D = rbinom(n(), 1, tu_vec[as.character(arm)]))
      agg <- df_soc |>
        group_by(community_id) |>
        summarize(D_bar = mean(D), .groups = "drop")
      df_env <- df_env |> left_join(agg, by = "community_id")
      df_env$D <- df_env$D_bar
    }

    # ------------------- Participation outcome generation ------------------ #
    if (config_run$experiment_type == "community_level") {
      # AR(1) innovation variance for socio-economic outcome, defaulting to 1 if not provided
      var_e_p <- if (!is.null(config_run$soc_outcome_AR1_var)) config_run$soc_outcome_AR1_var else 1
      # Between-community random intercept variance: decoupled from AR1 variance.
      # Use ICC alone on its own scale: var_u = ICC / (1 - ICC)
      if (config_run$soc_outcome_ICC >= 1 || config_run$soc_outcome_ICC < 0) stop("soc_outcome_ICC must be in [0,1)")
      var_u_p <- if (config_run$soc_outcome_ICC == 0) 0 else (config_run$soc_outcome_ICC / (1 - config_run$soc_outcome_ICC))
      u_comm_p <- rnorm(config_run$n_communities, 0, sqrt(pmax(var_u_p, 0)))
      sigma_e_p <- if (var_e_p <= 0) 0 else sqrt(var_e_p)
      df_soc <- df_soc |>
        left_join(tibble(community_id = 1:config_run$n_communities, u_comm_p), by = "community_id") |>
        group_by(community_id) |>
        arrange(time, .by_group = TRUE) |>
        mutate(e_ar = ar1_series(time, rho = config_run$soc_outcome_AR1_rho, sigma = sigma_e_p)) |>
        ungroup()

      log_mu_base <- log(pmax(config_run$soc_outcome_base_mean, 1e-6)) + df_soc$u_comm_p + df_soc$e_ar
    # Per-arm multiplicative effects, realized only when D==1
  mrr  <- get_soc_ate_pct(config_run)
      mult <- rep(1, nrow(df_soc))
      if (!is.null(mrr)) for (nm in names(mrr)) mult[df_soc$arm == nm] <- mult[df_soc$arm == nm] * mrr[[nm]]
  log_mu <- log_mu_base + log(mult) * df_soc$D + df_soc$tau_strata * df_soc$D
      mu <- pmax(exp(log_mu), 1e-6)
      if (identical(tolower(config_run$soc_outcome_dist), "none")) {
        # Deterministic: use expected mean directly
        df_soc$Y <- mu
      } else {
        df_soc$Y <- r_count(length(mu), mu, family = config_run$soc_outcome_dist, theta = config_run$soc_outcome_theta)
      }
    } else {
      # individual_within_community
      var_e_p <- if (!is.null(config_run$soc_outcome_AR1_var)) config_run$soc_outcome_AR1_var else 1
      if (config_run$soc_outcome_ICC >= 1 || config_run$soc_outcome_ICC < 0) stop("soc_outcome_ICC must be in [0,1)")
      var_u_p <- if (config_run$soc_outcome_ICC == 0) 0 else (config_run$soc_outcome_ICC / (1 - config_run$soc_outcome_ICC))
      u_comm_p <- rnorm(config_run$n_communities, 0, sqrt(pmax(var_u_p, 0)))
      sigma_e_p <- if (var_e_p <= 0) 0 else sqrt(var_e_p)
      df_soc <- df_soc |>
        left_join(tibble(community_id = 1:config_run$n_communities, u_comm_p), by = "community_id") |>
        group_by(community_id, indiv) |>
        arrange(time, .by_group = TRUE) |>
        mutate(e_ar = ar1_series(time, rho = config_run$soc_outcome_AR1_rho, sigma = sigma_e_p)) |>
        ungroup()

      log_mu_base <- log(pmax(config_run$soc_outcome_base_mean, 1e-6)) + df_soc$u_comm_p + df_soc$e_ar
  mrr  <- get_soc_ate_pct(config_run)
      mult <- rep(1, nrow(df_soc))
      if (!is.null(mrr)) for (nm in names(mrr)) mult[df_soc$arm == nm] <- mult[df_soc$arm == nm] * mrr[[nm]]
  log_mu <- log_mu_base + log(mult) * df_soc$D + df_soc$tau_strata * df_soc$D
      mu <- pmax(exp(log_mu), 1e-6)
      if (identical(tolower(config_run$soc_outcome_dist), "none")) {
        df_soc$Y <- mu
      } else {
        df_soc$Y <- r_count(length(mu), mu, family = config_run$soc_outcome_dist, theta = config_run$soc_outcome_theta)
      }
    }

    # ----------------------------- Participation OLS ----------------------------- #
    df_soc$arm <- factor(df_soc$arm, levels = config_run$arms)
    # Optional stratifier controls
    strata_terms <- NULL
    if (isTRUE(config_run$use_strata_controls)) {
      st <- c("year_in_program","ngo_id","tribe_id")
      st <- st[st %in% names(df_soc)]
      if (length(st) > 0) {
        strata_terms <- paste0("factor(", st, ")")
      }
    }
    base_terms <- c("arm", "factor(time)", strata_terms)
    rhs_no_fe <- paste(base_terms[!is.na(base_terms) & base_terms != ""], collapse = " + ")
    form_p <- if (isTRUE(config_run$cluster_fe_yn)) {
      as.formula(paste("Y ~", paste(c(rhs_no_fe, "factor(community_id)"), collapse = " + ")))
    } else {
      as.formula(paste("Y ~", rhs_no_fe))
    }
    fit_itt_p <- lm(form_p, data = df_soc)
    if (isTRUE(config_run$cluster_se)) {
      vc_p <- tryCatch(sandwich::vcovCL(fit_itt_p, cluster = ~ community_id), error = function(e) sandwich::vcovHC(fit_itt_p, type = config_run$hc_type))
      ct_p <- lmtest::coeftest(fit_itt_p, vcov. = vc_p)
      tt_itt_p <- tibble(term = rownames(ct_p), estimate = ct_p[,1], std_error = ct_p[,2], statistic = ct_p[,3], p_value = ct_p[,4])
    } else {
      tt_itt_p  <- robust_test(fit_itt_p, hc_type = config_run$hc_type)
    }

    # Collect ITT & TOT p-values per treatment arm
    trt_arms <- setdiff(config_run$arms, "control")
    p_itt_p_vec <- setNames(rep(NA_real_, length(trt_arms)), trt_arms)
    p_tot_p_vec <- setNames(rep(NA_real_, length(trt_arms)), trt_arms)

    # ITT (OLS) p-values per arm
    for (a in trt_arms) {
      term_a <- paste0("arm", a)
      p_itt_p_vec[a] <- tt_itt_p$p_value[tt_itt_p$term == term_a]
    }

    # TOT (IV) per arm: run separate IV using only control + that arm for clarity
    for (a in trt_arms) {
      df_sub <- df_soc[df_soc$arm %in% c("control", a), , drop = FALSE]
      # Re-factor for subset to avoid singularities
      df_sub$arm <- droplevels(df_sub$arm)
      df_sub$Z_eval <- as.integer(df_sub$arm == a)
      use_fe_iv_p <- isTRUE(config_run$cluster_fe_yn) && !(config_run$experiment_type == "community_level")
      if (isTRUE(config_run$cluster_fe_yn) && config_run$experiment_type == "community_level") use_fe_iv_p <- FALSE
      strata_rhs <- if (!is.null(strata_terms) && length(strata_terms) > 0) paste("+", paste(strata_terms, collapse = " + ")) else ""
      if (use_fe_iv_p) {
        iv_form_sub <- as.formula(paste(
          "Y ~ D + factor(time)", strata_rhs, "+ factor(community_id) | Z_eval + factor(time)", strata_rhs, "+ factor(community_id)"
        ))
      } else {
        iv_form_sub <- as.formula(paste(
          "Y ~ D + factor(time)", strata_rhs, "| Z_eval + factor(time)", strata_rhs
        ))
      }
      iv_fit_sub <- tryCatch(AER::ivreg(iv_form_sub, data = df_sub), error = function(e) NULL)
      if (!is.null(iv_fit_sub)) {
        if (isTRUE(config_run$cluster_se)) {
          vc_iv_sub <- tryCatch(sandwich::vcovCL(iv_fit_sub, cluster = ~ community_id), error = function(e) sandwich::vcovHC(iv_fit_sub, type = config_run$hc_type))
          ct_iv_sub <- lmtest::coeftest(iv_fit_sub, vcov. = vc_iv_sub)
        } else {
          ct_iv_sub <- lmtest::coeftest(iv_fit_sub, vcov. = sandwich::vcovHC(iv_fit_sub, type = config_run$hc_type))
        }
        p_tot_p_vec[a] <- if ("D" %in% rownames(ct_iv_sub)) ct_iv_sub["D","Pr(>|t|)"] else NA_real_
      } else {
        p_tot_p_vec[a] <- NA_real_
      }
    }

    # ---------------------------- Environmental generation ---------------------------- #
    # AR(1) innovation variance for environmental outcome. If not provided, default to base_sd^2
    var_e_v <- if (!is.null(config_run$env_outcome_AR1_var)) config_run$env_outcome_AR1_var else (config_run$env_outcome_base_sd^2)
    # Between-community variance decoupled from AR1 variance; use ICC only
    if (config_run$env_outcome_ICC >= 1 || config_run$env_outcome_ICC < 0) stop("env_outcome_ICC must be in [0,1)")
    var_u_v <- if (config_run$env_outcome_ICC == 0) 0 else (config_run$env_outcome_ICC / (1 - config_run$env_outcome_ICC))
    u_comm_v <- rnorm(config_run$n_communities, 0, sqrt(pmax(var_u_v, 0)))
    sigma_e_v <- if (var_e_v <= 0) 0 else sqrt(var_e_v)
    df_env <- df_env |>
      left_join(tibble(community_id = 1:config_run$n_communities, u_comm_v), by = "community_id") |>
      group_by(community_id) |>
      arrange(time, .by_group = TRUE) |>
  mutate(e_ar = ar1_series(time, rho = config_run$env_outcome_AR1_rho, sigma = sigma_e_v)) |>
      ungroup()

    # Ensure D exists from compliance block
  if (!"D" %in% names(df_env)) stop("Column D missing in df_env after compliance")
  df_env$arm <- factor(df_env$arm, levels = config_run$arms)

  vdm <- get_env_add(config_run)
    add <- rep(0, nrow(df_env))
    if (!is.null(vdm)) for (nm in names(vdm)) add[df_env$arm == nm] <- add[df_env$arm == nm] + vdm[[nm]]
  df_env$Y <- config_run$env_outcome_base_mean + df_env$u_comm_v + df_env$e_ar + add * df_env$D + df_env$tau_strata * df_env$D

    # ------------------------------- Environmental OLS -------------------------------- #
    # Environmental ITT with optional stratifier controls
    strata_terms_v <- NULL
    if (isTRUE(config_run$use_strata_controls)) {
      st <- c("year_in_program","ngo_id","tribe_id")
      st <- st[st %in% names(df_env)]
      if (length(st) > 0) strata_terms_v <- paste0("factor(", st, ")")
    }
    base_terms_v <- c("arm", "factor(time)", strata_terms_v)
    rhs_no_fe_v <- paste(base_terms_v[!is.na(base_terms_v) & base_terms_v != ""], collapse = " + ")
    form_v <- if (isTRUE(config_run$cluster_fe_yn)) {
      as.formula(paste("Y ~", paste(c(rhs_no_fe_v, "factor(community_id)"), collapse = " + ")))
    } else {
      as.formula(paste("Y ~", rhs_no_fe_v))
    }
    fit_itt_v <- lm(form_v, data = df_env)
    if (isTRUE(config_run$cluster_se)) {
      vc_v <- tryCatch(sandwich::vcovCL(fit_itt_v, cluster = ~ community_id), error = function(e) sandwich::vcovHC(fit_itt_v, type = config_run$hc_type))
      ct_v <- lmtest::coeftest(fit_itt_v, vcov. = vc_v)
      tt_itt_v <- tibble(term = rownames(ct_v), estimate = ct_v[,1], std_error = ct_v[,2], statistic = ct_v[,3], p_value = ct_v[,4])
    } else {
      tt_itt_v  <- robust_test(fit_itt_v, hc_type = config_run$hc_type)
    }
  # ITT per arm (environmental)
    p_itt_v_vec <- setNames(rep(NA_real_, length(trt_arms)), trt_arms)
    for (a in trt_arms) {
      term_a <- paste0("arm", a)
      p_itt_v_vec[a] <- tt_itt_v$p_value[tt_itt_v$term == term_a]
    }
  # TOT per arm (environmental) via separate IV on control + that arm
    p_tot_v_vec <- setNames(rep(NA_real_, length(trt_arms)), trt_arms)
    for (a in trt_arms) {
      df_sub_v <- df_env[df_env$arm %in% c("control", a), , drop = FALSE]
      df_sub_v$arm <- droplevels(df_sub_v$arm)
      df_sub_v$Z_eval <- as.integer(df_sub_v$arm == a)
      use_fe_iv_v <- isTRUE(config_run$cluster_fe_yn) && !(config_run$experiment_type == "community_level")
      if (isTRUE(config_run$cluster_fe_yn) && config_run$experiment_type == "community_level") use_fe_iv_v <- FALSE
      strata_rhs_v <- if (!is.null(strata_terms_v) && length(strata_terms_v) > 0) paste("+", paste(strata_terms_v, collapse = " + ")) else ""
      if (use_fe_iv_v) {
        iv_form_sub_v <- as.formula(paste(
          "Y ~ D + factor(time)", strata_rhs_v, "+ factor(community_id) | Z_eval + factor(time)", strata_rhs_v, "+ factor(community_id)"
        ))
      } else {
        iv_form_sub_v <- as.formula(paste(
          "Y ~ D + factor(time)", strata_rhs_v, "| Z_eval + factor(time)", strata_rhs_v
        ))
      }
      iv_fit_sub_v <- tryCatch(AER::ivreg(iv_form_sub_v, data = df_sub_v), error = function(e) NULL)
      if (!is.null(iv_fit_sub_v)) {
        if (isTRUE(config_run$cluster_se)) {
          vc_iv_sub_v <- tryCatch(sandwich::vcovCL(iv_fit_sub_v, cluster = ~ community_id), error = function(e) sandwich::vcovHC(iv_fit_sub_v, type = config_run$hc_type))
          ct_iv_sub_v <- lmtest::coeftest(iv_fit_sub_v, vcov. = vc_iv_sub_v)
        } else {
          ct_iv_sub_v <- lmtest::coeftest(iv_fit_sub_v, vcov. = sandwich::vcovHC(iv_fit_sub_v, type = config_run$hc_type))
        }
        p_tot_v_vec[a] <- if ("D" %in% rownames(ct_iv_sub_v)) ct_iv_sub_v["D","Pr(>|t|)"] else NA_real_
      } else {
        p_tot_v_vec[a] <- NA_real_
      }
    }

    # Return list elements with per-arm naming
    out_list <- list()
    for (a in trt_arms) {
      out_list[[paste0("p_itt_p_", a)]]  <- p_itt_p_vec[a]
      out_list[[paste0("p_tot_p_", a)]]  <- p_tot_p_vec[a]
      out_list[[paste0("p_itt_v_", a)]]  <- p_itt_v_vec[a]
      out_list[[paste0("p_tot_v_", a)]]  <- p_tot_v_vec[a]
    }
    if (isTRUE(capture_data)) {
      out_list$last_df_soc <- df_soc
      out_list$last_df_env <- df_env
      if (!is.null(diag_strata_counts)) out_list$diag_strata_counts <- diag_strata_counts
      if (!is.null(diag_te_draws)) out_list$diag_te_draws <- diag_te_draws
    }
    out_list
  }

  # Run sims for a given sweep value; override_arm allows independent sweeps per treatment arm
  run_for_value <- function(val, override_arm = NULL) {
    config_sweep <- config
    # Socio-economic ATE sweep: use config-native names only
    if (sweep_param %in% c("soc_outcome_ate_pct_single","soc_outcome_ate_pct_multi")) {
      # Choose which container to update based on presence/preference
      mrr_container <- if (!is.null(config_sweep$soc_outcome_ate_pct_single) || sweep_param == "soc_outcome_ate_pct_single") {
        "soc_outcome_ate_pct_single"
      } else if (!is.null(config_sweep$soc_outcome_ate_pct_multi) || sweep_param == "soc_outcome_ate_pct_multi") {
        "soc_outcome_ate_pct_multi"
      } else stop("soc_outcome_ate_pct_* must be defined in config")
      local_arm <- if (!is.null(override_arm)) override_arm else sweep_arm
      if (!is.null(local_arm)) {
        if (is.null(config_sweep[[mrr_container]])) stop(mrr_container, " must be defined in config")
        config_sweep[[mrr_container]][[local_arm]] <- val
      } else {
        # If a scalar val provided without specifying arm, recycle across all treatment arms
        treatment_arms <- setdiff(config_sweep$arms, "control")
        if (length(val) == 1 && is.null(names(val))) {
          config_sweep[[mrr_container]] <- setNames(rep(val, length(treatment_arms)), treatment_arms)
        } else {
          # Expect a named vector; validate names
          if (is.null(names(val))) stop("When providing multiple ", mrr_container, " values, they must be named by treatment arms")
          missing_names <- setdiff(treatment_arms, names(val))
          if (length(missing_names) > 0) stop(mrr_container, " sweep value missing names for: ", paste(missing_names, collapse=","))
          config_sweep[[mrr_container]] <- val[treatment_arms]
        }
      }
    }
    # Environmental ATE sweep: support config-native names and canonical
    if (sweep_param %in% c("env_outcome_ate_single","env_outcome_ate_multi")) {
      vdm_container <- if (!is.null(config_sweep$env_outcome_ate_single) || sweep_param == "env_outcome_ate_single") {
        "env_outcome_ate_single"
      } else if (!is.null(config_sweep$env_outcome_ate_multi) || sweep_param == "env_outcome_ate_multi") {
        "env_outcome_ate_multi"
      } else stop("env_outcome_ate_* must be defined in config")
      local_arm <- if (!is.null(override_arm)) override_arm else sweep_arm
      if (!is.null(local_arm)) {
        if (is.null(config_sweep[[vdm_container]])) stop(vdm_container, " must be defined in config")
        config_sweep[[vdm_container]][[local_arm]] <- val
      } else {
        treatment_arms <- setdiff(config_sweep$arms, "control")
        if (length(val) == 1 && is.null(names(val))) {
          config_sweep[[vdm_container]] <- setNames(rep(val, length(treatment_arms)), treatment_arms)
        } else {
          if (is.null(names(val))) stop("When providing multiple ", vdm_container, " values, they must be named by treatment arms")
          missing_names <- setdiff(treatment_arms, names(val))
          if (length(missing_names) > 0) stop(vdm_container, " sweep value missing names for: ", paste(missing_names, collapse=","))
          config_sweep[[vdm_container]] <- val[treatment_arms]
        }
      }
    }
    if (sweep_param == "n_communities")          config_sweep$n_communities        <- as.integer(val)
    if (sweep_param == "avg_ind_obs_per_comm")   config_sweep$avg_ind_obs_per_comm <- as.integer(val)
    if (sweep_param == "soc_outcome_T") {
        config_sweep$soc_outcome_T        <- as.integer(val)
        config_sweep$soc_outcome_T_months <- seq_len(config_sweep$soc_outcome_T)
    }
    if (sweep_param == "soc_outcome_ICC")      config_sweep$soc_outcome_ICC <- val
    if (sweep_param == "soc_outcome_AR1_rho")      config_sweep$soc_outcome_AR1_rho <- val
    if (sweep_param == "env_outcome_T") {
        config_sweep$env_outcome_T        <- as.integer(val)
        config_sweep$env_outcome_T_months <- seq_len(config_sweep$env_outcome_T)
    }
    if (sweep_param == "env_outcome_ICC")      config_sweep$env_outcome_ICC <- val
  if (sweep_param == "env_outcome_AR1_rho")      config_sweep$env_outcome_AR1_rho <- val
  if (sweep_param == "soc_outcome_AR1_var")  config_sweep$soc_outcome_AR1_var <- val
  if (sweep_param == "env_outcome_AR1_var")  config_sweep$env_outcome_AR1_var <- val
    if (sweep_param == "alloc_ratio")          config_sweep$alloc_ratio <- val  # binary fallback path only

  # Monte Carlo over sims (parallelizable layer: inner)
  if (isTRUE(parallel_inner)) {
    res <- future_map(1:config_sweep$sims, ~ one_run(config_sweep), .options = furrr_options(seed = TRUE))
  } else {
    res <- map(1:config_sweep$sims, ~ one_run(config_sweep))
  }

    get_rate <- function(xvec) mean(xvec < config_sweep$alpha, na.rm = TRUE)
    trt_arms <- setdiff(config_sweep$arms, "control")
    # Per-arm power calculations
    pow_list <- list()
    for (a in trt_arms) {
      p_itt_p_vec <- sapply(res, `[[`, paste0("p_itt_p_", a))
      p_tot_p_vec <- sapply(res, `[[`, paste0("p_tot_p_", a))
      p_itt_v_vec <- sapply(res, `[[`, paste0("p_itt_v_", a))
      p_tot_v_vec <- sapply(res, `[[`, paste0("p_tot_v_", a))
      # Socio-economic outcome (generic count outcome)
      pow_list[[paste0("power_ITT_soc_outcome_", a)]] <- get_rate(p_itt_p_vec)
      pow_list[[paste0("power_TOT_soc_outcome_", a)]] <- get_rate(p_tot_p_vec)
  # Environmental outcome
      pow_list[[paste0("power_ITT_env_outcome_", a)]] <- get_rate(p_itt_v_vec)
      pow_list[[paste0("power_TOT_env_outcome_", a)]] <- get_rate(p_tot_v_vec)
  # Backward compatibility aliases (not used downstream here)
  # pow_list[[paste0("power_ITT_outcome_", a)]] <- pow_list[[paste0("power_ITT_soc_outcome_", a)]]
  # pow_list[[paste0("power_TOT_outcome_", a)]] <- pow_list[[paste0("power_TOT_soc_outcome_", a)]]
    }
    # ---------------- Add per-arm effect size columns to output ---------------- #
    trt_arms <- setdiff(config_sweep$arms, "control")
    # Meeting rate ratios (multiplicative)
    mrr_vals <- get_soc_ate_pct(config_sweep)
    if (length(mrr_vals) == 0 && length(trt_arms) > 0) {
      mrr_vals <- setNames(rep(1, length(trt_arms)), trt_arms)
    }
    # Environmental additive effects
    vdm_vals <- get_env_add(config_sweep)
    if (length(vdm_vals) == 0 && length(trt_arms) > 0) {
      vdm_vals <- setNames(rep(0, length(trt_arms)), trt_arms)
    }

    out_row <- tibble(
      sweep_param    = sweep_param,
      sweep_value    = val,
      experiment_type = config_sweep$experiment_type,
      outcome_dist    = config_sweep$soc_outcome_dist
    )
    # Bind power columns
    for (nm in names(pow_list)) out_row[[nm]] <- pow_list[[nm]]
    if (length(trt_arms) == 1) {
      # duplicate unsuffixed for single-arm convenience
      a <- trt_arms[1]
      out_row$power_ITT_soc_outcome <- out_row[[paste0("power_ITT_soc_outcome_", a)]]
      out_row$power_TOT_soc_outcome <- out_row[[paste0("power_TOT_soc_outcome_", a)]]
      out_row$power_ITT_env_outcome <- out_row[[paste0("power_ITT_env_outcome_", a)]]
      out_row$power_TOT_env_outcome <- out_row[[paste0("power_TOT_env_outcome_", a)]]
  # Backward compatibility (old names)
  out_row$power_ITT_outcome <- out_row$power_ITT_soc_outcome
  out_row$power_TOT_outcome <- out_row$power_TOT_soc_outcome
    }
    # Append effect size columns explicitly (renamed: arm_soc_outcome_ate_pct_*)
    for (a in trt_arms) {
      col_mrr <- paste0("arm_soc_outcome_ate_pct_", a)
  col_vdm_new <- paste0("arm_env_outcome_ate_", a)
      out_row[[col_mrr]]    <- unname(mrr_vals[a])
      out_row[[col_vdm_new]] <- unname(vdm_vals[a])
    }
    if (!is.null(override_arm)) out_row$swept_arm <- override_arm
    out_row
  }

  # Run across sweep with progress bar
  if (sweep_each_arm && !(sweep_param %in% c("soc_outcome_ate_pct_single","soc_outcome_ate_pct_multi","env_outcome_ate_single","env_outcome_ate_multi"))) {
    warning("sweep_each_arm currently only implemented for sweep_param in {'soc_outcome_ate_pct_single','soc_outcome_ate_pct_multi','env_outcome_ate_single','env_outcome_ate_multi'}; ignoring flag.")
    sweep_each_arm <- FALSE
  }

  if (sweep_each_arm) {
    trt_arms_all <- setdiff(config$arms, "control")
    if (length(trt_arms_all) < 1) stop("No treatment arms to sweep.")
    # For each treatment arm, run independent sweep varying only that arm's effect while holding others at baseline
    arm_results_list <- list()
  p <- progressr::progressor(along = trt_arms_all)
    for (a in trt_arms_all) {
      p(message = paste0("Sweeping arm ", a))
      if (isTRUE(parallel_outer)) {
        results_a <- future_map_dfr(
          sweep_values,
          ~ run_for_value(.x, override_arm = a),
          .progress = FALSE,
          .options = furrr_options(seed = TRUE)  # reproducible per-arm parallel sweep
        )
      } else {
        results_a <- map_dfr(
          sweep_values,
          ~ run_for_value(.x, override_arm = a)
        )
      }
      arm_results_list[[a]] <- results_a
    }
    results <- bind_rows(arm_results_list)
  } else {
    if (isTRUE(parallel_outer)) {
      with_progress({
          results <- future_map_dfr(
            sweep_values,
            run_for_value,
            .progress = TRUE,
            .options = furrr_options(seed = TRUE)  # reproducible parallel RNG
          )
      })
    } else {
      results <- map_dfr(sweep_values, run_for_value)
    }
  }
  # Note: env_outcome_ate_* sweep handling occurs inside run_for_value; no additional handling needed here.

  # Reorder: place arm_soc_outcome_ate_pct_* columns first (user request)
  effect_ate_cols <- grep("^arm_soc_outcome_ate_pct_", names(results), value = TRUE)
  preferred_order <- c("arm_soc_outcome_ate_pct_T1", "arm_soc_outcome_ate_pct_T2")
  effect_ate_cols <- unique(c(preferred_order[preferred_order %in% effect_ate_cols], setdiff(effect_ate_cols, preferred_order)))
  if (length(effect_ate_cols) > 0) {
    remaining <- setdiff(names(results), effect_ate_cols)
    results <- results[, c(effect_ate_cols, remaining)]
  }

  cat("Simulation complete. Saving CSVs to tabs/ and figures to figs/...\n")

  # Set output directories (tabs for CSV/tables; figs for images)
  tabs_dir <- "/Users/flavioamalagutti/Documents/Work/GitHub/Research_Projects/Biltong/Simulations/Biltong_power_simulations/output/tabs"
  figs_dir <- "/Users/flavioamalagutti/Documents/Work/GitHub/Research_Projects/Biltong/Simulations/Biltong_power_simulations/output/figs"

  if (!dir.exists(tabs_dir)) dir.create(tabs_dir, recursive = TRUE)
  if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)
  # Save CSV
  csv_path  <- file.path(tabs_dir, glue::glue("{outfile_stem}_{sweep_param}.csv"))
  write_csv(results, csv_path)

  # Figure
  # Reshape for plotting; handle per-arm columns
  long_res <- results |>
    pivot_longer(
      cols = matches("^power_(ITT|TOT)_(soc_outcome|env_outcome)(_.*)?$"),
      names_to = "metric",
      values_to = "power"
    ) |>
    mutate(
      estimator = sub("^power_(ITT|TOT).*", "\\1", metric),
      outcome = case_when(
        grepl("soc_outcome", metric) ~ "soc_outcome",
        grepl("env_outcome", metric) ~ "env_outcome",
        TRUE ~ NA_character_
      ),
      arm = sub("^power_(ITT|TOT)_(soc_outcome|env_outcome)_?", "", metric),
      arm = ifelse(arm %in% c("soc_outcome","env_outcome",""), NA_character_, arm)
    )

  # If independent sweeps per arm, restrict each row's varying power series to the swept arm only for plotting x vs its own effect.
  if (sweep_each_arm && "swept_arm" %in% names(long_res)) {
    # Keep only power rows where the arm matches the arm that was actively swept so curves reflect own effect trajectory
    long_res <- long_res |> filter(is.na(arm) | arm == swept_arm)
  }

  # Determine if multi-arm
  multi_arm <- any(!is.na(long_res$arm))
  long_res <- long_res |>
    mutate(series = if (multi_arm) {
      paste(estimator, arm, outcome, sep = " - ")
    } else {
      paste(estimator, outcome, sep = " - ")
    })

  # ---------------- Tidy long table (one row per arm x estimator x outcome x sweep value) --------------- #
  # Keep only rows with a specific arm
  # Build selection vector safely (avoid non-standard evaluation error if swept_arm absent)
  base_cols <- c("sweep_param","sweep_value","arm","estimator","outcome","power")
  if ("swept_arm" %in% names(long_res)) base_cols <- c(base_cols, "swept_arm")
  power_long <- long_res |> filter(!is.na(arm)) |> select(all_of(base_cols))
  # If independent sweeps, keep only the own-arm trajectories for tidy output (avoid misleading cross-arm rows)
  if (sweep_each_arm && "swept_arm" %in% names(power_long)) {
    power_long <- power_long |> filter(arm == swept_arm)
  }
  # Pivot effects for socio-economic ATE, preserving swept_arm as an identifier if present
  include_swept <- ("swept_arm" %in% names(results))
  id_cols <- c("sweep_value", if (include_swept) "swept_arm")
  effect_mrr_long <- results |>
    select(all_of(c(id_cols)), starts_with("arm_soc_outcome_ate_pct_")) |>
    pivot_longer(cols = -all_of(id_cols), names_to = "effect_var", values_to = "arm_soc_outcome_ate_pct") |>
    mutate(arm = sub("^arm_soc_outcome_ate_pct_", "", effect_var)) |>
    select(-effect_var)
  effect_env_long <- results |>
    select(all_of(c(id_cols)), starts_with("arm_env_outcome_ate_")) |>
    pivot_longer(cols = -all_of(id_cols), names_to = "effect_var", values_to = "arm_env_outcome_ate") |>
    mutate(arm = sub("^arm_env_outcome_ate_", "", effect_var)) |>
    select(-effect_var)
  effects_join <- full_join(effect_mrr_long, effect_env_long, by = c(id_cols, "arm"))
  long_table <- power_long |>
    left_join(effects_join, by = c(id_cols, "arm")) |>
    mutate(
      scenario_soc_outcome_ate_pct = sweep_value,
      is_swept = ifelse(!is.na(arm) & "swept_arm" %in% names(power_long), arm == swept_arm, TRUE)
    ) |>
    arrange(outcome, estimator, arm, sweep_value)

  # Legacy cleanups removed: meeting_rate_ratio no longer supported

  # Column ordering for long_table: sweep_param, outcome, swept_arm, arm, estimator, sweep_value, then the rest.
  # Ensure the renamed effect columns appear early but after the primary identifiers.
  desired_first <- c("sweep_param","outcome","swept_arm","arm","estimator","sweep_value")
  existing_first <- desired_first[desired_first %in% names(long_table)]
  remaining_cols <- setdiff(names(long_table), existing_first)
  long_table <- long_table[, c(existing_first, remaining_cols)]
  # long_table already contains swept_arm where applicable
  # Write tidy table CSV
  long_csv_path  <- file.path(tabs_dir, glue::glue("{outfile_stem}_{sweep_param}_long.csv"))
  write_csv(long_table, long_csv_path)

  # -------------------- Export example of the last simulated dataset -------------------- #
  # Build the config corresponding to the last sweep value to simulate one final dataset
  last_val <- sweep_values[length(sweep_values)]
  config_last <- config
  # Apply the last sweep to config_last (mirror logic from run_for_value)
  if (sweep_param %in% c("soc_outcome_ate_pct_single","soc_outcome_ate_pct_multi")) {
    trt_arms_all <- setdiff(config_last$arms, "control")
    mrr_container <- if (sweep_param == "soc_outcome_ate_pct_single" || !is.null(config_last$soc_outcome_ate_pct_single)) {
      "soc_outcome_ate_pct_single"
    } else if (sweep_param == "soc_outcome_ate_pct_multi" || !is.null(config_last$soc_outcome_ate_pct_multi)) {
      "soc_outcome_ate_pct_multi"
    } else stop("soc_outcome_ate_pct_* must be defined in config")
    if (sweep_each_arm) {
      last_arm <- tail(trt_arms_all, 1)
      if (is.null(config_last[[mrr_container]])) stop(mrr_container, " must be defined in config")
      config_last[[mrr_container]][[last_arm]] <- last_val
    } else if (!is.null(sweep_arm)) {
      if (is.null(config_last[[mrr_container]])) stop(mrr_container, " must be defined in config")
      config_last[[mrr_container]][[sweep_arm]] <- last_val
    } else {
      config_last[[mrr_container]] <- setNames(rep(last_val, length(trt_arms_all)), trt_arms_all)
    }
  }
  if (sweep_param == "n_communities")          config_last$n_communities        <- as.integer(last_val)
  if (sweep_param == "avg_ind_obs_per_comm")   config_last$avg_ind_obs_per_comm <- as.integer(last_val)
  if (sweep_param == "soc_outcome_T") {
    config_last$soc_outcome_T        <- as.integer(last_val)
    config_last$soc_outcome_T_months <- seq_len(config_last$soc_outcome_T)
  }
  if (sweep_param == "soc_outcome_ICC")      config_last$soc_outcome_ICC <- last_val
  if (sweep_param == "soc_outcome_AR1_rho")  config_last$soc_outcome_AR1_rho <- last_val
  if (sweep_param == "env_outcome_T") {
    config_last$env_outcome_T        <- as.integer(last_val)
    config_last$env_outcome_T_months <- seq_len(config_last$env_outcome_T)
  }
  if (sweep_param == "env_outcome_ICC")      config_last$env_outcome_ICC <- last_val
  if (sweep_param == "env_outcome_AR1_rho")  config_last$env_outcome_AR1_rho <- last_val
  if (sweep_param == "soc_outcome_AR1_var")  config_last$soc_outcome_AR1_var <- last_val
  if (sweep_param == "env_outcome_AR1_var")  config_last$env_outcome_AR1_var <- last_val
  if (sweep_param == "alloc_ratio")          config_last$alloc_ratio <- last_val
  if (sweep_param %in% c("env_outcome_ate_single","env_outcome_ate_multi")) {
    trt_arms_all <- setdiff(config_last$arms, "control")
    vdm_container <- if (sweep_param == "env_outcome_ate_single" || !is.null(config_last$env_outcome_ate_single)) {
      "env_outcome_ate_single"
    } else if (sweep_param == "env_outcome_ate_multi" || !is.null(config_last$env_outcome_ate_multi)) {
      "env_outcome_ate_multi"
    } else stop("env_outcome_ate_* must be defined in config")
    if (!is.null(sweep_arm)) {
      if (is.null(config_last[[vdm_container]])) stop(vdm_container, " must be defined in config")
      config_last[[vdm_container]][[sweep_arm]] <- last_val
    } else {
      config_last[[vdm_container]] <- setNames(rep(last_val, length(trt_arms_all)), trt_arms_all)
    }
  }

  # Simulate one dataset with these settings and export
  last_run <- one_run(config_last, capture_data = TRUE)
  last_soc <- last_run$last_df_soc |> mutate(outcome = "soc_outcome")
  last_env <- last_run$last_df_env |> mutate(outcome = "env_outcome")
  last_combined <- bind_rows(last_soc, last_env)
  last_data_csv <- file.path(tabs_dir, glue::glue("{outfile_stem}_{sweep_param}_last_sim_data.csv"))
  write_csv(last_combined, last_data_csv)
  # Optional diagnostics export for the last run
  last_diag_files <- list()
  if (isTRUE(config$print_strata_diag)) {
    if (!is.null(last_run$diag_strata_counts)) {
      arm_counts_csv <- file.path(tabs_dir, glue::glue("{outfile_stem}_{sweep_param}_last_arm_counts.csv"))
      write_csv(last_run$diag_strata_counts, arm_counts_csv)
      last_diag_files$arm_counts_csv <- arm_counts_csv
    }
    if (!is.null(last_run$diag_te_draws)) {
      te_draws_csv <- file.path(tabs_dir, glue::glue("{outfile_stem}_{sweep_param}_last_te_draws.csv"))
      write_csv(last_run$diag_te_draws, te_draws_csv)
      last_diag_files$te_draws_csv <- te_draws_csv
    }
  }

  base_theme <- theme_minimal(base_size = 12) +
    theme(
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.direction = "horizontal",
      legend.title     = element_blank(),
      legend.key.width = unit(1.2, "cm"),
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 10),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(hjust = 0, size = 8)
      )

  # Helper formatters for readable caption text
  percent_str <- function(x, accuracy = 1) scales::percent(x, accuracy = accuracy)
  number_str <- function(x, accuracy = 0.1) scales::number(x, accuracy = accuracy, big.mark = ",")
  format_named <- function(x, formatter, fallback = "not specified") {
    if (is.null(x) || length(x) == 0) return(fallback)
    if (is.null(names(x)) || any(names(x) == "")) names(x) <- paste0("V", seq_along(x))
    paste(paste0(names(x), "=", formatter(x)), collapse = ", ")
  }

  soc_ar1_var   <- if (is.null(config$soc_outcome_AR1_var)) 1 else config$soc_outcome_AR1_var
  env_ar1_var   <- if (is.null(config$env_outcome_AR1_var)) (config$env_outcome_base_sd^2) else config$env_outcome_AR1_var
  strata_range  <- paste(range(config$year_in_program), collapse = " to ")
  se_text <- if (isTRUE(config$cluster_se)) {
    glue::glue("Cluster-robust ({config$hc_type})")
  } else {
    glue::glue("Eicker-White {config$hc_type}")
  }

  caption_text <- glue::glue_collapse(
    c(
      glue::glue("Arms: {paste(config$arms, collapse = ', ')}. Allocation: {format_named(config$alloc_ratios, function(v) percent_str(v, accuracy = 1))}"),
      glue::glue("Design: {config$n_communities} communities; avg individuals={number_str(config$avg_ind_obs_per_comm, accuracy = 1)} (sd={number_str(config$sd_indiv_per_comm, accuracy = 1)})."),
      glue::glue("Observations: Socio-economic T={config$soc_outcome_T} at months [{paste(config$soc_outcome_T_months, collapse = ', ')}]; Environmental T={config$env_outcome_T} at months [{paste(config$env_outcome_T_months, collapse = ', ')}]."),
      glue::glue("Stratification: on year, NGO, and tribe with outcome distributions N(0, var: year={number_str(config$year_in_program_ate_var, accuracy = 0.1)}, NGO={number_str(config$ngo_id_ate_var, accuracy = 0.1)}, tribe={number_str(config$tribe_id_ate_var, accuracy = 0.1)})."),
      glue::glue("Take-up (ITT vs. TOT): {format_named(config$take_up, function(v) percent_str(v, accuracy = 1))}"),
      glue::glue("Soc outcome: dist={config$soc_outcome_dist}, theta={number_str(config$soc_outcome_theta, accuracy = 0.01)}, baseline={number_str(config$soc_outcome_base_mean, accuracy = 1)}, ICC={number_str(config$soc_outcome_ICC, accuracy = 0.001)}, AR1 rho={number_str(config$soc_outcome_AR1_rho, accuracy = 0.01)}, AR1 var={number_str(soc_ar1_var, accuracy = 0.01)}."),
      glue::glue("Env outcome: baseline mean={number_str(config$env_outcome_base_mean, accuracy = 1)}, sd={number_str(config$env_outcome_base_sd, accuracy = 1)}, ICC={number_str(config$env_outcome_ICC, accuracy = 0.001)}, AR1 rho={number_str(config$env_outcome_AR1_rho, accuracy = 0.01)}, AR1 var={number_str(env_ar1_var, accuracy = 0.01)}."),
      glue::glue("Inference: alpha={percent_str(config$alpha, accuracy = 1)}, sims={config$sims}, seed={seed}, SE={se_text}")
    ),
    sep = "\n"
  )

  # Socio-economic outcome plot
  part_df <- long_res |> filter(outcome == "soc_outcome")
  # Dynamic palette for participation (distinct colors per series)
  part_series_all <- sort(unique(part_df$series))
  if (multi_arm) {
    trt_arms <- setdiff(config$arms, "control")
  # Series keys must match what's in data: we used outcome tokens (soc_outcome)
  itt_order <- paste("ITT", trt_arms, "soc_outcome", sep = " - ")
  tot_order <- paste("TOT", trt_arms, "soc_outcome", sep = " - ")
    part_order <- c(itt_order, tot_order)
    # Some safety: keep only those that exist
    part_order <- part_order[part_order %in% part_series_all]
    part_cols_vec <- get_palette(length(part_order), "Dark2")
    part_cols <- setNames(part_cols_vec, part_order)
  part_labels <- c(paste0("ITT: ", trt_arms, " - Socio-economic"),
           paste0("TOT: ", trt_arms, " - Socio-economic"))
    # Align labels to existing order (after filtering existence)
    lab_map <- setNames(part_labels, c(itt_order, tot_order))
    part_labels_final <- lab_map[part_order]
    part_scale <- scale_color_manual(
      values = part_cols,
      breaks = part_order,
      labels = part_labels_final,
      guide = guide_legend(nrow = 1, byrow = TRUE, direction = "horizontal", label.position = "right")
    )
  } else {
  # Single arm: keep previous coloring (two series max)
    part_order <- part_series_all
    part_cols <- setNames(get_palette(length(part_order), "Dark2"), part_order)
    # Friendly labels with colon style
  part_labels_final <- sub("^(ITT) - (Socio-economic)$", "ITT: Socio-economic", part_order)
  part_labels_final <- sub("^(TOT) - (Socio-economic)$", "TOT: Socio-economic", part_labels_final)
  part_labels_final <- sub("^(ITT) - ([^-]+) - Socio-economic$", "ITT: \\2 - Socio-economic", part_labels_final)
  part_labels_final <- sub("^(TOT) - ([^-]+) - Socio-economic$", "TOT: \\2 - Socio-economic", part_labels_final)
    part_scale <- scale_color_manual(
      values = part_cols,
      breaks = part_order,
      labels = part_labels_final,
      guide = guide_legend(direction = "horizontal", label.position = "right")
    )
  }
  # Simplified legend: only two entries (ITT solid, TOT dashed); color & linetype both map to estimator.
  # NOTE: Arms are still distinguished in grouping but share aesthetics; if per-arm distinction is needed later,
  # we can add facetting or minor alpha jitter. Assumption: user prioritizes compact legend over per-arm color.
  part_df <- part_df |> filter(!is.na(arm)) |> mutate(group_id = interaction(estimator, arm, drop = TRUE))
  plt_part <- ggplot(part_df, aes(x = sweep_value, y = power,
                                  group = group_id,
                                  color = estimator, linetype = estimator)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.0, show.legend = FALSE) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#666666") +
    scale_color_manual(values = c(ITT = "#1b9e77", TOT = "#d95f02"), name = NULL) +
    scale_linetype_manual(values = c(ITT = "solid", TOT = "dashed"), name = NULL) +
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1)) +
    scale_x_continuous(breaks = unique(results$sweep_value)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1), expand = c(0, 0)) +
    labs(
      x = sweep_param,
      y = "Power",
      title = glue::glue("Power vs {sweep_param} (Socio-economic){if (multi_arm) ' by Arm' else ''}"),
      subtitle = glue::glue("Design: {unique(results$experiment_type)} | Dist: {unique(results$outcome_dist)}"),
      caption = caption_text
    ) + base_theme

  # Environmental outcome plot
  env_df <- long_res |> filter(outcome == "env_outcome")
  env_series_all <- sort(unique(env_df$series))
  if (multi_arm) {
    trt_arms <- setdiff(config$arms, "control")
  # Series keys must match what's in data: we used outcome tokens (env_outcome)
  itt_order_v <- paste("ITT", trt_arms, "env_outcome", sep = " - ")
  tot_order_v <- paste("TOT", trt_arms, "env_outcome", sep = " - ")
    env_order <- c(itt_order_v, tot_order_v)
    env_order <- env_order[env_order %in% env_series_all]
    env_cols_vec <- get_palette(length(env_order), "Set1")
    env_cols <- setNames(env_cols_vec, env_order)
  env_labels <- c(paste0("ITT: ", trt_arms, " - Environmental"),
          paste0("TOT: ", trt_arms, " - Environmental"))
    lab_map_v <- setNames(env_labels, c(itt_order_v, tot_order_v))
    env_labels_final <- lab_map_v[env_order]
    env_scale <- scale_color_manual(
      values = env_cols,
      breaks = env_order,
      labels = env_labels_final,
      guide = guide_legend(nrow = 1, byrow = TRUE, direction = "horizontal", label.position = "right")
    )
  } else {
    env_order <- env_series_all
    env_cols <- setNames(get_palette(length(env_order), "Set1"), env_order)
    env_labels_final <- sub("^(ITT) - (Environmental)$", "ITT: Environmental", env_order)
    env_labels_final <- sub("^(TOT) - (Environmental)$", "TOT: Environmental", env_labels_final)
    env_labels_final <- sub("^(ITT) - ([^-]+) - Environmental$", "ITT: \\2 - Environmental", env_labels_final)
    env_labels_final <- sub("^(TOT) - ([^-]+) - Environmental$", "TOT: \\2 - Environmental", env_labels_final)
    env_scale <- scale_color_manual(
      values = env_cols,
      breaks = env_order,
      labels = env_labels_final,
      guide = guide_legend(direction = "horizontal", label.position = "right")
    )
  }
  env_df <- env_df |> filter(!is.na(arm)) |> mutate(group_id = interaction(estimator, arm, drop = TRUE))
  plt_env <- ggplot(env_df, aes(x = sweep_value, y = power,
                                group = group_id,
                                color = estimator, linetype = estimator)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.0, show.legend = FALSE) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#666666") +
    scale_color_manual(values = c(ITT = "#1b9e77", TOT = "#d95f02"), name = NULL) +
    scale_linetype_manual(values = c(ITT = "solid", TOT = "dashed"), name = NULL) +
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 1)) +
    scale_x_continuous(breaks = unique(results$sweep_value)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1), expand = c(0, 0)) +
    labs(
      x = sweep_param,
      y = "Power",
      title = glue::glue("Power vs {sweep_param} (Environmental){if (multi_arm) ' by Arm' else ''}"),
      subtitle = glue::glue("Design: {unique(results$experiment_type)}"),
      caption = caption_text
    ) + base_theme

  # File paths for separate plots
  part_png  <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}_soc_outcome.png"))
  part_pdf  <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}_soc_outcome.pdf"))
  env_png   <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}_env_outcome.png"))
  env_pdf   <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}_env_outcome.pdf"))

  ggsave(part_png,  plt_part, width = 8, height = 5, dpi = 300)
  ggsave(part_pdf,  plt_part, width = 8, height = 5)
  ggsave(env_png,   plt_env,  width = 8, height = 5, dpi = 300)
  ggsave(env_pdf,   plt_env,  width = 8, height = 5)
  message("Exported results to:\n", csv_path, "\n", long_csv_path, "\n",
  part_png, "\n", part_pdf, "\n", env_png,  "\n", env_pdf,
    "\n",
    last_data_csv,
    if (length(last_diag_files) > 0) paste0("\n", paste(unlist(last_diag_files), collapse = "\n")) else "")

  list(results = results,
    csv = csv_path,
    last_sim_data_csv = last_data_csv,
    last_diag = last_diag_files,
    # Preferred names
    soc_outcome_png = part_png,
    soc_outcome_pdf = part_pdf,
    env_outcome_png = env_png,
    env_outcome_pdf = env_pdf,
    plot_soc_outcome = plt_part,
    plot_env_outcome = plt_env,
    long_table = long_table,
    long_csv = long_csv_path,
    # Backward compatibility aliases (kept minimal; removed VCU-specific aliases)
    outcome_png = part_png,
    outcome_pdf = part_pdf,
    participation_png = part_png,
    participation_pdf = part_pdf,
    plot_outcome = plt_part,
    plot_participation = plt_part)
}