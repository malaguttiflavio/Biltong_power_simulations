# ------------------------------ Core Simulator ----------------------------- #

# Parallelized using {furrr} for multicore processing
# Uses 14 cores based on MacBook Pro M4 Pro
# Supports multi-arm designs: arms[1] is control; per-arm effects via ate_pct (formerly avg_treatment_effect_pct) and VCU_delta_mean
simulate_power <- function(
  config,
  sweep_param = c("ate_pct","n_communities","avg_indiv_per_comm",
                  "T_outcome","ICC_outcome","AR1_outcome","alloc_ratio"),
  sweep_arm = NULL,
  sweep_values = NULL,
  sweep_each_arm = FALSE,  # when TRUE and sweep_param is arm-level (e.g. ate_pct), run an independent sweep for each treatment arm
  parallel_layer = c("inner","outer","both","none"), # choose which layer to parallelize
  outfile_stem = "power_results",
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

    # ate_pct (formerly avg_treatment_effect_pct): multiplicative (>0); only for treatment arms; default 1
    if (!is.null(cfg$avg_treatment_effect_pct) && is.null(cfg$ate_pct)) {
      cfg$ate_pct <- cfg$avg_treatment_effect_pct
      message("[compat] Mapped deprecated 'avg_treatment_effect_pct' to 'ate_pct'.")
      cfg$avg_treatment_effect_pct <- NULL
    }
    if (is.null(cfg$ate_pct)) {
      cfg$ate_pct <- setNames(rep(1, length(treatment_arms)), treatment_arms)
    } else {
      if (is.null(names(cfg$ate_pct)) && length(cfg$ate_pct) == 1 && length(treatment_arms) == 1) {
        cfg$ate_pct <- setNames(cfg$ate_pct, treatment_arms)
      }
      if (is.null(names(cfg$ate_pct))) stop("ate_pct must be a named vector with treatment arm names (omit 'control')")
      missing_ate <- setdiff(treatment_arms, names(cfg$ate_pct))
      if (length(missing_ate) > 0) {
        cfg$ate_pct <- c(cfg$ate_pct, setNames(rep(1, length(missing_ate)), missing_ate))
      }
      cfg$ate_pct <- cfg$ate_pct[setdiff(names(cfg$ate_pct), "control")]
    }

    # VCU_delta_mean: additive effect; default 0
    if (is.null(cfg$VCU_delta_mean)) {
      cfg$VCU_delta_mean <- setNames(rep(0, length(treatment_arms)), treatment_arms)
    } else {
      if (is.null(names(cfg$VCU_delta_mean)) && length(cfg$VCU_delta_mean) == 1 && length(treatment_arms) == 1) {
        cfg$VCU_delta_mean <- setNames(cfg$VCU_delta_mean, treatment_arms)
      }
      if (is.null(names(cfg$VCU_delta_mean))) stop("VCU_delta_mean must be a named vector with treatment arm names (omit 'control')")
      missing_vdm <- setdiff(treatment_arms, names(cfg$VCU_delta_mean))
      if (length(missing_vdm) > 0) {
        cfg$VCU_delta_mean <- c(cfg$VCU_delta_mean, setNames(rep(0, length(missing_vdm)), missing_vdm))
      }
      cfg$VCU_delta_mean <- cfg$VCU_delta_mean[setdiff(names(cfg$VCU_delta_mean), "control")]
    }

    cfg
  }
  config <- normalize_config(config)
  # Defaults for new inference options
  if (is.null(config$cluster_se))     config$cluster_se     <- TRUE
  if (is.null(config$use_cluster_fe)) config$use_cluster_fe <- FALSE  # If TRUE keeps factor(community_id) in model

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
      n_i <- pmax(1L, round(rnorm(cfg_local$n_communities, cfg_local$avg_indiv_per_comm, cfg_local$sd_indiv_per_comm)))
    } else {
      n_i <- rep(cfg_local$avg_indiv_per_comm, cfg_local$n_communities)
    }
    comm$N_indiv <- n_i

    if (cfg_local$experiment_type == "community_level") {
      # Participation panel: community-time
      df_p <- comm |>
        rowwise() |>
  mutate(times = list(cfg_local$months_outcome)) |>
        unnest(times) |>
        mutate(time = times, .keep = "unused") |>
        ungroup()
      # VCU panel: community-time
      df_v <- comm |>
        rowwise() |>
  mutate(times = list(cfg_local$months_VCU)) |>
        unnest(times) |>
        mutate(time = times, .keep = "unused") |>
        ungroup()
    } else {
      # individual_within_community
      df_p <- comm |>
        rowwise() |>
  mutate(indiv = list(1:N_indiv), times = list(cfg_local$months_outcome)) |>
        unnest(c(indiv, times)) |>
        rename(time = times) |>
        ungroup()
      # VCU remains at community-time
      df_v <- comm |>
        rowwise() |>
  mutate(times = list(cfg_local$months_VCU)) |>
        unnest(times) |>
        mutate(time = times, .keep = "unused") |>
        ungroup()
    }

    list(comm = comm, df_p = df_p, df_v = df_v)
  }

  # One Monte Carlo run returning p-values for ITT and TOT for both outcomes, per treatment arm
  one_run <- function(config_run) {
    smp <- make_sample(config_run)

    # ----------------------------- Assignment ----------------------------- #
    # Prefer multi-arm assignment using arms + alloc_ratios; fallback to binary alloc_ratio
    if (!is.null(config_run$alloc_ratios) && !is.null(config_run$arms)) {
      stopifnot(abs(sum(as.numeric(config_run$alloc_ratios[config_run$arms])) - 1) < 1e-8)
      arm_comm <- tibble(
        community_id = smp$comm$community_id,
        arm = factor(
          sample(config_run$arms,
                 size = nrow(smp$comm),
                 replace = TRUE,
                 prob = as.numeric(config_run$alloc_ratios[config_run$arms])),
          levels = config_run$arms
        )
      )
    } else {
      prob <- if (!is.null(config_run$alloc_ratio)) config_run$alloc_ratio else 0.5
      arm_comm <- tibble(
        community_id = smp$comm$community_id,
        arm = factor(ifelse(rbinom(nrow(smp$comm), 1, prob) == 1, "T1", "control"),
                     levels = c("control","T1"))
      )
      if (is.null(config_run$arms)) config_run$arms <- levels(arm_comm$arm)
      if (is.null(config_run$alloc_ratios)) config_run$alloc_ratios <- c(control = 1 - prob, T1 = prob)
    }

    if (config_run$experiment_type == "community_level") {
      df_p <- smp$df_p |> left_join(arm_comm, by = "community_id")
      df_v <- smp$df_v |> left_join(arm_comm, by = "community_id")
    } else {
      # individual_within_community: inherit the community arm
      df_p <- smp$df_p |> left_join(arm_comm, by = "community_id")
      df_v <- smp$df_v |> left_join(arm_comm, by = "community_id")
    }

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
      df_p <- df_p |> left_join(Dc, by = "community_id")
      df_v <- df_v |> left_join(Dc, by = "community_id")
    } else {
      df_p <- df_p |> mutate(D = rbinom(n(), 1, tu_vec[as.character(arm)]))
      agg <- df_p |>
        group_by(community_id) |>
        summarize(D_bar = mean(D), .groups = "drop")
      df_v <- df_v |> left_join(agg, by = "community_id")
      df_v$D <- df_v$D_bar
    }

    # ------------------- Participation outcome generation ------------------ #
    if (config_run$experiment_type == "community_level") {
      var_e_p <- 1
      var_u_p <- ifelse(config_run$ICC_outcome >= 0.999, 1e6, config_run$ICC_outcome/(1-config_run$ICC_outcome))
      u_comm_p <- rnorm(config_run$n_communities, 0, sqrt(var_u_p))
      df_p <- df_p |>
        left_join(tibble(community_id = 1:config_run$n_communities, u_comm_p), by = "community_id") |>
        group_by(community_id) |>
        arrange(time, .by_group = TRUE) |>
        mutate(e_ar = ar1_series(time, rho = config_run$AR1_outcome, sigma = sqrt(var_e_p))) |>
        ungroup()

      log_mu_base <- log(pmax(config_run$outcome_base_mean, 1e-6)) + df_p$u_comm_p + df_p$e_ar
      # Per-arm multiplicative effects, realized only when D==1
  mrr  <- config_run$ate_pct
      mult <- rep(1, nrow(df_p))
      if (!is.null(mrr)) for (nm in names(mrr)) mult[df_p$arm == nm] <- mult[df_p$arm == nm] * mrr[[nm]]
      log_mu <- log_mu_base + log(mult) * df_p$D
      mu <- pmax(exp(log_mu), 1e-6)
      df_p$Y <- r_count(length(mu), mu, family = config_run$outcome_dist, theta = config_run$theta_outcome)
    } else {
      var_e_p <- 1
      var_u_p <- ifelse(config_run$ICC_outcome >= 0.999, 1e6, config_run$ICC_outcome/(1-config_run$ICC_outcome))
      u_comm_p <- rnorm(config_run$n_communities, 0, sqrt(var_u_p))
      df_p <- df_p |>
        left_join(tibble(community_id = 1:config_run$n_communities, u_comm_p), by = "community_id") |>
        group_by(community_id, indiv) |>
        arrange(time, .by_group = TRUE) |>
        mutate(e_ar = ar1_series(time, rho = config_run$AR1_outcome, sigma = sqrt(var_e_p))) |>
        ungroup()

      log_mu_base <- log(pmax(config_run$outcome_base_mean, 1e-6)) + df_p$u_comm_p + df_p$e_ar
  mrr  <- config_run$ate_pct
      mult <- rep(1, nrow(df_p))
      if (!is.null(mrr)) for (nm in names(mrr)) mult[df_p$arm == nm] <- mult[df_p$arm == nm] * mrr[[nm]]
      log_mu <- log_mu_base + log(mult) * df_p$D
      mu <- pmax(exp(log_mu), 1e-6)
      df_p$Y <- r_count(length(mu), mu, family = config_run$outcome_dist, theta = config_run$theta_outcome)
    }

    # ----------------------------- Participation OLS ----------------------------- #
    df_p$arm <- factor(df_p$arm, levels = config_run$arms)
    form_p <- if (isTRUE(config_run$use_cluster_fe)) {
      Y ~ arm + factor(time) + factor(community_id)
    } else {
      Y ~ arm + factor(time)
    }
    fit_itt_p <- lm(form_p, data = df_p)
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
      df_sub <- df_p[df_p$arm %in% c("control", a), , drop = FALSE]
      # Re-factor for subset to avoid singularities
      df_sub$arm <- droplevels(df_sub$arm)
      df_sub$Z_eval <- as.integer(df_sub$arm == a)
  use_fe_iv_p <- isTRUE(config_run$use_cluster_fe) && !(config_run$experiment_type == "community_level")
  if (isTRUE(config_run$use_cluster_fe) && config_run$experiment_type == "community_level") use_fe_iv_p <- FALSE
      iv_form_sub <- if (use_fe_iv_p) {
        formula(Y ~ D + factor(time) + factor(community_id) | Z_eval + factor(time) + factor(community_id))
      } else {
        Y ~ D + factor(time) | Z_eval + factor(time)
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

    # ---------------------------- VCU generation ---------------------------- #
    var_e_v <- config_run$VCU_base_sd^2
    var_u_v <- ifelse(config_run$ICC_VCU >= 0.999, 1e6, (var_e_v * config_run$ICC_VCU) / (1 - config_run$ICC_VCU))
    u_comm_v <- rnorm(config_run$n_communities, 0, sqrt(pmax(var_u_v, 1e-8)))
    df_v <- df_v |>
      left_join(tibble(community_id = 1:config_run$n_communities, u_comm_v), by = "community_id") |>
      group_by(community_id) |>
      arrange(time, .by_group = TRUE) |>
      mutate(e_ar = ar1_series(time, rho = config_run$AR1_VCU, sigma = sqrt(pmax(var_e_v, 1e-8)))) |>
      ungroup()

    # Ensure D exists from compliance block
  if (!"D" %in% names(df_v)) stop("Column D missing in df_v after compliance")
  df_v$arm <- factor(df_v$arm, levels = config_run$arms)

  vdm <- config_run$VCU_delta_mean
    add <- rep(0, nrow(df_v))
    if (!is.null(vdm)) for (nm in names(vdm)) add[df_v$arm == nm] <- add[df_v$arm == nm] + vdm[[nm]]
  df_v$Y <- config_run$VCU_base_mean + df_v$u_comm_v + df_v$e_ar + add * df_v$D

    # ------------------------------- VCU OLS -------------------------------- #
    form_v <- if (isTRUE(config_run$use_cluster_fe)) {
      Y ~ arm + factor(time) + factor(community_id)
    } else {
      Y ~ arm + factor(time)
    }
    fit_itt_v <- lm(form_v, data = df_v)
    if (isTRUE(config_run$cluster_se)) {
      vc_v <- tryCatch(sandwich::vcovCL(fit_itt_v, cluster = ~ community_id), error = function(e) sandwich::vcovHC(fit_itt_v, type = config_run$hc_type))
      ct_v <- lmtest::coeftest(fit_itt_v, vcov. = vc_v)
      tt_itt_v <- tibble(term = rownames(ct_v), estimate = ct_v[,1], std_error = ct_v[,2], statistic = ct_v[,3], p_value = ct_v[,4])
    } else {
      tt_itt_v  <- robust_test(fit_itt_v, hc_type = config_run$hc_type)
    }
    # ITT per arm (VCU)
    p_itt_v_vec <- setNames(rep(NA_real_, length(trt_arms)), trt_arms)
    for (a in trt_arms) {
      term_a <- paste0("arm", a)
      p_itt_v_vec[a] <- tt_itt_v$p_value[tt_itt_v$term == term_a]
    }
    # TOT per arm (VCU) via separate IV on control + that arm
    p_tot_v_vec <- setNames(rep(NA_real_, length(trt_arms)), trt_arms)
    for (a in trt_arms) {
      df_sub_v <- df_v[df_v$arm %in% c("control", a), , drop = FALSE]
      df_sub_v$arm <- droplevels(df_sub_v$arm)
      df_sub_v$Z_eval <- as.integer(df_sub_v$arm == a)
  use_fe_iv_v <- isTRUE(config_run$use_cluster_fe) && !(config_run$experiment_type == "community_level")
  if (isTRUE(config_run$use_cluster_fe) && config_run$experiment_type == "community_level") use_fe_iv_v <- FALSE
      iv_form_sub_v <- if (use_fe_iv_v) {
        formula(Y ~ D + factor(time) + factor(community_id) | Z_eval + factor(time) + factor(community_id))
      } else {
        Y ~ D + factor(time) | Z_eval + factor(time)
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
    out_list
  }

  # Run sims for a given sweep value; override_arm allows independent sweeps per treatment arm
  run_for_value <- function(val, override_arm = NULL) {
    config_sweep <- config
  if (sweep_param == "ate_pct") {
      local_arm <- if (!is.null(override_arm)) override_arm else sweep_arm
      if (!is.null(local_arm)) {
        if (is.null(config_sweep$ate_pct)) stop("ate_pct must be defined in config")
        config_sweep$ate_pct[[local_arm]] <- val
      } else {
        # If a scalar val provided without specifying arm, recycle across all treatment arms
        treatment_arms <- setdiff(config_sweep$arms, "control")
        if (length(val) == 1 && is.null(names(val))) {
          config_sweep$ate_pct <- setNames(rep(val, length(treatment_arms)), treatment_arms)
        } else {
          # Expect a named vector; validate names
            if (is.null(names(val))) stop("When providing multiple ate_pct values, they must be named by treatment arms")
            missing_names <- setdiff(treatment_arms, names(val))
            if (length(missing_names) > 0) stop("ate_pct sweep value missing names for: ", paste(missing_names, collapse=","))
            config_sweep$ate_pct <- val[treatment_arms]
        }
      }
    }
    if (sweep_param == "n_communities")        config_sweep$n_communities      <- as.integer(val)
    if (sweep_param == "avg_indiv_per_comm")   config_sweep$avg_indiv_per_comm <- as.integer(val)
    if (sweep_param == "T_outcome") {
      config_sweep$T_outcome     <- as.integer(val)
      config_sweep$months_outcome <- seq_len(config_sweep$T_outcome)
    }
    if (sweep_param == "ICC_outcome")          config_sweep$ICC_outcome <- val
    if (sweep_param == "AR1_outcome")          config_sweep$AR1_outcome <- val
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
      pow_list[[paste0("power_ITT_participation_", a)]] <- get_rate(p_itt_p_vec)
      pow_list[[paste0("power_TOT_participation_", a)]] <- get_rate(p_tot_p_vec)
      pow_list[[paste0("power_ITT_VCU_", a)]] <- get_rate(p_itt_v_vec)
      pow_list[[paste0("power_TOT_VCU_", a)]] <- get_rate(p_tot_v_vec)
      # Backward compatibility (single-arm only): also create unsuffixed columns
    }
    # ---------------- Add per-arm effect size columns to output ---------------- #
    trt_arms <- setdiff(config_sweep$arms, "control")
    # Meeting rate ratios (multiplicative)
  mrr_vals <- config_sweep$ate_pct
    if (length(mrr_vals) == 0 && length(trt_arms) > 0) {
      mrr_vals <- setNames(rep(1, length(trt_arms)), trt_arms)
    }
    # VCU additive effects
  vdm_vals <- config_sweep$VCU_delta_mean
    if (length(vdm_vals) == 0 && length(trt_arms) > 0) {
      vdm_vals <- setNames(rep(0, length(trt_arms)), trt_arms)
    }

    out_row <- tibble(
      sweep_param = sweep_param,
      sweep_value = val,
  experiment_type = config_sweep$experiment_type,
  outcome_dist = config_sweep$outcome_dist
    )
    # Bind power columns
    for (nm in names(pow_list)) out_row[[nm]] <- pow_list[[nm]]
    if (length(trt_arms) == 1) {
      # duplicate unsuffixed for single-arm convenience
      a <- trt_arms[1]
      out_row$power_ITT_participation <- out_row[[paste0("power_ITT_participation_", a)]]
      out_row$power_TOT_participation <- out_row[[paste0("power_TOT_participation_", a)]]
      out_row$power_ITT_VCU <- out_row[[paste0("power_ITT_VCU_", a)]]
      out_row$power_TOT_VCU <- out_row[[paste0("power_TOT_VCU_", a)]]
    }
    # Append effect size columns explicitly (renamed: arm_ate_pct_*)
    for (a in trt_arms) {
      col_mrr <- paste0("arm_ate_pct_", a)
      col_vdm <- paste0("arm_VCU_delta_mean_", a)
      out_row[[col_mrr]] <- unname(mrr_vals[a])
      out_row[[col_vdm]] <- unname(vdm_vals[a])
    }
    if (!is.null(override_arm)) out_row$swept_arm <- override_arm
    out_row
  }

  # Run across sweep with progress bar
  if (sweep_each_arm && sweep_param != "ate_pct") {
    warning("sweep_each_arm currently only implemented for sweep_param = 'ate_pct'; ignoring flag.")
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

  # Reorder: place arm_ate_pct_* columns first (user request)
  effect_ate_cols <- grep("^arm_ate_pct_", names(results), value = TRUE)
  preferred_order <- c("arm_ate_pct_T1", "arm_ate_pct_T2")
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
  alloc_caption <- if (!is.null(config$alloc_ratios)) {
    paste0("alloc={", paste(names(config$alloc_ratios), config$alloc_ratios, sep=":", collapse=", "), "}")
  } else {
    paste0("alloc_ratio=", config$alloc_ratio)
  }
  takeup_caption <- if (!is.null(config$take_up)) {
    paste0("take_up={", paste(names(config$take_up), config$take_up, sep=":", collapse=", "), "}")
  } else {
    "take_up not set"
  }

  # Reshape for plotting; handle per-arm columns
  long_res <- results |>
    pivot_longer(
      cols = matches("^power_(ITT|TOT)_(participation|VCU)(_.*)?$"),
      names_to = "metric",
      values_to = "power"
    ) |>
    mutate(
      estimator = sub("^power_(ITT|TOT).*", "\\1", metric),
      outcome = ifelse(grepl("participation", metric), "Participation", "VCU"),
      arm = sub("^power_(ITT|TOT)_(participation|VCU)_?", "", metric),
      arm = ifelse(arm %in% c("participation","VCU",""), NA_character_, arm)
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
  # Keep only rows with a specific arm (avoid duplicated unsuffixed single-arm legacy columns)
  # Build selection vector safely (avoid non-standard evaluation error if swept_arm absent)
  base_cols <- c("sweep_param","sweep_value","arm","estimator","outcome","power")
  if ("swept_arm" %in% names(long_res)) base_cols <- c(base_cols, "swept_arm")
  power_long <- long_res |> filter(!is.na(arm)) |> select(all_of(base_cols))
  # If independent sweeps, keep only the own-arm trajectories for tidy output (avoid misleading cross-arm rows)
  if (sweep_each_arm && "swept_arm" %in% names(power_long)) {
    power_long <- power_long |> filter(arm == swept_arm)
  }
  # Pivot effects for ate_pct, preserving swept_arm as an identifier if present
  include_swept <- ("swept_arm" %in% names(results))
  id_cols <- c("sweep_value", if (include_swept) "swept_arm")
  effect_mrr_long <- results |>
    select(all_of(c(id_cols)), starts_with("arm_ate_pct_")) |>
    pivot_longer(cols = -all_of(id_cols), names_to = "effect_var", values_to = "arm_ate_pct") |>
    mutate(arm = sub("^arm_ate_pct_", "", effect_var)) |>
    select(-effect_var)
  effect_vcu_long <- results |>
    select(all_of(c(id_cols)), starts_with("arm_VCU_delta_mean_")) |>
    pivot_longer(cols = -all_of(id_cols), names_to = "effect_var", values_to = "arm_VCU_delta_mean") |>
    mutate(arm = sub("^arm_VCU_delta_mean_", "", effect_var)) |>
    select(-effect_var)
  effects_join <- full_join(effect_mrr_long, effect_vcu_long, by = c(id_cols, "arm"))
  long_table <- power_long |>
    left_join(effects_join, by = c(id_cols, "arm")) |>
    mutate(
      scenario_ate_pct = sweep_value,
      # arm_ate_pct already present from effects_join
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

  base_theme <- theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.key.width = unit(1.2, "cm"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0, size = 8)
    )

  caption_text <- glue::glue(
    "{alloc_caption}, n_communities={config$n_communities}, avg_indiv_per_comm={config$avg_indiv_per_comm},T_outcome={config$T_outcome}, \n",
    "year strata range={paste(range(config$year_in_program), collapse='-')}, NGOs up to 7, tribes up to 5, \n",
    "alpha={config$alpha}, sims={config$sims}, seed={seed}, arms={paste(config$arms, collapse=',')}, {takeup_caption}, ICC_outcome={config$ICC_outcome};\n",
    "outcome_dist={config$outcome_dist}, theta={config$theta_outcome}, base_mean_outcome={config$outcome_base_mean};\n",
    "VCU_base_mean={config$VCU_base_mean}, VCU_base_sd={config$VCU_base_sd}, T_VCU={config$T_VCU}, ICC_VCU={config$ICC_VCU}, AR1_outcome={config$AR1_outcome}, AR1_VCU={config$AR1_VCU}, SE=Eicker-White {config$hc_type}"
  )

  # Participation plot
  part_df <- long_res |> filter(outcome == "Participation")
  # Dynamic palette for participation (distinct colors per series)
  part_series_all <- sort(unique(part_df$series))
  if (multi_arm) {
    trt_arms <- setdiff(config$arms, "control")
    itt_order <- paste("ITT", trt_arms, "Participation", sep = " - ")
    tot_order <- paste("TOT", trt_arms, "Participation", sep = " - ")
    part_order <- c(itt_order, tot_order)
    # Some safety: keep only those that exist
    part_order <- part_order[part_order %in% part_series_all]
    part_cols_vec <- get_palette(length(part_order), "Dark2")
    part_cols <- setNames(part_cols_vec, part_order)
    part_labels <- c(paste0("ITT: ", trt_arms, " - Participation"),
                     paste0("TOT: ", trt_arms, " - Participation"))
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
    part_labels_final <- sub("^(ITT) - (Participation)$", "ITT: Participation", part_order)
    part_labels_final <- sub("^(TOT) - (Participation)$", "TOT: Participation", part_labels_final)
    part_labels_final <- sub("^(ITT) - ([^-]+) - Participation$", "ITT: \2 - Participation", part_labels_final)
    part_labels_final <- sub("^(TOT) - ([^-]+) - Participation$", "TOT: \2 - Participation", part_labels_final)
    part_scale <- scale_color_manual(
      values = part_cols,
      breaks = part_order,
      labels = part_labels_final,
      guide = guide_legend(direction = "horizontal", label.position = "right")
    )
  }
  plt_part <- ggplot(part_df, aes(x = sweep_value, y = power, group = series, color = series)) +
    geom_line(size = 0.7) +
    geom_point(size = 2.2) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#666666") +
    part_scale +
    scale_x_continuous(breaks = unique(results$sweep_value)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1), expand = c(0, 0)) +
    labs(
      x = sweep_param,
      y = "Power",
  title = glue::glue("Power vs {sweep_param} (Participation){if (multi_arm) ' by Arm' else ''}"),
  subtitle = glue::glue("Design: {unique(results$experiment_type)} | Dist: {unique(results$outcome_dist)}"),
      caption = caption_text
    ) + base_theme

  # VCU plot
  vcu_df <- long_res |> filter(outcome == "VCU")
  vcu_series_all <- sort(unique(vcu_df$series))
  if (multi_arm) {
    trt_arms <- setdiff(config$arms, "control")
    itt_order_v <- paste("ITT", trt_arms, "VCU", sep = " - ")
    tot_order_v <- paste("TOT", trt_arms, "VCU", sep = " - ")
    vcu_order <- c(itt_order_v, tot_order_v)
    vcu_order <- vcu_order[vcu_order %in% vcu_series_all]
    vcu_cols_vec <- get_palette(length(vcu_order), "Set1")
    vcu_cols <- setNames(vcu_cols_vec, vcu_order)
    vcu_labels <- c(paste0("ITT: ", trt_arms, " - VCU"),
                    paste0("TOT: ", trt_arms, " - VCU"))
    lab_map_v <- setNames(vcu_labels, c(itt_order_v, tot_order_v))
    vcu_labels_final <- lab_map_v[vcu_order]
    vcu_scale <- scale_color_manual(
      values = vcu_cols,
      breaks = vcu_order,
      labels = vcu_labels_final,
      guide = guide_legend(nrow = 1, byrow = TRUE, direction = "horizontal", label.position = "right")
    )
  } else {
    vcu_order <- vcu_series_all
    vcu_cols <- setNames(get_palette(length(vcu_order), "Set1"), vcu_order)
    vcu_labels_final <- sub("^(ITT) - (VCU)$", "ITT: VCU", vcu_order)
    vcu_labels_final <- sub("^(TOT) - (VCU)$", "TOT: VCU", vcu_labels_final)
    vcu_labels_final <- sub("^(ITT) - ([^-]+) - VCU$", "ITT: \2 - VCU", vcu_labels_final)
    vcu_labels_final <- sub("^(TOT) - ([^-]+) - VCU$", "TOT: \2 - VCU", vcu_labels_final)
    vcu_scale <- scale_color_manual(
      values = vcu_cols,
      breaks = vcu_order,
      labels = vcu_labels_final,
      guide = guide_legend(direction = "horizontal", label.position = "right")
    )
  }
  plt_vcu <- ggplot(vcu_df, aes(x = sweep_value, y = power, group = series, color = series)) +
    geom_line(size = 0.7) +
    geom_point(size = 2.2) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#666666") +
    vcu_scale +
    scale_x_continuous(breaks = unique(results$sweep_value)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1), expand = c(0, 0)) +
    labs(
      x = sweep_param,
      y = "Power",
  title = glue::glue("Power vs {sweep_param} (VCU){if (multi_arm) ' by Arm' else ''}"),
  subtitle = glue::glue("Design: {unique(results$experiment_type)}"),
      caption = caption_text
    ) + base_theme

  # File paths for separate plots
  part_png  <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}_participation.png"))
  part_pdf  <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}_participation.pdf"))
  vcu_png   <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}_VCU.png"))
  vcu_pdf   <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}_VCU.pdf"))

  ggsave(part_png,  plt_part, width = 8, height = 5, dpi = 300)
  ggsave(part_pdf,  plt_part, width = 8, height = 5)
  ggsave(vcu_png,   plt_vcu,  width = 8, height = 5, dpi = 300)
  ggsave(vcu_pdf,   plt_vcu,  width = 8, height = 5)
  message("Exported results to:\n", csv_path, "\n", long_csv_path, "\n",
    part_png, "\n", part_pdf, "\n", vcu_png,  "\n", vcu_pdf)

  list(results = results,
       csv = csv_path,
       participation_png = part_png,
       participation_pdf = part_pdf,
       VCU_png = vcu_png,
       VCU_pdf = vcu_pdf,
       plot_participation = plt_part,
    plot_VCU = plt_vcu,
    long_table = long_table,
    long_csv = long_csv_path)
}