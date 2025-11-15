simulate_dynamic_power <- function(
  config,
  # Accept historical name "effect_p2" for compatibility
  sweep_param = c("effect_last_period","effect_p2"),
  sweep_values,
  outfile_stem = "dynamic_power",
  parallel_layer = c("inner","outer","both","none"),
  seed = 1,
  output_dir_base = NULL
) {
  set.seed(seed)
  sweep_param <- match.arg(sweep_param)
  parallel_layer <- match.arg(parallel_layer)
  parallel_inner <- parallel_layer %in% c("inner","both")
  parallel_outer <- parallel_layer %in% c("outer","both")

  # ---------- Normalize config and set defaults ---------- #
  normalize <- function(cfg) {
    # Required basics
    if (is.null(cfg$n_communities)) stop("n_communities required")
    if (is.null(cfg$avg_ind_obs_per_comm)) stop("avg_ind_obs_per_comm required")
    if (is.null(cfg$sims)) cfg$sims <- 500
    if (is.null(cfg$alpha)) cfg$alpha <- 0.05
    if (is.null(cfg$arms)) cfg$arms <- c("control","T1","T2")
    stopifnot(cfg$arms[1] == "control")

  # Periods and adaptation controls
    if (is.null(cfg$T_periods)) cfg$T_periods <- 2L
  if (is.null(cfg$fixed_control)) cfg$fixed_control <- TRUE
  # Randomization unit: "individual" (default, legacy) or "community"
  if (is.null(cfg$randomize_at)) cfg$randomize_at <- "individual"
  cfg$randomize_at <- match.arg(cfg$randomize_at, c("individual","community"))
    if (is.null(cfg$control_anchor_share)) cfg$control_anchor_share <- 0.3   # share of units locked into never-treated if fixed_control=TRUE
    if (is.null(cfg$min_arm_prob)) cfg$min_arm_prob <- 0.05                  # lower bound per arm each period
    if (is.null(cfg$min_control_prob)) cfg$min_control_prob <- 0.10          # lower bound for control if not fixed
    if (is.null(cfg$adapt_rule)) cfg$adapt_rule <- "softmax"                 # "softmax" or "epsilon_greedy"
    if (is.null(cfg$temperature)) cfg$temperature <- 1.0                     # softmax temperature
    if (is.null(cfg$epsilon)) cfg$epsilon <- 0.10                             # epsilon for epsilon-greedy
    if (is.null(cfg$adapt_scope)) cfg$adapt_scope <- "global"                # "global" or "by_community"
    if (is.null(cfg$use_ipw)) cfg$use_ipw <- FALSE                           # store/use inverse prob. weights in final-period regression

    # Period 1 randomization ratios
    if (is.null(cfg$alloc_ratios_p1)) cfg$alloc_ratios_p1 <- c(control = 1/3, T1 = 1/3, T2 = 1/3)

    # Eligibility option (optional Caria-style targeting by top performers each period)
    if (is.null(cfg$use_eligibility)) cfg$use_eligibility <- FALSE
    if (is.null(cfg$target_top_share)) cfg$target_top_share <- 0.2
    if (is.null(cfg$alloc_ratios_eligible)) cfg$alloc_ratios_eligible <- c(control = 0.10, T1 = 0.45, T2 = 0.45)
    if (is.null(cfg$alloc_ratios_ineligible)) cfg$alloc_ratios_ineligible <- c(control = 1.00, T1 = 0.00, T2 = 0.00)

    # Take-up (compliance)
    if (is.null(cfg$take_up_by_period)) {
      # default: full take-up on treated arms, zero on control for all periods
      K <- length(cfg$arms); names_default <- cfg$arms
      tpl <- setNames(as.list(rep(NA, cfg$T_periods)), paste0("t", seq_len(cfg$T_periods)))
      for (tt in seq_len(cfg$T_periods)) {
        v <- rep(1, K); names(v) <- names_default; v["control"] <- 0
        tpl[[tt]] <- v
      }
      cfg$take_up_by_period <- tpl
    }

    # Outcome DGP parameters
    if (is.null(cfg$y_base_mean)) cfg$y_base_mean <- 10
    if (is.null(cfg$y_sd)) cfg$y_sd <- 5
    if (is.null(cfg$ability_sd)) cfg$ability_sd <- 3
    if (is.null(cfg$icc_y)) cfg$icc_y <- 0.05
    if (cfg$icc_y < 0 || cfg$icc_y >= 1) stop("icc_y must be in [0,1)")
    # Community characteristic (time-invariant) and interaction with treatment
    if (is.null(cfg$community_char_sd)) cfg$community_char_sd <- 1.0
    if (is.null(cfg$y_char_coef)) cfg$y_char_coef <- 0.0
    if (is.null(cfg$trt_char_interaction)) {
      cfg$trt_char_interaction <- setNames(rep(0, length(setdiff(cfg$arms, "control"))), setdiff(cfg$arms, "control"))
    } else {
      # ensure names and fill missing with 0
      v <- cfg$trt_char_interaction
      nm <- setdiff(cfg$arms, "control"); tmp <- rep(0, length(nm)); names(tmp) <- nm
      tmp[names(v)] <- v
      cfg$trt_char_interaction <- tmp
    }

    # Effects: allow per-period vectors; if not supplied, fall back to effect_p1/effect_p2 and recycle
    make_effects <- function() {
      eff_list <- setNames(vector("list", cfg$T_periods), paste0("t", seq_len(cfg$T_periods)))
      # back-compat with effect_p1/effect_p2
      eff_p1 <- if (!is.null(cfg$effect_p1)) cfg$effect_p1 else setNames(rep(0, length(cfg$arms)-1), setdiff(cfg$arms, "control"))
      eff_p2 <- if (!is.null(cfg$effect_p2)) cfg$effect_p2 else eff_p1
      for (tt in seq_len(cfg$T_periods)) {
        if (!is.null(cfg$effects_by_period) && !is.null(cfg$effects_by_period[[tt]])) {
          eff_list[[tt]] <- cfg$effects_by_period[[tt]]
        } else if (tt == 1) {
          eff_list[[tt]] <- eff_p1
        } else if (tt == 2) {
          eff_list[[tt]] <- eff_p2
        } else {
          eff_list[[tt]] <- eff_p2
        }
      }
      eff_list
    }
    cfg$effects_by_period <- make_effects()

    # Normalize probability vectors
    norm <- function(x, nm) {
      if (is.null(names(x))) names(x) <- cfg$arms
      x <- x[cfg$arms]
      if (any(is.na(x))) stop(nm, " must include all arms")
      x / sum(x)
    }
    cfg$alloc_ratios_p1 <- norm(cfg$alloc_ratios_p1, "alloc_ratios_p1")
    cfg$alloc_ratios_eligible <- norm(cfg$alloc_ratios_eligible, "alloc_ratios_eligible")
    cfg$alloc_ratios_ineligible <- norm(cfg$alloc_ratios_ineligible, "alloc_ratios_ineligible")

    # Derived
    cfg$trt_arms <- setdiff(cfg$arms, "control")
    cfg
  }
  config <- normalize(config)

  # ---------- Helper: build assignment probabilities from past outcomes ---------- #
  build_probs_from_history <- function(df_hist, cfg) {
    # df_hist contains last-period outcomes and assignments
    # Return named numeric vector of length |arms| with probabilities for next period
    means <- sapply(cfg$arms, function(a) {
      sel <- df_hist$arm == a & df_hist$D == 1
      if (!any(sel)) return(NA_real_)
      mean(df_hist$Y[sel], na.rm = TRUE)
    })
    # Replace NA with grand mean to avoid collapse
    grand <- mean(df_hist$Y, na.rm = TRUE)
    means[is.na(means)] <- grand

    if (cfg$adapt_rule == "softmax") {
      z <- (means - max(means)) / max(1e-8, cfg$temperature)
      probs <- exp(z)
      probs <- pmax(probs / sum(probs), cfg$min_arm_prob)
      probs <- probs / sum(probs)
    } else if (cfg$adapt_rule == "epsilon_greedy") {
      best <- which.max(means)
      probs <- rep(cfg$epsilon / length(cfg$arms), length(cfg$arms))
      probs[best] <- 1 - cfg$epsilon + cfg$epsilon / length(cfg$arms)
    } else {
      stop("Unknown adapt_rule: ", cfg$adapt_rule)
    }
    # Enforce minimum control probability if control not fixed
    if (!config$fixed_control) {
      probs["control"] <- max(probs["control"], cfg$min_control_prob)
      probs <- probs / sum(probs)
    }
    probs
  }

  # ---------- One Monte Carlo run ---------- #
  one_run <- function(cfg) {
    nC <- cfg$n_communities
    N_i <- rep(cfg$avg_ind_obs_per_comm, nC)
    comm <- tibble(community_id = 1:nC, N = N_i)
    df <- comm |>
      rowwise() |>
      mutate(indiv = list(seq_len(N))) |>
      unnest(indiv) |>
      ungroup()

    # Random effects: community intercept
    var_u <- if (cfg$icc_y == 0) 0 else cfg$icc_y/(1-cfg$icc_y)
    u_c <- rnorm(nC, 0, sqrt(var_u))
    df <- df |> left_join(tibble(community_id = 1:nC, u_c), by = "community_id")

  # Community-level, time-invariant characteristic
  c_char <- rnorm(nC, 0, cfg$community_char_sd)
  df <- df |> left_join(tibble(community_id = 1:nC, c_char), by = "community_id")

    # Individual ability
    df$ability <- rnorm(nrow(df), 0, cfg$ability_sd)

    # Fixed control anchor (optional; applies only to individual-level randomization)
    df$control_anchor <- 0L
    if (cfg$fixed_control && cfg$randomize_at == "individual") {
      # Lock a fraction of individuals into never-treated status
      df$control_anchor <- as.integer(runif(nrow(df)) < cfg$control_anchor_share)
      df$control_anchor[is.na(df$control_anchor)] <- 0L
    }

    # Containers to store per-period variables
    Tt <- cfg$T_periods
    for (tt in seq_len(Tt)) {
      df[[paste0("arm_p", tt)]] <- NA_character_
      df[[paste0("D", tt)]] <- 0L
      df[[paste0("Y", tt)]] <- NA_real_
      df[[paste0("p_assign_p", tt)]] <- NA_real_
    }

    # -------- Period 1 assignment (unconstrained) -------- #
    if (cfg$randomize_at == "community") {
      # Draw one arm per community; all individuals inherit it
      probs <- cfg$alloc_ratios_p1[cfg$arms]
      comm_draws <- tibble(community_id = 1:nC) |>
        rowwise() |>
        mutate(arm = sample(cfg$arms, size = 1, replace = TRUE, prob = as.numeric(probs))) |>
        ungroup()
      df <- df |>
        left_join(comm_draws, by = "community_id") |>
        mutate(arm_p1 = factor(arm, levels = cfg$arms)) |>
        select(-arm)
      # p_assign is the probability of the chosen arm under the per-community draw
      df$p_assign_p1 <- as.numeric(cfg$alloc_ratios_p1[as.character(df$arm_p1)])
    } else {
      p1_probs <- as.numeric(cfg$alloc_ratios_p1[cfg$arms])
      df$arm_p1 <- factor(sample(cfg$arms, size = nrow(df), replace = TRUE, prob = p1_probs), levels = cfg$arms)
      df$p_assign_p1 <- p1_probs[as.integer(df$arm_p1)]
      # Enforce fixed control anchors only for individual-level randomization
      if (cfg$fixed_control) {
        df$arm_p1[df$control_anchor == 1] <- factor("control", levels = cfg$arms)
        df$p_assign_p1[df$control_anchor == 1] <- cfg$alloc_ratios_p1["control"]
      }
    }
    # Compliance -> D1
    take1 <- cfg$take_up_by_period[[1]]; names(take1) <- cfg$arms
    df$D1 <- 0L
    for (a in cfg$trt_arms) df$D1[df$arm_p1 == a] <- rbinom(sum(df$arm_p1 == a), 1, take1[[a]])
    # Outcomes
  eff1 <- setNames(rep(0, length(cfg$arms)), cfg$arms); eff1[names(cfg$effects_by_period[[1]])] <- cfg$effects_by_period[[1]]
  # Treatment-char interaction applies only to treated arms
  gamma <- setNames(rep(0, length(cfg$arms)), cfg$arms);
  gamma[names(cfg$trt_char_interaction)] <- cfg$trt_char_interaction
  trt_effect1 <- eff1[as.character(df$arm_p1)] + gamma[as.character(df$arm_p1)] * df$c_char
  mu1 <- cfg$y_base_mean + df$ability + df$u_c + cfg$y_char_coef * df$c_char + df$D1 * trt_effect1
    df$Y1 <- rnorm(nrow(df), mean = mu1, sd = cfg$y_sd)

    # -------- Subsequent periods t = 2..T -------- #
    if (Tt >= 2) {
      for (tt in 2:Tt) {
        # Build probabilities for this period
        if (cfg$randomize_at == "community") {
          # Community-level bandit: pick one arm per community based on period t-1 performance
          # Global prior means by arm from last period
          arm_prev <- df[[paste0("arm_p", tt-1)]]
          y_prev <- df[[paste0("Y", tt-1)]]
          global_means <- tapply(y_prev, arm_prev, function(z) mean(z, na.rm = TRUE))
          # Fill missing arms with grand mean
          grand <- mean(y_prev, na.rm = TRUE)
          global_means <- setNames(rep(grand, length(cfg$arms)), cfg$arms)
          gm_obs <- tapply(y_prev, arm_prev, function(z) mean(z, na.rm = TRUE))
          global_means[names(gm_obs)] <- gm_obs

          # Function to produce probs from Q values
          probs_from_Q <- function(Q, cfg) {
            if (cfg$adapt_rule == "softmax") {
              z <- (Q - max(Q, na.rm = TRUE)) / max(1e-8, cfg$temperature)
              p <- exp(z)
              p <- pmax(p / sum(p), cfg$min_arm_prob)
              if (!cfg$fixed_control) p["control"] <- max(p["control"], cfg$min_control_prob)
              p / sum(p)
            } else if (cfg$adapt_rule == "epsilon_greedy") {
              best <- which.max(Q)
              p <- rep(cfg$epsilon / length(Q), length(Q)); names(p) <- names(Q)
              p[best] <- 1 - cfg$epsilon + cfg$epsilon / length(Q)
              if (!cfg$fixed_control) p["control"] <- max(p["control"], cfg$min_control_prob)
              p / sum(p)
            } else stop("Unknown adapt_rule: ", cfg$adapt_rule)
          }

          # Choose per-community arm
          df[[paste0("arm_p", tt)]] <- NA_character_
          df[[paste0("p_assign_p", tt)]] <- NA_real_
          for (cc in 1:nC) {
            idx <- which(df$community_id == cc)
            if (length(idx) == 0) next
            prev_arm <- unique(as.character(df[[paste0("arm_p", tt-1)]][idx]))
            # Community reward: mean prev outcome in this community
            r_c <- mean(df[[paste0("Y", tt-1)]][idx], na.rm = TRUE)
            Q <- global_means
            if (length(prev_arm) == 1 && prev_arm %in% names(Q)) Q[prev_arm] <- r_c
            p_vec <- probs_from_Q(Q, cfg)
            chosen <- sample(cfg$arms, size = 1, replace = TRUE, prob = as.numeric(p_vec[cfg$arms]))
            df[[paste0("arm_p", tt)]][idx] <- chosen
            df[[paste0("p_assign_p", tt)]][idx] <- as.numeric(p_vec[chosen])
          }
          df[[paste0("arm_p", tt)]] <- factor(df[[paste0("arm_p", tt)]], levels = cfg$arms)

        } else if (cfg$use_eligibility) {
          # Caria-style eligibility by top performers within community based on Y_{t-1}
          q <- cfg$target_top_share
          df <- df |>
            group_by(community_id) |>
            mutate(cut_prev = quantile(.data[[paste0("Y", tt-1)]], probs = 1 - q, type = 7),
                   eligible = as.integer(.data[[paste0("Y", tt-1)]] >= cut_prev)) |>
            ungroup() |>
            select(-cut_prev)
          # Build named, aligned probability maps for eligible/ineligible
          prob_map_e <- cfg$alloc_ratios_eligible[cfg$arms]; names(prob_map_e) <- cfg$arms
          prob_map_i <- cfg$alloc_ratios_ineligible[cfg$arms]; names(prob_map_i) <- cfg$arms
          # Assign
          draw_arm <- function(n, prob_map) {
            probs <- as.numeric(prob_map[cfg$arms])
            factor(sample(cfg$arms, size = n, replace = TRUE, prob = probs), levels = cfg$arms)
          }
          arm_vec <- character(nrow(df))
          idx_e <- which(df$eligible == 1)
          idx_i <- which(df$eligible == 0)
          if (length(idx_e) > 0) arm_vec[idx_e] <- as.character(draw_arm(length(idx_e), prob_map_e))
          if (length(idx_i) > 0) arm_vec[idx_i] <- as.character(draw_arm(length(idx_i), prob_map_i))
          # Enforce fixed control anchors
          if (cfg$fixed_control) {
            arm_vec[df$control_anchor == 1] <- "control"
          }
          df[[paste0("arm_p", tt)]] <- factor(arm_vec, levels = cfg$arms)
          # Save assignment probabilities used
          df[[paste0("p_assign_p", tt)]] <- ifelse(df$eligible == 1,
            prob_map_e[as.character(df[[paste0("arm_p", tt)]])],
            prob_map_i[as.character(df[[paste0("arm_p", tt)]])]
          )
        } else {
          # Bandit-style global adaptation using last-period outcomes (softmax/epsilon-greedy)
          hist <- tibble(arm = df[[paste0("arm_p", tt-1)]],
                         D = df[[paste0("D", tt-1)]],
                         Y = df[[paste0("Y", tt-1)]])
          probs_vec <- build_probs_from_history(hist, cfg)
          # Assign with these probs
          arm_vec <- as.character(sample(cfg$arms, size = nrow(df), replace = TRUE, prob = probs_vec))
          # Enforce fixed control anchors
          if (cfg$fixed_control) {
            arm_vec[df$control_anchor == 1] <- "control"
          }
          df[[paste0("arm_p", tt)]] <- factor(arm_vec, levels = cfg$arms)
          df[[paste0("p_assign_p", tt)]] <- probs_vec[as.integer(df[[paste0("arm_p", tt)]])]
        }

        # Compliance at period t
        take_t <- cfg$take_up_by_period[[tt]]; names(take_t) <- cfg$arms
        Dname <- paste0("D", tt); ArmName <- paste0("arm_p", tt)
        df[[Dname]] <- 0L
        for (a in cfg$trt_arms) {
          idx <- which(df[[ArmName]] == a)
          if (length(idx) > 0) df[[Dname]][idx] <- rbinom(length(idx), 1, take_t[[a]])
        }

        # Outcomes at period t
        eff_t <- setNames(rep(0, length(cfg$arms)), cfg$arms); eff_t[names(cfg$effects_by_period[[tt]])] <- cfg$effects_by_period[[tt]]
        gamma <- setNames(rep(0, length(cfg$arms)), cfg$arms); gamma[names(cfg$trt_char_interaction)] <- cfg$trt_char_interaction
        trt_effect_t <- eff_t[as.character(df[[ArmName]])] + gamma[as.character(df[[ArmName]])] * df$c_char
        mu_t <- cfg$y_base_mean + df$ability + df$u_c + cfg$y_char_coef * df$c_char + df[[Dname]] * trt_effect_t
        df[[paste0("Y", tt)]] <- rnorm(nrow(df), mean = mu_t, sd = cfg$y_sd)
      }
    }

    # -------- Inference on last period (ITT) -------- #
    last <- cfg$T_periods
    df[[paste0("arm_p", last)]] <- factor(df[[paste0("arm_p", last)]], levels = cfg$arms)
    # Base ANCOVA with cluster-robust SEs
    form_str <- paste0("Y", last, " ~ ", "arm_p", last, " + ", "Y", last-1)
    if (last == 1) form_str <- paste0("Y1 ~ arm_p1")
    base_form <- as.formula(form_str)
    fit <- lm(base_form, data = df)

    if (cfg$use_ipw && last >= 2) {
      # simple IPW via frequency weights from assignment probabilities of the last period
      w <- 1 / pmax(1e-6, df[[paste0("p_assign_p", last)]])
      fit <- lm(base_form, data = df, weights = w)
    }

    ct_tbl <- robust_test(fit, cluster = ~ community_id, hc_type = "HC1")
    get_p <- function(term) {
      row <- ct_tbl[ct_tbl$term == term, , drop = FALSE]
      if (nrow(row) == 1) row$p_value else NA_real_
    }
    out <- list()
    for (a in config$trt_arms) {
      out[[paste0("p_", a)]] <- get_p(paste0("arm_p", last, a))
    }
    out
  }

  # ---------- Sweep wrapper ---------- #
  run_for_value <- function(val) {
    cfg <- config
    if (sweep_param %in% c("effect_last_period","effect_p2")) {
      trt <- cfg$trt_arms
      # set the last-period effect for all trt arms to val (others unchanged)
      last <- cfg$T_periods
      eff_list <- cfg$effects_by_period
      eff_list[[last]] <- setNames(rep(val, length(trt)), trt)
      cfg$effects_by_period <- eff_list
    }

    # Monte Carlo
    if (parallel_inner) {
      res <- future_map(1:cfg$sims, ~ one_run(cfg), .options = furrr_options(seed = TRUE))
    } else {
      res <- map(1:cfg$sims, ~ one_run(cfg))
    }
    # aggregate
    out_df <- tibble(sweep_value = val)
    for (a in config$trt_arms) {
      pvals <- sapply(res, `[[`, paste0("p_", a))
      out_df[[paste0("power_", a)]] <- mean(pvals < config$alpha, na.rm = TRUE)
    }
    out_df
  }

  # ---------- Run across sweep ---------- #
  if (parallel_outer) {
    with_progress({
      results <- future_map_dfr(sweep_values, run_for_value, .progress = TRUE, .options = furrr_options(seed = TRUE))
    })
  } else {
    results <- map_dfr(sweep_values, run_for_value)
  }

  # ---------- Save outputs ---------- #
  if (is.null(output_dir_base)) {
    output_dir_base <- file.path(getwd(), "output", "dynamic")
  }
  tabs_dir <- file.path(output_dir_base, "tabs")
  figs_dir <- file.path(output_dir_base, "figs")
  if (!dir.exists(tabs_dir)) dir.create(tabs_dir, recursive = TRUE)
  if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)
  csv_path <- file.path(tabs_dir, glue::glue("{outfile_stem}_{sweep_param}.csv"))
  readr::write_csv(results, csv_path)

  # Plot
  melt <- results |>
    pivot_longer(cols = starts_with("power_"), names_to = "arm", values_to = "power") |>
    mutate(arm = sub("^power_", "", arm))

  # Labels per request
  arms_label <- if (length(config$trt_arms) > 1) "multi" else "single"
  design_label <- paste0(config$randomize_at, ", T=", config$T_periods)
  dist_label <- "normal"

  plt <- ggplot(melt, aes(x = sweep_value, y = power, color = arm, linetype = arm)) +
    geom_line(linewidth = 1) +
    geom_point() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    labs(
      x = sweep_param,
      y = "Power",
      title = glue::glue("Power for different {sweep_param}"),
      subtitle = glue::glue("Design: {design_label} | Distribution: {dist_label} | Arms: {arms_label}")
    ) +
    theme_minimal(base_size = 12)
  png_path <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}.png"))
  pdf_path <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}.pdf"))
  ggsave(png_path, plt, width = 8, height = 5, dpi = 300)
  ggsave(pdf_path, plt, width = 8, height = 5)

  list(results = results, csv = csv_path, png = png_path, pdf = pdf_path)
}
