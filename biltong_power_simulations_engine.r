# ------------------------------ Core Simulator ----------------------------- #

# Parallelized using {furrr} for multicore processing
# Uses 14 cores based on MacBook Pro M4 Pro
simulate_power <- function(
  config,
  sweep_param = c("effect_meeting","n_communities","avg_indiv_per_comm",
                  "T_meeting","take_up_T","ICC_meeting","AR1_meeting",
                  "alloc_ratio"),
  sweep_values = NULL,
  outfile_stem = "power_results",
  seed = 1
) {
  set.seed(seed)
  sweep_param <- match.arg(sweep_param)
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
    comm <- tibble(
      community_id = 1:cfg_local$n_communities,
      year_in_program = cfg_local$year_in_program,
      ngo_id = cfg_local$ngo_id,
      tribe_id = cfg_local$tribe_id
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
        mutate(times = list(cfg_local$times_meeting)) |>
        unnest(times) |>
        mutate(time = times, .keep = "unused") |>
        ungroup()
      # VCU panel: community-time
      df_v <- comm |>
        rowwise() |>
        mutate(times = list(cfg_local$times_VCU)) |>
        unnest(times) |>
        mutate(time = times, .keep = "unused") |>
        ungroup()
    } else {
      # individual_within_community
      df_p <- comm |>
        rowwise() |>
        mutate(indiv = list(1:N_indiv), times = list(cfg_local$times_meeting)) |>
        unnest(c(indiv, times)) |>
        rename(time = times) |>
        ungroup()
      # VCU remains at community-time
      df_v <- comm |>
        rowwise() |>
        mutate(times = list(cfg_local$times_VCU)) |>
        unnest(times) |>
        mutate(time = times, .keep = "unused") |>
        ungroup()
    }

    list(comm = comm, df_p = df_p, df_v = df_v)
  }

  # One Monte Carlo run returning p-values for ITT and TOT for both outcomes
  one_run <- function(cfg_local) {
    smp <- make_sample(cfg_local)

    # Assignment
    if (cfg_local$experiment_type == "community_level") {
      # Create treatment assignment Z at the community level:
      # Z == 1 means the community was assigned to treatment, Z == 0 means control.
      # Assignment is balanced within strata (year_in_program, ngo_id, tribe_id).
      Zc <- stratified_assign(smp$comm, treat_share = cfg_local$alloc_ratio) |>
        transmute(community_id, Z)
      df_p <- smp$df_p |> left_join(Zc, by = "community_id")
      df_v <- smp$df_v |> left_join(Zc, by = "community_id")
    } else {
      ind <- smp$df_p |>
        distinct(community_id, indiv) |>
        left_join(smp$comm, by = "community_id")
      # Create individual-level treatment assignment Z within strata (1/0).
      # For the individual-within-community design we assign at the individual level
      # and later aggregate to community-level shares (Z_share) for VCU outcomes.
      ind <- stratified_assign(ind, treat_share = cfg_local$alloc_ratio)
      df_p <- smp$df_p |> left_join(ind[,c("community_id","indiv","Z")], by = c("community_id","indiv"))
      # VCU intensity vars
      Zshare <- ind |>
        group_by(community_id) |>
        summarize(Z_share = mean(Z), .groups = "drop")
      df_v <- smp$df_v |> left_join(Zshare, by = "community_id")
    }

    # Compliance
    if (cfg_local$experiment_type == "community_level") {
      # Create compliance indicator D at the community level (actual uptake).
      # D is drawn conditionally on assignment Z using take_up_T / take_up_C.
      Dc <- smp$comm |>
        left_join(Zc, by = "community_id") |>
        mutate(D = ifelse(Z==1, rbinom(n(), 1, cfg_local$take_up_T), rbinom(n(), 1, cfg_local$take_up_C))) |>
        select(community_id, D)
      df_p <- df_p |> left_join(Dc, by = "community_id")
      df_v <- df_v |> left_join(Dc, by = "community_id")
    } else {
      # Individual-level compliance: create D for each individual conditional on Z.
      # Then aggregate to community-level average compliance D_bar and Z_share for VCU.
      df_p <- df_p |>
        mutate(D = ifelse(Z==1, rbinom(n(), 1, cfg_local$take_up_T), rbinom(n(), 1, cfg_local$take_up_C)))
      agg <- df_p |>
        group_by(community_id) |>
        summarize(D_bar = mean(D), Z_share = mean(Z), .groups = "drop")
      df_v <- df_v |> left_join(agg, by = "community_id")
      # For VCU we use the aggregated intensity measures
      df_v$Z <- df_v$Z_share
      df_v$D <- df_v$D_bar
    }

    # Participation outcome generation (counts)
    if (cfg_local$experiment_type == "community_level") {
      var_e_p <- 1
      var_u_p <- ifelse(cfg_local$ICC_meeting >= 0.999, 1e6, cfg_local$ICC_meeting/(1-cfg_local$ICC_meeting))
      u_comm_p <- rnorm(cfg_local$n_communities, 0, sqrt(var_u_p))
      df_p <- df_p |>
        left_join(tibble(community_id = 1:cfg_local$n_communities, u_comm_p), by = "community_id") |>
        group_by(community_id) |>
        arrange(time, .by_group = TRUE) |>
        mutate(e_ar = ar1_series(time, rho = cfg_local$AR1_meeting, sigma = sqrt(var_e_p))) |>
        ungroup()

      log_mu_base <- log(pmax(cfg_local$meeting_base_mean, 1e-6)) + df_p$u_comm_p + df_p$e_ar
      log_mu <- log_mu_base + log(cfg_local$meeting_rate_ratio) * df_p$D  # effect realized if D==1
      mu <- pmax(exp(log_mu), 1e-6)
      df_p$Y <- r_count(length(mu), mu, family = cfg_local$participation_dist, theta = cfg_local$theta_meeting)
    } else {
      var_e_p <- 1
      var_u_p <- ifelse(cfg_local$ICC_meeting >= 0.999, 1e6, cfg_local$ICC_meeting/(1-cfg_local$ICC_meeting))
      u_comm_p <- rnorm(cfg_local$n_communities, 0, sqrt(var_u_p))
      df_p <- df_p |>
        left_join(tibble(community_id = 1:cfg_local$n_communities, u_comm_p), by = "community_id") |>
        group_by(community_id, indiv) |>
        arrange(time, .by_group = TRUE) |>
        mutate(e_ar = ar1_series(time, rho = cfg_local$AR1_meeting, sigma = sqrt(var_e_p))) |>
        ungroup()

      log_mu_base <- log(pmax(cfg_local$meeting_base_mean, 1e-6)) + df_p$u_comm_p + df_p$e_ar
      log_mu <- log_mu_base + log(cfg_local$meeting_rate_ratio) * df_p$D
      mu <- pmax(exp(log_mu), 1e-6)
      df_p$Y <- r_count(length(mu), mu, family = cfg_local$participation_dist, theta = cfg_local$theta_meeting)
    }

    # Participation OLS with FE
    fit_itt_p <- lm(Y ~ Z + factor(time) + factor(community_id), data = df_p)
    fit_tot_p <- AER::ivreg(Y ~ D + factor(time) + factor(community_id) |
                               Z + factor(time) + factor(community_id), data = df_p)
    tt_itt_p <- robust_test(fit_itt_p, hc_type = cfg_local$hc_type)
    tt_tot_p <- robust_test(fit_tot_p, hc_type = cfg_local$hc_type)
    p_itt_p <- tt_itt_p$p_value[tt_itt_p$term == "Z"]
    p_tot_p <- tt_tot_p$p_value[tt_tot_p$term == "D"]

    # VCU generation: continuous at community-time
    var_e_v <- cfg_local$VCU_base_sd^2
    var_u_v <- ifelse(cfg_local$ICC_VCU >= 0.999, 1e6, (var_e_v * cfg_local$ICC_VCU) / (1 - cfg_local$ICC_VCU))
    u_comm_v <- rnorm(cfg_local$n_communities, 0, sqrt(pmax(var_u_v, 1e-8)))
    df_v <- df_v |>
      left_join(tibble(community_id = 1:cfg_local$n_communities, u_comm_v), by = "community_id") |>
      group_by(community_id) |>
      arrange(time, .by_group = TRUE) |>
      mutate(e_ar = ar1_series(time, rho = cfg_local$AR1_VCU, sigma = sqrt(pmax(var_e_v, 1e-8)))) |>
      ungroup()

    # Check if D column already exists from compliance step above
    if (!"D" %in% names(df_v)) {
      if (cfg_local$experiment_type == "community_level") {
        # Re-add compliance for VCU if missing
        df_v <- df_v |> left_join(Dc, by = "community_id")
      } else {
        # For individual-within-community, create aggregated compliance
        agg_vcu <- df_p |>
          group_by(community_id) |>
          summarize(D_bar = mean(D, na.rm = TRUE), Z_share = mean(Z, na.rm = TRUE), .groups = "drop")
        df_v <- df_v |> left_join(agg_vcu, by = "community_id")
        df_v$Z <- df_v$Z_share
        df_v$D <- df_v$D_bar
      }
    }
    
    # Final check for required columns
    if (!"D" %in% names(df_v)) stop("Column D still missing in df_v after compliance correction")
    if (!"u_comm_v" %in% names(df_v)) stop("Column u_comm_v missing in df_v")
    if (!"e_ar" %in% names(df_v)) stop("Column e_ar missing in df_v")
    
    df_v$Y <- cfg_local$VCU_base_mean + df_v$u_comm_v + df_v$e_ar + cfg_local$VCU_delta_mean * df_v$D

    # VCU OLS with FE
    fit_itt_v <- lm(Y ~ Z + factor(time) + factor(community_id), data = df_v)
    fit_tot_v <- AER::ivreg(Y ~ D + factor(time) + factor(community_id) |
                               Z + factor(time) + factor(community_id), data = df_v)
    tt_itt_v <- robust_test(fit_itt_v, hc_type = cfg_local$hc_type)
    tt_tot_v <- robust_test(fit_tot_v, hc_type = cfg_local$hc_type)
    p_itt_v <- tt_itt_v$p_value[tt_itt_v$term == "Z"]
    p_tot_v <- tt_tot_v$p_value[tt_tot_v$term == "D"]

    list(p_itt_p = p_itt_p, p_tot_p = p_tot_p, p_itt_v = p_itt_v, p_tot_v = p_tot_v)
  }

  # Run sims for a given sweep value
  run_for_value <- function(val) {
    cfg2 <- config
    if (sweep_param == "effect_meeting") cfg2$meeting_rate_ratio <- val
    if (sweep_param == "n_communities")  cfg2$n_communities <- as.integer(val)
    if (sweep_param == "avg_indiv_per_comm") cfg2$avg_indiv_per_comm <- as.integer(val)
    if (sweep_param == "T_meeting") {
      cfg2$T_meeting <- as.integer(val)
      cfg2$times_meeting <- seq_len(cfg2$T_meeting)
    }
    if (sweep_param == "take_up_T") cfg2$take_up_T <- val
    if (sweep_param == "ICC_meeting") cfg2$ICC_meeting <- val
    if (sweep_param == "AR1_meeting") cfg2$AR1_meeting <- val
    if (sweep_param == "alloc_ratio") cfg2$alloc_ratio <- val

    res <- replicate(cfg2$sims, one_run(cfg2), simplify = FALSE)

    get_rate <- function(xvec) mean(xvec < cfg2$alpha, na.rm = TRUE)
    p_itt_p <- sapply(res, `[[`, "p_itt_p"); pow_itt_p <- get_rate(p_itt_p)
    p_tot_p <- sapply(res, `[[`, "p_tot_p"); pow_tot_p <- get_rate(p_tot_p)
    p_itt_v <- sapply(res, `[[`, "p_itt_v"); pow_itt_v <- get_rate(p_itt_v)
    p_tot_v <- sapply(res, `[[`, "p_tot_v"); pow_tot_v <- get_rate(p_tot_v)

    tibble(
      sweep_param = sweep_param,
      sweep_value = val,
      experiment_type = cfg2$experiment_type,
      participation_dist = cfg2$participation_dist,
      power_ITT_participation = pow_itt_p,
      power_TOT_participation = pow_tot_p,
      power_ITT_VCU = pow_itt_v,
      power_TOT_VCU = pow_tot_v
    )
  }

  # Run across sweep with progress bar
  with_progress({
    results <- future_map_dfr(sweep_values, run_for_value, .progress = TRUE)
  })

  cat("Simulation complete. Saving CSV and figure files to both output locations...\n")

  # Set output directory
  output_dir <- "/Users/flavioamalagutti/Library/CloudStorage/GoogleDrive-fmalagutti@ucsb.edu/Shared drives/emlab/projects/current-projects/arnhold-biltong/biltong_rct_design/Power_simulations"
  secondary_output_dir <- "/Users/flavioamalagutti/Documents/Work/GitHub/Research_Projects/Biltong/Simulations/Biltong_power_simulations"
  
  # Create output directories if they don't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  if (!dir.exists(secondary_output_dir)) {
    dir.create(secondary_output_dir, recursive = TRUE)
  }
  
  # Save CSV
  csv_path <- file.path(output_dir, glue("{outfile_stem}_{sweep_param}.csv"))
  csv_path2 <- file.path(secondary_output_dir, glue("{outfile_stem}_{sweep_param}.csv"))
  write_csv(results, csv_path)
  
  # Figure
  plt <- results |>
    pivot_longer(cols = starts_with("power_"), names_to = "metric", values_to = "power") |>
    mutate(series = case_when(
      metric == "power_ITT_participation" ~ "ITT - Participation",
      metric == "power_TOT_participation" ~ "TOT - Participation",
      metric == "power_ITT_VCU" ~ "ITT - VCU",
      metric == "power_TOT_VCU" ~ "TOT - VCU",
      TRUE ~ metric
    )) |>
    ggplot(aes(x = sweep_value, y = power, group = series, color = series)) +
    geom_line(size = 0.7) +
    geom_point(size = 2.2) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "#666666") +
    scale_color_manual(values = c(
      `ITT - Participation` = "#1b9e77",
      `TOT - Participation` = "#d95f02",
      `ITT - VCU` = "#7570b3",
      `TOT - VCU` = "#e7298a"
    )) +
    scale_x_continuous(breaks = unique(results$sweep_value)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1), expand = c(0, 0)) +
    labs(
      x = sweep_param,
      y = "Power",
      title = glue("Power vs {sweep_param}"),
      subtitle = glue("Design: {unique(results$experiment_type)} | Participation: {unique(results$participation_dist)}"),
      caption = glue(
        "alpha={config$alpha}, sims={config$sims}, seed={seed}, alloc_ratio={config$alloc_ratio}, n_communities={config$n_communities}, avg_indiv_per_comm={config$avg_indiv_per_comm}, T_meeting={config$T_meeting}, T_VCU={config$T_VCU};\n",
        "year strata range={paste(range(config$year_in_program), collapse='-')}, NGOs up to 7, tribes up to 5, participation_dist={config$participation_dist}, theta={config$theta_meeting}, base_mean_meeting={config$meeting_base_mean};\n",
        "VCU_base_mean={config$VCU_base_mean}, VCU_base_sd={config$VCU_base_sd}, ICC_meeting={config$ICC_meeting}, ICC_VCU={config$ICC_VCU}, AR1_meeting={config$AR1_meeting}, AR1_VCU={config$AR1_VCU}, take_up_T={config$take_up_T}, take_up_C={config$take_up_C}, SE=Eicker-White {config$hc_type}" # nolint: line_length_linter.
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.width = unit(1.2, "cm"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(hjust = 0, size = 8)
    )

  png_path <- file.path(output_dir, glue("{outfile_stem}_{sweep_param}.png"))
  pdf_path <- file.path(output_dir, glue("{outfile_stem}_{sweep_param}.pdf"))
  png_path2 <- file.path(secondary_output_dir, glue("{outfile_stem}_{sweep_param}.png"))
  pdf_path2 <- file.path(secondary_output_dir, glue("{outfile_stem}_{sweep_param}.pdf"))
  ggsave(png_path, plt, width = 8, height = 5, dpi = 300)
  ggsave(pdf_path, plt, width = 8, height = 5)
  # Duplicate exports to secondary location
  write_csv(results, csv_path2)
  ggsave(png_path2, plt, width = 8, height = 5, dpi = 300)
  ggsave(pdf_path2, plt, width = 8, height = 5)
  # Optional message
  message("Exported results to:\n", csv_path, "\n", csv_path2)

  list(results = results, csv = csv_path, png = png_path, pdf = pdf_path, plot = plt)
}