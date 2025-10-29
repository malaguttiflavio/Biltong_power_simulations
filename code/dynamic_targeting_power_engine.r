# ---------------- Dynamic Targeting Power Simulator ---------------- #

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(furrr)
  library(progressr)
  library(lmtest)
  library(sandwich)
})

simulate_dynamic_power <- function(
  config,
  sweep_param = c("effect_p2"),
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

    # Period 1 randomization ratios
    if (is.null(cfg$alloc_ratios_p1)) cfg$alloc_ratios_p1 <- c(control = 1/3, T1 = 1/3, T2 = 1/3)
    # Period 2 randomization ratios for eligible and ineligible
    if (is.null(cfg$alloc_ratios_p2_eligible)) cfg$alloc_ratios_p2_eligible <- c(control = 0.1, T1 = 0.45, T2 = 0.45)
    if (is.null(cfg$alloc_ratios_p2_ineligible)) cfg$alloc_ratios_p2_ineligible <- c(control = 1.0, T1 = 0.0, T2 = 0.0)

    # Eligibility share (top performers by Y1 within community)
    if (is.null(cfg$target_top_share)) cfg$target_top_share <- 0.2

    # Take-up (compliance) default full
    if (is.null(cfg$take_up_p1)) cfg$take_up_p1 <- c(control = 0, T1 = 1, T2 = 1)
    if (is.null(cfg$take_up_p2)) cfg$take_up_p2 <- c(control = 0, T1 = 1, T2 = 1)

    # Outcome DGP parameters
    if (is.null(cfg$y_base_mean)) cfg$y_base_mean <- 10
    if (is.null(cfg$y_sd)) cfg$y_sd <- 5
    if (is.null(cfg$ability_sd)) cfg$ability_sd <- 3
    if (is.null(cfg$icc_y)) cfg$icc_y <- 0.05
    if (cfg$icc_y < 0 || cfg$icc_y >= 1) stop("icc_y must be in [0,1)")
    if (is.null(cfg$effect_p1)) cfg$effect_p1 <- c(T1 = 0, T2 = 0)
    if (is.null(cfg$effect_p2)) cfg$effect_p2 <- c(T1 = 2, T2 = 3)

    # Convert alloc ratios to probs normalized
    norm <- function(x, nm) {
      if (is.null(names(x))) names(x) <- cfg$arms
      x <- x[cfg$arms]
      if (any(is.na(x))) stop(nm, " must include all arms")
      x / sum(x)
    }
    cfg$alloc_ratios_p1 <- norm(cfg$alloc_ratios_p1, "alloc_ratios_p1")
    cfg$alloc_ratios_p2_eligible <- norm(cfg$alloc_ratios_p2_eligible, "alloc_ratios_p2_eligible")
    cfg$alloc_ratios_p2_ineligible <- norm(cfg$alloc_ratios_p2_ineligible, "alloc_ratios_p2_ineligible")

    # Derived
    cfg$trt_arms <- setdiff(cfg$arms, "control")

    cfg
  }
  config <- normalize(config)

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

    # Individual ability
    df$ability <- rnorm(nrow(df), 0, cfg$ability_sd)

    # -------- Period 1 assignment (unconstrained) -------- #
    p1_probs <- as.numeric(cfg$alloc_ratios_p1[cfg$arms])
    df$arm_p1 <- factor(sample(cfg$arms, size = nrow(df), replace = TRUE, prob = p1_probs), levels = cfg$arms)
    # Compliance -> D1
    df$D1 <- 0
    for (a in cfg$trt_arms) df$D1[df$arm_p1 == a] <- rbinom(sum(df$arm_p1 == a), 1, cfg$take_up_p1[[a]])

    # Period 1 outcomes
    eff1 <- setNames(rep(0, length(cfg$arms)), cfg$arms); eff1[names(cfg$effect_p1)] <- cfg$effect_p1
    mu1 <- cfg$y_base_mean + df$ability + df$u_c + eff1[as.character(df$arm_p1)]*df$D1
    df$Y1 <- rnorm(nrow(df), mean = mu1, sd = cfg$y_sd)

    # -------- Determine eligibility (top performers by Y1 within community) -------- #
    q <- cfg$target_top_share
    df <- df |>
      group_by(community_id) |>
      mutate(cut = quantile(Y1, probs = 1 - q, type = 7), eligible = as.integer(Y1 >= cut)) |>
      ungroup() |>
      select(-cut)

    # -------- Period 2 assignment targeted -------- #
    # Eligible get alloc_ratios_p2_eligible; ineligible get alloc_ratios_p2_ineligible
    probs_elig <- as.numeric(cfg$alloc_ratios_p2_eligible[cfg$arms])
    probs_inel <- as.numeric(cfg$alloc_ratios_p2_ineligible[cfg$arms])
    draw_arm <- function(n, probs) factor(sample(cfg$arms, size = n, replace = TRUE, prob = probs), levels = cfg$arms)
    df$arm_p2 <- NA
    df$arm_p2[df$eligible == 1] <- as.character(draw_arm(sum(df$eligible == 1), probs_elig))
    df$arm_p2[df$eligible == 0] <- as.character(draw_arm(sum(df$eligible == 0), probs_inel))
    df$arm_p2 <- factor(df$arm_p2, levels = cfg$arms)
    # Compliance -> D2
    df$D2 <- 0
    for (a in cfg$trt_arms) df$D2[df$arm_p2 == a] <- rbinom(sum(df$arm_p2 == a), 1, cfg$take_up_p2[[a]])

    # Period 2 outcomes
    eff2 <- setNames(rep(0, length(cfg$arms)), cfg$arms); eff2[names(cfg$effect_p2)] <- cfg$effect_p2
    mu2 <- cfg$y_base_mean + df$ability + df$u_c + eff2[as.character(df$arm_p2)]*df$D2
    df$Y2 <- rnorm(nrow(df), mean = mu2, sd = cfg$y_sd)

    # -------- Inference (ITT) on period 2 -------- #
    # OLS: Y2 ~ arm_p2 + Y1 + eligible, cluster-robust by community
    df$arm_p2 <- factor(df$arm_p2, levels = cfg$arms)
    base_form <- as.formula("Y2 ~ arm_p2 + Y1 + eligible")
    fit <- lm(base_form, data = df)
    vc <- tryCatch(sandwich::vcovCL(fit, cluster = ~ community_id), error = function(e) sandwich::vcovHC(fit, type = "HC1"))
    ct <- lmtest::coeftest(fit, vcov. = vc)
    get_p <- function(term) if (term %in% rownames(ct)) ct[term, "Pr(>|t|)"] else NA_real_
    out <- list(
      p_T1 = get_p("arm_p2T1"),
      p_T2 = get_p("arm_p2T2")
    )
    out
  }

  # ---------- Sweep wrapper ---------- #
  run_for_value <- function(val) {
    cfg <- config
    if (sweep_param == "effect_p2") {
      trt <- cfg$trt_arms
      cfg$effect_p2 <- setNames(rep(val, length(trt)), trt)
    }

    # Monte Carlo
    if (parallel_inner) {
      res <- future_map(1:cfg$sims, ~ one_run(cfg), .options = furrr_options(seed = TRUE))
    } else {
      res <- map(1:cfg$sims, ~ one_run(cfg))
    }
    p_T1 <- sapply(res, `[[`, "p_T1")
    p_T2 <- sapply(res, `[[`, "p_T2")
    tibble(
      sweep_value = val,
      power_T1 = mean(p_T1 < config$alpha, na.rm = TRUE),
      power_T2 = mean(p_T2 < config$alpha, na.rm = TRUE)
    )
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
  # Determine base output directory
  if (is.null(output_dir_base)) {
    # Fallback: write under working directory's output/dynamic
    output_dir_base <- file.path(getwd(), "output", "dynamic")
  }
  tabs_dir <- file.path(output_dir_base, "tabs")
  figs_dir <- file.path(output_dir_base, "figs")
  if (!dir.exists(tabs_dir)) dir.create(tabs_dir, recursive = TRUE)
  if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)
  csv_path <- file.path(tabs_dir, glue::glue("{outfile_stem}_{sweep_param}.csv"))
  readr::write_csv(results, csv_path)

  plt <- ggplot(results, aes(x = sweep_value)) +
    geom_line(aes(y = power_T1, color = "T1"), linewidth = 1) +
    geom_point(aes(y = power_T1, color = "T1")) +
    geom_line(aes(y = power_T2, color = "T2"), linewidth = 1, linetype = "dashed") +
    geom_point(aes(y = power_T2, color = "T2")) +
    scale_color_manual(values = c(T1 = "#1b9e77", T2 = "#d95f02"), name = "Arm") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    labs(x = sweep_param, y = "Power", title = "Dynamic targeting: period-2 power vs effect size") +
    theme_minimal(base_size = 12)
  png_path <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}.png"))
  pdf_path <- file.path(figs_dir, glue::glue("{outfile_stem}_{sweep_param}.pdf"))
  ggsave(png_path, plt, width = 8, height = 5, dpi = 300)
  ggsave(pdf_path, plt, width = 8, height = 5)

  list(results = results, csv = csv_path, png = png_path, pdf = pdf_path)
}
