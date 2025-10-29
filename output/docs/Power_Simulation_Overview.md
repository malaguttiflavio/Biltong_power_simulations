# Biltong/GRASS Power Simulation — Practical Guide

This guide explains what the simulator does, what to edit, and how to read its outputs. It’s concise and geared for day-to-day use.

## What the program does

It estimates statistical power for GRASS/Meat Naturally RCTs by:
- Randomizing communities to treatment arms (with stratification),
- Simulating compliance (take-up),
- Generating outcomes over time with clustering and AR(1) correlation,
- Estimating ITT (OLS) and TOT (IV) effects with robust/cluster-robust SEs,
- Repeating many times and computing power = share of significant results.

Outputs are saved to:
- Tables: `output/tabs`
- Figures: `output/figs`

## Files

- `Biltong_power_simulations_master.r` — sets the scenario and runs sweeps.
- `biltong_power_simulations_engine.r` — simulates data, runs estimators, writes outputs.
- `biltong_power_simulations_utils.r` — helpers (AR(1), Poisson/NegBin draws, robust tests).

## What you edit in the master script

- Sample and design
  - `n_communities`, `avg_ind_obs_per_comm`, `sd_indiv_per_comm`
  - `arms_*` and `alloc_ratios_*` (single vs multi-arm presets)
  - `take_up_*` (per-arm compliance probabilities)
  - `experiment_type` = `"community_level"` or `"individual_within_community"`
  - Stratification vectors: `year_in_program`, `ngo_id`, `tribe_id`
- Outcomes and noise
  - Socio-economic counts: `soc_outcome_dist` ("negbin" | "poisson" | "none"), `soc_outcome_theta`, `soc_outcome_ICC`, `soc_outcome_AR1_rho`, `soc_outcome_AR1_var`
  - Environmental continuous: `env_outcome_base_mean`, `env_outcome_base_sd`, `env_outcome_ICC`, `env_outcome_AR1_rho`, `env_outcome_AR1_var`
  - Time grids: `soc_outcome_T`, `soc_outcome_T_months`; `env_outcome_T`, `env_outcome_T_months`
- Effects by arm (omit "control")
  - Soc (multiplicative “rate ratio”): `soc_outcome_ate_pct_*` (e.g., `c(T1 = 1.20)`)
  - Env (additive): `env_outcome_ate_*` (e.g., `c(T1 = 10)`)
- Inference and repetitions
  - `sims`, `alpha`, `cluster_se`, `hc_type`, `cluster_fe_yn`

## How one dataset is generated

### Randomization (with stratification)
Communities are assigned to arms using `alloc_ratios`, balanced within strata (`year_in_program`, `ngo_id`, `tribe_id`).

### Compliance (take-up)
Realized treatment is a Bernoulli draw per observation with probability equal to the arm’s `take_up`. ITT uses assignment. TOT instruments realized treatment with assignment.

### Outcome models

Socio-economic (counts; Negative Binomial or Poisson)

- Log-mean model for community c at time t:

$$
\log \mu_{ct} = \log(\text{base}) + u_c + e_{ct} + [\log(m_a) + \tau_{\text{strata},c}]\, D_{ct}
$$

- $u_c$ is the cluster random intercept from ICC; $e_{ct}$ is AR(1): $e_{ct} = \rho e_{c,t-1} + \varepsilon_{ct}$.
- $m_a$ is the treatment arm’s multiplicative rate ratio; $\tau_{\text{strata},c}$ adds stratum-specific TE when treated.
- Draws:
  - Negative binomial if `soc_outcome_dist == "negbin"`:  
  $$Y_{ct} \sim \text{NegBin}(\mu_{ct}, \theta),\quad p = \frac{\theta}{\theta + \mu_{ct}}$$
  - Poisson if `"poisson"`: $Y_{ct} \sim \text{Poisson}(\mu_{ct})$
  - Deterministic if `"none"`: $Y_{ct} = \mu_{ct}$

Environmental (continuous)

$$
Y_{ct} = \text{base} + u_c + e_{ct} + [\delta_a + \tau_{\text{strata},c}]\, D_{ct}
$$

- Same ICC/AR(1) structure. $\delta_a$ is the additive arm effect.

ICC and AR(1)

- Cluster intercept variance:  
$$\sigma_u^2 = \frac{\text{ICC}}{1 - \text{ICC}}$$
- AR(1): $e_{ct} = \rho e_{c,t-1} + \varepsilon_{ct}$ with innovation variance `*_AR1_var`.

## Estimation per dataset

- ITT (OLS):

$$
Y \sim \text{arm indicators} + \text{time FE} + [\text{optional strata FE}] + [\text{optional community FE}]
$$

- TOT (IV, per arm using only control + that arm): first stage $D \sim Z + \cdots$, second stage $Y \sim \hat D + \cdots$.
- SEs: cluster-robust when `cluster_se=TRUE`, else HC-type (`hc_type`).

## Power and sweeps

Repeat the above `sims` times. For each arm and outcome:

- Power(ITT) = share of $p$-values < `alpha` in the ITT model.
- Power(TOT) = share of $p$-values < `alpha` in the IV model.

Typical sweeps vary effect size (rate ratio or additive effect) and plot power vs. effect; you can also sweep ICC, AR(1), T, or sample size.

## Outputs you’ll see

- Tables (wide and long) in `output/tabs`, including the last simulated dataset for a quick realism check.
- Figures in `output/figs` with captions listing arms, allocations, number/timing of observations, take-up, distributions, and inference settings.

## Quick checks

- Negative binomial in use: set `soc_outcome_dist = "negbin"` and confirm in the caption; the last simulated socio-economic data should have variance > mean.
- Avoid community FE in ITT for pure community-level RCTs (assignment is time-invariant within community).
- Watch cores when parallelizing.

## Minimal math summary

- AR(1): $e_{ct} = \rho e_{c,t-1} + \varepsilon_{ct}$
- Cluster intercept variance: $\sigma_u^2 = \dfrac{\text{ICC}}{1-\text{ICC}}$
- Socio-economic mean: $\log \mu_{ct} = \log(\text{base}) + u_c + e_{ct} + [\log(m_a)+\tau_{\text{strata},c}]D_{ct}$
- Environmental level: $Y_{ct} = \text{base} + u_c + e_{ct} + [\delta_a+\tau_{\text{strata},c}]D_{ct}$
- Power: $\widehat{\text{Power}} = \frac{1}{\text{sims}} \sum_{s=1}^{\text{sims}} \mathbf{1}\{p_s < \alpha\}$
