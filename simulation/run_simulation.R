# ============================================================
# Monte Carlo Simulation for Adjusted DID Inference
# ============================================================
#
# This script reproduces the simulation tables reported in the manuscript.
#
# Outputs:
#   output/simulation_results_adid.csv
#   output/table1_tau0_se_ratio_rejection.csv
#   output/table2_tau05_se_ratio_rejection.csv
#
# Notes:
#   - tau = 0.0: rejection_rate is Type I error rate.
#   - tau = 0.5: rejection_rate is power.
# ============================================================

source("adid/R/adid.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(future)
  library(furrr)
})

# ------------------------------------------------------------
# Global simulation configuration
# ------------------------------------------------------------

SIM_CONFIG <- list(
  R = 1000,
  B_boot = 1000,
  seed = 2026,
  n_values = c(500, 1000, 2000),
  J = 50,
  tau_values = c(0.0, 0.5),
  school_sd_values = c(0.0, 1.0, 2.0),
  treat_assignment_values = c("individual", "mixed", "school"),
  beta_P = 0.8,
  beta_Y = 1.2,
  beta_C = 1.0,
  treat_strength = 0.6,
  school_treat_sd = 1.0
)

OUTPUT_DIR <- "output"

num_cores <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = num_cores)

message("Parallel processing activated. Workers: ", num_cores)

# ------------------------------------------------------------
# Data generation
# ------------------------------------------------------------

generate_data <- function(
    n = 1000,
    J = 50,
    tau = 0.0,
    beta_P = 0.8,
    beta_Y = 1.2,
    beta_C = 1.0,
    school_sd = 0.0,
    treat_strength = 0.6,
    treat_assignment = c("individual", "mixed", "school"),
    school_treat_sd = 1.0,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  treat_assignment <- match.arg(treat_assignment)
  
  school <- sample(seq_len(J), size = n, replace = TRUE)
  
  u_school <- rnorm(J, mean = 0, sd = school_sd)
  S <- u_school[school]
  
  U_school_effect <- rnorm(J, mean = 0, sd = 1)
  U <- U_school_effect[school] + rnorm(n, mean = 0, sd = 1)
  
  if (treat_assignment == "individual") {
    pr_g <- plogis(treat_strength * U)
    G <- rbinom(n, size = 1, prob = pr_g)
  }
  
  if (treat_assignment == "mixed") {
    school_propensity <- rnorm(J, mean = 0, sd = school_treat_sd)
    pr_g <- plogis(treat_strength * U + school_propensity[school])
    G <- rbinom(n, size = 1, prob = pr_g)
  }
  
  if (treat_assignment == "school") {
    school_treat <- rbinom(J, size = 1, prob = 0.5)
    G <- school_treat[school]
  }
  
  if (length(unique(G)) < 2) {
    stop("Only one treatment group was generated.")
  }
  
  C <- beta_C * U + rnorm(n, sd = 1)
  P <- beta_P * U + 0.6 * S + rnorm(n, sd = 1)
  Y <- tau * G + beta_Y * U + 0.6 * S + rnorm(n, sd = 1)
  
  data.frame(
    id = seq_len(n),
    school = school,
    U = U,
    S = S,
    G = G,
    C = C,
    P = P,
    Y = Y
  )
}

# ------------------------------------------------------------
# One replication
# ------------------------------------------------------------

one_rep <- function(
    n,
    J,
    tau,
    beta_P,
    beta_Y,
    beta_C,
    school_sd,
    treat_strength,
    treat_assignment,
    school_treat_sd,
    B_boot,
    seed
) {
  dat <- generate_data(
    n = n,
    J = J,
    tau = tau,
    beta_P = beta_P,
    beta_Y = beta_Y,
    beta_C = beta_C,
    school_sd = school_sd,
    treat_strength = treat_strength,
    treat_assignment = treat_assignment,
    school_treat_sd = school_treat_sd,
    seed = seed
  )
  
  # Fixed-coefficient method for paper comparison only.
  # This method is not included in the user-facing adid package API.
  X <- cbind(Intercept = 1, C = dat$C, G = dat$G)
  fit_post <- lm.fit(x = X, y = dat$Y)
  fit_pre <- lm.fit(x = X, y = dat$P)
  
  delta_hat <- fit_post$coefficients["C"] / fit_pre$coefficients["C"]
  adjusted_change <- dat$Y - delta_hat * dat$P
  fixed_fit <- summary(lm(adjusted_change ~ dat$G))
  
  se_fixed <- coef(fixed_fit)[2, 2]
  
  fit_delta <- adid(
    data = dat,
    pre = "P",
    post = "Y",
    treat = "G",
    compass = "C",
    cluster = "school",
    method = c("delta", "delta_cluster")
  )
  
  fit_boot <- adid(
    data = dat,
    pre = "P",
    post = "Y",
    treat = "G",
    compass = "C",
    cluster = "school",
    method = c("bootstrap", "bootstrap_cluster"),
    B = B_boot,
    seed = seed
  )
  
  data.frame(
    beta_hat = fit_delta$estimate[fit_delta$method == "delta"],
    se_fixed = se_fixed,
    se_boot_iid = fit_boot$se[fit_boot$method == "bootstrap"],
    se_delta_iid = fit_delta$se[fit_delta$method == "delta"],
    se_boot_cluster = fit_boot$se[fit_boot$method == "bootstrap_cluster"],
    se_delta_cluster = fit_delta$se[fit_delta$method == "delta_cluster"]
  )
}

# ------------------------------------------------------------
# Repeated simulation under one condition
# ------------------------------------------------------------

summarize_method <- function(se, beta, emp_sd) {
  valid <- is.finite(se) & se > 0 & is.finite(beta)
  
  if (sum(valid) == 0) {
    return(data.frame(
      se_ratio = NA_real_,
      rejection_rate = NA_real_,
      n_valid = 0
    ))
  }
  
  se_valid <- se[valid]
  beta_valid <- beta[valid]
  z_stat <- beta_valid / se_valid
  
  data.frame(
    se_ratio = sqrt(mean(se_valid^2, na.rm = TRUE)) / emp_sd,
    rejection_rate = mean(abs(z_stat) > 1.96, na.rm = TRUE),
    n_valid = sum(valid)
  )
}

run_condition <- function(
    R,
    n,
    J,
    tau,
    beta_P,
    beta_Y,
    beta_C,
    school_sd,
    treat_strength,
    treat_assignment,
    school_treat_sd,
    B_boot,
    seed
) {
  set.seed(seed)
  seeds <- sample.int(1e8, R)
  
  reps <- furrr::future_map_dfr(
    seq_len(R),
    function(r) {
      tryCatch(
        {
          one_rep(
            n = n,
            J = J,
            tau = tau,
            beta_P = beta_P,
            beta_Y = beta_Y,
            beta_C = beta_C,
            school_sd = school_sd,
            treat_strength = treat_strength,
            treat_assignment = treat_assignment,
            school_treat_sd = school_treat_sd,
            B_boot = B_boot,
            seed = seeds[r]
          )
        },
        error = function(e) {
          data.frame(
            beta_hat = NA_real_,
            se_fixed = NA_real_,
            se_boot_iid = NA_real_,
            se_delta_iid = NA_real_,
            se_boot_cluster = NA_real_,
            se_delta_cluster = NA_real_
          )
        }
      )
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
  
  reps <- reps[is.finite(reps$beta_hat), ]
  
  if (nrow(reps) == 0) {
    stop("All simulation replications failed.")
  }
  
  emp_sd <- stats::sd(reps$beta_hat)
  
  dplyr::bind_rows(
    summarize_method(reps$se_fixed, reps$beta_hat, emp_sd) |>
      dplyr::mutate(method = "fixed"),
    summarize_method(reps$se_boot_iid, reps$beta_hat, emp_sd) |>
      dplyr::mutate(method = "bootstrap_iid"),
    summarize_method(reps$se_delta_iid, reps$beta_hat, emp_sd) |>
      dplyr::mutate(method = "delta_iid"),
    summarize_method(reps$se_boot_cluster, reps$beta_hat, emp_sd) |>
      dplyr::mutate(method = "bootstrap_cluster"),
    summarize_method(reps$se_delta_cluster, reps$beta_hat, emp_sd) |>
      dplyr::mutate(method = "delta_cluster")
  ) |>
    dplyr::mutate(
      R_success = nrow(reps),
      R_requested = R,
      n = n,
      J = J,
      tau = tau,
      beta_P = beta_P,
      beta_Y = beta_Y,
      beta_C = beta_C,
      common_trend_violation = beta_Y - beta_P,
      school_sd = school_sd,
      treat_assignment = treat_assignment,
      school_treat_sd = school_treat_sd,
      rejection_type = ifelse(abs(tau) < 1e-12, "type1_error", "power")
    ) |>
    dplyr::relocate(
      n, J, treat_assignment, beta_C,
      beta_P, beta_Y, common_trend_violation,
      tau, school_sd, method,
      se_ratio, rejection_rate, rejection_type
    )
}

# ------------------------------------------------------------
# Main simulation grid
# ------------------------------------------------------------

run_grid_main <- function(config = SIM_CONFIG) {
  grid <- tidyr::expand_grid(
    n = config$n_values,
    J = config$J,
    tau = config$tau_values,
    school_sd = config$school_sd_values,
    treat_assignment = config$treat_assignment_values
  )
  
  purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    p <- grid[i, ]
    
    message(sprintf(
      "Running condition %d/%d: n=%d, J=%d, assignment=%s, tau=%.1f, school_sd=%.1f",
      i, nrow(grid), p$n, p$J, p$treat_assignment, p$tau, p$school_sd
    ))
    
    run_condition(
      R = config$R,
      n = p$n,
      J = p$J,
      tau = p$tau,
      beta_P = config$beta_P,
      beta_Y = config$beta_Y,
      beta_C = config$beta_C,
      school_sd = p$school_sd,
      treat_strength = config$treat_strength,
      treat_assignment = p$treat_assignment,
      school_treat_sd = config$school_treat_sd,
      B_boot = config$B_boot,
      seed = config$seed + i
    )
  })
}

# ------------------------------------------------------------
# Paper table formatting
# ------------------------------------------------------------

make_rejection_table <- function(grid_res, tau_value) {
  grid_res |>
    dplyr::filter(abs(tau - tau_value) < 1e-12) |>
    dplyr::group_by(treat_assignment, school_sd, method) |>
    dplyr::summarise(
      se_ratio = mean(se_ratio, na.rm = TRUE),
      rejection_rate = mean(rejection_rate, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      method_label = dplyr::recode(
        method,
        "fixed" = "Fixed",
        "bootstrap_iid" = "Individual BS",
        "delta_iid" = "Delta",
        "bootstrap_cluster" = "Cluster BS",
        "delta_cluster" = "Cluster Delta"
      ),
      method_label = factor(
        method_label,
        levels = c("Fixed", "Individual BS", "Delta", "Cluster BS", "Cluster Delta")
      ),
      treat_label = dplyr::recode(
        treat_assignment,
        "individual" = "Individual",
        "mixed" = "Mixed",
        "school" = "School"
      ),
      treat_label = factor(treat_label, levels = c("Individual", "Mixed", "School")),
      cell = sprintf("%.3f (%.3f)", se_ratio, rejection_rate)
    ) |>
    dplyr::select(treat_label, school_sd, method_label, cell) |>
    tidyr::pivot_wider(names_from = method_label, values_from = cell) |>
    dplyr::arrange(treat_label, school_sd)
}

# ------------------------------------------------------------
# Run full simulation
# ------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

full_time <- system.time({
  grid_res <- run_grid_main(SIM_CONFIG)
})

print(full_time)

table1_tau0 <- make_rejection_table(grid_res, tau_value = 0.0)
table2_tau05 <- make_rejection_table(grid_res, tau_value = 0.5)

write.csv(
  grid_res,
  file.path(OUTPUT_DIR, "simulation_results_adid.csv"),
  row.names = FALSE
)

write.csv(
  table1_tau0,
  file.path(OUTPUT_DIR, "table1_tau0_se_ratio_rejection.csv"),
  row.names = FALSE
)

write.csv(
  table2_tau05,
  file.path(OUTPUT_DIR, "table2_tau05_se_ratio_rejection.csv"),
  row.names = FALSE
)

message("Simulation completed.")
