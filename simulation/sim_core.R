# ============================================================
# simulation/sim_core.R
# ============================================================
# 조정된 이중차분(aDID) 추론 논문의 Monte Carlo 모의실험 공통 함수
#
# 이 파일은 다음 실행 스크립트에서 불러온다.
#   - simulation/run_chunk.R
#   - simulation/run_all_chunks.R
#   - simulation/build_tables.R
#
# 설계 원칙
#   전체 모의실험 격자를 한 번에 실행하지 않고,
#   tau x 처치배정 방식 단위로 나누어 실행한다.
#
# 예시 chunk
#   - tau = 0.0, assignment = individual
#   - tau = 0.0, assignment = mixed
#   - tau = 0.0, assignment = school
#   - tau = 0.5, assignment = individual
#   - tau = 0.5, assignment = mixed
#   - tau = 0.5, assignment = school
#
# 각 chunk 결과는 output/chunks/에 따로 저장하고,
# build_tables.R에서 이를 모아 논문 표 형태로 재구성한다.
# ============================================================


# ------------------------------------------------------------
# 경로 설정
# ------------------------------------------------------------

.infer_sim_dir <- function() {
  if (exists("SIM_DIR", envir = .GlobalEnv, inherits = FALSE)) {
    return(get("SIM_DIR", envir = .GlobalEnv))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  cwd <- normalizePath(getwd(), mustWork = FALSE)

  if (basename(cwd) == "simulation" && file.exists(file.path(cwd, "sim_core.R"))) {
    return(cwd)
  }

  if (file.exists(file.path(cwd, "simulation", "sim_core.R"))) {
    return(file.path(cwd, "simulation"))
  }

  cwd
}

SIM_DIR <- normalizePath(.infer_sim_dir(), mustWork = FALSE)
PROJECT_DIR <- normalizePath(file.path(SIM_DIR, ".."), mustWork = FALSE)

ADID_FILE <- file.path(PROJECT_DIR, "adid", "R", "adid.R")

if (!file.exists(ADID_FILE)) {
  stop(
    "adid/R/adid.R 파일을 찾을 수 없습니다. 예상 경로: ",
    ADID_FILE,
    call. = FALSE
  )
}

source(ADID_FILE)


# ------------------------------------------------------------
# 필요한 패키지
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(future)
  library(furrr)
})


# ------------------------------------------------------------
# 작은 보조 함수
# ------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

default_workers <- function() {
  cores <- parallel::detectCores()

  if (is.na(cores)) {
    return(1L)
  }

  max(1L, cores - 1L)
}

set_future_plan <- function(workers = default_workers()) {
  workers <- as.integer(workers)

  if (is.na(workers) || workers < 1) {
    workers <- 1L
  }

  if (.Platform$OS.type == "windows" || !future::supportsMulticore()) {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::multicore, workers = workers)
  }

  invisible(workers)
}

safe_label <- function(x) {
  x <- as.character(x)
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

num_label <- function(x) {
  out <- format(x, scientific = FALSE, trim = TRUE)
  out <- gsub("-", "m", out)
  out <- gsub("[.]", "p", out)
  out
}

parse_key_value_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()

  if (length(args) == 0) {
    return(out)
  }

  for (arg in args) {
    kv <- strsplit(arg, "=", fixed = TRUE)[[1]]

    if (length(kv) >= 2) {
      key <- kv[1]
      value <- paste(kv[-1], collapse = "=")
      out[[key]] <- value
    }
  }

  out
}

parse_numeric_vector <- function(x) {
  as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
}

parse_character_vector <- function(x) {
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

parse_logical_scalar <- function(x) {
  val <- as.logical(x)
  isTRUE(val)
}


# ------------------------------------------------------------
# 기본 모의실험 설정
# ------------------------------------------------------------

SIM_CONFIG <- list(
  run_label = "main",

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

OUTPUT_DIR <- file.path(PROJECT_DIR, "output")
CHUNK_DIR <- file.path(OUTPUT_DIR, "chunks")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CHUNK_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)


# ------------------------------------------------------------
# 자료 생성 함수
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
  if (!is.null(seed)) {
    set.seed(seed)
  }

  treat_assignment <- match.arg(treat_assignment)

  school <- sample(seq_len(J), size = n, replace = TRUE)

  u_school <- rnorm(J, mean = 0, sd = school_sd)
  S <- u_school[school]

  U_school_effect <- rnorm(J, mean = 0, sd = 1)
  U <- U_school_effect[school] + rnorm(n, mean = 0, sd = 1)

  if (treat_assignment == "individual") {
    pr_g <- stats::plogis(treat_strength * U)
    G <- stats::rbinom(n, size = 1, prob = pr_g)
  }

  if (treat_assignment == "mixed") {
    school_propensity <- rnorm(J, mean = 0, sd = school_treat_sd)
    pr_g <- stats::plogis(treat_strength * U + school_propensity[school])
    G <- stats::rbinom(n, size = 1, prob = pr_g)
  }

  if (treat_assignment == "school") {
    school_treat <- stats::rbinom(J, size = 1, prob = 0.5)
    G <- school_treat[school]
  }

  if (length(unique(G)) < 2) {
    stop("처치집단이 하나만 생성되었습니다.", call. = FALSE)
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
# Monte Carlo 1회 반복
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

  # 고정계수 방식은 논문 내 비교 목적으로만 계산한다.
  # 이 방식은 조정계수 추정의 불확실성을 반영하지 않으므로
  # 패키지 공개 함수 adid()에는 포함하지 않는다.
  X <- cbind(
    Intercept = 1,
    C = dat$C,
    G = dat$G
  )

  fit_post <- stats::lm.fit(x = X, y = dat$Y)
  fit_pre <- stats::lm.fit(x = X, y = dat$P)

  delta_hat <- fit_post$coefficients["C"] / fit_pre$coefficients["C"]

  if (!is.finite(delta_hat)) {
    stop("고정계수 방식의 delta_hat이 유한하지 않습니다.", call. = FALSE)
  }

  adjusted_change <- dat$Y - delta_hat * dat$P

  fixed_dat <- data.frame(
    adjusted_change = adjusted_change,
    G = dat$G
  )

  fixed_fit <- summary(stats::lm(adjusted_change ~ G, data = fixed_dat))
  se_fixed <- coef(fixed_fit)["G", "Std. Error"]

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
# 표준오차 방법별 요약
# ------------------------------------------------------------

summarize_method <- function(se, beta, emp_sd, true_tau) {
  valid <- is.finite(se) & se > 0 & is.finite(beta)

  if (sum(valid) == 0 || !is.finite(emp_sd) || emp_sd <= 0) {
    return(data.frame(
      se_ratio = NA_real_,
      rejection_rate = NA_real_,
      coverage_rate = NA_real_,
      n_valid = sum(valid)
    ))
  }

  se_valid <- se[valid]
  beta_valid <- beta[valid]

  # 검정은 H0: beta = 0 기준으로 수행한다.
  # tau = 0.0이면 rejection_rate는 제1종 오류율이다.
  # tau = 0.5이면 rejection_rate는 검정력이다.
  z_stat <- beta_valid / se_valid

  lower_bound <- beta_valid - 1.96 * se_valid
  upper_bound <- beta_valid + 1.96 * se_valid

  data.frame(
    se_ratio = sqrt(mean(se_valid^2, na.rm = TRUE)) / emp_sd,
    rejection_rate = mean(abs(z_stat) > 1.96, na.rm = TRUE),
    coverage_rate = mean(
      lower_bound <= true_tau & upper_bound >= true_tau,
      na.rm = TRUE
    ),
    n_valid = sum(valid)
  )
}


# ------------------------------------------------------------
# 하나의 조건에서 반복 모의실험 실행
# ------------------------------------------------------------

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

  reps <- reps[is.finite(reps$beta_hat), , drop = FALSE]

  if (nrow(reps) == 0) {
    stop("모든 모의실험 반복이 실패했습니다.", call. = FALSE)
  }

  emp_sd <- stats::sd(reps$beta_hat)

  dplyr::bind_rows(
    summarize_method(reps$se_fixed, reps$beta_hat, emp_sd, tau) |>
      dplyr::mutate(method = "fixed"),

    summarize_method(reps$se_boot_iid, reps$beta_hat, emp_sd, tau) |>
      dplyr::mutate(method = "bootstrap_iid"),

    summarize_method(reps$se_delta_iid, reps$beta_hat, emp_sd, tau) |>
      dplyr::mutate(method = "delta_iid"),

    summarize_method(reps$se_boot_cluster, reps$beta_hat, emp_sd, tau) |>
      dplyr::mutate(method = "bootstrap_cluster"),

    summarize_method(reps$se_delta_cluster, reps$beta_hat, emp_sd, tau) |>
      dplyr::mutate(method = "delta_cluster")
  ) |>
    dplyr::mutate(
      R_success = nrow(reps),
      R_requested = R,
      B_boot = B_boot,
      emp_sd = emp_sd,

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
      n,
      J,
      treat_assignment,
      beta_C,
      beta_P,
      beta_Y,
      common_trend_violation,
      tau,
      school_sd,
      method,
      se_ratio,
      rejection_rate,
      coverage_rate,
      rejection_type,
      emp_sd,
      n_valid,
      R_success,
      R_requested,
      B_boot
    )
}


# ------------------------------------------------------------
# chunk 결과 파일명 관리
# ------------------------------------------------------------

make_chunk_file <- function(tau, treat_assignment, run_label = "main") {
  file.path(
    CHUNK_DIR,
    paste0(
      "chunk_",
      safe_label(run_label),
      "_tau",
      num_label(tau),
      "_",
      treat_assignment,
      ".csv"
    )
  )
}


# ------------------------------------------------------------
# 하나의 chunk 실행: tau x 처치배정 방식
# ------------------------------------------------------------

run_chunk <- function(
    tau,
    treat_assignment,
    config = SIM_CONFIG,
    run_label = config$run_label %||% "main",
    overwrite = FALSE,
    workers = default_workers()
) {
  if (!treat_assignment %in% config$treat_assignment_values) {
    stop(
      "알 수 없는 처치배정 방식입니다: ",
      treat_assignment,
      call. = FALSE
    )
  }

  outfile <- make_chunk_file(
    tau = tau,
    treat_assignment = treat_assignment,
    run_label = run_label
  )

  if (file.exists(outfile) && !overwrite) {
    message("이미 존재하는 chunk를 건너뜁니다: ", outfile)
    return(read.csv(outfile, stringsAsFactors = FALSE))
  }

  workers <- set_future_plan(workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  grid <- tidyr::expand_grid(
    n = config$n_values,
    J = config$J,
    school_sd = config$school_sd_values
  )

  assignment_offset <- match(
    treat_assignment,
    config$treat_assignment_values
  )

  tau_offset <- as.integer(round(tau * 1000))

  out <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
    p <- grid[i, ]

    message(sprintf(
      "chunk [%s] 실행 중: tau=%.3f, assignment=%s, condition %d/%d, n=%d, J=%d, school_sd=%.3f",
      run_label,
      tau,
      treat_assignment,
      i,
      nrow(grid),
      p$n,
      p$J,
      p$school_sd
    ))

    condition_seed <- config$seed +
      100000L * assignment_offset +
      1000L * tau_offset +
      i

    run_condition(
      R = config$R,
      n = p$n,
      J = p$J,
      tau = tau,
      beta_P = config$beta_P,
      beta_Y = config$beta_Y,
      beta_C = config$beta_C,
      school_sd = p$school_sd,
      treat_strength = config$treat_strength,
      treat_assignment = treat_assignment,
      school_treat_sd = config$school_treat_sd,
      B_boot = config$B_boot,
      seed = condition_seed
    )
  })

  out <- out |>
    dplyr::mutate(
      run_label = run_label,
      .before = 1
    )

  write.csv(out, outfile, row.names = FALSE)

  message("chunk 저장 완료: ", outfile)

  out
}


# ------------------------------------------------------------
# 저장된 chunk 결과 모으기
# ------------------------------------------------------------

collect_chunks <- function(run_label = "main") {
  pattern <- paste0(
    "^chunk_",
    safe_label(run_label),
    "_.*[.]csv$"
  )

  files <- list.files(
    CHUNK_DIR,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop(
      "run_label = `",
      run_label,
      "`에 해당하는 chunk 파일을 찾을 수 없습니다. 검색 위치: ",
      CHUNK_DIR,
      call. = FALSE
    )
  }

  purrr::map_dfr(
    files,
    read.csv,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 논문 표 형식으로 변환
# ------------------------------------------------------------

make_rejection_table <- function(grid_res, tau_value, treat_assignment = NULL) {
  out <- grid_res[abs(grid_res$tau - tau_value) < 1e-12, , drop = FALSE]

  if (!is.null(treat_assignment)) {
    out <- out[out$treat_assignment == treat_assignment, , drop = FALSE]
  }

  if (nrow(out) == 0) {
    stop(
      "tau = ",
      tau_value,
      if (!is.null(treat_assignment)) paste0(", assignment = ", treat_assignment) else "",
      "에 해당하는 결과가 없습니다.",
      call. = FALSE
    )
  }

  out |>
    dplyr::group_by(treat_assignment, n, school_sd, method) |>
    dplyr::summarise(
      se_ratio = mean(se_ratio, na.rm = TRUE),
      rejection_rate = mean(rejection_rate, na.rm = TRUE),
      coverage_rate = mean(coverage_rate, na.rm = TRUE),
      n_valid = sum(n_valid, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      method_label = dplyr::recode(
        method,
        "fixed" = "고정계수",
        "bootstrap_iid" = "개인 부트스트랩",
        "delta_iid" = "델타방법",
        "bootstrap_cluster" = "군집 부트스트랩",
        "delta_cluster" = "군집 델타방법"
      ),
      method_label = factor(
        method_label,
        levels = c(
          "고정계수",
          "개인 부트스트랩",
          "델타방법",
          "군집 부트스트랩",
          "군집 델타방법"
        )
      ),
      treat_label = dplyr::recode(
        treat_assignment,
        "individual" = "개인 배정",
        "mixed" = "혼합 배정",
        "school" = "학교 배정"
      ),
      treat_label = factor(
        treat_label,
        levels = c("개인 배정", "혼합 배정", "학교 배정")
      ),
      cell = sprintf(
        "%.3f (%.3f) [%.3f]",
        se_ratio,
        rejection_rate,
        coverage_rate
      )
    ) |>
    dplyr::select(
      treat_label,
      n,
      school_sd,
      method_label,
      cell
    ) |>
    tidyr::pivot_wider(
      names_from = method_label,
      values_from = cell
    ) |>
    dplyr::arrange(
      treat_label,
      n,
      school_sd
    )
}
