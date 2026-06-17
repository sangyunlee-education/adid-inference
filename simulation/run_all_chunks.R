# ============================================================
# simulation/run_all_chunks.R
# ============================================================
# 전체 모의실험 chunk를 순차적으로 실행한다.
#
# 내부적으로 run_chunk()를 반복 호출한다.
# 파일을 여러 개 만들 필요 없이, 이 스크립트 하나가 모든 tau x 처치배정 방식
# 조합을 차례로 실행한다.
#
# 실행 예시
#   Rscript simulation/run_all_chunks.R label=pilot R=50 B_boot=50
#   Rscript simulation/run_all_chunks.R label=main R=1000 B_boot=1000
#
# 일부 조건만 실행하고 싶은 경우
#   Rscript simulation/run_all_chunks.R label=pilot R=50 B_boot=50 tau_values=0 treat_assignment_values=individual,mixed
# ============================================================


# ------------------------------------------------------------
# 현재 스크립트 위치 확인 및 공통 함수 불러오기
# ------------------------------------------------------------

.script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  normalizePath(getwd(), mustWork = FALSE)
}

SIM_DIR <- .script_dir()

source(file.path(SIM_DIR, "sim_core.R"))


# ------------------------------------------------------------
# 명령행 인자 해석
# ------------------------------------------------------------

args <- parse_key_value_args()

config <- SIM_CONFIG

if (!is.null(args$label)) {
  config$run_label <- args$label
}

if (!is.null(args$R)) {
  config$R <- as.integer(args$R)
}

if (!is.null(args$B_boot)) {
  config$B_boot <- as.integer(args$B_boot)
}

if (!is.null(args$seed)) {
  config$seed <- as.integer(args$seed)
}

if (!is.null(args$n_values)) {
  config$n_values <- parse_numeric_vector(args$n_values)
}

if (!is.null(args$school_sd_values)) {
  config$school_sd_values <- parse_numeric_vector(args$school_sd_values)
}

if (!is.null(args$tau_values)) {
  config$tau_values <- parse_numeric_vector(args$tau_values)
}

if (!is.null(args$treat_assignment_values)) {
  config$treat_assignment_values <- parse_character_vector(args$treat_assignment_values)
}

if (!is.null(args$J)) {
  config$J <- as.integer(args$J)
}

if (!is.null(args$beta_P)) {
  config$beta_P <- as.numeric(args$beta_P)
}

if (!is.null(args$beta_Y)) {
  config$beta_Y <- as.numeric(args$beta_Y)
}

if (!is.null(args$beta_C)) {
  config$beta_C <- as.numeric(args$beta_C)
}

if (!is.null(args$treat_strength)) {
  config$treat_strength <- as.numeric(args$treat_strength)
}

if (!is.null(args$school_treat_sd)) {
  config$school_treat_sd <- as.numeric(args$school_treat_sd)
}

workers <- if (!is.null(args$workers)) {
  as.integer(args$workers)
} else {
  default_workers()
}

overwrite <- parse_logical_scalar(args$overwrite %||% "FALSE")
run_label <- config$run_label %||% "main"

chunk_grid <- tidyr::expand_grid(
  tau = config$tau_values,
  treat_assignment = config$treat_assignment_values
)


# ------------------------------------------------------------
# 전체 chunk 실행
# ------------------------------------------------------------

message("============================================")
message("aDID 전체 모의실험 chunk 실행")
message("============================================")
message("run_label              = ", run_label)
message("R                      = ", config$R)
message("B_boot                 = ", config$B_boot)
message("tau_values             = ", paste(config$tau_values, collapse = ", "))
message("treat_assignment_values= ", paste(config$treat_assignment_values, collapse = ", "))
message("n_values               = ", paste(config$n_values, collapse = ", "))
message("school_sd_values       = ", paste(config$school_sd_values, collapse = ", "))
message("workers                = ", workers)
message("overwrite              = ", overwrite)
message("chunk 수               = ", nrow(chunk_grid))
message("============================================")

for (i in seq_len(nrow(chunk_grid))) {
  tau_i <- chunk_grid$tau[i]
  assignment_i <- chunk_grid$treat_assignment[i]

  message("")
  message("--------------------------------------------")
  message("chunk ", i, "/", nrow(chunk_grid))
  message("tau = ", tau_i)
  message("assignment = ", assignment_i)
  message("--------------------------------------------")

  run_chunk(
    tau = tau_i,
    treat_assignment = assignment_i,
    config = config,
    run_label = run_label,
    overwrite = overwrite,
    workers = workers
  )
}

message("")
message("모든 chunk 실행 완료")
