# ============================================================
# simulation/run_chunk.R
# ============================================================
# 하나의 모의실험 chunk를 실행한다.
#
# chunk 단위
#   tau x 처치배정 방식
#
# 저장 위치
#   output/chunks/
#
# 실행 예시
#   Rscript simulation/run_chunk.R tau=0 assignment=individual
#   Rscript simulation/run_chunk.R tau=0 assignment=mixed
#   Rscript simulation/run_chunk.R tau=0 assignment=school
#   Rscript simulation/run_chunk.R tau=0.5 assignment=individual
#   Rscript simulation/run_chunk.R tau=0.5 assignment=mixed
#   Rscript simulation/run_chunk.R tau=0.5 assignment=school
#
# 빠른 점검 실행
#   Rscript simulation/run_chunk.R tau=0 assignment=individual label=pilot R=50 B_boot=50
#
# 최종 실행
#   Rscript simulation/run_chunk.R tau=0 assignment=individual label=main R=1000 B_boot=1000
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

tau <- as.numeric(args$tau %||% "0")
assignment <- args$assignment %||% "individual"

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


# ------------------------------------------------------------
# 실행
# ------------------------------------------------------------

message("============================================")
message("aDID 모의실험 chunk 실행")
message("============================================")
message("run_label          = ", run_label)
message("tau                = ", tau)
message("assignment         = ", assignment)
message("R                  = ", config$R)
message("B_boot             = ", config$B_boot)
message("n_values           = ", paste(config$n_values, collapse = ", "))
message("J                  = ", config$J)
message("school_sd_values   = ", paste(config$school_sd_values, collapse = ", "))
message("workers            = ", workers)
message("overwrite          = ", overwrite)
message("============================================")

run_chunk(
  tau = tau,
  treat_assignment = assignment,
  config = config,
  run_label = run_label,
  overwrite = overwrite,
  workers = workers
)

message("chunk 실행 완료")
