# ============================================================
# simulation/build_tables.R
# ============================================================
# 저장된 chunk 결과를 모아 논문 표를 생성한다.
#
# 실행 예시
#   Rscript simulation/build_tables.R
#   Rscript simulation/build_tables.R label=pilot
#   Rscript simulation/build_tables.R label=main
#
# 입력
#   output/chunks/chunk_<label>_*.csv
#
# 출력
#   output/simulation_results_adid_<label>.csv
#   output/tables/table_<label>_tau0_all_assignments.csv
#   output/tables/table_<label>_tau0_individual.csv
#   output/tables/table_<label>_tau0_mixed.csv
#   output/tables/table_<label>_tau0_school.csv
#   output/tables/table_<label>_tau0p5_all_assignments.csv
#   output/tables/table_<label>_tau0p5_individual.csv
#   output/tables/table_<label>_tau0p5_mixed.csv
#   output/tables/table_<label>_tau0p5_school.csv
#
# 표 셀 형식
#   표준오차비율 (기각률) [포함률]
#
# tau = 0.0일 때 기각률은 제1종 오류율을 의미한다.
# tau = 0.5일 때 기각률은 검정력을 의미한다.
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

run_label <- args$label %||% SIM_CONFIG$run_label %||% "main"
run_label_safe <- safe_label(run_label)


# ------------------------------------------------------------
# chunk 결과 모으기
# ------------------------------------------------------------

message("============================================")
message("aDID 모의실험 표 생성")
message("============================================")
message("run_label = ", run_label)
message("============================================")

grid_res <- collect_chunks(run_label = run_label)

combined_file <- file.path(
  OUTPUT_DIR,
  paste0("simulation_results_adid_", run_label_safe, ".csv")
)

write.csv(
  grid_res,
  combined_file,
  row.names = FALSE
)

message("통합 결과 저장 완료: ", combined_file)


# ------------------------------------------------------------
# tau별 전체표와 처치배정 방식별 세부표 생성
# ------------------------------------------------------------

tau_values <- sort(unique(grid_res$tau))
assignments <- intersect(
  SIM_CONFIG$treat_assignment_values,
  unique(grid_res$treat_assignment)
)

for (tau_value in tau_values) {
  tau_safe <- num_label(tau_value)

  table_all <- make_rejection_table(
    grid_res = grid_res,
    tau_value = tau_value
  )

  outfile_all <- file.path(
    TABLE_DIR,
    paste0(
      "table_",
      run_label_safe,
      "_tau",
      tau_safe,
      "_all_assignments.csv"
    )
  )

  write.csv(
    table_all,
    outfile_all,
    row.names = FALSE
  )

  message("전체표 저장 완료: ", outfile_all)

  for (assignment in assignments) {
    has_rows <- nrow(
      grid_res[
        abs(grid_res$tau - tau_value) < 1e-12 &
          grid_res$treat_assignment == assignment,
        ,
        drop = FALSE
      ]
    ) > 0

    if (!has_rows) {
      next
    }

    tab <- make_rejection_table(
      grid_res = grid_res,
      tau_value = tau_value,
      treat_assignment = assignment
    )

    outfile <- file.path(
      TABLE_DIR,
      paste0(
        "table_",
        run_label_safe,
        "_tau",
        tau_safe,
        "_",
        assignment,
        ".csv"
      )
    )

    write.csv(
      tab,
      outfile,
      row.names = FALSE
    )

    message("세부표 저장 완료: ", outfile)
  }
}

message("============================================")
message("표 생성 완료")
message("저장 위치: ", TABLE_DIR)
message("============================================")
