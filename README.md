# 조정된 이중차분(aDID) 추론 코드

이 저장소는 조정된 이중차분(adjusted difference-in-differences; aDID) 추정량의 표준오차 계산과 Monte Carlo 모의실험 재현을 위한 R 코드입니다.


> 조정된 이중차분 추정량의 통계적 추론: 델타방법 표준오차와 다층자료로의 확장

## 저장소 구성

```text
adid-inference/
├── adid/
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   └── R/
│       └── adid.R
├── simulation/
│   ├── sim_core.R
│   ├── run_chunk.R
│   ├── run_all_chunks.R
│   ├── build_tables.R
│   └── README.md
├── examples/
│   └── example_basic.R
├── output/
│   ├── chunks/
│   └── tables/
├── README.md
├── LICENSE
└── .gitignore
```

## 핵심 기능

`adid/R/adid.R`는 다음 표준오차 방법을 지원합니다.

| 방법 | 설명 |
|---|---|
| `delta` | 개인 단위 독립성을 가정한 델타방법 표준오차 |
| `delta_cluster` | 군집 강건 델타방법 표준오차 |
| `bootstrap` | 개인 단위 부트스트랩 표준오차 |
| `bootstrap_cluster` | 군집 부트스트랩 표준오차 |
| `all` | 위 네 가지 표준오차를 모두 계산 |

고정계수 방식은 조정계수 추정의 불확실성을 반영하지 않기 때문에 공개 함수 `adid()`에는 포함하지 않았습니다. 단, 논문 내 비교를 위해 모의실험 코드에서는 고정계수 방식의 표준오차를 함께 계산합니다.

## 설치 및 준비

이 저장소는 연구 재현용 코드 저장소입니다. R 패키지로 설치하지 않고도 바로 실행할 수 있습니다.

```r
source("adid/R/adid.R")
```

모의실험 실행에는 다음 패키지가 필요합니다.

```r
install.packages(c("dplyr", "tidyr", "purrr", "future", "furrr"))
```

## 기본 사용 예시

```r
source("adid/R/adid.R")

set.seed(2026)
n <- 500
J <- 50
school <- sample(seq_len(J), size = n, replace = TRUE)
U <- rnorm(n)
G <- rbinom(n, size = 1, prob = plogis(0.6 * U))
C <- U + rnorm(n)
P <- 0.8 * U + rnorm(n)
Y <- 0.5 * G + 1.2 * U + rnorm(n)

dat <- data.frame(
  school = school,
  G = G,
  C = C,
  P = P,
  Y = Y
)

fit <- adid(
  data = dat,
  pre = "P",
  post = "Y",
  treat = "G",
  compass = "C",
  cluster = "school",
  method = c("delta", "delta_cluster"),
  apa = TRUE
)

fit
attr(fit, "apa")
```

## 모의실험 실행

전체 모의실험은 시간이 오래 걸리므로 먼저 작은 반복 수로 점검한 뒤 최종 반복 수로 실행하는 것을 권장합니다.

### 1. 빠른 점검 실행

```bash
Rscript simulation/run_all_chunks.R label=pilot R=50 B_boot=50
Rscript simulation/build_tables.R label=pilot
```

### 2. 최종 실행

```bash
Rscript simulation/run_all_chunks.R label=main R=1000 B_boot=1000
Rscript simulation/build_tables.R label=main
```

### 3. 특정 chunk만 실행

```bash
Rscript simulation/run_chunk.R tau=0 assignment=individual label=main R=1000 B_boot=1000
Rscript simulation/run_chunk.R tau=0 assignment=mixed label=main R=1000 B_boot=1000
Rscript simulation/run_chunk.R tau=0 assignment=school label=main R=1000 B_boot=1000
Rscript simulation/run_chunk.R tau=0.5 assignment=individual label=main R=1000 B_boot=1000
Rscript simulation/run_chunk.R tau=0.5 assignment=mixed label=main R=1000 B_boot=1000
Rscript simulation/run_chunk.R tau=0.5 assignment=school label=main R=1000 B_boot=1000
```

## 출력 파일

모의실험 결과는 다음 위치에 생성됩니다.

```text
output/chunks/   # tau x 처치배정 방식별 중간 결과
output/tables/   # 논문 표 형태로 정리한 결과
```

대표 출력 파일은 다음과 같습니다.

```text
output/simulation_results_adid_main.csv
output/tables/table_main_tau0_all_assignments.csv
output/tables/table_main_tau0_individual.csv
output/tables/table_main_tau0_mixed.csv
output/tables/table_main_tau0_school.csv
output/tables/table_main_tau0p5_all_assignments.csv
output/tables/table_main_tau0p5_individual.csv
output/tables/table_main_tau0p5_mixed.csv
output/tables/table_main_tau0p5_school.csv
```

표의 각 셀은 다음 형식입니다.

```text
표준오차비율 (기각률) [포함률]
```

- `tau = 0.0`: 기각률은 제1종 오류율입니다.
- `tau = 0.5`: 기각률은 검정력입니다.

## 재현성 메모

- `run_chunk.R`와 `run_all_chunks.R`는 `label` 인자를 받아 결과 파일명을 구분합니다.
- 이미 생성된 chunk 파일은 기본적으로 다시 계산하지 않습니다.
- 같은 조건을 다시 계산하려면 `overwrite=TRUE`를 지정하십시오.

```bash
Rscript simulation/run_chunk.R tau=0 assignment=individual label=main overwrite=TRUE
```

## 인용

```text
이상윤, 김용남. (2026). 조정된 이중차분 추정량의 통계적 추론: 델타방법 표준오차와 다층자료로의 확장. 교육평가연구, 39(2),
```

## 라이선스

이 저장소의 코드는 MIT 라이선스를 따릅니다.
