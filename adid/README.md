# adid R 함수

이 폴더는 조정된 이중차분(adjusted DID; aDID) 추정량과 표준오차 계산 함수를 담고 있습니다.

## 주요 함수

- `adid()`: aDID 추정량과 표준오차 계산
- `apa_adid()`: `adid()` 결과를 APA 형식 문장으로 변환
- `print_apa()`: APA 형식 문장 출력

## 지원 표준오차

- `delta`: 개인 단위 독립성을 가정한 델타방법 표준오차
- `delta_cluster`: 군집 강건 델타방법 표준오차
- `bootstrap`: 개인 단위 부트스트랩 표준오차
- `bootstrap_cluster`: 군집 부트스트랩 표준오차
- `all`: 네 가지 표준오차를 모두 계산

## 사용 예시

```r
source("adid/R/adid.R")

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
