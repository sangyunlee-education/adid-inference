# Monte Carlo 모의실험

이 폴더는 논문 표를 재현하기 위한 Monte Carlo 모의실험 코드를 담고 있습니다.

## 파일 구성

- `sim_core.R`: 자료 생성, 1회 반복, 조건별 반복 실행, 표 생성 함수
- `run_chunk.R`: 하나의 `tau x 처치배정 방식` chunk 실행
- `run_all_chunks.R`: 모든 chunk를 순차 실행
- `build_tables.R`: 저장된 chunk 결과를 모아 논문 표 생성

## 빠른 점검 실행

```bash
Rscript simulation/run_all_chunks.R label=pilot R=50 B_boot=50
Rscript simulation/build_tables.R label=pilot
```

## 최종 실행

```bash
Rscript simulation/run_all_chunks.R label=main R=1000 B_boot=1000
Rscript simulation/build_tables.R label=main
```

## 일부 조건만 실행

```bash
Rscript simulation/run_chunk.R tau=0 assignment=individual label=pilot R=50 B_boot=50
```

## 표 해석

표의 각 셀은 다음 형식입니다.

```text
표준오차비율 (기각률) [포함률]
```

- `tau = 0.0`: 기각률은 제1종 오류율입니다.
- `tau = 0.5`: 기각률은 검정력입니다.
