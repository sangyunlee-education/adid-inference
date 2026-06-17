# ============================================================
# examples/example_basic.R
# ============================================================
# aDID 함수 기본 사용 예시
# ============================================================

source("adid/R/adid.R")

set.seed(2026)

n <- 500
J <- 50

school <- sample(seq_len(J), size = n, replace = TRUE)
U <- rnorm(n)
G <- rbinom(n, size = 1, prob = stats::plogis(0.6 * U))
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

print(fit)
print_apa(fit)
