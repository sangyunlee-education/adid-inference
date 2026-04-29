# ============================================================
# Basic Example for the adid Package
# ============================================================

source("../adid/R/adid.R")

set.seed(2026)

n <- 1000
J <- 50

school <- sample(seq_len(J), size = n, replace = TRUE)

school_effect <- rnorm(J, mean = 0, sd = 1)
S <- school_effect[school]

U_school <- rnorm(J, mean = 0, sd = 1)
U <- U_school[school] + rnorm(n, mean = 0, sd = 1)

G <- rbinom(n, size = 1, prob = plogis(0.6 * U))

C <- 1.0 * U + rnorm(n, sd = 1)

P <- 0.8 * U + 0.6 * S + rnorm(n, sd = 1)
Y <- 0.5 * G + 1.2 * U + 0.6 * S + rnorm(n, sd = 1)

dat <- data.frame(
  id = seq_len(n),
  school = school,
  G = G,
  C = C,
  P = P,
  Y = Y
)

fit_all <- adid(
  data = dat,
  pre = "P",
  post = "Y",
  treat = "G",
  compass = "C",
  cluster = "school",
  method = "all",
  B = 100,
  seed = 2026,
  apa = TRUE
)

print(fit_all)
print_apa(fit_all)
