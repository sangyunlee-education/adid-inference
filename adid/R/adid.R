# ============================================================
# adid/R/adid.R
# ============================================================
# 조정된 이중차분(adjusted DID; aDID) 추정량과 표준오차 계산 함수
#
# 주요 공개 함수
#   - adid()      : aDID 추정량과 표준오차 계산
#   - apa_adid()  : 결과를 APA 형식 문장으로 변환
#   - print_apa() : APA 형식 문장 출력
#
# 지원하는 표준오차 방법
#   - "delta"             : 개인 단위 독립성을 가정한 델타방법 표준오차
#   - "delta_cluster"     : 군집 강건 델타방법 표준오차
#   - "bootstrap"         : 개인 단위 부트스트랩 표준오차
#   - "bootstrap_cluster" : 군집 부트스트랩 표준오차
#   - "all"               : 위 네 가지 방법을 모두 계산
#
# 고정계수 방식은 패키지 공개 API에서 제외하였다.
# 해당 방식은 조정계수 추정의 불확실성을 표준오차에 반영하지 않기 때문이다.
# ============================================================


# ------------------------------------------------------------
# 내부 점검 함수
# ------------------------------------------------------------

.check_required_vars <- function(data, vars) {
  vars <- unlist(vars)
  vars <- vars[!is.na(vars) & nzchar(vars)]
  missing_vars <- setdiff(vars, names(data))

  if (length(missing_vars) > 0) {
    stop(
      "자료에 다음 변수가 없습니다: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.check_numeric_vector <- function(x, name) {
  if (!is.numeric(x) && !is.integer(x)) {
    stop(
      "변수 `", name, "`는 숫자형이어야 합니다.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.check_binary_treat <- function(g) {
  if (is.factor(g)) {
    g <- as.character(g)
  }

  if (is.character(g)) {
    if (!all(g %in% c("0", "1"))) {
      stop("처치 변수는 0/1로 코딩되어야 합니다.", call. = FALSE)
    }
    g <- as.numeric(g)
  }

  if (is.logical(g)) {
    g <- as.integer(g)
  }

  if (!is.numeric(g) && !is.integer(g)) {
    stop("처치 변수는 0/1로 코딩되어야 합니다.", call. = FALSE)
  }

  if (!all(g %in% c(0, 1))) {
    stop("처치 변수는 0/1로 코딩되어야 합니다.", call. = FALSE)
  }

  as.numeric(g)
}

.check_no_missing <- function(data, vars) {
  cc <- stats::complete.cases(data[, vars, drop = FALSE])

  if (!all(cc)) {
    stop(
      "현재 함수는 결측값을 자동 처리하지 않습니다. ",
      "adid()를 실행하기 전에 결측값을 제거하거나 대체하십시오.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.check_b_pre <- function(b_pre) {
  if (!is.finite(b_pre) || abs(b_pre) < sqrt(.Machine$double.eps)) {
    stop(
      "사전검사 변수에 대한 컴패스 계수가 0이거나 0에 매우 가깝습니다.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------
# 영향함수 계산 보조 함수
# ------------------------------------------------------------

#' OLS 계수의 영향함수 계산
#'
#' @param y 숫자형 종속변수 벡터.
#' @param X 설계행렬.
#' @return OLS 계수, 잔차, 영향함수를 포함한 리스트.
#' @keywords internal
ols_if <- function(y, X) {
  fit <- stats::lm.fit(x = X, y = y)

  beta_hat <- fit$coefficients

  if (any(!is.finite(beta_hat))) {
    stop("OLS 계수 추정에 실패했습니다. 설계행렬의 선형종속성을 확인하십시오.", call. = FALSE)
  }

  resid <- as.numeric(y - X %*% beta_hat)
  A <- crossprod(X) / nrow(X)

  Ainv <- tryCatch(
    solve(A),
    error = function(e) {
      stop("OLS 영향함수 계산 중 행렬 역산에 실패했습니다.", call. = FALSE)
    }
  )

  IF <- sweep(X, 1, resid, "*") %*% t(Ainv)
  colnames(IF) <- colnames(X)

  list(
    beta = beta_hat,
    resid = resid,
    IF = IF
  )
}


#' 집단 평균의 영향함수 계산
#'
#' @param y 숫자형 결과변수 벡터.
#' @param g 0/1 처치집단 지시변수.
#' @param target 평균을 계산할 집단값. 0 또는 1.
#' @return 집단 평균과 영향함수를 포함한 리스트.
#' @keywords internal
group_mean_if <- function(y, g, target) {
  ind <- as.numeric(g == target)
  p <- mean(ind)

  if (p <= 0 || p >= 1) {
    stop("처치집단 또는 비교집단 중 하나가 비어 있습니다.", call. = FALSE)
  }

  mu <- mean(y[ind == 1])
  IF <- ind / p * (y - mu)

  list(
    mu = mu,
    IF = IF
  )
}


# ------------------------------------------------------------
# aDID 추정량과 기울기
# ------------------------------------------------------------

#' 기본 모수 벡터에서 aDID 추정량 계산
#'
#' @param theta 이름이 있는 숫자형 벡터. b_post, b_pre, mu_post_1, mu_post_0, mu_pre_1, mu_pre_0를 포함해야 함.
#' @return 조정된 이중차분 추정량.
#' @keywords internal
adid_from_theta <- function(theta) {
  b_post <- unname(theta["b_post"])
  b_pre <- unname(theta["b_pre"])
  mu_post_1 <- unname(theta["mu_post_1"])
  mu_post_0 <- unname(theta["mu_post_0"])
  mu_pre_1 <- unname(theta["mu_pre_1"])
  mu_pre_0 <- unname(theta["mu_pre_0"])

  .check_b_pre(b_pre)

  delta_hat <- b_post / b_pre

  unname(
    (mu_post_1 - mu_post_0) -
      delta_hat * (mu_pre_1 - mu_pre_0)
  )
}


#' aDID 추정량의 델타방법 기울기 계산
#'
#' @param theta 이름이 있는 숫자형 벡터. b_post, b_pre, mu_post_1, mu_post_0, mu_pre_1, mu_pre_0를 포함해야 함.
#' @return 각 기본 모수에 대한 편미분 벡터.
#' @keywords internal
get_adid_gradient <- function(theta) {
  b_post <- unname(theta["b_post"])
  b_pre <- unname(theta["b_pre"])
  mu_pre_1 <- unname(theta["mu_pre_1"])
  mu_pre_0 <- unname(theta["mu_pre_0"])

  .check_b_pre(b_pre)

  d_pre <- mu_pre_1 - mu_pre_0

  c(
    b_post = -d_pre / b_pre,
    b_pre = b_post * d_pre / (b_pre^2),
    mu_post_1 = 1,
    mu_post_0 = -1,
    mu_pre_1 = -b_post / b_pre,
    mu_pre_0 = b_post / b_pre
  )
}


# ------------------------------------------------------------
# 빠른 점추정 함수
# ------------------------------------------------------------

#' aDID 점추정량만 빠르게 계산
#'
#' 이 함수는 표준오차를 계산하지 않고 점추정량만 계산한다.
#' 부트스트랩 반복 안에서 영향함수와 분산-공분산 행렬을 매번 계산하지 않기 위해 사용한다.
#'
#' @param data 자료프레임.
#' @param pre 사전검사 변수명.
#' @param post 사후검사 변수명.
#' @param treat 0/1 처치 변수명.
#' @param compass 컴패스 변수명.
#' @return 점추정량, 조정계수, 기본 모수 벡터를 포함한 리스트.
#' @keywords internal
adid_point <- function(data, pre, post, treat, compass) {
  required_vars <- c(pre, post, treat, compass)

  .check_required_vars(data, required_vars)
  .check_no_missing(data, required_vars)
  .check_numeric_vector(data[[pre]], pre)
  .check_numeric_vector(data[[post]], post)
  .check_numeric_vector(data[[compass]], compass)

  y_pre <- as.numeric(data[[pre]])
  y_post <- as.numeric(data[[post]])
  g <- .check_binary_treat(data[[treat]])
  cvar <- as.numeric(data[[compass]])

  X <- cbind(
    Intercept = 1,
    compass = cvar,
    treat = g
  )

  fit_post <- stats::lm.fit(x = X, y = y_post)
  fit_pre <- stats::lm.fit(x = X, y = y_pre)

  theta_hat <- c(
    b_post = unname(fit_post$coefficients["compass"]),
    b_pre = unname(fit_pre$coefficients["compass"]),
    mu_post_1 = mean(y_post[g == 1]),
    mu_post_0 = mean(y_post[g == 0]),
    mu_pre_1 = mean(y_pre[g == 1]),
    mu_pre_0 = mean(y_pre[g == 0])
  )

  estimate <- adid_from_theta(theta_hat)
  adjustment_coef <- unname(theta_hat["b_post"] / theta_hat["b_pre"])

  list(
    estimate = estimate,
    adjustment_coef = adjustment_coef,
    theta = theta_hat
  )
}


# ------------------------------------------------------------
# 기본 모수와 분산-공분산 행렬
# ------------------------------------------------------------

#' 기본 모수와 분산-공분산 행렬 계산
#'
#' @param data 자료프레임.
#' @param pre 사전검사 변수명.
#' @param post 사후검사 변수명.
#' @param treat 0/1 처치 변수명.
#' @param compass 컴패스 변수명.
#' @param cluster 선택 사항. 군집 변수명.
#' @return 기본 모수 추정치, 분산-공분산 행렬, 영향함수를 포함한 리스트.
#' @keywords internal
get_theta_vcov <- function(data, pre, post, treat, compass, cluster = NULL) {
  required_vars <- c(pre, post, treat, compass)

  if (!is.null(cluster)) {
    required_vars <- c(required_vars, cluster)
  }

  .check_required_vars(data, required_vars)
  .check_no_missing(data, required_vars)
  .check_numeric_vector(data[[pre]], pre)
  .check_numeric_vector(data[[post]], post)
  .check_numeric_vector(data[[compass]], compass)

  y_pre <- as.numeric(data[[pre]])
  y_post <- as.numeric(data[[post]])
  g <- .check_binary_treat(data[[treat]])
  cvar <- as.numeric(data[[compass]])

  cl <- if (!is.null(cluster)) data[[cluster]] else NULL
  n <- nrow(data)

  if (!is.null(cl) && length(unique(cl)) < 2) {
    stop("군집 강건 추론을 위해서는 최소 두 개 이상의 군집이 필요합니다.", call. = FALSE)
  }

  X <- cbind(
    Intercept = 1,
    compass = cvar,
    treat = g
  )

  fit_post <- ols_if(y = y_post, X = X)
  fit_pre <- ols_if(y = y_pre, X = X)

  mu_post_1 <- group_mean_if(y_post, g, target = 1)
  mu_post_0 <- group_mean_if(y_post, g, target = 0)
  mu_pre_1 <- group_mean_if(y_pre, g, target = 1)
  mu_pre_0 <- group_mean_if(y_pre, g, target = 0)

  theta_hat <- c(
    b_post = unname(fit_post$beta["compass"]),
    b_pre = unname(fit_pre$beta["compass"]),
    mu_post_1 = unname(mu_post_1$mu),
    mu_post_0 = unname(mu_post_0$mu),
    mu_pre_1 = unname(mu_pre_1$mu),
    mu_pre_0 = unname(mu_pre_0$mu)
  )

  IF <- cbind(
    b_post = fit_post$IF[, "compass"],
    b_pre = fit_pre$IF[, "compass"],
    mu_post_1 = mu_post_1$IF,
    mu_post_0 = mu_post_0$IF,
    mu_pre_1 = mu_pre_1$IF,
    mu_pre_0 = mu_pre_0$IF
  )

  if (is.null(cl)) {
    vcov_theta <- crossprod(IF) / (n^2)
  } else {
    IF_cluster <- rowsum(IF, group = cl, reorder = FALSE)
    vcov_theta <- crossprod(as.matrix(IF_cluster)) / (n^2)
  }

  list(
    theta_hat = theta_hat,
    vcov_theta = vcov_theta,
    IF = IF
  )
}


# ------------------------------------------------------------
# 델타방법 표준오차
# ------------------------------------------------------------

#' aDID 델타방법 표준오차 계산
#'
#' @param data 자료프레임.
#' @param pre 사전검사 변수명.
#' @param post 사후검사 변수명.
#' @param treat 0/1 처치 변수명.
#' @param compass 컴패스 변수명.
#' @param cluster 선택 사항. 군집 변수명.
#' @return 추정량, 표준오차, 조정계수, 기본 모수, 기울기, 분산-공분산 행렬을 포함한 리스트.
#' @keywords internal
adid_delta <- function(data, pre, post, treat, compass, cluster = NULL) {
  tv <- get_theta_vcov(
    data = data,
    pre = pre,
    post = post,
    treat = treat,
    compass = compass,
    cluster = cluster
  )

  theta_hat <- tv$theta_hat
  grad <- get_adid_gradient(theta_hat)

  estimate <- adid_from_theta(theta_hat)
  variance <- as.numeric(t(grad) %*% tv$vcov_theta %*% grad)
  se <- sqrt(max(variance, 0))

  delta_hat <- unname(theta_hat["b_post"] / theta_hat["b_pre"])

  list(
    estimate = estimate,
    se = se,
    adjustment_coef = delta_hat,
    theta = theta_hat,
    gradient = grad,
    vcov_theta = tv$vcov_theta
  )
}


# ------------------------------------------------------------
# 부트스트랩 표준오차
# ------------------------------------------------------------

#' aDID 부트스트랩 표준오차 계산
#'
#' @param data 자료프레임.
#' @param pre 사전검사 변수명.
#' @param post 사후검사 변수명.
#' @param treat 0/1 처치 변수명.
#' @param compass 컴패스 변수명.
#' @param cluster 선택 사항. 군집 변수명. cluster_boot = TRUE이면 반드시 필요함.
#' @param B 부트스트랩 반복 횟수.
#' @param cluster_boot TRUE이면 군집 부트스트랩을 수행함.
#' @param seed 선택 사항. 난수 시드.
#' @param max_attempts 부트스트랩 재시도 최대 횟수.
#' @return 추정량, 부트스트랩 표준오차, 부트스트랩 추정값, 성공한 반복 횟수를 포함한 리스트.
#' @keywords internal
adid_bootstrap <- function(data,
                           pre,
                           post,
                           treat,
                           compass,
                           cluster = NULL,
                           B = 1000,
                           cluster_boot = FALSE,
                           seed = NULL,
                           max_attempts = B * 10) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    stop("B는 양의 정수여야 합니다.", call. = FALSE)
  }

  B <- as.integer(B)
  max_attempts <- as.integer(max_attempts)

  if (cluster_boot && is.null(cluster)) {
    stop("cluster_boot = TRUE이면 cluster를 지정해야 합니다.", call. = FALSE)
  }

  n <- nrow(data)

  point <- adid_point(
    data = data,
    pre = pre,
    post = post,
    treat = treat,
    compass = compass
  )

  if (cluster_boot) {
    cl <- data[[cluster]]
    idx_list <- split(seq_len(n), cl, drop = TRUE)
    clusters <- names(idx_list)

    if (length(clusters) < 2) {
      stop("군집 부트스트랩을 위해서는 최소 두 개 이상의 군집이 필요합니다.", call. = FALSE)
    }
  }

  draw_once <- function() {
    if (!cluster_boot) {
      idx <- sample.int(n, size = n, replace = TRUE)
      boot_dat <- data[idx, , drop = FALSE]
    } else {
      sampled_clusters <- sample(
        clusters,
        size = length(clusters),
        replace = TRUE
      )

      boot_idx <- unlist(
        idx_list[sampled_clusters],
        use.names = FALSE
      )

      boot_dat <- data[boot_idx, , drop = FALSE]
    }

    adid_point(
      data = boot_dat,
      pre = pre,
      post = post,
      treat = treat,
      compass = compass
    )$estimate
  }

  boot_estimates <- numeric(B)
  success_count <- 0L
  attempts <- 0L

  while (success_count < B && attempts < max_attempts) {
    attempts <- attempts + 1L

    val <- tryCatch(
      draw_once(),
      error = function(e) NA_real_
    )

    if (is.finite(val)) {
      success_count <- success_count + 1L
      boot_estimates[success_count] <- val
    }
  }

  boot_estimates <- boot_estimates[seq_len(success_count)]

  if (success_count < B) {
    warning(
      "요청한 B = ", B, "회 중 ", success_count,
      "회의 부트스트랩 반복만 성공했습니다.",
      call. = FALSE
    )
  }

  list(
    estimate = point$estimate,
    se = if (length(boot_estimates) >= 2) stats::sd(boot_estimates) else NA_real_,
    adjustment_coef = point$adjustment_coef,
    boot_estimates = boot_estimates,
    B_success = length(boot_estimates),
    B_requested = B
  )
}


# ------------------------------------------------------------
# APA 형식 출력 보조 함수
# ------------------------------------------------------------

#' APA 형식 p값 표기
#'
#' @param p p값.
#' @return APA 형식 p값 문자열.
#' @keywords internal
apa_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < .001,
      "p < .001",
      paste0("p = ", sub("^0", "", sprintf("%.3f", p)))
    )
  )
}


#' APA 형식 숫자 표기
#'
#' @param x 숫자값.
#' @param digits 소수점 이하 자릿수.
#' @param omit_leading_zero TRUE이면 절댓값이 1보다 작은 값에서 앞의 0을 생략함.
#' @return 숫자 문자열.
#' @keywords internal
apa_num <- function(x, digits = 3, omit_leading_zero = TRUE) {
  out <- sprintf(paste0("%.", digits, "f"), x)

  if (omit_leading_zero) {
    out <- sub("^-0[.]", "-.", out)
    out <- sub("^0[.]", ".", out)
  }

  out
}


#' aDID 결과를 APA 형식 문장으로 변환
#'
#' @param result adid()가 반환한 자료프레임.
#' @param digits 소수점 이하 자릿수.
#' @return APA 형식 결과 문장 벡터.
#' @export
apa_adid <- function(result, digits = 3) {
  apply(result, 1, function(row) {
    method <- row[["method"]]
    estimate <- as.numeric(row[["estimate"]])
    se <- as.numeric(row[["se"]])
    z <- as.numeric(row[["z"]])
    p <- as.numeric(row[["p_value"]])
    ci_low <- as.numeric(row[["conf_low"]])
    ci_high <- as.numeric(row[["conf_high"]])

    paste0(
      method, " 표준오차를 사용했을 때, 조정된 이중차분 추정량은 ",
      "b = ", apa_num(estimate, digits = digits), ", ",
      "SE = ", apa_num(se, digits = digits), ", ",
      "z = ", apa_num(z, digits = digits), ", ",
      apa_p(p), ", ",
      "95% CI [", apa_num(ci_low, digits = digits), ", ",
      apa_num(ci_high, digits = digits), "]였다."
    )
  })
}


# ------------------------------------------------------------
# 공개 함수: aDID 추정
# ------------------------------------------------------------

#' 조정된 이중차분 추정량
#'
#' @description
#' 조정된 이중차분(adjusted difference-in-differences; aDID) 추정량과
#' 선택한 표준오차를 계산한다. 고정계수 방식은 조정계수 추정의 불확실성을
#' 반영하지 않기 때문에 공개 함수에는 포함하지 않았다.
#'
#' @param data 자료프레임.
#' @param pre 사전검사 변수명.
#' @param post 사후검사 변수명.
#' @param treat 0/1로 코딩된 처치 변수명.
#' @param compass 컴패스 변수명.
#' @param cluster 선택 사항. 군집 변수명.
#' @param method 사용할 표준오차 방법. "delta", "delta_cluster", "bootstrap", "bootstrap_cluster", "all" 중 하나 이상.
#' @param B 부트스트랩 반복 횟수.
#' @param seed 선택 사항. 난수 시드.
#' @param conf_level 신뢰수준. 기본값은 0.95.
#' @param apa TRUE이면 APA 형식 결과 문장을 속성으로 함께 저장함.
#'
#' @return 추정량, 표준오차, z 통계량, p값, 신뢰구간, 조정계수를 담은 자료프레임.
#'
#' @examples
#' # fit <- adid(
#' #   data = dat,
#' #   pre = "P",
#' #   post = "Y",
#' #   treat = "G",
#' #   compass = "C",
#' #   cluster = "school",
#' #   method = c("delta", "delta_cluster"),
#' #   apa = TRUE
#' # )
#' # attr(fit, "apa")
#'
#' @export
adid <- function(data,
                 pre,
                 post,
                 treat,
                 compass,
                 cluster = NULL,
                 method = "delta",
                 B = 1000,
                 seed = NULL,
                 conf_level = 0.95,
                 apa = FALSE) {
  if (!is.data.frame(data)) {
    stop("data는 자료프레임이어야 합니다.", call. = FALSE)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      conf_level <= 0 || conf_level >= 1) {
    stop("conf_level은 0과 1 사이의 숫자여야 합니다.", call. = FALSE)
  }

  method <- match.arg(
    method,
    choices = c("delta", "delta_cluster", "bootstrap", "bootstrap_cluster", "all"),
    several.ok = TRUE
  )

  if ("all" %in% method) {
    method <- c("delta", "delta_cluster", "bootstrap", "bootstrap_cluster")
  }

  method <- unique(method)

  if (any(c("delta_cluster", "bootstrap_cluster") %in% method) && is.null(cluster)) {
    stop("군집 표준오차 방법을 사용하려면 cluster를 지정해야 합니다.", call. = FALSE)
  }

  required_vars <- c(pre, post, treat, compass)

  if (!is.null(cluster)) {
    required_vars <- c(required_vars, cluster)
  }

  .check_required_vars(data, required_vars)
  .check_no_missing(data, required_vars)

  alpha <- 1 - conf_level
  zcrit <- stats::qnorm(1 - alpha / 2)

  results <- lapply(method, function(m) {
    if (m == "delta") {
      fit <- adid_delta(
        data = data,
        pre = pre,
        post = post,
        treat = treat,
        compass = compass
      )
    }

    if (m == "delta_cluster") {
      fit <- adid_delta(
        data = data,
        pre = pre,
        post = post,
        treat = treat,
        compass = compass,
        cluster = cluster
      )
    }

    if (m == "bootstrap") {
      fit <- adid_bootstrap(
        data = data,
        pre = pre,
        post = post,
        treat = treat,
        compass = compass,
        B = B,
        cluster_boot = FALSE,
        seed = seed
      )
    }

    if (m == "bootstrap_cluster") {
      fit <- adid_bootstrap(
        data = data,
        pre = pre,
        post = post,
        treat = treat,
        compass = compass,
        cluster = cluster,
        B = B,
        cluster_boot = TRUE,
        seed = seed
      )
    }

    z <- fit$estimate / fit$se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

    data.frame(
      method = m,
      estimate = fit$estimate,
      se = fit$se,
      z = z,
      p_value = p,
      conf_low = fit$estimate - zcrit * fit$se,
      conf_high = fit$estimate + zcrit * fit$se,
      adjustment_coef = fit$adjustment_coef,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, results)
  rownames(out) <- NULL

  if (apa) {
    attr(out, "apa") <- apa_adid(out)
  }

  out
}


#' APA 형식 결과 출력
#'
#' @param result adid()가 반환한 자료프레임.
#' @param digits 소수점 이하 자릿수.
#' @return APA 형식 문장을 보이지 않게 반환함.
#' @export
print_apa <- function(result, digits = 3) {
  sentences <- apa_adid(result, digits = digits)
  cat(paste(sentences, collapse = "\n"), "\n")
  invisible(sentences)
}
