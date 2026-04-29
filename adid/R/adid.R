# R/adid.R
# ============================================================
# Core functions for an R package implementing the adjusted DID
# estimator and its standard errors.
#
# Package name: adid
#
# Main user-facing function:
#   adid()
#
# Supported standard error methods:
#   - "delta"             : delta-method SE under individual-level independence
#   - "delta_cluster"     : cluster-robust delta-method SE
#   - "bootstrap"         : individual-level bootstrap SE
#   - "bootstrap_cluster" : cluster bootstrap SE
#   - "all"               : all of the above
#
# The fixed-coefficient method is intentionally excluded from the package API
# because it ignores uncertainty in the estimated adjustment coefficient.
# ============================================================

#' Influence Function for OLS Coefficients
#'
#' @param y Numeric outcome vector.
#' @param X Design matrix.
#' @return A list containing OLS coefficients, residuals, and influence functions.
#' @keywords internal
ols_if <- function(y, X) {
  fit <- lm.fit(x = X, y = y)
  
  beta_hat <- fit$coefficients
  resid <- as.numeric(y - X %*% beta_hat)
  
  Ainv <- solve(crossprod(X) / nrow(X))
  IF <- (X * resid) %*% t(Ainv)
  
  list(
    beta = beta_hat,
    resid = resid,
    IF = IF
  )
}

#' Influence Function for a Group Mean
#'
#' @param y Numeric outcome vector.
#' @param g Binary treatment/group indicator.
#' @param target Target group value, either 0 or 1.
#' @return A list containing the group mean and its influence function.
#' @keywords internal
group_mean_if <- function(y, g, target) {
  ind <- as.numeric(g == target)
  p <- mean(ind)
  
  if (p <= 0 || p >= 1) {
    stop("One of the treatment groups is empty.", call. = FALSE)
  }
  
  mu <- mean(y[ind == 1])
  IF <- ind / p * (y - mu)
  
  list(mu = mu, IF = IF)
}

#' Compute the aDID Estimator from the Basic Parameter Vector
#'
#' @param theta Named numeric vector: b_post, b_pre, mu_post_1, mu_post_0, mu_pre_1, mu_pre_0.
#' @return Numeric scalar, the adjusted DID estimate.
#' @keywords internal
adid_from_theta <- function(theta) {
  b_post <- theta["b_post"]
  b_pre <- theta["b_pre"]
  mu_post_1 <- theta["mu_post_1"]
  mu_post_0 <- theta["mu_post_0"]
  mu_pre_1 <- theta["mu_pre_1"]
  mu_pre_0 <- theta["mu_pre_0"]
  
  delta_hat <- b_post / b_pre
  
  (mu_post_1 - mu_post_0) - delta_hat * (mu_pre_1 - mu_pre_0)
}

#' Gradient of the aDID Estimator
#'
#' @param theta Named numeric vector: b_post, b_pre, mu_post_1, mu_post_0, mu_pre_1, mu_pre_0.
#' @return Numeric vector of partial derivatives.
#' @keywords internal
get_adid_gradient <- function(theta) {
  b_post <- theta["b_post"]
  b_pre <- theta["b_pre"]
  mu_pre_1 <- theta["mu_pre_1"]
  mu_pre_0 <- theta["mu_pre_0"]
  
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

#' Compute Basic Parameters and Their Variance-Covariance Matrix
#'
#' @param data Data frame.
#' @param pre Name of the pretest variable.
#' @param post Name of the posttest variable.
#' @param treat Name of the binary treatment indicator.
#' @param compass Name of the compass variable.
#' @param cluster Optional name of the cluster variable.
#' @return A list with theta_hat, vcov_theta, and influence functions.
#' @keywords internal
get_theta_vcov <- function(data, pre, post, treat, compass, cluster = NULL) {
  dat <- data
  
  y_pre <- dat[[pre]]
  y_post <- dat[[post]]
  g <- dat[[treat]]
  cvar <- dat[[compass]]
  
  if (!all(g %in% c(0, 1))) {
    stop("The treatment indicator must be coded as 0/1.", call. = FALSE)
  }
  
  cl <- if (!is.null(cluster)) dat[[cluster]] else NULL
  n <- nrow(dat)
  
  X <- cbind(Intercept = 1, compass = cvar, treat = g)
  
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

#' Delta-Method Standard Error for Adjusted DID
#'
#' @param data Data frame.
#' @param pre Name of the pretest variable.
#' @param post Name of the posttest variable.
#' @param treat Name of the binary treatment indicator.
#' @param compass Name of the compass variable.
#' @param cluster Optional name of the cluster variable.
#' @return A list containing estimate, SE, adjustment coefficient, theta, gradient, and vcov.
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
  
  delta_hat <- theta_hat["b_post"] / theta_hat["b_pre"]
  
  list(
    estimate = estimate,
    se = se,
    adjustment_coef = delta_hat,
    theta = theta_hat,
    gradient = grad,
    vcov_theta = tv$vcov_theta
  )
}

#' Bootstrap Standard Error for Adjusted DID
#'
#' @param data Data frame.
#' @param pre Name of the pretest variable.
#' @param post Name of the posttest variable.
#' @param treat Name of the binary treatment indicator.
#' @param compass Name of the compass variable.
#' @param cluster Optional name of the cluster variable. Required if cluster_boot = TRUE.
#' @param B Number of bootstrap replications.
#' @param cluster_boot Logical. If TRUE, perform cluster bootstrap.
#' @param seed Optional random seed.
#' @param max_attempts Maximum number of bootstrap attempts.
#' @return A list containing estimate, bootstrap SE, bootstrap estimates, and successful bootstrap count.
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
  if (!is.null(seed)) set.seed(seed)
  
  if (cluster_boot && is.null(cluster)) {
    stop("cluster must be supplied when cluster_boot = TRUE.", call. = FALSE)
  }
  
  n <- nrow(data)
  
  point <- adid_delta(
    data = data,
    pre = pre,
    post = post,
    treat = treat,
    compass = compass
  )
  
  draw_once <- function() {
    if (!cluster_boot) {
      idx <- sample.int(n, size = n, replace = TRUE)
      boot_dat <- data[idx, , drop = FALSE]
    } else {
      cl <- data[[cluster]]
      clusters <- unique(cl)
      sampled_clusters <- sample(clusters, size = length(clusters), replace = TRUE)
      
      boot_dat <- do.call(
        rbind,
        lapply(seq_along(sampled_clusters), function(k) {
          tmp <- data[cl == sampled_clusters[k], , drop = FALSE]
          tmp[[cluster]] <- k
          tmp
        })
      )
    }
    
    adid_delta(
      data = boot_dat,
      pre = pre,
      post = post,
      treat = treat,
      compass = compass
    )$estimate
  }
  
  boot_estimates <- numeric(0)
  attempts <- 0
  
  while (length(boot_estimates) < B && attempts < max_attempts) {
    attempts <- attempts + 1
    
    val <- tryCatch(
      draw_once(),
      error = function(e) NA_real_
    )
    
    if (is.finite(val)) {
      boot_estimates <- c(boot_estimates, val)
    }
  }
  
  list(
    estimate = point$estimate,
    se = ifelse(length(boot_estimates) >= 2, stats::sd(boot_estimates), NA_real_),
    adjustment_coef = point$adjustment_coef,
    boot_estimates = boot_estimates,
    B_success = length(boot_estimates),
    B_requested = B
  )
}

#' APA-Style p Value
#'
#' @param p Numeric p value.
#' @return Character string formatted in APA style.
#' @keywords internal
apa_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < .001, "p < .001", paste0("p = ", sub("^0", "", sprintf("%.3f", p))))
  )
}

#' APA-Style Number
#'
#' @param x Numeric value.
#' @param digits Number of digits after the decimal point.
#' @param omit_leading_zero Logical. If TRUE, omit the leading zero for values between -1 and 1.
#' @return Character string.
#' @keywords internal
apa_num <- function(x, digits = 3, omit_leading_zero = TRUE) {
  out <- sprintf(paste0("%.", digits, "f"), x)
  if (omit_leading_zero) {
    out <- sub("^-0[.]", "-.", out)
    out <- sub("^0[.]", ".", out)
  }
  out
}

#' APA-Style Result Sentence for Adjusted DID
#'
#' @param result A data frame returned by adid().
#' @param digits Number of digits after the decimal point.
#' @return Character vector of APA-style result sentences.
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
      "Using the ", method, " standard error, the adjusted DID estimate was ",
      "b = ", apa_num(estimate, digits = digits), ", ",
      "SE = ", apa_num(se, digits = digits), ", ",
      "z = ", apa_num(z, digits = digits), ", ",
      apa_p(p), ", ",
      "95% CI [", apa_num(ci_low, digits = digits), ", ",
      apa_num(ci_high, digits = digits), "]."
    )
  })
}

#' Adjusted Difference-in-Differences Estimator
#'
#' @description
#' Computes the adjusted difference-in-differences estimator and selected
#' standard errors. The fixed-coefficient method is not included because it
#' ignores uncertainty in the estimated adjustment coefficient.
#'
#' @param data Data frame.
#' @param pre Character. Name of the pretest variable.
#' @param post Character. Name of the posttest variable.
#' @param treat Character. Name of the binary treatment indicator coded 0/1.
#' @param compass Character. Name of the compass variable.
#' @param cluster Optional character. Name of the cluster variable.
#' @param method Character vector. One or more of "delta", "delta_cluster", "bootstrap", "bootstrap_cluster", or "all".
#' @param B Number of bootstrap replications for bootstrap methods.
#' @param seed Optional random seed.
#' @param conf_level Confidence level. Default is 0.95.
#' @param apa Logical. If TRUE, include APA-style result sentences as an attribute.
#'
#' @return A data frame with estimate, standard error, z statistic, p value,
#' confidence interval, and adjustment coefficient. If apa = TRUE, APA-style
#' result sentences are attached as the "apa" attribute.
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
                 method = c("delta", "delta_cluster", "bootstrap", "bootstrap_cluster"),
                 B = 1000,
                 seed = NULL,
                 conf_level = 0.95,
                 apa = FALSE) {
  method <- match.arg(
    method,
    choices = c("delta", "delta_cluster", "bootstrap", "bootstrap_cluster", "all"),
    several.ok = TRUE
  )
  
  if ("all" %in% method) {
    method <- c("delta", "delta_cluster", "bootstrap", "bootstrap_cluster")
  }
  
  if (any(c("delta_cluster", "bootstrap_cluster") %in% method) && is.null(cluster)) {
    stop("cluster must be supplied for cluster methods.", call. = FALSE)
  }
  
  required_vars <- c(pre, post, treat, compass, cluster)
  required_vars <- required_vars[!is.null(required_vars)]
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(
      "The following variables are missing from data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  
  alpha <- 1 - conf_level
  zcrit <- stats::qnorm(1 - alpha / 2)
  
  results <- lapply(method, function(m) {
    if (m == "delta") {
      fit <- adid_delta(data, pre, post, treat, compass)
    }
    
    if (m == "delta_cluster") {
      fit <- adid_delta(data, pre, post, treat, compass, cluster = cluster)
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

#' Print APA-Style Results
#'
#' @param result A data frame returned by adid().
#' @param digits Number of digits after the decimal point.
#' @return Invisibly returns the APA-style sentences.
#' @export
print_apa <- function(result, digits = 3) {
  sentences <- apa_adid(result, digits = digits)
  cat(paste(sentences, collapse = "\n"), "\n")
  invisible(sentences)
}
