# Independent mean calculation for finite-difference checks of effects.
share_means <- function(x, beta) {
  eta <- tcrossprod(x, beta)
  numerator <- exp(sweep(eta, 1L, apply(eta, 1L, max), "-"))
  numerator / rowSums(numerator)
}

simulated_margins <- function(N = 300L, J = 3L, seed = 1L) {
  set.seed(seed)
  x <- cbind("(Intercept)" = 1, income = rnorm(N), indicator = rbinom(N, 1, 0.5))
  beta <- rbind(reference = c(0, 0, 0),
                second = c(0.4, 0.8, -0.6),
                third = c(-0.3, -0.5, 1))[seq_len(J), , drop = FALSE]
  colnames(beta) <- colnames(x)
  p <- share_means(x, beta)
  y <- t(vapply(seq_len(N), function(i) {
    as.numeric(rmultinom(1L, size = 30L, prob = p[i, ])) / 30
  }, numeric(J)))
  colnames(y) <- colnames(p)
  fmlogit(y, x, method = "BFGS", print.level = 0L)
}

test_that("Marginal effects equal numerical derivatives of predicted shares", {
  object <- simulated_margins()
  step <- 1e-5

  # Check the mean design and a specified point, keeping other covariates fixed.
  for (at in list(NULL, c(1.25, 0.2))) {
    result <- effects.fmlogit(object, at = at)
    x <- matrix(if (is.null(at)) colMeans(object$x) else c(1, at), nrow = 1L)

    for (k in 2:ncol(x)) {
      x_low <- x_high <- x
      x_low[, k] <- x_low[, k] - step
      x_high[, k] <- x_high[, k] + step
      expected <- (share_means(x_high, object$coefficient) -
                     share_means(x_low, object$coefficient)) / (2 * step)
      expect_equal(as.numeric(result$effects[, k - 1L]),
                   as.numeric(expected), tolerance = 1e-7)
    }
  }
})

test_that("Discrete effects equal minimum-to-maximum changes in predicted shares", {
  object <- simulated_margins()
  at <- c(0.5, 0.25)

  for (marg.type in c("atmean", "aveacr")) {
    x <- if (marg.type == "atmean") matrix(c(1, at), nrow = 1L) else object$x
    result <- effects.fmlogit(object, effect = "discrete", marg.type = marg.type,
                             at = if (marg.type == "atmean") at else NULL)

    # For the binary indicator, this is the change from zero to one.
    for (k in 2:ncol(x)) {
      x_low <- x_high <- x
      x_low[, k] <- min(object$x[, k])
      x_high[, k] <- max(object$x[, k])
      expected <- share_means(x_high, object$coefficient) -
        share_means(x_low, object$coefficient)
      expect_equal(as.numeric(result$effects[, k - 1L]),
                   as.numeric(colMeans(expected)), tolerance = 1e-8)
    }

    # A covariate redistributes shares; it cannot change their total.
    expect_equal(as.numeric(colSums(result$effects)), c(0, 0), tolerance = 1e-8)
  }
})

test_that("Average marginal effects average individual expected-share derivatives", {
  object <- simulated_margins()
  result <- effects.fmlogit(object, marg.type = "aveacr")
  step <- 1e-5

  for (k in 2:ncol(object$x)) {
    x_low <- x_high <- object$x
    x_low[, k] <- x_low[, k] - step
    x_high[, k] <- x_high[, k] + step
    individual <- (share_means(x_high, object$coefficient) -
                     share_means(x_low, object$coefficient)) / (2 * step)

    # Check both the heterogeneous effects and their sample average.
    expect_equal(result$marg.list[[k - 1L]], individual, tolerance = 1e-7)
    expect_equal(as.numeric(result$effects[, k - 1L]),
                 as.numeric(colMeans(individual)), tolerance = 1e-7)
  }
  expect_equal(as.numeric(colSums(result$effects)), c(0, 0), tolerance = 1e-12)
})

test_that("Two-outcome margins reduce to the binary-logit formula", {
  object <- simulated_margins(J = 2L)
  p_other <- share_means(object$x, object$coefficient)[, 2L]
  slopes <- object$coefficient[2L, -1L]

  # dp_2/dx_k = p_2 (1 - p_2) beta_2k, and dp_1/dx_k = -dp_2/dx_k.
  expected <- colMeans(outer(p_other * (1 - p_other), slopes))
  result <- effects.fmlogit(object, marg.type = "aveacr")
  expect_equal(unname(result$effects), unname(rbind(-expected, expected)),
               tolerance = 1e-8)

  p_mean <- share_means(matrix(colMeans(object$x), nrow = 1L),
                        object$coefficient)[1L, 2L]
  expected <- p_mean * (1 - p_mean) * slopes
  result <- effects.fmlogit(object)
  expect_equal(unname(result$effects), unname(rbind(-expected, expected)),
               tolerance = 1e-8)
})

test_that("Baseline effect uncertainty includes covariance between outcome slopes", {
  object <- simulated_margins()
  p_at <- as.numeric(share_means(matrix(c(1, 0, 0), nrow = 1L),
                                 object$coefficient))
  slopes <- c("second:income", "third:income")
  sigma <- 0.01

  # At zero covariates, mean shares depend only on the fixed intercepts.
  # The baseline income effect is -p_1 (p_2 beta_2 + p_3 beta_3),
  # so its variance is known exactly when only the two slopes are uncertain.
  for (rho in c(0, 0.85)) {
    object$vcov_full[,] <- 0
    object$vcov_full[slopes, slopes] <- sigma^2 * matrix(c(1, rho, rho, 1), 2L)
    expected <- sigma * p_at[1L] * sqrt(p_at[2L]^2 + p_at[3L]^2 + 2 * rho * p_at[2L] * p_at[3L])

    set.seed(123)
    result <- effects.fmlogit(object, varlist = "income", at = c(0, 0),
                              se = TRUE, R = 2000L)
    # Allow Monte Carlo error, while checking independent and correlated slopes.
    expect_lt(abs(result$se[1L, 1L] / expected - 1), 0.1)
  }
})

test_that("Krinsky-Robb standard errors agree with the delta method", {
  object <- simulated_margins()
  # Small coefficient uncertainty makes the delta approximation accurate.
  object$vcov_full <- object$vcov_full / 100
  theta <- object$raw.estimates
  J <- ncol(object$y)
  step <- 1e-5

  for (effect in c("marginal", "discrete")) {
    for (marg.type in c("atmean", "aveacr")) {
      at <- if (marg.type == "atmean") c(1.25, 0.2) else NULL
      effect_at <- function(theta) {
        object$raw.estimates <- theta
        as.numeric(effects.fmlogit(object, effect = effect, marg.type = marg.type,
                                   varlist = "income", at = at)$effects)
      }

      # Numerically differentiate the requested statistic with respect to theta.
      gradient <- vapply(seq_along(theta), function(q) {
        lower <- upper <- theta
        lower[q] <- lower[q] - step
        upper[q] <- upper[q] + step
        (effect_at(upper) - effect_at(lower)) / (2 * step)
      }, numeric(J))

      # Delta-method variance: gradient %*% vcov(theta) %*% t(gradient).
      expected <- sqrt(diag(gradient %*% object$vcov_full %*% t(gradient)))
      set.seed(123)
      result <- effects.fmlogit(object, effect = effect, marg.type = marg.type,
                               varlist = "income", at = at, se = TRUE, R = 1000L)
      expect_lt(max(abs(as.numeric(result$se) / expected - 1)), 0.1)
    }
  }
})
