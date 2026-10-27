# Calculate mean shares independently to check effects by finite differences.
shareMeans <- function(x, beta) {
  eta <- tcrossprod(x, beta)
  numerator <- exp(sweep(eta, 1L, apply(eta, 1L, max), "-"))
  numerator/rowSums(numerator)
}

simulatedMargins <- function(N = 300L, J = 3L, seed = 1L) {
  set.seed(seed)
  x <- cbind("(Intercept)" = 1, income = rnorm(N), indicator = rbinom(N, 1, 0.5))
  beta <- rbind(reference = c(0, 0, 0),
                second = c(0.4, 0.8, -0.6),
                third = c(-0.3, -0.5, 1))[seq_len(J), , drop = FALSE]
  colnames(beta) <- colnames(x)
  p <- shareMeans(x, beta)
  y <- t(vapply(seq_len(N), function(i) {
    as.numeric(rmultinom(1L, size = 25L, prob = p[i, ]))/25
  }, numeric(J)))
  colnames(y) <- colnames(p)
  fmlogit(y, x, method = "BHHH", print.level = 0L)
}

test_that("Marginal effects equal numerical derivatives of predicted shares", {
  object <- simulatedMargins(N = 300L, seed = 1)
  step <- 1e-5

  # Check the mean design and a specified point, keeping other covariates fixed.
  for (at in list(NULL, c(1.25, 0.2))) {
    result <- effects.fmlogit(object, at = at)
    x <- matrix(if (is.null(at)) colMeans(object$x) else c(1, at), nrow = 1L)

    for (k in 2:ncol(x)) {
      x.low <- x.high <- x
      x.low[, k] <- x.low[, k] - step
      x.high[, k] <- x.high[, k] + step
      expected <- (shareMeans(x.high, object$coefficient) -
                   shareMeans(x.low, object$coefficient)) / (2*step)
      expect_equal(as.numeric(result$effects[, k-1L]),
                   as.numeric(expected), tolerance = 1e-7)
    }
  }
})

test_that("Discrete effects equal minimum-to-maximum changes in predicted shares", {
  object <- simulatedMargins(N = 300L, seed = 1)
  at <- c(0.5, 0.25)

  for (marg.type in c("atmean", "aveacr")) {
    x <- if (marg.type == "atmean") matrix(c(1, at), nrow = 1L) else object$x
    result <- effects.fmlogit(object, effect = "discrete", marg.type = marg.type,
                              at = if (marg.type == "atmean") at else NULL)

    # For the binary indicator, this is the change from zero to one.
    for (k in 2:ncol(x)) {
      x.low <- x.high <- x
      x.low[, k] <- min(object$x[, k])
      x.high[, k] <- max(object$x[, k])
      expected <- shareMeans(x.high, object$coefficient) -
        shareMeans(x.low, object$coefficient)
      expect_equal(as.numeric(result$effects[, k-1L]),
                   as.numeric(colMeans(expected)), tolerance = 1e-8)
    }

    # A covariate redistributes shares; it cannot change their total.
    expect_equal(as.numeric(colSums(result$effects)), c(0, 0), tolerance = 1e-8)
  }
})

test_that("Average marginal effects average individual expected-share derivatives", {
  object <- simulatedMargins(N = 300L, seed = 1)
  result <- effects.fmlogit(object, marg.type = "aveacr")
  step <- 1e-5

  for (k in 2:ncol(object$x)) {
    x.low <- x.high <- object$x
    x.low[, k] <- x.low[, k] - step
    x.high[, k] <- x.high[, k] + step
    individual <- (shareMeans(x.high, object$coefficient) -
                   shareMeans(x.low, object$coefficient)) / (2*step)

    # Check both the heterogeneous effects and their sample average.
    expect_equal(result$marg.list[[k-1L]], individual, tolerance = 1e-7)
    expect_equal(as.numeric(result$effects[, k-1L]),
                 as.numeric(colMeans(individual)), tolerance = 1e-7)
  }
  expect_equal(as.numeric(colSums(result$effects)), c(0, 0), tolerance = 1e-12)
})

test_that("Two-outcome margins reduce to the binary-logit formula", {
  object <- simulatedMargins(N = 300, J = 2L, seed = 1)
  p.other <- shareMeans(object$x, object$coefficient)[, 2L]
  slopes <- object$coefficient[2L, -1L]

  # dp_2/dx_k = p_2 (1 - p_2) beta_2k, and dp_1/dx_k = -dp_2/dx_k.
  expected <- colMeans(outer(p.other * (1 - p.other), slopes))
  result <- effects.fmlogit(object, marg.type = "aveacr")
  expect_equal(unname(result$effects), unname(rbind(-expected, expected)), tolerance = 1e-8)

  p.mean <- shareMeans(matrix(colMeans(object$x), nrow = 1L), object$coefficient)[1L, 2L]
  expected <- p.mean * (1 - p.mean)*slopes
  result <- effects.fmlogit(object)
  expect_equal(unname(result$effects), unname(rbind(-expected, expected)), tolerance = 1e-8)
})

test_that("Baseline effect uncertainty includes covariance between outcome slopes", {
  object <- simulatedMargins(N = 300L, seed = 1)
  p.at <- as.numeric(shareMeans(matrix(c(1, 0, 0), nrow = 1L), object$coefficient))
  slopes <- c("second:income", "third:income")
  sigma <- 0.01

  # At zero covariates, mean shares depend only on the fixed intercepts.
  # The baseline income effect is -p_1 (p_2 beta_2 + p_3 beta_3),
  # so its variance is known exactly when only the two slopes are uncertain.
  for (rho in c(0, 0.85)) {
    object$vcov.full[, ] <- 0
    object$vcov.full[slopes, slopes] <- sigma^2*matrix(c(1, rho, rho, 1), nrow = 2L)
    expected <- sigma*p.at[1L]*sqrt(p.at[2L]^2 + p.at[3L]^2 + 2*rho*p.at[2L]*p.at[3L])

    set.seed(123)
    result <- effects.fmlogit(object, varlist = "income", at = c(0, 0), se = TRUE, R = 1000L)
    # Allow Monte Carlo error, while checking independent and correlated slopes.
    expect_lt(abs(result$se[1L, 1L]/expected - 1), 0.1)
  }
})

test_that("Krinsky-Robb standard errors agree with the delta method", {
  object <- simulatedMargins(N = 300L, seed = 1)
  # Small coefficient uncertainty makes the delta approximation accurate.
  object$vcov.full <- object$vcov.full/100
  theta <- object$raw.estimates
  J <- ncol(object$y)
  step <- 1e-5

  for (effect in c("marginal", "discrete")) {
    for (marg.type in c("atmean", "aveacr")) {
      at <- if (marg.type == "atmean") c(1.25, 0.2) else NULL
      effectAt <- function(theta) {
        object$raw.estimates <- theta
        as.numeric(effects.fmlogit(object, effect = effect, marg.type = marg.type,
                                   varlist = "income", at = at)$effects)
      }

      # Numerically differentiate the requested statistic with respect to theta.
      gradient <- vapply(seq_along(theta), function(q) {
        lower <- upper <- theta
        lower[q] <- lower[q] - step
        upper[q] <- upper[q] + step
        (effectAt(upper) - effectAt(lower)) / (2*step)
      }, numeric(J))

      # Delta-method variance: gradient %*% vcov(theta) %*% t(gradient).
      expected <- sqrt(diag(gradient %*% object$vcov.full %*% t(gradient)))
      set.seed(123)
      result <- effects.fmlogit(object, effect = effect, marg.type = marg.type,
                                varlist = "income", at = at, se = TRUE, R = 1000L)
      expect_lt(max(abs(as.numeric(result$se)/expected - 1)), 0.1)
    }
  }
})
