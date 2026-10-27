weighted_example <- function() {
  s <- simulated_shares(N = 200L, seed = 12L)
  weights <- ifelse(s$x[, "x1"] > 0, 5, 1)
  fmlogit(s$y_dirichlet, s$x, weights = weights,
          method = "BFGS", print.level = 0L)
}

test_that("Population effects use weighted averages and a weighted mean profile", {
  fit <- weighted_example()
  weights <- fit$weights / sum(fit$weights)

  for (effect in c("marginal", "discrete")) {
    result <- effects.fmlogit(fit, effect = effect, marg.type = "aveacr", varlist = "x1")
    expected <- colSums(result$marg.list$x1 * weights)
    expect_equal(as.numeric(result$effects), as.numeric(expected))
  }

  at <- colSums(fit$x[, -1L, drop = FALSE] * weights)
  expect_equal(effects.fmlogit(fit)$effects, effects.fmlogit(fit, at = at)$effects)

  # A supplied profile describes that point, regardless of the averaging weights.
  other <- fit
  other$weights <- rev(fit$weights)
  expect_equal(effects.fmlogit(fit, at = at)$effects,
               effects.fmlogit(other, at = at)$effects)
})

test_that("Weight units do not change effects or their simulated uncertainty", {
  s <- simulated_shares(N = 200L, seed = 12L)
  weights <- ifelse(s$x[, "x1"] > 0, 5, 1)
  # Both default BHHH fits start independently: large raw weights must not stall it.
  fit <- fmlogit(s$y_dirichlet, s$x, weights = weights, print.level = 0L)
  rescaled <- fmlogit(s$y_dirichlet, s$x, weights = 1e4 * weights, print.level = 0L)
  expect_equal(fit$raw.estimates, rescaled$raw.estimates, tolerance = 1e-7)
  expect_equal(vcov(fit), vcov(rescaled), tolerance = 1e-7)
  expect_equal(rescaled$likelihood, 1e4 * fit$likelihood, tolerance = 1e-7)
  expect_equal(rescaled$weights, 1e4 * weights)
  set.seed(123)
  result <- effects.fmlogit(fit, marg.type = "aveacr", varlist = "x1", se = TRUE, R = 400L)
  set.seed(123)
  comparison <- effects.fmlogit(rescaled, marg.type = "aveacr", varlist = "x1", se = TRUE, R = 400L)
  expect_equal(result$effects, comparison$effects, tolerance = 1e-7)
  expect_equal(result$se, comparison$se, tolerance = 1e-7)
  expect_equal(result$cov.list, comparison$cov.list, tolerance = 1e-7)

  # Equal averaging weights recover the ordinary sample average.
  fit$weights[] <- 1
  equal <- effects.fmlogit(fit, marg.type = "aveacr", varlist = "x1")
  expect_equal(as.numeric(equal$effects), as.numeric(colMeans(equal$marg.list$x1)))
})

test_that("A zero-weight outlier does not define the reference profile or endpoints", {
  fit <- weighted_example()
  outlier <- fit$x[1L, , drop = FALSE]
  outlier[, "x1"] <- 100
  augmented <- fmlogit(rbind(fit$y, fit$y[1L, ]), rbind(fit$x, outlier),
                      weights = c(fit$weights, 0), beta0 = fit$raw.estimates,
                      method = "BFGS", print.level = 0L)

  for (effect in c("marginal", "discrete")) {
    for (marg.type in c("atmean", "aveacr")) {
      expected <- effects.fmlogit(fit, effect = effect, marg.type = marg.type, varlist = "x1")
      result <- effects.fmlogit(augmented, effect = effect, marg.type = marg.type, varlist = "x1")
      expect_equal(result$effects, expected$effects, tolerance = 1e-4)
    }
  }
})

test_that("Weighted Krinsky-Robb SEs agree with the conditional delta method", {
  fit <- weighted_example()
  fit$vcov_full <- fit$vcov_full / 100  # Keep the linear approximation accurate.
  theta <- fit$raw.estimates
  step <- 1e-5
  effect_at <- function(theta) {
    fit$raw.estimates <- theta
    as.numeric(effects.fmlogit(fit, marg.type = "aveacr", varlist = "x1")$effects)
  }
  gradient <- vapply(seq_along(theta), function(q) {
    lower <- upper <- theta
    lower[q] <- lower[q] - step
    upper[q] <- upper[q] + step
    (effect_at(upper) - effect_at(lower)) / (2 * step)
  }, numeric(ncol(fit$y)))

  expected <- sqrt(diag(gradient %*% fit$vcov_full %*% t(gradient)))
  set.seed(123)
  result <- effects.fmlogit(fit, marg.type = "aveacr", varlist = "x1", se = TRUE, R = 1000L)
  expect_lt(max(abs(as.numeric(result$se) / expected - 1)), 0.1)
})
