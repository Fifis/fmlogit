test_that("Starting coefficients match the processed design", {
  s <- simulatedShares(N = 300L)
  theta.init <- initEst(s$y.multinomial, s$x)
  expect_length(theta.init, 2L*ncol(s$x))
  expect_true(all(is.finite(theta.init)))
  
  fit <- fmlogit(s$y.multinomial, s$x, beta0 = theta.init, method = "BHHH")
  expect_equal(colnames(fit$x), colnames(s$x))
  expect_true(all(fit$coefficient[1L, ] == 0))
  expect_equal(dim(fit$coefficient), c(3L, 5L))
  expect_equal(dim(fit$vcov.full), c(10L, 10L))
  expect_equal(fit$vcov.full, t(fit$vcov.full), tolerance = 1e-10)
  expect_equal(stats::vcov(fit), fit$vcov.full)
  expect_equal(rowSums(fit$fitted.values), rep(1, nrow(s$x)), tolerance = 1e-12)
  expect_true(fit$conv.code %in% c(0L, 1L, 2L, 8L))
  expect_true(all(is.finite(fit$raw.estimates)))
})

test_that("Intercept-only starting values reproduce the mean-share ratio", {
  y <- cbind(base = c(0.2, 0.3, 0.4), other = c(0.8, 0.7, 0.6))
  theta.init <- initEst(y, matrix(1, nrow = nrow(y), ncol = 1L))
  expect_equal(unname(theta.init), qlogis(mean(y[, "other"])))
  expect_identical(names(theta.init), "other:(Intercept)")
})

test_that("Both multinomial and Dirichlet DGPs are admitted", {
  s <- simulatedShares(N = 300L)
  for (y in list(s$y.dirichlet, s$y.multinomial)) {
    fit <- fmlogit(y, s$x, method = "BFGS")
    expect_true(all(is.finite(fit$coefficient)))
    expect_equal(rowSums(fit$fitted.values), rep(1, nrow(y)), tolerance = 1e-12)
    expect_equal(length(fit$vcov), ncol(y))
    expect_equal(names(fit$estimates), colnames(y)[-1L])
  }
})

test_that("Rescaling relative observation weights leaves HC0 inference unchanged", {
  s <- simulatedShares(N = 300L)
  set.seed(123)
  weights <- runif(nrow(s$y.dirichlet), 0.5, 2)  # Relative observation weights
  fit1 <- fmlogit(s$y.dirichlet, s$x, weights = weights)
  fit2 <- fmlogit(s$y.dirichlet, s$x, weights = 10*weights, beta0 = fit1$raw.estimates)
  expect_equal(fit1$raw.estimates, fit2$raw.estimates, tolerance = 1e-3)
  expect_equal(fit1$vcov.full, fit2$vcov.full, tolerance = 1e-3)
})

test_that("Missing rows and weights stay aligned", {
  s <- simulatedShares(N = 300L)
  y <- s$y.multinomial
  x <- s$x
  y[3L, 2L] <- NA_real_
  x[5L, 2L] <- NA_real_
  weights <- seq_len(nrow(y))/nrow(y)
  fit <- fmlogit(y, x, weights = weights)
  expect_equal(fit$row.no, setdiff(seq_len(nrow(y)), c(3L, 5L)))
  expect_equal(fit$weights, weights[fit$row.no])
})

test_that("Rank-deficient designs and invalid shares are rejected", {
  s <- simulatedShares(N = 300L)
  x.bad <- cbind(s$x, duplicate = s$x[, 2L])
  expect_error(fmlogit(s$y.multinomial, x.bad), "rank-deficient")
  y.bad <- s$y.multinomial
  y.bad[1L, 1L] <- 2
  expect_error(fmlogit(y.bad, s$x), "shares")
})

test_that("The joint HC0 covariance matches a numerical Hessian", {
  s <- simulatedShares(N = 300L)
  fit <- fmlogit(s$y.dirichlet, s$x, method = "BHHH")
  x <- fit$x
  y <- fit$y
  weights <- fit$weights
  J <- ncol(y)
  K <- ncol(x)
  p <- fit$fitted.values
  
  scores <- do.call(cbind, lapply(seq.int(2L, J), function(j) {
    x*(weights*(y[, j] - rowSums(y)*p[, j]))
  }))
  
  # An independent objective checks the analytic joint covariance.
  loglik <- function(theta) {
    beta <- rbind(rep(0, K), matrix(theta, nrow = J-1L, ncol = K, byrow = TRUE))
    eta <- tcrossprod(x, beta)
    shift <- apply(eta, 1L, max)
    log.denom <- shift + log(rowSums(exp(sweep(eta, 1L, shift, FUN = "-"))))
    sum(weights*(rowSums(y*eta) - rowSums(y)*log.denom))
  }
  bread.num <- -stats::optimHess(fit$raw.estimates, fn = loglik)
  bread.inv <- solve(bread.num)
  vcov.num <- bread.inv %*% crossprod(scores) %*% bread.inv
  expect_equal(unname(fit$vcov.full), unname(vcov.num), tolerance = 1e-4)
})

test_that("Fitted shares and response residuals are consistent", {
  s <- simulatedShares(N = 300L)
  fit <- fmlogit(s$y.dirichlet, s$x, print.level = 0L)
  p <- as.matrix(fitted(fit))
  expect_equal(p, fit$fitted.values)
  expect_equal(as.matrix(predict(fit)), p)
  expect_equal(as.matrix(residuals(fit)), fit$y - p)
  expect_equal(as.numeric(rowSums(residuals(fit))), rep(0, nrow(fit$y)),
               tolerance = 1e-12)
})

test_that("Predictions recover known multinomial shares and honor new coefficients", {
  s <- simulatedShares(N = 300L)
  fit <- fmlogit(s$y.multinomial, s$x, print.level = 0L)
  newdata <- s$x[1:3, -1L, drop = FALSE]
  p <- predict(fit, newdata = newdata, newbeta = s$beta)
  expect_equal(unname(as.matrix(p)), unname(s$p[1:3, ]), tolerance = 1e-12)

  p <- predict(fit, newdata = newdata[1L, ], newbeta = s$beta)
  expect_equal(as.numeric(p[1L, ]), as.numeric(s$p[1L, ]), tolerance = 1e-12)

  # All-zero coefficients imply equal shares, including without newdata.
  p <- predict(fit, newbeta = 0*fit$coefficient)
  expect_equal(unname(as.matrix(p)), matrix(1/3, nrow = nrow(fit$y), ncol = 3L),
               tolerance = 1e-12)
})
