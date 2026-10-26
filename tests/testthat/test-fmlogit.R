simulatedShares <- function(n = 500L, phi = 10, m = 10L, seed = 1L) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * 4L), ncol = 4L))
  colnames(X) <- c("(Intercept)", paste0("X", 1:4))
  B <- rbind(
    rep(0, 5L),
    c(1, 1, 1, 1, 1),
    c(1, -1, 1, -1, 1)
  )
  eta <- X %*% t(B)
  e <- exp(sweep(eta, 1L, apply(eta, 1L, max), FUN = "-"))
  P <- e / rowSums(e)
  gammas <- matrix(rgamma(n * 3L, shape = phi * as.vector(P)), nrow = n)
  Y_dirichlet <- gammas / rowSums(gammas)
  Y_multinom <- t(vapply(seq_len(n), function(i) {
    as.vector(rmultinom(1L, size = m, prob = P[i, ])) / m
  }, numeric(3L)))
  colnames(Y_dirichlet) <- colnames(Y_multinom) <- paste0("share", 1:3)
  list(x = X, beta = B, p = P, y1 = Y_dirichlet, y2 = Y_multinom)
}

test_that("Starting coefficients match the processed design", {
  s <- simulatedShares(n = 300L)
  b0 <- initEst(s$y2, s$x)
  expect_length(b0, 2L * ncol(s$x))
  expect_true(all(is.finite(b0)))
  
  fit <- fmlogit(s$y2, s$x, beta0 = b0, method = "BHHH")
  expect_equal(colnames(fit$x), colnames(s$x))
  expect_true(all(fit$coefficient[1L, ] == 0))
  expect_equal(dim(fit$coefficient), c(3L, 5L))
  expect_equal(dim(fit$vcov_full), c(10L, 10L))
  expect_equal(fit$vcov_full, t(fit$vcov_full), tolerance = 1e-10)
  expect_equal(stats::vcov(fit), fit$vcov_full)
  expect_equal(rowSums(fit$fitted.values), rep(1, nrow(s$x)), tolerance = 1e-12)
  expect_true(fit$conv_code %in% c(0L, 1L, 2L, 8L))
  expect_true(all(is.finite(fit$raw.estimates)))
})

test_that("Both multi-nomial and Dirichlet DGPs are admitted", {
  s <- simulatedShares(n = 300L)
  for (Y in list(s$y1, s$y2)) {
    fit <- fmlogit(Y, s$x, method = "BFGS")
    expect_true(all(is.finite(fit$coefficient)))
    expect_equal(rowSums(fit$fitted.values), rep(1, nrow(Y)), tolerance = 1e-12)
    expect_equal(length(fit$vcov), ncol(Y))
    expect_equal(names(fit$estimates), colnames(Y)[-1L])
  }
})

test_that("Rescaling probability weights leaves HC0 inference unchanged", {
  s <- simulatedShares(n = 300L)
  set.seed(123)
  w <- runif(nrow(s$y1), 0.5, 2)  # Importance weights
  fit1 <- fmlogit(s$y1, s$x, weights = w)
  fit2 <- fmlogit(s$y1, s$x, weights = 10 * w, beta0 = fit1$raw.estimates)
  expect_equal(fit1$raw.estimates, fit2$raw.estimates, tolerance = 1e-3)
  expect_equal(fit1$vcov_full, fit2$vcov_full, tolerance = 1e-3)
})

test_that("Missing rows and weights stay aligned", {
  s <- simulatedShares(n = 300L)
  Y <- s$y2
  X <- s$x
  Y[3L, 2L] <- NA_real_
  X[5L, 2L] <- NA_real_
  w <- seq_len(nrow(Y)) / nrow(Y)
  fit <- fmlogit(Y, X, weights = w)
  expect_equal(fit$rowNo, setdiff(seq_len(nrow(Y)), c(3L, 5L)))
  expect_equal(fit$weights, w[fit$rowNo])
})

test_that("Rank-deficient designs and invalid shares are rejected", {
  s <- simulatedShares(n = 300L)
  X_bad <- cbind(s$x, duplicate = s$x[, 2L])
  expect_error(fmlogit(s$y2, X_bad), "rank-deficient")
  Y_bad <- s$y2
  Y_bad[1L, 1L] <- 2
  expect_error(fmlogit(Y_bad, s$x), "shares")
})

test_that("The joint HC0 covariance matches a numerical Hessian", {
  s <- simulatedShares(n = 300L)
  fit <- fmlogit(s$y1, s$x, method = "BHHH")
  X <- fit$x
  Y <- fit$y
  w <- fit$weights
  J <- ncol(Y)
  K <- ncol(X)
  P <- fit$fitted.values
  
  scores <- do.call(cbind, lapply(2:J, function(j) {
    X * (w * (Y[, j] - rowSums(Y) * P[, j]))
  }))
  
  loglik <- function(theta) {
    B <- rbind(rep(0, K), matrix(theta, J - 1L, K, byrow = TRUE))
    eta <- tcrossprod(X, B)
    shift <- apply(eta, 1L, max)
    log_denom <- shift + log(rowSums(exp(sweep(eta, 1L, shift, FUN = "-"))))
    sum(w * (rowSums(Y * eta) - rowSums(Y) * log_denom))
  }
  H_num <- -stats::optimHess(fit$raw.estimates, fn = loglik)
  H_inv <- solve(H_num)
  V_num <- H_inv %*% crossprod(scores) %*% H_inv
  expect_equal(unname(fit$vcov_full), unname(V_num), tolerance = 1e-4)
})
