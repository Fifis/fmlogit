test_that("fmlogit goes through on simulated data", {
  n <- 200 # Observations
  j <- 3   # Alternatives
  k <- 4   # Explanatory variables
  set.seed(1)
  x <- matrix(rnorm(n * k), ncol = k)
  x1 <- cbind(1, x) # Covariates
  colnames(x) <- paste0("X", 1:k)
  colnames(x1) <- c("(Intercept)", colnames(x))
  
  # True coefficients: j x (k + 1)
  # The first alternative is the baseline (hence the normalisation to 0)
  # Last column contains the intercepts.
  betas <- rbind(
    0,
    c(1, 1, 1, 1, 1),
    c( 1, -1,  1,  -1,  1)
  )
  rownames(betas) <- paste0("share", 1:j)
  colnames(betas) <- colnames(x1)
  
  # Multinomial-logit conditional means
  eta <- x1 %*% t(betas)
  # Numerically stable way of computing the exponential terms in the
  # multinomial logit model. The subtraction does not change the resulting
  # probabilities (or expected shares). It only prevents numerical overflow
  # when exponentiating large numbers.
  e <- exp(eta - apply(eta, 1, max))
  p <- e / rowSums(e)  # p_ij
  
  # Dirichlet DGP with concentration parameter phi
  phi <- 10
  g <- matrix(rgamma(n * j, shape = phi * as.vector(p)), nrow = n)
  y1 <- g / rowSums(g)
  colnames(y1) <- paste0("share", 1:j)
  
  # Observed fractional shares: 0, 1/m 2/m, ..., (m-1)/m, 1:
  m <- 10
  y2 <- t(vapply(seq_len(n), function(i) {
    as.vector(rmultinom(1, size = m, prob = p[i, ])) / m
  }, numeric(j)))
  colnames(y2) <- colnames(y1)
  
  # Estimate
  mod1 <- fmlogit(y = y1, x = x, method = "BHHH")
  mod2 <- fmlogit(y = y2, x = x, method = "BHHH")
  
  expect_identical(mod1$conv_code, 0L)
  expect_identical(mod2$conv_code, 0L)
  
  s1 <- summary.fmlogit(mod1)
  s2 <- summary.fmlogit(mod2)
  expect_identical(length(s1), 5L)
  expect_identical(length(s2), 5L)

  # Optimisation from good initial values must work
  mod2bfgs <- fmlogit(y = y2, x = x, beta0 = mod2$raw.estimates, method = "BFGS")
  expect_lt(max(abs(mod2$raw.estimates - mod2bfgs$raw.estimates)), 0.001)
})
