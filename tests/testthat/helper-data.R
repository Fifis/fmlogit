simulatedShares <- function(N = 500L, concentration = 10, trials = 10L, seed = 1L) {
  set.seed(seed)
  x <- cbind(1, matrix(rnorm(N*4L), ncol = 4L))
  colnames(x) <- c("(Intercept)", paste0("x", 1:4))
  beta <- rbind(
    rep(0, 5L),
    c(1, 1, 1, 1, 1),
    c(1, -1, 1, -1, 1)
  )
  eta <- tcrossprod(x, beta)
  numerator <- exp(sweep(eta, 1L, apply(eta, 1L, max), FUN = "-"))
  p <- numerator/rowSums(numerator)
  gamma.draws <- matrix(rgamma(N*3L, shape = concentration*as.vector(p)), nrow = N)
  y.dirichlet <- gamma.draws/rowSums(gamma.draws)
  y.multinomial <- t(vapply(seq_len(N), function(i) {
    as.vector(rmultinom(1L, size = trials, prob = p[i, ]))/trials
  }, numeric(3L)))
  colnames(y.dirichlet) <- colnames(y.multinomial) <- paste0("share", 1:3)
  list(x = x, beta = beta, p = p, y.dirichlet = y.dirichlet, y.multinomial = y.multinomial)
}
