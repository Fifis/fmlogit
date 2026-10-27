#' @importFrom Rdpack reprompt
#' @importFrom maxLik maxLik
#' @importFrom stats effects fitted predict residuals vcov

asterisk <- function(p.value, sigcode = c(0.10, 0.05, 0.01, 0.001)) {
  stopifnot(length(sigcode) > 0, all(diff(sigcode) < 0), !anyNA(sigcode))

  n.codes <- length(sigcode)
  code.index <- 0:n.codes

  symbols <- if (n.codes <= 3) {
    strrep("*", code.index)
  } else if (n.codes == 4) {
    c("", ".", strrep("*", seq_len(n.codes - 1)))
  } else {
    paste0(strrep("*", code.index %/% 2), ifelse(code.index %% 2 == 1, ".", ""))
  }

  level <- n.codes - findInterval(p.value, rev(sigcode), left.open = TRUE)
  symbols[level + 1]
}

# Shared notation: y is N x J, x is N x K with the intercept first, and beta
# is J x K with beta_1 = 0. theta stacks the K-vectors for outcomes 2, ..., J.
# eta_ij = x_i' beta_j is the linear index; p_ij is the conditional mean share,
# called xi_k(x; beta) in Mullahy (2015), with his outcome index k renamed j.
# xi does not denote the linear index.
softmax <- function(theta, x, J) {
  K <- ncol(x)
  beta <- rbind(rep(0, K), matrix(theta, nrow = J-1L, ncol = K, byrow = TRUE))

  eta <- tcrossprod(x, beta)       # N x J
  shift <- apply(eta, 1L, max)     # N
  shifted <- sweep(eta, 1L, shift, FUN = "-")
  numerator <- exp(shifted)
  denominator <- rowSums(numerator)

  list(
    p = numerator/denominator,
    beta = beta,
    eta = eta,
    log.denom = shift + log(denominator)
  )
}

# The same two-sided normal approximation is used for coefficients and effects.
waldTable <- function(est, se) {
  z <- est/se
  cbind(Estimate = est, `Std. Error` = se, `z value` = z,
        `Pr(>|z|)` = 2*stats::pnorm(abs(z), lower.tail = FALSE))
}

# Post-estimation inherits the retained observation weights from the fit.
fitWeights <- function(object) {
  N <- nrow(object$x)
  weights <- object$weights
  if (is.null(weights)) weights <- rep(1, N)
  if (!is.numeric(weights) || length(weights) != N ||
      any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
    stop("'object$weights' must be finite, nonnegative and aligned with 'object$x', with at least one positive weight.")
  }
  weights
}

# Each column estimates a mean over the distribution represented by the weights.
weightedColMeans <- function(x, weights) {
  weights <- weights/max(weights)  # Keep the normalising sum finite.
  colSums(x*(weights/sum(weights)))
}
