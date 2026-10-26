#' Obtain pairwise fractional-logit starting values
#'
#' @description
#' Fits a separate, share-weighted fractional binary logit for each non-baseline
#' outcome against the first outcome. These estimators target the multinomial
#' conditional-mean coefficients under correct specification, but are generally
#' not identical to the joint fractional multinomial QMLE in finite samples.
#'
#' @param y Numeric N-by-J matrix of observed shares, with J >= 2.
#'   Each row must sum to one.
#' @param x Numeric N-by-K *design matrix*, including a column of ones first.
#'   The columns must have the same order used during joint estimation.
#' @param weights Optional vector of N nonnegative observation weights.
#' @return Named numeric vector of length (J - 1) * K. Coefficients are grouped
#'   by non-baseline alternative; the baseline coefficients are omitted.
#'   
#' @details
#' For alternative j > 1, fit a fractional binary logit with response
#' \eqn{y_j / (y_1 + y_j)}{y_j / (y_1 + y_j)}, and prior weights
#' \eqn{w (y_1 + y_j)}{w * (y_1 + y_j)}.
#' Observations with \eqn{y_1 + y_j = 0}{y_1 + y_j = 0} make no contribution.
#' In the event of failed or non-finite pairwise estimates, the corresponding
#' coefficient block is initialized to zero (with a warning).
#' 
#' The reasonable initial value comes from `J-1` separate weighted binary-logit
#' regressions, each comparing one alternative with the baseline alternative.
#' Since the normalisation is such that the first `K` elements are zeros (and
#' do not enter the optimisation), we note that
#' \eqn{E(y_{ij}\mid x_i) = p_{ij}}{E(y_ij | x_i) = p_ij} and
#' \eqn{p_{ij}/p_{i1} = \exp (x_i' \beta_j)}{p_ij / p_i1 = exp (x_i' beta_j)}.
#' Now define \eqn{\pi_{ij} = \frac{p_{ij}}{p_{i1} + p_{ij}}}{pi_ij = p_ij / (p_i1 + p_ij)}.
#' Substituting the multinomial-logit probabilities gives
#' \eqn{\pi_{ij} = \frac{\exp(x_i' \beta_j)}{1 + \exp(x_i' \beta_j)} = 
#' \Lambda(x_i' \beta_j)}{pi_ij = exp(x_i' beta_j) / (1 + exp(x_i' beta_j)) =
#' Lambda(x_i' beta_j)}.
#' This is a pairwise application of the fractional quasi-likelihood logic in
#' \insertCite{papke1996econometric}{fmlogit} and
#' \insertCite{mullahy2015multivariate}{fmlogit}.
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso [fmlogit()]
#' 
#' 
#' @export
initEst <- function(y, x, weights = NULL) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  
  if (!is.numeric(y) || !is.numeric(x) ||
      length(dim(y)) != 2L || length(dim(x)) != 2L ||
      nrow(y) != nrow(x) || ncol(y) < 2L || ncol(x) < 1L ||
      any(!is.finite(y)) || any(!is.finite(x))) {
    stop("'y' and 'x' must be conformable, finite numeric matrices.")
  }
  if (any(abs(x[, 1L] - 1) > 1e-12)) {
    stop("'x' must have an intercept (a column of ones) in its first column.")
  }
  if (any(y < 0 | y > 1) || any(abs(rowSums(y) - 1) > 1e-8)) {
    stop("'y' must contain shares in [0, 1], with rows summing to one.")
  }
  
  N <- nrow(y)
  J <- ncol(y)
  K <- ncol(x)
  if (is.null(weights)) weights <- rep(1, N)
  if (!is.numeric(weights) || length(weights) != N ||
      any(!is.finite(weights)) || any(weights < 0) ||
      !any(weights > 0)) {
    stop("'weights' must be nonnegative, finite, aligned with 'y', and not all zero.")
  }
  
  if (is.null(colnames(y))) colnames(y) <- paste0("share", seq_len(J))
  if (is.null(colnames(x))) {
    colnames(x) <- c("(Intercept)", paste0("X", seq_len(K - 1L)))
  }
  
  B <- matrix(0, nrow = J - 1L, ncol = K, dimnames = list(colnames(y)[-1L], colnames(x)))
  
  for (h in seq.int(2L, J)) {
    pair_total <- y[, 1L] + y[, h]
    prior_weights <- pair_total * weights
    keep <- prior_weights > 0
    
    if (!any(keep)) {
      warning(sprintf("No positive pairwise weight for '%s'; using zeros.", colnames(y)[h]))
      next
    }
    
    response <- y[keep, h] / pair_total[keep]
    fit <- tryCatch(
      stats::glm.fit(
        x = x[keep, , drop = FALSE], y = response, weights = prior_weights[keep],
        family = stats::quasibinomial(link = "logit"), intercept = TRUE
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit) || !isTRUE(fit$converged) ||
        length(fit$coefficients) != K ||
        any(!is.finite(fit$coefficients))) {
      warning(sprintf("Pairwise GLM for '%s' failed; using zeros.", colnames(y)[h]))
      next
    }
    
    B[h - 1L, ] <- fit$coefficients
  }
  
  ans <- as.vector(t(B))
  names(ans) <- paste(rep(rownames(B), each = K),
                      rep(colnames(B), times = J - 1L), sep = ":")
  ans
}




#' Estimate a fractional multinomial logit model
#'
#' @description
#' Estimates a multinomial-logit conditional-mean model for compositional
#' outcomes by maximizing the weighted multinomial *quasi*-log-likelihood.
#' The first outcome is the reference category. Robust standard errors are
#' computed using the full joint, observation-level sandwich covariance matrix.
#'
#' @param y Numeric N-by-J matrix or data frame of observed shares, `J >= 2`.
#'   Values must lie in [0, 1] and rows must sum to one.
#'   The first column is automatically treated as the baseline.
#' @param x `N`-row matrix or data frame of regressors. Numeric, factor, and
#'   character columns are allowed in data frames. An intercept is added
#'   automatically as the first column. Any supplied constant columns are
#'   removed before creating the design matrix. No variable selection is done.
#' @param weights Optional non-negative observation weights, length `N`. The
#'   default is equal weights. These are relative sampling weights that multiply
#'   individual pseudo-likelihood values. Their behaviour is similar to
#'   the `weights` argument in `lm()`.
#' @param beta0 Optional starting coefficients, length `(J-1) * K`, where `K` is
#'   the number of columns in the *processed* design matrix including the
#'   intercept. Ordered by alternative, with the intercept first in each block.
#'   If NULL, use pairwise fractional binary logits from [initEst()].
#' @param method Optimiser for [maxLik::maxLik()]; default "BHHH". Use "CG"
#'   if the estimation is very slow.
#' @param iterlim Positive iteration limit passed to `maxLik()`.
#' @param reltol Positive relative convergence tolerance passed to `maxLik()`.
#' @param print.level Optimiser verbosity; default 1 (short output).
#' @param ... Additional arguments forwarded to [maxLik::maxLik()].
#'
#' @details
#' For the first, reference category, \eqn{\beta_1 = 0}{beta_1 = 0}.
#' The conditional means are
#' \eqn{p_{ij} = \exp(x_i' \beta_j) / \sum_h \exp(x_i' \beta_h)}{p_ij = exp(x_i' beta_j) / sum_h exp(x_i' beta_h)}.
#' The objective is \eqn{\sum_i w_i \sum_j y_{ij} \log p_ij}{sum_i w_i sum_j y_ij log(p_ij)},
#' as per \insertCite{mullahy2015multivariate}{fmlogit}.
#' The joint HC0 covariance uses the inverse weighted negative Hessian as
#' bread and the outer product of the weighted *full* scores as meat.
#' This is an independent-observation sandwich estimator: weights alone do not
#' account for survey clusters, strata, or finite-population corrections
#' \insertCite{papke1996econometric}{fmlogit}.
#' Rows with missing responses/covariates are removed before estimation;
#' the input weights are subset using exactly the same row indices.
#'
#' @return An object of class "fmlogit". The existing fields `raw.estimates`,
#'   `estimates`, `baseline`, `likelihood`, `conv_code`, `convergence`,
#'   `count`, `y`, `x`, `weights`, `rowNo`, and `coefficient` are retained.
#'   `coefficient` is J-by-K, with an all-zero first row.
#'   `vcov_full` is the complete `K(J-1)`-by-`K(J-1)` joint HC0 covariance matrix.
#'   `vcov` is a legacy list of J diagonal blocks of `vcov_full`, with a zero
#'   block for the normalized baseline (cross-category blocks are only in
#'   `vcov_full`). `fitted.values` holds the `N`-by-`J` predicted shares.
#'
#' @references
#' \insertAllCited{}
#' 
#' @export
fmlogit <- function(y, x, weights = NULL, beta0 = NULL,
                    method = "BHHH", iterlim = 2000L,
                    reltol = 1e-8, print.level = 1L, ...) {
  y <- as.matrix(y)
  if (!is.numeric(y) || ncol(y) < 2L || nrow(y) < 1L) {
    stop("'y' must have at least two numeric share columns and one row.")
  }
  N0 <- nrow(y)
  x_data <- as.data.frame(x, stringsAsFactors = FALSE)
  if (nrow(x_data) != N0) {
    stop("'x' and 'y' must have the same number of rows.")
  }
  if (is.null(colnames(y))) colnames(y) <- paste0("share", seq_len(ncol(y)))
  if (anyNA(colnames(y)) || anyDuplicated(colnames(y)) || any(!nzchar(colnames(y)))) {
    stop("'y' must have unique, nonempty column names.")
  }
  if (ncol(x_data) && (anyNA(names(x_data)) || anyDuplicated(names(x_data)) || any(!nzchar(names(x_data))))) {
    stop("'x' must have unique, nonempty column names.")
  }
  
  if (is.null(weights)) weights <- rep(1, N0)
  if (!is.numeric(weights) || length(weights) != N0 || any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must be a finite, nonnegative numeric vector of length nrow(y).")
  }
  
  keep_rows <- stats::complete.cases(y, x_data)
  rowNo <- which(keep_rows)
  if (!length(rowNo)) stop("No complete observations remain.")
  y <- y[keep_rows, , drop = FALSE]
  x_data <- x_data[keep_rows, , drop = FALSE]
  weights <- weights[keep_rows]
  if (!any(weights > 0)) stop("All retained observation weights are zero.")
  
  if (any(!is.finite(y)) || any(y < 0 | y > 1) || any(abs(rowSums(y) - 1) > 1e-8)) {
    stop("'y' must contain finite shares in [0, 1], each row summing to one.")
  }
  
  # Construct a single design matrix for starting values and estimation.
  # The intercept is deliberately the first column.
  if (ncol(x_data) > 0L) {
    good_type <- vapply(x_data, function(v) {
      is.numeric(v) || is.factor(v) || is.character(v) || is.logical(v)
    }, logical(1))
    if (!all(good_type)) {
      stop("'x' contains unsupported column types in columns ", which(!good_type)-1, ".")
    }
    x_data <- lapply(x_data, function(v) {
      if (is.factor(v)) droplevels(v) else v
    })
    x_data <- as.data.frame(x_data, check.names = FALSE)
    is_constant <- vapply(x_data, function(v) length(unique(v)) <= 1L, logical(1))
    x_data <- x_data[, !is_constant, drop = FALSE]
  }
  x <- if (ncol(x_data) > 0L) {
    stats::model.matrix(~ ., data = x_data)
  } else {
    matrix(1, nrow(y), 1L, dimnames = list(NULL, "(Intercept)"))
  }
  if (any(!is.finite(x))) stop("'x' contains non-finite design values.")
  if (colnames(x)[1L] != "(Intercept)" || any(x[, 1L] != 1)) {
    stop("Internal error: the design matrix must have an intercept first.")
  }
  weighted_x <- x[weights > 0, , drop = FALSE]
  rank_x <- qr(weighted_x)$rank
  if (rank_x != ncol(x)) {
    stop("The X design matrix is rank-deficient. Remove redundant regressors before fitting.")
  }
  
  N <- nrow(y)
  J <- ncol(y)
  K <- ncol(x)
  d <- K * (J - 1L)  # Number of coefficients
  totals <- rowSums(y)
  nonbase <- seq.int(2L, J)
  blocks <- lapply(seq_len(J - 1L), function(h) ((h - 1L) * K + 1L):(h * K))
  
  # Quasi-log-likelihood contributions: N-vector.
  logLik_obs <- function(theta) {
    z <- softmax(theta, x, J)
    weights * (rowSums(y * z$eta) - totals * z$log_denom)
  }
  
  # Observation-wise weighted scores: N x d, in optimizer parameter order.
  score_obs <- function(theta) {
    P <- softmax(theta, x, J)$probabilities
    residual <- y[, nonbase, drop = FALSE] - sweep(P[, nonbase, drop = FALSE], 1L, totals, FUN = "*")
    do.call(cbind, lapply(seq_len(J - 1L), function(h) x * (weights * residual[, h])))
  }
  
  # Exact Hessian of the weighted objective, d x d (negative semidefinite).
  hessian <- function(theta) {
    P <- softmax(theta, x, J)$probabilities
    H <- matrix(0, d, d)
    for (a in seq_len(J - 1L)) {
      for (b in seq_len(J - 1L)) {
        curvature <- weights * totals * P[, a + 1L] * (as.numeric(a == b) - P[, b + 1L])
        H[blocks[[a]], blocks[[b]]] <- crossprod(x, x * curvature)
      }
    }
    -H
  }
  
  if (is.null(beta0)) {
    beta0 <- rep(0, (J-1)*K)
    beta1 <- initEst(y = y, x = x, weights = weights)
    if (sum(logLik_obs(beta1)) > sum(logLik_obs(beta0))) beta0 <- beta1
  } else if (!is.numeric(beta0) || length(beta0) != d || any(!is.finite(beta0))) {
    stop("'beta0' must contain (J - 1) * K finite numeric coefficients  in alternative-major order.")
  }
  
  param_names <- paste(rep(colnames(y)[nonbase], each = K),
                       rep(colnames(x), times = J - 1L), sep = ":")
  names(beta0) <- param_names
  
  opt <- maxLik::maxLik(logLik = logLik_obs, grad = score_obs, hess = hessian,
                        start = beta0, method = method,
                        control = list(iterlim = iterlim, reltol = reltol,
                                       printLevel = print.level), ...)
  theta_hat <- as.numeric(opt$estimate)
  names(theta_hat) <- param_names
  if (any(!is.finite(theta_hat))) stop("The optimiser returned non-finite coefficients.")
  if (!as.integer(opt$code) %in% c(0L, 1L, 2L, 8L)) {
    warning("The optimizer did not report ordinary convergence: ", opt$message)
  }
  
  # Robust joint sandwich covariance: (H^{-1}) (S'S) (H^{-1}).
  H <- -hessian(theta_hat)
  if (any(!is.finite(H)) || rcond(H) < 1e-12) {
    stop("The joint weighted Hessian is singular or poorly conditioned; robust covariance estimation is unavailable.")
  }
  H_inv <- solve(H)
  S <- score_obs(theta_hat)
  V <- H_inv %*% crossprod(S) %*% H_inv
  V <- (V + t(V)) / 2  # Remove floating-point asymmetry
  dimnames(V) <- list(param_names, param_names)
  
  B <- rbind(rep(0, K), matrix(theta_hat, nrow = J - 1L, ncol = K, byrow = TRUE))
  dimnames(B) <- list(colnames(y), colnames(x))
  
  vcov <- vector("list", J)
  names(vcov) <- colnames(y)
  vcov[[1L]] <- matrix(0, K, K, dimnames = list(colnames(x), colnames(x)))
  sigmat <- matrix(NA_real_, J - 1L, K)
  for (h in seq_len(J - 1L)) {
    vh <- V[blocks[[h]], blocks[[h]], drop = FALSE]
    dimnames(vh) <- list(colnames(x), colnames(x))
    vcov[[h + 1L]] <- vh
    sigmat[h, ] <- sqrt(pmax(0, diag(vh)))
  }
  
  tabout <- lapply(seq_len(J - 1L), function(h) {
    est <- B[h + 1L, ]
    se <- sigmat[h, ]
    z <- est / se
    result <- cbind(
      Estimate = est,
      `Std. Error` = se,
      `z value` = z,
      `Pr(>|z|)` = 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    )
    rownames(result) <- colnames(x)
    result
  })
  names(tabout) <- colnames(y)[nonbase]
  
  fitted_values <- softmax(theta_hat, x, J)$probabilities
  dimnames(fitted_values) <- list(rownames(y), colnames(y))
  
  structure(list(
    raw.estimates = theta_hat,
    estimates = tabout,
    init.val = beta0,
    baseline = colnames(y)[1L],
    likelihood = opt$maximum,
    conv_code = as.integer(opt$code),
    convergence = paste(opt$type, paste(opt$iterations, "iterations"), opt$message, sep = ", "),
    count = c(Obs = N, Explanatories = K - 1L, Choices = J),
    y = y,
    x = x,
    weights = weights,
    rowNo = rowNo,
    coefficient = B,
    vcov = vcov,
    vcov_full = V,
    fitted.values = fitted_values
  ), class = "fmlogit")
}

#' Extract the full joint covariance matrix of an fmlogit fit
#'
#' @param object An object of class `"fmlogit"`.
#' @param ... Unused.
#' @return Full robust joint covariance matrix, excluding normalized baseline
#'   coefficients and retaining cross-outcome covariance blocks.
#' @export
vcov.fmlogit <- function(object, ...) {
  object$vcov_full
}
