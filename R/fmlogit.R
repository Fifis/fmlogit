#' Obtain starting values from pairwise fractional logits
#'
#' @description
#' Fits a separate fractional binary logit for each non-baseline outcome
#' against the first outcome, weighted by their combined share. Under correct
#' specification, they estimate the same conditional-mean
#' coefficients as the joint fractional multinomial model, though the
#' estimates generally differ in finite samples.
#'
#' @param y Numeric `N`-by-`J` matrix of observed shares, with `J >= 2`.
#'   Each row must sum to one.
#' @param x Numeric `N`-by-`K` *design matrix*, including a column of ones first.
#'   The columns must follow the same order as in joint estimation.
#' @param weights Optional vector of `N` non-negative observation weights.
#' @return Named numeric vector of length `(J-1)*K`. Coefficients are grouped
#'   by non-baseline outcome; the baseline coefficients are omitted.
#'
#' @details
#' For outcome `j > 1`, fit a fractional binary logit with response
#' \eqn{y_j / (y_1 + y_j)}{y_j / (y_1 + y_j)} and prior weights
#' \eqn{w (y_1 + y_j)}{w * (y_1 + y_j)}.
#' Observations with \eqn{y_1 + y_j = 0}{y_1 + y_j = 0} make no contribution.
#' In the event of failure or non-finite pairwise estimates, the corresponding
#' coefficients are initialised to zero, and a warning is issued.
#'
#' The first outcome has a zero coefficient vector, omitted from optimisation.
#' With this normalisation,
#' \eqn{\mathbb{E}(y_{ij}\mid x_i) = p_{ij}}{E(y_ij | x_i) = p_ij} and
#' \eqn{p_{ij}/p_{i1} = \exp (x_i' \beta_j)}{p_ij / p_i1 = exp(x_i' beta_j)}.
#' Now define \eqn{\pi_{ij} = \frac{p_{ij}}{p_{i1} + p_{ij}}}{pi_ij = p_ij / (p_i1 + p_ij)}.
#' Substituting the multinomial-logit conditional means gives
#' \eqn{\pi_{ij} = \frac{\exp(x_i' \beta_j)}{1 + \exp(x_i' \beta_j)} =
#' \Lambda(x_i' \beta_j)}{pi_ij = exp(x_i' beta_j) / (1 + exp(x_i' beta_j)) =
#' Lambda(x_i' beta_j)}.
#' This applies the fractional quasi-likelihood argument of
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
    stop("The first column of 'x' must contain ones for the intercept.")
  }
  if (any(y < 0 | y > 1) || any(abs(rowSums(y) - 1) > 1e-8)) {
    stop("'y' must contain shares in [0, 1], with rows summing to one.")
  }

  N <- nrow(y)
  J <- ncol(y)
  K <- ncol(x)
  if (is.null(weights)) weights <- rep(1, N)
  if (!is.numeric(weights) || length(weights) != N ||
      any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
    stop("'weights' must be non-negative, finite, aligned with 'y', and not all zero.")
  }
  # Rescaling the weights improves numerical stability without changing the estimates.
  weights <- weights / mean(weights)

  if (is.null(colnames(y))) colnames(y) <- paste0("share", seq_len(J))
  if (is.null(colnames(x))) {
    colnames(x) <- c("(Intercept)", if (K > 1L) paste0("X", seq_len(K-1L)))
  }

  xnames <- colnames(x)
  ynames <- colnames(y)
  beta <- matrix(0, nrow = J, ncol = K, dimnames = list(ynames, xnames))

  for (j in seq.int(2L, J)) {
    pair.total <- y[, 1L] + y[, j]
    prior.weights <- pair.total*weights
    keep <- prior.weights > 0

    if (!any(keep)) {
      warning(sprintf("No positive pairwise weight for '%s'; using zeros.", ynames[j]))
      next
    }

    response <- y[keep, j]/pair.total[keep]
    fit <- tryCatch(
      stats::glm.fit(
        x = x[keep, , drop = FALSE], y = response, weights = prior.weights[keep],
        family = stats::quasibinomial(link = "logit"), intercept = TRUE
      ),
      error = function(e) NULL
    )

    if (is.null(fit) || !isTRUE(fit$converged) ||
        length(fit$coefficients) != K || any(!is.finite(fit$coefficients))) {
      warning(sprintf("Pairwise GLM for '%s' failed; using zeros.", ynames[j]))
      next
    }

    beta[j, ] <- fit$coefficients
  }

  theta <- as.vector(t(beta[-1L, , drop = FALSE]))
  names(theta) <- paste(rep(ynames[-1L], each = K),
                        rep(xnames, times = J-1L), sep = ":")
  theta
}
#' Estimate a fractional multinomial logit model
#'
#' @description
#' Estimates a multinomial-logit conditional-mean model for compositional
#' outcomes by maximising the weighted multinomial *quasi*-log-likelihood.
#' The first outcome is the reference category. Robust standard errors are
#' computed using the full joint sandwich covariance matrix, with each
#' observation contributing a separate score.
#'
#' @param y Numeric `N`-by-`J` matrix or data frame of observed shares,
#'   `J >= 2`. Values must lie in [0, 1] and rows must sum to one.
#'   The first column is automatically treated as the baseline.
#' @param x `N`-row matrix or data frame of regressors. Numeric, factor, and
#'   character columns are allowed in data frames. An intercept is added
#'   automatically as the first column. Any supplied constant columns are
#'   removed before creating the design matrix. No variable selection is made.
#' @param weights Optional non-negative observation weights, length `N`. The
#'   default is equal weights. These relative weights multiply the individual
#'   quasi-log-likelihood contributions; rescaling them leaves estimates and
#'   the HC0 covariance unchanged. They are inherited by mean profiles and
#'   average partial effects.
#'   The weights are normalised to mean one for numerical stability;
#'   the object retains the original weights and quasi-log-likelihood scale.
#' @param beta0 Optional starting coefficients, length `(J-1)*K`, where `K` is
#'   the number of columns in the *processed* design matrix including the
#'   intercept. Coefficients are grouped by non-baseline outcome, with the
#'   intercept first in each block. If NULL, compare the pairwise fractional
#'   binary logits from [initEst()] with zero coefficients and use whichever
#'   gives the larger quasi-log-likelihood.
#' @param method Optimiser for [maxLik::maxLik()]; default `"BHHH"`. Use `"CG"`
#'   if estimation is slow.
#' @param iterlim Positive iteration limit passed to `maxLik()`.
#' @param reltol Positive relative convergence tolerance passed to `maxLik()`.
#' @param print.level Optimiser verbosity; default 1 (short output).
#' @param ... Additional arguments forwarded to [maxLik::maxLik()].
#'
#' @details
#' For the first, reference category, \eqn{\beta_1 = 0}{beta_1 = 0}.
#' The linear index is \eqn{\eta_{ij} = x_i'\beta_j}{eta_ij = x_i' beta_j}.
#' The conditional mean shares are
#' \eqn{p_{ij} = \exp(x_i' \beta_j) / \sum_h \exp(x_i' \beta_h)}{p_ij = exp(x_i' beta_j) / sum_h exp(x_i' beta_h)}.
#' Here, \eqn{p_{ij}} corresponds to \eqn{\xi_k(x;\beta)} in
#' \insertCite{mullahy2015multivariate}{fmlogit},
#' with his outcome index \eqn{k} relabelled \eqn{j}. He normalises the last
#' outcome instead; fitted shares and partial effects are invariant to the
#' choice of baseline.
#' The objective is \eqn{\sum_i w_i \sum_j y_{ij} \log p_{ij}}{sum_i w_i sum_j y_ij log(p_ij)},
#' as described by \insertCite{mullahy2015multivariate}{fmlogit}.
#' The joint HC0 covariance places the sum of outer products of the weighted
#' *full* scores between two copies of the inverse weighted negative Hessian.
#' For unweighted scores \eqn{s_i}, this is
#' \eqn{H^{-1}(\sum_i w_i^2 s_i s_i')H^{-1}}, where
#' \eqn{H = -\sum_i w_i\partial s_i/\partial\theta'}.
#' The squared weights arise because each score already contains its weight.
#' This sandwich estimator treats observations as independent: weights alone
#' do not account for survey clusters, strata, or finite-population corrections
#' \insertCite{papke1996econometric}{fmlogit}.
#' Rows with missing responses or covariates are removed before estimation;
#' the corresponding weights are removed with them.
#'
#' @return An object of class `"fmlogit"`. `raw.estimates` contains the stacked
#'   non-baseline coefficients; `estimates` contains their robust Wald tables.
#'   `baseline` names the first outcome, and `likelihood` is the maximised
#'   quasi-log-likelihood. `conv.code` and `convergence` describe optimisation.
#'   `count`, `y`, `x`, `weights`, and `row.no` describe the retained sample.
#'   `coefficient` is `J`-by-`K`, with an all-zero first row. `vcov.full` is
#'   the complete `K(J-1)`-by-`K(J-1)` joint HC0 covariance matrix.
#'   `vcov` is a legacy list of `J` diagonal blocks of `vcov.full`, with a zero
#'   block for the normalised baseline (cross-category blocks are only in
#'   `vcov.full`). `fitted.values` holds the `N`-by-`J` fitted shares.
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
  x.data <- as.data.frame(x, stringsAsFactors = FALSE)
  if (nrow(x.data) != N0) {
    stop("'x' and 'y' must have the same number of rows.")
  }
  if (is.null(colnames(y))) colnames(y) <- paste0("share", seq_len(ncol(y)))
  if (anyNA(colnames(y)) || anyDuplicated(colnames(y)) || any(!nzchar(colnames(y)))) {
    stop("'y' must have unique, non-empty column names.")
  }
  if (ncol(x.data) && (anyNA(names(x.data)) || anyDuplicated(names(x.data)) ||
                       any(!nzchar(names(x.data))))) {
    stop("'x' must have unique, non-empty column names.")
  }

  if (is.null(weights)) weights <- rep(1, N0)
  if (!is.numeric(weights) || length(weights) != N0 ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("'weights' must be a finite, non-negative numeric vector of length nrow(y).")
  }

  keep.rows <- stats::complete.cases(y, x.data)
  row.no <- which(keep.rows)
  if (!length(row.no)) stop("No complete observations remain.")
  y <- y[keep.rows, , drop = FALSE]
  x.data <- x.data[keep.rows, , drop = FALSE]
  weights <- weights[keep.rows]
  if (!any(weights > 0)) stop("All retained observation weights are zero.")
  weight.scale <- mean(weights)
  optim.weights <- weights / weight.scale

  if (any(!is.finite(y)) || any(y < 0 | y > 1) || any(abs(rowSums(y) - 1) > 1e-8)) {
    stop("'y' must contain finite shares in [0, 1], each row summing to one.")
  }

  # Construct a single design matrix for starting values and estimation.
  # The intercept occupies the first column.
  if (ncol(x.data) > 0L) {
    good.type <- vapply(x.data, function(v) {
      is.numeric(v) || is.factor(v) || is.character(v) || is.logical(v)
    }, logical(1))
    if (!all(good.type)) {
      stop("'x' contains unsupported column types in columns ", which(!good.type), ".")
    }
    x.data <- lapply(x.data, function(v) {
      if (is.factor(v)) droplevels(v) else v
    })
    x.data <- as.data.frame(x.data, check.names = FALSE)
    is.constant <- vapply(x.data, function(v) length(unique(v)) <= 1L, logical(1))
    x.data <- x.data[, !is.constant, drop = FALSE]
  }
  x <- if (ncol(x.data) > 0L) {
    stats::model.matrix(~ ., data = x.data)
  } else {
    matrix(1, nrow = nrow(y), ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  }
  if (any(!is.finite(x))) stop("'x' contains non-finite design values.")
  if (colnames(x)[1L] != "(Intercept)" || any(x[, 1L] != 1)) {
    stop("Internal error: the design matrix must have an intercept first.")
  }
  rank.x <- qr(x[weights > 0, , drop = FALSE])$rank
  if (rank.x != ncol(x)) {
    stop("The x design matrix is rank-deficient. Remove redundant regressors before fitting.")
  }

  N <- nrow(y)
  J <- ncol(y)
  K <- ncol(x)
  d <- K * (J-1L)  # Number of free coefficients, excluding the first, baseline outcome.
  xnames <- colnames(x)
  ynames <- colnames(y)
  totals <- rowSums(y)
  nonbase <- seq.int(2L, J)
  blocks <- lapply(nonbase, function(j) seq_len(K) + (j-2L)*K)

  # Each observation contributes one quasi-log-likelihood value: an N-vector.
  loglikObs <- function(theta) {
    model <- softmax(theta, x, J)
    optim.weights * (rowSums(y*model$eta) - totals*model$log.denom)
  }

  # Weighted scores psi_i = w_i s_i: N x d, in the optimiser's parameter order.
  scoreObs <- function(theta) {
    p <- softmax(theta, x, J)$p
    residual <- y - p*totals
    do.call(cbind, lapply(nonbase, function(j) x * (optim.weights*residual[, j])))
  }

  # Exact Hessian of the weighted objective: d x d, negative semidefinite.
  hessian <- function(theta) {
    p <- softmax(theta, x, J)$p
    hess <- matrix(0, nrow = d, ncol = d)
    for (j in nonbase) {
      for (h in nonbase) {
        curvature <- optim.weights*totals*p[, j] * (as.numeric(j == h) - p[, h])
        hess[blocks[[j-1L]], blocks[[h-1L]]] <- crossprod(x, x*curvature)
      }
    }
    -hess
  }

  if (is.null(beta0)) {
    theta.init <- rep(0, d)
    theta.pairwise <- initEst(y = y, x = x, weights = weights)
    if (sum(loglikObs(theta.pairwise)) > sum(loglikObs(theta.init))) theta.init <- theta.pairwise
  } else {
    if (!is.numeric(beta0) || length(beta0) != d || any(!is.finite(beta0))) {
      stop("'beta0' must contain (J-1)*K finite numeric coefficients, grouped by outcome.")
    }
    theta.init <- beta0
  }

  param.names <- paste(rep(ynames[nonbase], each = K),
                       rep(xnames, times = J-1L), sep = ":")
  names(theta.init) <- param.names

  opt <- maxLik::maxLik(logLik = loglikObs, grad = scoreObs, hess = hessian,
                        start = theta.init, method = method,
                        control = list(iterlim = iterlim, reltol = reltol,
                                       printLevel = print.level), ...)
  theta.hat <- as.numeric(opt$estimate)
  names(theta.hat) <- param.names
  if (any(!is.finite(theta.hat))) stop("The optimiser returned non-finite coefficients.")
  if (!as.integer(opt$code) %in% c(0L, 1L, 2L, 8L)) {
    warning("The optimiser did not report ordinary convergence: ", opt$message)
  }

  # Joint HC0: H^{-1} sum_i(w_i^2 s_i s_i') H^{-1}.
  # H uses w_i once; scores already include w_i, so their products use w_i^2.
  # Using optim.weights gives the same covariance: the common scale cancels.
  bread <- -hessian(theta.hat)
  if (any(!is.finite(bread)) || rcond(bread) < 1e-12) {
    stop("The joint weighted Hessian is singular or poorly conditioned; robust covariance estimation is unavailable.")
  }
  bread.inv <- solve(bread)
  scores <- scoreObs(theta.hat)
  vcov.full <- bread.inv %*% crossprod(scores) %*% bread.inv
  vcov.full <- (vcov.full + t(vcov.full))/2  # Remove asymmetry due to rounding.
  dimnames(vcov.full) <- list(param.names, param.names)

  model <- softmax(theta.hat, x, J)
  beta <- model$beta
  dimnames(beta) <- list(ynames, xnames)

  vcov <- vector("list", J)
  names(vcov) <- ynames
  vcov[[1L]] <- matrix(0, nrow = K, ncol = K, dimnames = list(xnames, xnames))
  se.mat <- matrix(0, nrow = J, ncol = K, dimnames = list(ynames, xnames))
  for (j in nonbase) {
    vcov.j <- vcov.full[blocks[[j-1L]], blocks[[j-1L]], drop = FALSE]
    dimnames(vcov.j) <- list(xnames, xnames)
    vcov[[j]] <- vcov.j
    se.mat[j, ] <- sqrt(pmax(0, diag(vcov.j)))
  }

  tabout <- lapply(nonbase, function(j) {
    result <- waldTable(beta[j, ], se.mat[j, ])
    colnames(result) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(result) <- xnames
    result
  })
  names(tabout) <- ynames[nonbase]

  p <- model$p
  dimnames(p) <- list(rownames(y), ynames)

  structure(list(
    raw.estimates = theta.hat,
    estimates = tabout,
    init.val = theta.init,
    baseline = ynames[1L],
    likelihood = opt$maximum*weight.scale,
    conv.code = as.integer(opt$code),
    convergence = paste(opt$type, paste(opt$iterations, "iterations"), opt$message, sep = ", "),
    count = c(Obs = N, Explanatories = K-1L, Choices = J),
    y = y,
    x = x,
    weights = weights,
    row.no = row.no,
    coefficient = beta,
    vcov = vcov,
    vcov.full = vcov.full,
    fitted.values = p
  ), class = "fmlogit")
}

#' Extract the full joint covariance matrix of an fmlogit fit
#'
#' @param object An object of class `"fmlogit"`.
#' @param ... Unused.
#' @return The full robust joint covariance matrix, excluding the normalised
#'   baseline coefficients and retaining covariance blocks between outcomes.
#' @export
vcov.fmlogit <- function(object, ...) {
  object$vcov.full
}
