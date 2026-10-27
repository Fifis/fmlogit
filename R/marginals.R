#' Partial effects of the covariates
#'
#' @description
#' Calculates marginal derivatives or discrete changes in the fitted shares of
#' a fractional multinomial logit model, either at the weighted mean of the
#' design values or averaged over the retained observations.
#'
#' @param object An object of class \code{"fmlogit"}.
#' @param effect Either \code{"marginal"} for derivatives or \code{"discrete"}
#'   for changes from the minimum to the maximum among positive-weight rows.
#' @param marg.type Either \code{"atmean"} (the default) to evaluate effects at
#'   the weighted mean of the design values, or \code{"aveacr"} for weighted
#'   average effects.
#' @param se Logical; compute simulation-based Krinsky-Robb standard errors?
#' @param varlist Character vector of column names in \code{object$x}, in the
#'   desired output order. By default, all columns except the first, intercept
#'   column are selected. Use \code{"(Intercept)"} to select the intercept.
#' @param at Optional finite numeric vector of length \code{K - 1}, giving the
#'   non-intercept design values in the order of \code{colnames(object$x)[-1]}.
#'   Available only with \code{marg.type = "atmean"}.
#' @param R Number of joint coefficient draws for Krinsky-Robb standard errors;
#'   an integer of at least two when \code{se = TRUE}. Defaults to 1000.
#' @param ... Unused additional arguments.
#'
#' @details
#' With the first outcome normalised to \eqn{\beta_1 = 0}{beta_1 = 0}, the
#' conditional mean shares and their marginal derivatives are
#' \deqn{p_{ij} = \frac{\exp(x_i' \beta_j)}{\sum_h \exp(x_i' \beta_h)},}
#' \deqn{\frac{\partial p_{ij}}{\partial x_{ik}} =
#'   p_{ij}\left(\beta_{jk} - \sum_h p_{ih}\beta_{hk}\right).}
#' Discrete effects replace design column \eqn{k} by its observed minimum and
#' maximum, keeping all other columns fixed, and take the difference between
#' the resulting shares.
#' Effects sum to zero across outcomes because the fitted shares sum to one.
#'
#' Both summaries inherit \code{object$weights}, normalised to sum to one.
#' For \code{"atmean"}, the reference profile uses weighted design means,
#' unless \code{at} supplies their values. For \code{"aveacr"}, the effects are
#' computed for each retained observation and then averaged with these weights.
#' Discrete endpoints are the minimum and maximum among positive-weight rows.
#' Individual effects remain unweighted, including those for zero-weight rows.
#' Effects refer to processed design columns. For factor indicators and
#' interactions, changing one column does not reconstruct the original data
#' or update related design columns.
#'
#' The simulation-based Krinsky-Robb method draws the entire non-baseline
#' coefficient vector from its estimated joint normal distribution,
#' \eqn{\theta^{(r)} \sim N(\hat\theta, \widehat{\mathrm{Var}}(\hat\theta))}.
#' Each draw preserves the baseline normalisation and recomputes the requested
#' effect using the same design values and observation weights. The ordinary,
#' unweighted sample covariance across draws estimates the uncertainty in the
#' effects; simulation draws do not carry survey weights. This uses
#' \code{object$vcov.full}, including covariance between outcomes and covariates,
#' and treats the observed design and weights as fixed. This measures
#' uncertainty arising from coefficient estimation, rather than the weighted
#' dispersion of individual effects. It does not include uncertainty in the
#' covariate distribution or account for survey
#' clusters, strata, or finite-population corrections. Use \code{set.seed()}
#' for reproducible simulations. The reported z statistics and two-sided
#' p-values use a normal approximation.
#'
#' @return An object of class \code{"fmlogit.margins"} with:
#' \describe{
#' \item{effects}{A \code{J}-by-selected-covariates matrix of effects.}
#' \item{se}{A matching matrix of standard errors, present if \code{se = TRUE}.}
#' \item{ztable}{A named list of estimate, standard error, z statistic and p-value
#'   tables, present if \code{se = TRUE}.}
#' \item{marg.list}{A named list of \code{N}-by-\code{J} observation-level effect
#'   matrices, present for \code{marg.type = "aveacr"}.}
#' \item{weights}{The retained observation weights inherited from the fit.}
#' \item{cov.list}{A named list of \code{J}-by-\code{J} covariance matrices
#'   between outcome effects for each selected covariate, present if
#'   \code{se = TRUE}. These are used for inference on aggregate effects in
#'   \code{wtp()}.}
#' \item{R}{The number of simulation draws, or zero if \code{se = FALSE}.}
#' \item{expl}{A brief description of the effects and inference.}
#' }
#'
#' @examples
#' x <- data.frame(income = seq(-1, 1, length.out = 30))
#' y <- cbind(base = 0.3 + 0.05 * x$income,
#'            other = 0.7 - 0.05 * x$income)
#' fit <- fmlogit(y, x, print.level = 0)
#' effects.fmlogit(fit)
#' effects.fmlogit(fit, effect = "discrete", varlist = "income")
#'
#' @export effects.fmlogit
#' @exportS3Method
effects.fmlogit <- function(object, effect = c("marginal", "discrete"),
                           marg.type = "atmean", se = FALSE,
                           varlist = NULL, at = NULL, R = 1000, ...) {
  if (!inherits(object, "fmlogit")) stop("'object' must be an fmlogit fit.")
  effect <- match.arg(effect)
  marg.type <- match.arg(marg.type, c("atmean", "aveacr"))
  if (!is.logical(se) || length(se) != 1L || is.na(se)) {
    stop("'se' must be TRUE or FALSE.")
  }
  if (se && (!is.numeric(R) || length(R) != 1L || !is.finite(R) ||
             R < 2 || R != floor(R) || R > .Machine$integer.max)) {
    stop("'R' must be an integer of at least two when 'se = TRUE'.")
  }

  # Match fmlogit(): x is N x K with the intercept first; y is N x J.
  # theta.hat stacks K coefficients for each non-baseline outcome (j = 2, ..., J).
  x <- object$x
  y <- object$y
  N <- nrow(y)
  J <- ncol(y)
  K <- ncol(x)
  d <- K*(J-1L)
  theta.hat <- object$raw.estimates
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != N || K < 1L ||
      any(!is.finite(x)) || any(x[, 1L] != 1) ||
      !is.numeric(theta.hat) || length(theta.hat) != d ||
      any(!is.finite(theta.hat))) {
    stop("'object' must contain the current fmlogit design and coefficients.")
  }
  weights <- fitWeights(object)
  positive <- weights > 0

  xnames <- colnames(x)
  ynames <- colnames(y)
  if (is.null(varlist)) varlist <- xnames[-1L]
  if (!is.character(varlist) || !length(varlist) ||
      anyNA(varlist) || anyDuplicated(varlist)) {
    stop("'varlist' must select distinct covariate names; use '(Intercept)' for an intercept-only model.")
  }
  var.col.no <- match(varlist, xnames)
  if (anyNA(var.col.no)) stop("Unrecognised covariate names in 'varlist'.")
  n.vars <- length(var.col.no)

  if (!is.null(at) && marg.type != "atmean") {
    stop("'at' is only supported for 'marg.type = \"atmean\"'.")
  }
  marg.type.out <- if (marg.type == "aveacr") {
    "weighted average across observations,"
  } else if (is.null(at)) {
    "at the weighted mean,"
  } else {
    "at the evaluation point,"
  }
  x.eval <- x
  if (marg.type == "atmean") {
    if (is.null(at)) at <- weightedColMeans(x[, -1L, drop = FALSE], weights)
    if (!is.numeric(at) || !is.null(dim(at)) || length(at) != K-1L ||
        any(!is.finite(at))) {
      stop("'at' must be a finite numeric vector of length K - 1, excluding the intercept.")
    }
    x.eval <- matrix(c(1, at), nrow = 1L, dimnames = list(NULL, xnames))
  }
  eval.weights <- if (marg.type == "aveacr") weights else 1

  # Evaluate the same statistic for the fitted coefficients and each simulated draw.
  evaluate <- function(theta, keep = FALSE) {
    effects.mat <- matrix(NA_real_, nrow = J, ncol = n.vars,
                          dimnames = list(ynames, varlist))
    marg.list <- if (keep) stats::setNames(vector("list", n.vars), varlist) else NULL
    if (effect == "marginal") {
      # p_ij = exp(x_i' beta_j) / sum_h exp(x_i' beta_h), beta_1 = 0.
      # Use the estimator's softmax calculation to avoid numerical overflow.
      model <- softmax(theta, x.eval, J)
      p <- model$p
      beta <- model$beta
    }
    for (k in var.col.no) {
      if (effect == "marginal") {
        # dp_ij/dx_ik = p_ij (beta_jk - sum_h p_ih beta_hk).
        beta.k <- beta[, k]
        beta.bar <- as.vector(p %*% beta.k)
        marg.mat <- sweep(p, 2L, beta.k, FUN = "*") - p*beta.bar
      } else {
        # Use the minimum and maximum among observations with positive weights.
        # Delta_ijk = p_j(x_i with x_ik = max x_k) - p_j(x_i with x_ik = min x_k).
        x.min <- x.max <- x.eval
        x.min[, k] <- min(x[positive, k])
        x.max[, k] <- max(x[positive, k])
        marg.mat <- softmax(theta, x.max, J)$p - softmax(theta, x.min, J)$p
      }
      dimnames(marg.mat) <- list(rownames(x.eval), ynames)
      # AME_jk = sum_i w_i m_ijk / sum_i w_i; a single profile has one row.
      effects.mat[, xnames[k]] <- weightedColMeans(marg.mat, eval.weights)
      if (keep) marg.list[[xnames[k]]] <- marg.mat
    }
    list(effects = effects.mat, marg.list = marg.list)
  }

  result <- evaluate(theta.hat, keep = marg.type == "aveacr")
  outlist <- list(effects = result$effects, weights = weights)
  if (se) {
    vcov.full <- object$vcov.full
    if (!is.matrix(vcov.full) || !is.numeric(vcov.full) ||
        !identical(dim(vcov.full), c(d, d)) || any(!is.finite(vcov.full)) ||
        !isSymmetric(vcov.full, tol = 1e-8)) {
      stop("'object$vcov.full' must be a finite, symmetric d-by-d joint covariance matrix.")
    }
    # An eigenvalue decomposition also permits singular HC0 covariance matrices.
    # Set negative eigenvalues to zero only when they are small enough to reflect rounding.
    vcov.full <- (vcov.full + t(vcov.full))/2
    eigen.vcov <- eigen(vcov.full, symmetric = TRUE)
    tol <- 1e-8*max(abs(eigen.vcov$values))
    if (any(eigen.vcov$values < -tol)) {
      stop("The joint covariance matrix must be positive semidefinite.")
    }
    root <- sweep(eigen.vcov$vectors, 2L, sqrt(pmax(eigen.vcov$values, 0)), FUN = "*")

    # Simulation-based Krinsky-Robb inference:
    # theta^(r) = theta.hat + root z^(r), z^(r) ~ N(0, I_d).
    # Joint draws retain covariance between all coefficients and keep beta_1 = 0.
    theta.draws <- sweep(matrix(stats::rnorm(R*d), nrow = R) %*% t(root),
                         2L, theta.hat, FUN = "+")
    marg.draws <- matrix(NA_real_, nrow = R, ncol = length(outlist$effects))
    for (r in seq_len(R)) {
      marg.draws[r, ] <- evaluate(theta.draws[r, ])$effects
    }
    # Each draw has already used the observation weights. The draws are independent
    # and identically distributed, so use ordinary covariance (denominator R - 1).
    cov.list <- lapply(seq_len(n.vars), function(k) {
      block <- seq_len(J) + (k-1L)*J
      vcov.k <- stats::cov(marg.draws[, block, drop = FALSE])
      dimnames(vcov.k) <- list(ynames, ynames)
      vcov.k
    })
    names(cov.list) <- varlist
    se.mat <- vapply(cov.list, function(vcov.k) sqrt(pmax(0, diag(vcov.k))), numeric(J))
    dimnames(se.mat) <- dimnames(outlist$effects)
    tabout <- lapply(varlist, function(variable) {
      waldTable(outlist$effects[, variable], se.mat[, variable])
    })
    names(tabout) <- varlist
    outlist$se <- se.mat
    outlist$ztable <- tabout
    outlist$cov.list <- cov.list
  }
  if (marg.type == "aveacr") outlist$marg.list <- result$marg.list
  outlist$R <- if (se) R else 0
  outlist$expl <- paste(
    effect, "effect", marg.type.out,
    if (se) "Krinsky-Robb standard error calculated" else "standard error not computed"
  )
  structure(outlist, class = "fmlogit.margins")
}
