#' Partial effects of the covariates
#'
#' @description
#' Calculates marginal derivatives or discrete changes in the fitted shares of
#' a fractional multinomial logit model, at the mean design or averaged over
#' the retained observations.
#'
#' @param object An object of class \code{"fmlogit"}.
#' @param effect Either \code{"marginal"} for derivatives or \code{"discrete"}
#'   for changes from the minimum to the maximum among positive-weight rows.
#' @param marg.type Either \code{"atmean"} (the default) to evaluate effects at
#'   the weighted mean design, or \code{"aveacr"} for weighted average effects.
#' @param se Logical; compute simulation-based Krinsky-Robb standard errors?
#' @param varlist Character vector of column names in \code{object$x}, in the
#'   desired output order. By default, all columns except the first, intercept
#'   column are selected. Use \code{"(Intercept)"} to select the intercept.
#' @param at Optional finite numeric vector of length \code{K - 1}, giving the
#'   non-intercept design values in the order of \code{colnames(object$x)[-1]}.
#'   Only supported with \code{marg.type = "atmean"}.
#' @param R Number of joint coefficient draws for Krinsky-Robb standard errors;
#'   an integer of at least two when \code{se = TRUE}. Defaults to 1000.
#' @param ... Unused additional arguments.
#'
#' @details
#' With the first outcome normalized to \eqn{\beta_1 = 0}{beta_1 = 0}, the
#' conditional mean shares and their marginal derivatives are
#' \deqn{p_{ij} = \frac{\exp(x_i' \beta_j)}{\sum_h \exp(x_i' \beta_h)},}
#' \deqn{\frac{\partial p_{ij}}{\partial x_{ik}} =
#'   p_{ij}\left(\beta_{jk} - \sum_h p_{ih}\beta_{hk}\right).}
#' Discrete effects replace design column \eqn{k} by its observed minimum and
#' maximum, keeping all other columns fixed, and subtract the resulting shares.
#' Effects sum to zero across outcomes because the fitted shares sum to one.
#'
#' Both summaries inherit \code{object$weights}, normalised to sum to one.
#' For \code{"atmean"}, the reference profile uses weighted design means,
#' unless \code{at} supplies their values. For \code{"aveacr"}, the effects are
#' computed for each retained observation and then averaged with these weights.
#' Discrete endpoints are the minimum and maximum among positive-weight rows.
#' Individual effects remain unweighted values, including on zero-weight rows.
#' Effects refer to processed design columns. For factor indicators and
#' interactions, changing one column does not reconstruct the original data
#' or update related design columns.
#'
#' The simulation-based Krinsky-Robb method draws the entire non-baseline
#' coefficient vector from its estimated joint normal distribution,
#' \eqn{\theta^{(r)} \sim N(\hat\theta, \widehat{\mathrm{Var}}(\hat\theta))}.
#' Each draw preserves the baseline normalization and recomputes the requested
#' effect using the same design values and observation weights. The ordinary,
#' unweighted sample covariance across draws estimates effect uncertainty;
#' simulation draws do not carry survey weights. This uses
#' \code{object$vcov_full}, including covariance between outcomes and covariates,
#' and conditions on the observed design and weights. This is coefficient
#' uncertainty, not the weighted dispersion of individual effects. It does not
#' include uncertainty in the covariate distribution or account for survey
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
#'   \code{se = TRUE}. Used for joint aggregate inference in \code{wtp()}.}
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
  # theta_hat stacks K coefficients per non-baseline outcome (j = 2, ..., J).
  x <- object$x
  y <- object$y
  N <- nrow(y)
  J <- ncol(y)
  K <- ncol(x)
  d <- K * (J - 1L)
  theta_hat <- object$raw.estimates
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != N || K < 1L ||
      any(!is.finite(x)) || any(x[, 1L] != 1) ||
      !is.numeric(theta_hat) || length(theta_hat) != d ||
      any(!is.finite(theta_hat))) {
    stop("'object' must contain the current fmlogit design and coefficients.")
  }
  weights <- fit_weights(object)
  positive <- weights > 0

  xnames <- colnames(x)
  ynames <- colnames(y)
  if (is.null(varlist)) varlist <- xnames[-1L]
  if (!is.character(varlist) || !length(varlist) ||
      anyNA(varlist) || anyDuplicated(varlist)) {
    stop("'varlist' must select distinct covariate names; use '(Intercept)' for an intercept-only model.")
  }
  var_colNo <- match(varlist, xnames)
  if (anyNA(var_colNo)) stop("Unrecognized covariate names in 'varlist'.")
  n_vars <- length(var_colNo)

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
  x_eval <- x
  if (marg.type == "atmean") {
    if (is.null(at)) at <- weighted_colMeans(x[, -1L, drop = FALSE], weights)
    if (!is.numeric(at) || !is.null(dim(at)) || length(at) != K - 1L ||
        any(!is.finite(at))) {
      stop("'at' must be a finite numeric vector of length K - 1, excluding the intercept.")
    }
    x_eval <- matrix(c(1, at), nrow = 1L, dimnames = list(NULL, xnames))
  }
  eval_weights <- if (marg.type == "aveacr") weights else 1

  # Evaluate the same statistic for the fitted coefficients and every KR draw.
  evaluate <- function(theta, keep = FALSE) {
    effects_mat <- matrix(NA_real_, J, n_vars, dimnames = list(ynames, varlist))
    marg_list <- if (keep) stats::setNames(vector("list", n_vars), varlist) else NULL
    if (effect == "marginal") {
      # p_ij = exp(x_i' beta_j) / sum_h exp(x_i' beta_h), beta_1 = 0.
      # Reuse the estimator's stable softmax rather than exponentiating directly.
      model <- softmax(theta, x_eval, J)
      p <- model$p
      beta <- model$beta
    }
    for (k in var_colNo) {
      if (effect == "marginal") {
        # dp_ij/dx_ik = p_ij (beta_jk - sum_h p_ih beta_hk).
        beta_k <- beta[, k]
        beta_bar <- as.vector(p %*% beta_k)
        marg_mat <- sweep(p, 2L, beta_k, FUN = "*") - p * beta_bar
      } else {
        # Use endpoints of the support with positive observation weight.
        # Delta_ijk = p_j(x_i with x_ik = max x_k) - p_j(x_i with x_ik = min x_k).
        x_min <- x_max <- x_eval
        x_min[, k] <- min(x[positive, k])
        x_max[, k] <- max(x[positive, k])
        marg_mat <- softmax(theta, x_max, J)$p - softmax(theta, x_min, J)$p
      }
      dimnames(marg_mat) <- list(rownames(x_eval), ynames)
      # AME_jk = sum_i w_i m_ijk / sum_i w_i; at a profile there is one row.
      effects_mat[, xnames[k]] <- weighted_colMeans(marg_mat, eval_weights)
      if (keep) marg_list[[xnames[k]]] <- marg_mat
    }
    list(effects = effects_mat, marg.list = marg_list)
  }

  result <- evaluate(theta_hat, keep = marg.type == "aveacr")
  outlist <- list(effects = result$effects, weights = weights)
  if (se) {
    vcov_full <- object$vcov_full
    if (!is.matrix(vcov_full) || !is.numeric(vcov_full) ||
        !identical(dim(vcov_full), c(d, d)) || any(!is.finite(vcov_full)) ||
        !isSymmetric(vcov_full, tol = 1e-8)) {
      stop("'object$vcov_full' must be a finite, symmetric d-by-d joint covariance matrix.")
    }
    # A symmetric eigen square root also permits singular HC0 covariance.
    # Clip only round-off negatives; materially negative eigenvalues are invalid.
    vcov_full <- (vcov_full + t(vcov_full)) / 2
    eigen_vcov <- eigen(vcov_full, symmetric = TRUE)
    tol <- 1e-8 * max(abs(eigen_vcov$values))
    if (any(eigen_vcov$values < -tol)) {
      stop("The joint covariance matrix must be positive semidefinite.")
    }
    root <- sweep(eigen_vcov$vectors, 2L, sqrt(pmax(eigen_vcov$values, 0)), FUN = "*")

    # Simulation-based Krinsky-Robb inference:
    # theta^(r) = theta_hat + root z^(r), z^(r) ~ N(0, I_d).
    # Joint draws retain all cross-coefficient covariance and beta_1 = 0.
    theta_draws <- sweep(matrix(stats::rnorm(R * d), nrow = R) %*% t(root),
                         2L, theta_hat, FUN = "+")
    marg_draws <- matrix(NA_real_, R, length(outlist$effects))
    for (r in seq_len(R)) {
      marg_draws[r, ] <- evaluate(theta_draws[r, ])$effects
    }
    # Each draw has already used the observation weights. Across IID draws,
    # use ordinary covariance (R - 1 denominator), retaining outcome covariance.
    cov_list <- lapply(seq_len(n_vars), function(k) {
      block <- seq_len(J) + (k - 1L) * J
      vcov_k <- stats::cov(marg_draws[, block, drop = FALSE])
      dimnames(vcov_k) <- list(ynames, ynames)
      vcov_k
    })
    names(cov_list) <- varlist
    se_mat <- vapply(cov_list, function(vcov_k) sqrt(pmax(0, diag(vcov_k))), numeric(J))
    dimnames(se_mat) <- dimnames(outlist$effects)
    tabout <- lapply(varlist, function(variable) {
      wald_table(outlist$effects[, variable], se_mat[, variable])
    })
    names(tabout) <- varlist
    outlist$se <- se_mat
    outlist$ztable <- tabout
    outlist$cov.list <- cov_list
  }
  if (marg.type == "aveacr") outlist$marg.list <- result$marg.list
  outlist$R <- if (se) R else 0
  outlist$expl <- paste(
    effect, "effect", marg.type.out,
    if (se) "Krinsky-Robb standard error calculated" else "standard error not computed"
  )
  structure(outlist, class = "fmlogit.margins")
}
