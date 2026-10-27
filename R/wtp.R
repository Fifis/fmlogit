#' Aggregates of partial effects using supplied values
#'
#' Combine outcome-specific partial effects using supplied values, such as
#' category costs or returns. The function name is retained for compatibility;
#' it does not estimate willingness to pay from preferences.
#'
#' @param object An object of class \code{"fmlogit.margins"}.
#' @param wtp.vec Finite numeric vector of length \code{J}, in the outcome order
#'   of \code{rownames(object$effects)}, with the baseline first.
#' @param varlist Covariate names to select, in the desired output order.
#'   Defaults to all columns of \code{object$effects}.
#' @param indv.obs Also return effects for each observation? Requires margins
#'   computed with \code{marg.type = "aveacr"}.
#' @return An object of class \code{"fmlogit.wtp"} containing \code{wtp}, the
#'   aggregate effects with standard errors and tests where available; and,
#'   optionally, \code{wtp.obs}, the aggregate effects for each observation.
#' @details For outcome effects \eqn{m_{jk}} and supplied values \eqn{v_j},
#' the aggregate effect of covariate \eqn{k} is
#' \eqn{a_k = \sum_j v_j m_{jk}}. Effects can be derivatives or discrete
#' changes, evaluated at a reference profile or averaged across observations.
#' The supplied values are treated as fixed and determine the units of the aggregate.
#'
#' Average effects already include the fitted observation weights through
#' \code{effects.fmlogit()}; no further weighting is needed here. Individual
#' effects retain one row per observation, including zero-weight observations.
#'
#' When effect SEs are available, the aggregate variance is
#' \eqn{\mathrm{Var}(a_k) = v'\mathrm{Cov}(m_k)v}, using the joint outcome-effect
#' covariance from Krinsky-Robb draws in \code{object$cov.list}. This retains
#' covariance between outcomes while holding the observed design and weights fixed.
#' Older margins objects without this covariance must be recomputed with
#' \code{effects.fmlogit(..., se = TRUE)}. Plotting uses this same calculation.
#' Adding a common constant to the supplied values leaves estimates and SEs
#' unchanged because the outcome effects sum to zero.
#'
#' If the input object has no standard errors, only aggregate point effects are returned.
#' @examples
#' x <- data.frame(income = seq(-1, 1, length.out = 30))
#' y <- cbind(base = 0.3 + 0.05*x$income,
#'            other = 0.7 - 0.05*x$income)
#' fit <- fmlogit(y, x, print.level = 0)
#' margins <- effects.fmlogit(fit, marg.type = "aveacr")
#' wtp(margins, wtp.vec = c(0, 100), indv.obs = TRUE)
#' @export wtp

wtp <- function(object, wtp.vec, varlist = NULL, indv.obs = FALSE) {
  if (!inherits(object, "fmlogit.margins")) stop("'object' must be an fmlogit.margins object.")
  J <- nrow(object$effects)
  xnames <- colnames(object$effects)
  if (length(varlist) == 0L) varlist <- xnames
  var.col.no <- match(varlist, xnames)
  if (!is.character(varlist) || !length(varlist) ||
      anyNA(var.col.no) || anyDuplicated(varlist)) {
    stop("'varlist' must contain distinct names from 'object$effects'.")
  }
  n.vars <- length(var.col.no)
  if (!is.numeric(wtp.vec) || !is.null(dim(wtp.vec)) ||
      length(wtp.vec) != J || any(!is.finite(wtp.vec))) {
    stop("'wtp.vec' must be a finite numeric vector with one value per outcome.")
  }
  if (indv.obs && is.null(object$marg.list)) {
    stop("'indv.obs = TRUE' requires margins computed with 'marg.type = \"aveacr\"'.")
  }
  # Only value differences matter because the outcome effects sum to zero.
  values <- wtp.vec - wtp.vec[1L]
  effects.mat <- object$effects[, var.col.no, drop = FALSE]
  wtp.mean <- values %*% effects.mat

  if (object$R > 0) {
    if (is.null(object$cov.list)) {
      stop("Outcome-effect covariance is unavailable; recompute margins with effects.fmlogit(..., se = TRUE).")
    }
    wtp.se <- numeric(n.vars)
    for (k in seq_len(n.vars)) {
      cov.mat <- object$cov.list[[varlist[k]]]
      if (!is.matrix(cov.mat) || !is.numeric(cov.mat) ||
          !identical(dim(cov.mat), c(J, J)) || any(!is.finite(cov.mat))) {
        stop("'object$cov.list' must contain a finite J-by-J covariance matrix for each selected variable.")
      }
      wtp.se[k] <- sqrt(max(0, as.numeric(crossprod(values, cov.mat %*% values))))
    }
    tabout <- waldTable(as.vector(wtp.mean), wtp.se)
    rownames(tabout) <- varlist
  } else tabout <- wtp.mean
  if (indv.obs) {
    N <- nrow(object$marg.list[[1L]])
    wtp.mat <- matrix(NA_real_, nrow = N, ncol = n.vars,
                      dimnames = list(NULL, varlist))
    for (k in seq_len(n.vars)) {
      wtp.mat[, k] <- object$marg.list[[var.col.no[k]]] %*% values
    }
  }
  outlist <- list(wtp = tabout)
  if (indv.obs) outlist$wtp.obs <- wtp.mat
  return(structure(outlist, class = "fmlogit.wtp"))
}
