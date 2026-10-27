#' User-valued aggregates of partial effects
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
#' @param indv.obs Also return observation-level effects? Requires margins
#'   computed with \code{marg.type = "aveacr"}.
#' @return An object of class \code{"fmlogit.wtp"} containing \code{wtp}, the
#'   aggregate effects and, when available, approximate inference; and optionally
#'   \code{wtp.obs}, the observation-level aggregate effects.
#' @details For outcome effects \eqn{m_{jk}} and supplied values \eqn{v_j},
#' the aggregate effect of covariate \eqn{k} is
#' \eqn{a_k = \sum_j v_j m_{jk}}. Effects can be derivatives or discrete
#' changes, evaluated at a reference profile or averaged across observations.
#' Units depend on the supplied values. These values are treated as fixed.
#'
#' When effect SEs are available, a legacy independence approximation uses
#' \eqn{\mathrm{Var}(a_k) \approx \sum_j v_j^2\mathrm{Var}(m_{jk})}.
#' It omits covariance between outcome effects, which is not retained in
#' \code{fmlogit.margins} objects. The resulting SEs, z statistics and p-values
#' do not use the full joint uncertainty. Plotting uses this same calculation.
#'
#' If the input object has no standard errors, only aggregate point effects are returned.
#' @examples
#' x <- data.frame(income = seq(-1, 1, length.out = 30))
#' y <- cbind(base = 0.3 + 0.05 * x$income,
#'            other = 0.7 - 0.05 * x$income)
#' fit <- fmlogit(y, x, print.level = 0)
#' margins <- effects.fmlogit(fit, marg.type = "aveacr")
#' wtp(margins, wtp.vec = c(0, 100), indv.obs = TRUE)
#' @export wtp

wtp <- function(object, wtp.vec, varlist = NULL, indv.obs = FALSE) {
  if (!inherits(object, "fmlogit.margins")) stop("'object' must be an fmlogit.margins object.")
  J <- nrow(object$effects)
  xnames <- colnames(object$effects)
  if (length(varlist) == 0L) varlist <- xnames
  var_colNo <- match(varlist, xnames)
  if (!is.character(varlist) || !length(varlist) ||
      anyNA(var_colNo) || anyDuplicated(varlist)) {
    stop("'varlist' must contain distinct names from 'object$effects'.")
  }
  n_vars <- length(var_colNo)
  if (!is.numeric(wtp.vec) || !is.null(dim(wtp.vec)) ||
      length(wtp.vec) != J || any(!is.finite(wtp.vec))) {
    stop("'wtp.vec' must be a finite numeric vector with one value per outcome.")
  }
  if (indv.obs && is.null(object$marg.list)) {
    stop("'indv.obs = TRUE' requires margins computed with 'marg.type = \"aveacr\"'.")
  }
  effects_mat <- object$effects[, var_colNo, drop = FALSE]
  wtp_mean <- wtp.vec %*% effects_mat

  if (object$R > 0) {
    se_mat <- object$se[, var_colNo, drop = FALSE]
    # Legacy approximation: Var(sum_j v_j m_jk) ~= sum_j v_j^2 Var(m_jk).
    wtp_se <- sqrt(wtp.vec^2 %*% se_mat^2)
    tabout <- wald_table(as.vector(wtp_mean), as.vector(wtp_se))
    rownames(tabout) <- varlist
  } else tabout <- wtp_mean
  if (indv.obs) {
    N <- nrow(object$marg.list[[1L]])
    wtp_mat <- matrix(NA_real_, N, n_vars, dimnames = list(NULL, varlist))
    for (k in seq_len(n_vars)) {
      wtp_mat[, k] <- object$marg.list[[var_colNo[k]]] %*% wtp.vec
    }
  }
  outlist <- list(wtp = tabout)
  if (indv.obs) outlist$wtp.obs <- wtp_mat
  return(structure(outlist, class = "fmlogit.wtp"))
}
