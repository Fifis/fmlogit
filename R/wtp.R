#' Willingness to pay for fmlogit models
#'
#' Calculates the willingness to pay for fractional multinomial logit models.
#'
#' @param object An "fmlogit.margins" object.
#' @param wtp.vec A 1*J vector that contains the willingness to pay for each choice j.
#' @param varlist A string vector which provides the name of variables to calculate
#' the wtp. If missing, all variables in object will be calculated.
#' @param indv.obs Also return observation-level effects? Requires margins
#'   computed with \code{marg.type = "aveacr"}.
#' @return An object of class \code{"fmlogit.wtp"} containing \code{wtp}, the
#'   aggregate effects and, when available, approximate inference; and optionally
#'   \code{wtp.obs}, the observation-level aggregate effects.
#' @details This function calculates the aggregate effect of a variable on the
#' "willingness to pay" by linearly multiplying the average partial effect with ex-ante (arbitary)
#' willingness to pay numbers associated with each choice.
#'
#' Suppose there are three choices A,B,C, each with a willingness to pay (or cost, profit, budget),
#' of 100, 200, and 300. The discrete effect of variable x on A,B and C are 0.5, 0.5, and -1, with
#' standard error 0.2, 0.3 and 0.5. The aggregated discrete effect of x on the total willingness
#' to pay (or cost), is thus 100*0.5 + 200*0.5 + 300*(-1) = -150. And the standard error can be also
#' calculated to be 162.8 under the approximation that outcome effects are independent.
#' This legacy approximation omits covariance between outcome effects, which
#' is not retained in \code{fmlogit.margins} objects. The resulting standard errors,
#' z statistics and p-values do not use the full joint uncertainty.
#'
#' If the input object has no standard errors, only aggregate point effects are returned.
#' @examples
#' # results1 <- fmlogit(y, x)
#' # effects1 <- effects(results1,effect="marginal",se=T)
#' # assume that the WTP = 1,2,3,...J for each choice j.
#' wtp(effects1, seq(1:nrow(effects1$effects)))
#' @export wtp

wtp <- function(object, wtp.vec, varlist = NULL, indv.obs = FALSE) {
  if (!inherits(object, "fmlogit.margins")) stop("'object' must be an fmlogit.margins object.")
  j <- nrow(object$effects)
  xnames <- colnames(object$effects)
  if (length(varlist) == 0) varlist <- xnames
  var.col.no <- match(varlist, xnames)
  if (!is.character(varlist) || anyNA(var.col.no) || anyDuplicated(varlist)) {
    stop("'varlist' must contain distinct names from 'object$effects'.")
  }
  k <- length(var.col.no)
  if (!is.numeric(wtp.vec) || !is.null(dim(wtp.vec)) ||
      length(wtp.vec) != j || any(!is.finite(wtp.vec))) {
    stop("'wtp.vec' must be a finite numeric vector with one value per outcome.")
  }
  if (indv.obs && is.null(object$marg.list)) {
    stop("'indv.obs = TRUE' requires margins computed with 'marg.type = \"aveacr\"'.")
  }
  # wtp calcs
  betamat <- object$effects[, var.col.no, drop = FALSE]
  wtp_mean <- wtp.vec %*% betamat

  if (object$R > 0) {
    semat <- object$se[, var.col.no, drop = FALSE]
    # Legacy approximation: Var(sum_j w_j m_jk) ~= sum_j w_j^2 Var(m_jk).
    wtp_se <- sqrt(wtp.vec^2 %*% semat^2)
    # output tables
    tabout <- matrix(ncol = 4, nrow = k)
    tabout[, 1] <- wtp_mean
    tabout[, 2] <- wtp_se
    tabout[, 3] <- tabout[, 1] / tabout[, 2]
    tabout[, 4] <- 2 * stats::pnorm(abs(tabout[, 3]), lower.tail = FALSE)
    colnames(tabout) <- c("estimate", "std", "z", "p-value")
    rownames(tabout) <- varlist
  } else tabout <- wtp_mean
  if (indv.obs) {
    wtp_mat <- matrix(ncol = k, nrow = nrow(object$marg.list[[1]]))
    for (i in seq_len(k)) {
      wtp_mat[, i] <- object$marg.list[[var.col.no[i]]] %*% wtp.vec
    }
    colnames(wtp_mat) <- varlist
  }
  # output list
  outlist <- list()
  outlist$wtp <- tabout
  if (indv.obs) outlist$wtp.obs <- wtp_mat
  return(structure(outlist, class = "fmlogit.wtp"))
}
