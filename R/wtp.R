#' Effects on willingness to pay and other valued aggregates
#'
#' Combine outcome-specific partial effects using supplied values. When these
#' values represent willingness to pay (WTP), the result describes a change in
#' mean WTP. Costs, profits, or other values give the corresponding change in
#' those quantities.
#'
#' @param object An object of class \code{"fmlogit.margins"}.
#' @param wtp.vec Finite numeric vector of length \code{J}, in the outcome order
#'   of \code{rownames(object$effects)}, with the baseline first. Values are
#'   matched by position, even when the vector has names.
#' @param varlist Covariate names to select, in the desired output order.
#'   Defaults to all columns of \code{object$effects}.
#' @param indv.obs Also return effects for each observation? Requires margins
#'   computed with \code{marg.type = "aveacr"}.
#' @return An object of class \code{"fmlogit.wtp"} containing \code{wtp}, the
#'   aggregate effects with standard errors and tests where available; and,
#'   optionally, \code{wtp.obs}, the aggregate effects for each observation.
#' @details
#' Let \eqn{p_j(x)} be the fitted share of outcome \eqn{j}, and let \eqn{v_j}
#' be its supplied value. The mean value and the effect of covariate \eqn{k} are
#' \deqn{V(x) = \sum_j v_j p_j(x), \qquad a_k = \sum_j v_j m_{jk}.}
#' Here \eqn{m_{jk}} is a marginal derivative or a discrete change from
#' \code{effects.fmlogit()}. Thus \code{wtp()} returns an effect on \eqn{V},
#' rather than its level. The supplied values are held fixed: the calculation
#' measures changes arising from the redistribution of shares.
#'
#' Average effects already include the fitted observation weights through
#' \code{effects.fmlogit()}; no further weighting is needed here. Individual
#' effects retain one row per observation, including zero-weight observations.
#'
#' \subsection{When the result is an effect on willingness to pay}{
#' Suppose residents report the maximum annual amount they would pay for the
#' same park improvement, choosing one of three stated amounts: A, GBP100;
#' B, GBP200; or C, GBP300. Each observation records the proportions choosing
#' these amounts in a district. With \code{wtp.vec = c(100, 200, 300)},
#' \eqn{V(x)} is the fitted mean stated annual WTP per resident. An effect on
#' this mean is therefore an effect on WTP in a direct, literal sense.
#'
#' The original three-choice illustration considers share changes of
#' \eqn{(0.5, 0.5, -1)}. This is an idealised limiting case: complete
#' redistribution from C to an equal mixture of A and B. The aggregate change is
#' \deqn{100(0.5) + 200(0.5) + 300(-1) = -150.}
#' Mean annual WTP falls by GBP150 per resident; the result is a change, not
#' a negative WTP level. A smaller redistribution of \eqn{(0.05, 0.05, -0.10)}
#' gives a fall of GBP15. If the response categories are ranges rather than
#' exact amounts, assigning representative values gives an approximation to
#' mean WTP.
#'
#' Another valid application divides a fixed quantity among uses with constant
#' marginal WTP per unit. Under an additive monetary valuation, the same
#' calculation gives the change in WTP per unit of that quantity. The valuation
#' assumptions must be justified in the application. The function takes these
#' valuations as inputs; it does not infer them from price coefficients or
#' calculate the consumer-surplus measures of a utility-based choice model.
#' See Train (2009), Chapter 3, for those distinct calculations.
#' }
#'
#' \subsection{Another interpretation: expenditure per trip}{
#' Suppose the three outcomes are walking, bus travel, and taxi travel, with
#' fares of GBP0, GBP2, and GBP10. A change in their shares of
#' \eqn{(-0.10, 0.04, 0.06)} gives
#' \deqn{0(-0.10) + 2(0.04) + 10(0.06) = 0.68.}
#' Expected expenditure rises by GBP0.68 per trip. These supplied fares measure
#' expenditure; a WTP interpretation would require a separate valuation argument.
#' For a fixed total of \eqn{T} trips, multiply the result by \eqn{T} to obtain
#' the change in total expenditure. Changes in the number of trips or in the
#' fares themselves require additional calculations.
#' }
#'
#' \subsection{Reading and plotting the results}{
#' With standard errors, the returned \code{wtp} table has one row
#' per selected covariate. Its columns give the estimate, standard error,
#' z statistic, and two-sided p-value. Without standard errors, the point
#' estimates form one row, with covariates in columns.
#'
#' A marginal estimate is in value units per unit of the covariate; it is a
#' local derivative, not an exact finite change. A discrete estimate is the
#' value change for the contrast specified by \code{effects.fmlogit()}, such
#' as changing a binary indicator from zero to one. The result describes
#' either the weighted average effect or the effect at the chosen reference
#' profile, according to the input margins. A positive estimate means higher
#' mean WTP when values measure WTP, or greater expenditure when they measure
#' costs. The p-value tests whether the aggregate effect is zero; individual
#' outcome effects may be non-zero and cancel one another.
#'
#' \code{wtp.obs} contains effects for individual observations, not their WTP
#' levels. To obtain fitted levels of \eqn{V}, multiply the fitted shares by
#' the supplied value vector instead.
#'
#' For several covariates, show estimates as separate points with a zero
#' reference line and, where appropriate, approximate 95\% intervals formed
#' from the estimate plus or minus 1.96 standard errors. There is no natural
#' reason to join unrelated covariates by lines. To show how an effect varies
#' with a covariate, use \code{plot(fit, wtp.vec = ..., against = ...,
#' effect = "marginal", se = TRUE)}. This plots marginal effects at successive
#' profiles, holding other design columns at their weighted means; it does not
#' plot WTP levels or the average effect over the sample. Its bars extend one
#' standard error in each direction, rather than forming 95\% intervals.
#' }
#'
#' \subsection{Uncertainty}{
#' When effect SEs are available, the aggregate variance is
#' \eqn{\mathrm{Var}(a_k) = v'\mathrm{Cov}(m_k)v}, using the joint outcome-effect
#' covariance from Krinsky-Robb draws in \code{object$cov.list}. This retains
#' covariance between outcomes while holding the observed design and weights fixed.
#' Outcome effects sum to zero and should not be treated as independent when
#' combining their standard errors. Uncertainty in the supplied values is not
#' included in this calculation.
#' Older margins objects without this covariance must be recomputed with
#' \code{effects.fmlogit(..., se = TRUE)}. Plotting uses this same calculation.
#' Adding a common constant to the supplied values leaves estimates and SEs
#' unchanged because the outcome effects sum to zero.
#'
#' If the input object has no standard errors, only aggregate point effects are returned.
#' }
#' @seealso \code{\link{effects.fmlogit}}, \code{\link{plot.fmlogit}},
#'   \code{\link{summary.fmlogit}}.
#' @references
#' The original implementation and three-choice illustration:
#' \url{https://github.com/f1kidd/fmlogit/blob/master/R/wtp.R}.
#'
#' Train, K. E. (2009). \emph{Discrete Choice Methods with Simulation},
#' second edition, Chapter 3. Cambridge University Press.
#' \url{https://eml.berkeley.edu/books/choice2.html}.
#' @examples
#' # Artificial district shares choosing annual WTP amounts of GBP100, 200, 300.
#' x <- data.frame(information = seq(0, 1, length.out = 40))
#' y <- cbind(A = 0.2 + 0.05*x$information,
#'            B = 0.3 + 0.05*x$information,
#'            C = 0.5 - 0.10*x$information)
#' fit <- fmlogit(y, x, print.level = 0)
#' set.seed(1)
#' margins <- stats::effects(fit, marg.type = "aveacr", se = TRUE, R = 200)
#' wtp.values <- c(A = 100, B = 200, C = 300)
#' result <- wtp(margins, wtp.vec = wtp.values, indv.obs = TRUE)
#' summary(result)
#'
#' # Levels and effects are different quantities.
#' wtp.level <- drop(as.matrix(stats::fitted(fit)) %*% wtp.values)
#' wtp.level[1:6]
#' result$wtp.obs[1:6, , drop = FALSE]
#'
#' # Average effects with approximate 95% intervals.
#' tab <- result$wtp
#' k <- seq_len(nrow(tab))
#' lower <- tab[, "Estimate"] - 1.96*tab[, "Std. Error"]
#' upper <- tab[, "Estimate"] + 1.96*tab[, "Std. Error"]
#' graphics::plot(k, tab[, "Estimate"], ylim = range(0, lower, upper), xaxt = "n",
#'                xlab = "Covariate", ylab = "Effect on annual WTP (GBP)")
#' graphics::axis(1, at = k, labels = rownames(tab))
#' graphics::segments(k, lower, k, upper)
#' graphics::abline(h = 0, lty = 2)
#'
#' # A profile of marginal effects; these bars show plus/minus ONE standard error.
#' set.seed(2)
#' curves <- graphics::plot(fit, wtp.vec = wtp.values, varlist = "information",
#'                          against = "information", effect = "marginal",
#'                          t = 10, se = TRUE, R = 200,
#'                          xlab = "Information index",
#'                          ylab = "Marginal effect on annual WTP (GBP)")
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
