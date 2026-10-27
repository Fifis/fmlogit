#' Extract fitted shares, residuals, and predictions
#'
#' @name fitted.fmlogit
#' @aliases residuals.fmlogit predict.fmlogit
#' @description
#' Extract fitted conditional shares, response residuals, or predictions from
#' a fractional multinomial logit model.
#'
#' @param object An object of class \code{"fmlogit"}.
#' @param newdata Optional finite numeric matrix or data frame containing the
#'   processed design columns in \code{colnames(object$x)[-1]} order, excluding
#'   the intercept. A numeric vector supplies one observation. For several
#'   observations of a single regressor, use a one-column matrix or data frame.
#'   Factor indicators and interactions must already be encoded as in
#'   \code{object$x}. If \code{NULL}, use the retained estimation observations.
#' @param newbeta Optional finite \code{J}-by-\code{K} coefficient matrix, with
#'   outcomes and design columns in the same order as \code{object$coefficient},
#'   including the intercept column and baseline row. Applies even when
#'   \code{newdata = NULL}.
#' @param ... Unused additional arguments.
#'
#' @details
#' Predictions use \eqn{p_{ij} = \exp(x_i' \beta_j) /
#'   \sum_h \exp(x_i' \beta_h)}, evaluated with the estimator's stable softmax.
#' The intercept is added as the first design column. Response residuals are
#' \eqn{y_{ij} - p_{ij}}; their rows sum to zero.
#' The first outcome is the baseline. A common coefficient shift in
#' \code{newbeta} is removed so that its first row is zero before prediction.
#' Fitted values and residuals cover only the observations retained during
#' estimation; \code{object$row.no} identifies their original row positions.
#' Observation weights affect the fitted coefficients, but individual fitted
#' shares, predictions and residuals are not multiplied by those weights.
#' A population mean prediction requires a separate weighted average over the
#' relevant covariate distribution.
#'
#' @return A data frame with one column per outcome and one row per observation.
#'   \code{fitted()} returns the stored fitted shares, \code{residuals()} returns
#'   observed minus fitted shares, and \code{predict()} returns predicted shares.
#'
#' @examples
#' x <- data.frame(income = seq(-1, 1, length.out = 30))
#' y <- cbind(base = 0.3 + 0.05*x$income,
#'            other = 0.7 - 0.05*x$income)
#' fit <- fmlogit(y, x, print.level = 0)
#' fitted(fit)
#' residuals(fit)
#' predict(fit)
#'
#' # Predict at two new income values, or at a single point.
#' predict(fit, newdata = data.frame(income = c(-0.5, 0.5)))
#' predict(fit, newdata = 0)
#'
#' # Remove the income effect while retaining the estimated intercepts.
#' beta <- fit$coefficient
#' beta[, "income"] <- 0
#' predict(fit, newbeta = beta)
#'
#' @rdname fitted.fmlogit
#' @export fitted.fmlogit
#' @exportS3Method
fitted.fmlogit <- function(object, ...) {
  if (!inherits(object, "fmlogit")) stop("'object' must be an fmlogit fit.")
  as.data.frame(object$fitted.values)
}

#' @rdname fitted.fmlogit
#' @export residuals.fmlogit
#' @exportS3Method
residuals.fmlogit <- function(object, ...) {
  as.data.frame(object$y - as.matrix(fitted.fmlogit(object)))
}

#' @rdname fitted.fmlogit
#' @export predict.fmlogit
#' @exportS3Method
predict.fmlogit <- function(object, newdata = NULL, newbeta = NULL, ...) {
  if (!inherits(object, "fmlogit")) stop("'object' must be an fmlogit fit.")
  J <- ncol(object$y)
  K <- ncol(object$x)

  if (is.null(newdata)) {
    x <- object$x
    row.names <- rownames(object$y)
  } else {
    # A vector is one observation; matrices retain their observation rows.
    x <- if (is.null(dim(newdata))) matrix(newdata, nrow = 1L) else as.matrix(newdata)
    if (!is.numeric(x) || nrow(x) < 1L || ncol(x) != K-1L ||
        any(!is.finite(x))) {
      stop("'newdata' must contain finite numeric values in K - 1 design columns, excluding the intercept.")
    }
    row.names <- rownames(x)
    x <- cbind("(Intercept)" = 1, x)
  }

  beta <- if (is.null(newbeta)) object$coefficient else as.matrix(newbeta)
  if (!is.numeric(beta) || !identical(dim(beta), c(J, K)) ||
      any(!is.finite(beta))) {
    stop("'newbeta' must be a finite J-by-K coefficient matrix, including the baseline and intercept.")
  }
  # Common coefficient shifts cancel in expected-share ratios; normalize beta_1.
  beta <- sweep(beta, 2L, beta[1L, ], FUN = "-")
  theta <- as.vector(t(beta[-1L, , drop = FALSE]))
  p <- softmax(theta, x, J)$p
  dimnames(p) <- list(row.names, colnames(object$y))
  as.data.frame(p)
}
