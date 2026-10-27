#' Generate summary tables for fmlogit objects
#'
#' Generate tables of coefficient estimates, partial effects, and user-valued
#' aggregates from fractional multinomial logit results.
#'
#' @name summary.fmlogit
#' @aliases summary.fmlogit.margins
#' @aliases summary.fmlogit.wtp
#'
#' @param object An object of class \code{"fmlogit"}, \code{"fmlogit.margins"},
#'   or \code{"fmlogit.wtp"}.
#' @param varlist Covariate names to select, in the desired output order.
#'   Defaults to all available covariates. Unrecognized names are omitted
#'   with a warning.
#' @param sepline Use separate rows for estimates and SEs?
#' @param digits Number of significant digits; defaults to three.
#' @param add.info Include sample, baseline, or effect descriptions?
#' @param list Return a list instead of a data frame? Used by the fitted-model
#'   and margins methods.
#' @param sigcode Decreasing p-value cutoffs for significance symbols.
#' @param print Also print the returned table?
#' @param ... Unused additional arguments.
#' @return For fitted models and margins, a list with \code{estimates} and
#'   optional descriptive fields, or a data frame when \code{list = FALSE}.
#'   For aggregates without SEs, a matrix of point estimates; with SEs, a
#'   data frame or, when \code{sepline = TRUE}, a character vector.
#'
#' @details This module provides summary methods for three fmlogit objects:
#' \code{fmlogit}, \code{fmlogit.margins}, and \code{fmlogit.wtp}.
#'
#' Coefficient tables omit the first, baseline outcome because its coefficients
#' are fixed at zero. The intercept remains the first design column unless
#' \code{varlist} requests another output order. Partial-effect tables include
#' the baseline outcome: its share can change even though its coefficients
#' are zero.
#'
#' SEs appear only when available in the input object. Coefficient SEs use the
#' full joint HC0 covariance; partial-effect SEs use Krinsky-Robb simulation.
#' Aggregate SEs use the joint covariance of outcome effects from those draws.
#' The stored z statistics and two-sided normal p-values determine
#' the significance symbols. Fitted-model descriptions report the maximised
#' quasi-log-likelihood, not a fully specified probability likelihood.
#'
#' @examples
#' x <- data.frame(income = seq(-1, 1, length.out = 30))
#' y <- cbind(base = 0.3 + 0.05 * x$income,
#'            other = 0.7 - 0.05 * x$income)
#' fit <- fmlogit(y, x, print.level = 0)
#' summary(fit)
#' margins <- effects.fmlogit(fit)
#' summary(margins, list = FALSE)
#' summary(wtp(margins, wtp.vec = c(0, 100)))
#' @rdname summary.fmlogit
#' @export summary.fmlogit
#' @exportS3Method

summary.fmlogit <- function(object, varlist = NULL, sepline = FALSE, digits = 3,
                           add.info = TRUE, list = TRUE,
                           sigcode = c(0.1, 0.05, 0.01, 0.001), print = FALSE, ...) {
  if (!inherits(object, "fmlogit")) stop("Expect an fmlogit object. Wrong object type given.")
  ynames <- names(object$estimates)
  xnames <- rownames(object$estimates[[1L]])
  if (length(varlist) == 0) varlist <- xnames
  var_colNo <- match(varlist, xnames)
  if (anyNA(var_colNo)) warning("Some variables requested are not in the variable list. Those variables are omitted.")
  var_colNo <- var_colNo[!is.na(var_colNo)]
  J <- object$count[3]
  n_vars <- length(var_colNo)
  if (n_vars == 0L) stop("No matching variables in 'varlist'.")
  varlist <- xnames[var_colNo]
  est <- se_mat <- p_values <- matrix(NA_real_, n_vars, J - 1L,
                                     dimnames = list(varlist, ynames))
  for (j in seq.int(2L, J)) {
    tab <- object$estimates[[j - 1L]][var_colNo, , drop = FALSE]
    est[, j - 1L] <- tab[, 1L]
    se_mat[, j - 1L] <- tab[, 2L]
    p_values[, j - 1L] <- tab[, 4L]
  }
  store_mat <- format_summary(est, se_mat, p_values, sepline, digits, sigcode)
  if (sepline) rownames(store_mat)[seq.int(2L, 2L * n_vars, by = 2L)] <- " "
  sig_print <- paste("Significance code: ", paste0(sigcode, " ", asterisk(sigcode * 0.99999999, sigcode), collapse = "  "))
  if (add.info) {
    nc <- paste("N=", object$count[1], sep = "")
    llf <- paste("quasi-log-likelihood=", round(object$likelihood, digits = 2), sep = "")
    bl <- paste("Baseline outcome:", object$baseline)
  }
  if (list) {
    outlist <- list(estimates = store_mat)
    if (add.info) {
      outlist$N <- nc
      outlist$llf <- llf
      outlist$baseline <- bl
      outlist$sigcode <- sig_print
    }
    if (print) print(outlist)
    return(outlist)
  } else {
    if (add.info) {
      info <- matrix(ncol = J - 1, nrow = 4)
      info[, 1] <- c(nc, llf, bl, sig_print)
      store_mat <- rbind(store_mat, info)
    }
    if (print) print(store_mat)
    return(as.data.frame(store_mat))
  }
}

#' @rdname summary.fmlogit
#' @export summary.fmlogit.margins
#' @exportS3Method

summary.fmlogit.margins <- function(object, varlist = NULL, sepline = FALSE,
                                   digits = 3, add.info = TRUE, list = TRUE,
                                   sigcode = c(0.05, 0.01, 0.001),
                                   print = FALSE, ...) {
  if (!inherits(object, "fmlogit.margins")) stop("Expect an fmlogit.margins object. Wrong object type given.")
  ynames <- rownames(object$effects)
  xnames <- colnames(object$effects)
  if (length(varlist) == 0) varlist <- xnames
  var_colNo <- match(varlist, xnames)
  if (anyNA(var_colNo)) warning("Some variables requested are not in the variable list. Those variables are omitted.")
  var_colNo <- var_colNo[!is.na(var_colNo)]
  J <- length(ynames)
  n_vars <- length(var_colNo)
  if (n_vars == 0L) stop("No matching variables in 'varlist'.")
  varlist <- xnames[var_colNo]

  est <- t(object$effects[, var_colNo, drop = FALSE])
  se_mat <- p_values <- NULL
  if (object$R > 0) {
    se_mat <- p_values <- matrix(NA_real_, n_vars, J,
                                 dimnames = list(varlist, ynames))
    for (k in seq_len(n_vars)) {
      tab <- object$ztable[[var_colNo[k]]]
      se_mat[k, ] <- tab[, 2L]
      p_values[k, ] <- tab[, 4L]
    }
  }
  store_mat <- format_summary(est, se_mat, p_values, sepline, digits, sigcode)
  sig_print <- paste("Significance code: ", paste0(sigcode, " ", asterisk(sigcode * 0.99999999, sigcode), collapse = "  "))
  if (add.info) {
    expl <- object$expl
  }
  if (list) {
    outlist <- list(estimates = store_mat)
    if (add.info) {
      outlist$expl <- expl
      outlist$sigcode <- sig_print
    }
    if (print) print(outlist)
    return(outlist)
  } else {
    if (add.info) {
      info <- matrix(ncol = J, nrow = 2)
      info[, 1] <- c(expl, sig_print)
      store_mat <- rbind(store_mat, info)
    }
    if (print) print(store_mat)
    return(as.data.frame(store_mat))
  }
}

#' @rdname summary.fmlogit
#' @export summary.fmlogit.wtp
#' @exportS3Method

summary.fmlogit.wtp <- function(object, varlist = NULL, sepline = FALSE,
                                digits = 3, sigcode = c(0.05, 0.01, 0.001),
                                print = FALSE, ...) {
  if (!inherits(object, "fmlogit.wtp")) stop("Expect an fmlogit.wtp object. Wrong object type given.")
  has_se <- identical(colnames(object$wtp), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  xnames <- if (has_se) rownames(object$wtp) else colnames(object$wtp)
  if (length(varlist) == 0) varlist <- xnames
  var_colNo <- match(varlist, xnames)
  if (anyNA(var_colNo)) warning("Some variables requested are not in the variable list. Those variables are omitted.")
  var_colNo <- var_colNo[!is.na(var_colNo)]
  n_vars <- length(var_colNo)
  if (n_vars == 0L) stop("No matching variables in 'varlist'.")
  varlist <- xnames[var_colNo]
  if (!has_se) {
    store_mat <- object$wtp[, var_colNo, drop = FALSE]
    if (print) print(store_mat)
    return(store_mat)
  }
  sig_print <- paste("Significance code: ", paste0(sigcode, " ", asterisk(sigcode * 0.99999999, sigcode), collapse = "  "))
  tab <- object$wtp[var_colNo, , drop = FALSE]
  est <- matrix(tab[, 1L], ncol = 1L, dimnames = list(varlist, NULL))
  se_mat <- matrix(tab[, 2L], ncol = 1L)
  p_values <- matrix(tab[, 4L], ncol = 1L)
  store_mat <- format_summary(est, se_mat, p_values, sepline, digits, sigcode)
  if (!sepline) {
    store_mat <- as.data.frame(store_mat)
    colnames(store_mat) <- NULL
  } else {
    store_mat <- stats::setNames(as.vector(store_mat), rownames(store_mat))
  }
  if (print) {
    print(store_mat)
    print(sig_print)
  }
  return(store_mat)
}

# Format variables in rows and outcomes in columns, keeping p-values unrounded.
format_summary <- function(est, se_mat, p_values, sepline, digits, sigcode) {
  est <- signif(est, digits = digits)
  if (is.null(se_mat)) return(est)

  stars <- matrix(asterisk(p_values, sigcode), nrow(est), ncol(est))
  se_text <- matrix(paste0("(", signif(se_mat, digits = digits), ")"),
                    nrow(est), ncol(est))
  if (!sepline) {
    return(matrix(paste0(est, se_text, stars), nrow(est), ncol(est),
                  dimnames = dimnames(est)))
  }

  n_vars <- nrow(est)
  store_mat <- matrix("", 2L * n_vars, ncol(est),
                       dimnames = list(rep("", 2L * n_vars), colnames(est)))
  for (k in seq_len(n_vars)) {
    store_mat[2L * k - 1L, ] <- paste0(est[k, ], stars[k, ])
    store_mat[2L * k, ] <- se_text[k, ]
    rownames(store_mat)[2L * k - 1L] <- rownames(est)[k]
  }
  store_mat
}
