#' Generate summary tables for fmlogit objects
#'
#' Generate tables of coefficient estimates, partial effects, and aggregates
#' calculated from supplied outcome values in a fractional multinomial logit model.
#'
#' @name summary.fmlogit
#' @aliases summary.fmlogit.margins
#' @aliases summary.fmlogit.wtp
#'
#' @param object An object of class \code{"fmlogit"}, \code{"fmlogit.margins"},
#'   or \code{"fmlogit.wtp"}.
#' @param varlist Covariate names to select, in the desired output order.
#'   Defaults to all available covariates. Unrecognised names are omitted
#'   with a warning.
#' @param sepline Show estimates and their standard errors on separate rows?
#' @param digits Number of significant digits; defaults to three.
#' @param add.info Include details of the sample, baseline, or reported effects?
#' @param list Return a list instead of a data frame? Used by the fitted-model
#'   and margins methods.
#' @param sigcode Decreasing p-value thresholds for significance symbols.
#' @param print Also print the returned table?
#' @param ... Unused additional arguments.
#' @return For fitted models and margins, a list with \code{estimates} and
#'   optional descriptive fields, or a data frame when \code{list = FALSE}.
#'   For aggregates without SEs, a matrix of point estimates; with SEs, a
#'   data frame or, when \code{sepline = TRUE}, a character vector.
#'
#' @details These methods summarise objects of three classes:
#' \code{fmlogit}, \code{fmlogit.margins}, and \code{fmlogit.wtp}.
#'
#' Coefficient tables omit the first, baseline outcome because its coefficients
#' are fixed at zero. The intercept remains the first design column unless
#' \code{varlist} requests another output order. Partial-effect tables include
#' the baseline outcome: its share can change even though its coefficients
#' are zero.
#'
#' Standard errors are shown only when available in the input object.
#' Coefficient SEs use the full joint HC0 covariance; partial-effect SEs use
#' Krinsky-Robb simulation.
#' Aggregate SEs use the joint covariance of outcome effects from those draws.
#' The significance symbols follow the stored two-sided normal p-values.
#' Fitted-model summaries report the maximised quasi-log-likelihood, rather
#' than a likelihood for a fully specified distribution.
#'
#' @examples
#' x <- data.frame(income = seq(-1, 1, length.out = 30))
#' y <- cbind(base = 0.3 + 0.05*x$income,
#'            other = 0.7 - 0.05*x$income)
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
  if (!inherits(object, "fmlogit")) stop("'object' must be an fmlogit object.")
  ynames <- names(object$estimates)
  xnames <- rownames(object$estimates[[1L]])
  if (length(varlist) == 0) varlist <- xnames
  var.col.no <- match(varlist, xnames)
  if (anyNA(var.col.no)) warning("Some requested covariates were not found and have been omitted.")
  var.col.no <- var.col.no[!is.na(var.col.no)]
  J <- object$count[3]
  n.vars <- length(var.col.no)
  if (n.vars == 0L) stop("No matching variables in 'varlist'.")
  varlist <- xnames[var.col.no]
  est <- se.mat <- p.values <- matrix(NA_real_, nrow = n.vars, ncol = J-1L,
                                      dimnames = list(varlist, ynames))
  for (j in seq.int(2L, J)) {
    tab <- object$estimates[[j-1L]][var.col.no, , drop = FALSE]
    est[, j-1L] <- tab[, 1L]
    se.mat[, j-1L] <- tab[, 2L]
    p.values[, j-1L] <- tab[, 4L]
  }
  store.mat <- formatSummary(est, se.mat, p.values, sepline, digits, sigcode)
  if (sepline) rownames(store.mat)[seq.int(2L, 2L*n.vars, by = 2L)] <- " "
  sig.print <- paste("Significance code: ",
                     paste0(sigcode, " ", asterisk(sigcode*0.99999999, sigcode), collapse = "  "))
  if (add.info) {
    nc <- paste("N=", object$count[1], sep = "")
    llf <- paste("quasi-log-likelihood=", round(object$likelihood, digits = 2), sep = "")
    bl <- paste("Baseline outcome:", object$baseline)
  }
  if (list) {
    outlist <- list(estimates = store.mat)
    if (add.info) {
      outlist$N <- nc
      outlist$llf <- llf
      outlist$baseline <- bl
      outlist$sigcode <- sig.print
    }
    if (print) print(outlist)
    return(outlist)
  } else {
    if (add.info) {
      info <- matrix(ncol = J-1, nrow = 4)
      info[, 1] <- c(nc, llf, bl, sig.print)
      store.mat <- rbind(store.mat, info)
    }
    if (print) print(store.mat)
    return(as.data.frame(store.mat))
  }
}

#' @rdname summary.fmlogit
#' @export summary.fmlogit.margins
#' @exportS3Method

summary.fmlogit.margins <- function(object, varlist = NULL, sepline = FALSE,
                                    digits = 3, add.info = TRUE, list = TRUE,
                                    sigcode = c(0.05, 0.01, 0.001),
                                    print = FALSE, ...) {
  if (!inherits(object, "fmlogit.margins")) stop("'object' must be an fmlogit.margins object.")
  ynames <- rownames(object$effects)
  xnames <- colnames(object$effects)
  if (length(varlist) == 0) varlist <- xnames
  var.col.no <- match(varlist, xnames)
  if (anyNA(var.col.no)) warning("Some requested covariates were not found and have been omitted.")
  var.col.no <- var.col.no[!is.na(var.col.no)]
  J <- length(ynames)
  n.vars <- length(var.col.no)
  if (n.vars == 0L) stop("No matching variables in 'varlist'.")
  varlist <- xnames[var.col.no]

  est <- t(object$effects[, var.col.no, drop = FALSE])
  se.mat <- p.values <- NULL
  if (object$R > 0) {
    se.mat <- p.values <- matrix(NA_real_, nrow = n.vars, ncol = J,
                                 dimnames = list(varlist, ynames))
    for (k in seq_len(n.vars)) {
      tab <- object$ztable[[var.col.no[k]]]
      se.mat[k, ] <- tab[, 2L]
      p.values[k, ] <- tab[, 4L]
    }
  }
  store.mat <- formatSummary(est, se.mat, p.values, sepline, digits, sigcode)
  sig.print <- paste("Significance code: ",
                     paste0(sigcode, " ", asterisk(sigcode*0.99999999, sigcode), collapse = "  "))
  if (add.info) {
    expl <- object$expl
  }
  if (list) {
    outlist <- list(estimates = store.mat)
    if (add.info) {
      outlist$expl <- expl
      outlist$sigcode <- sig.print
    }
    if (print) print(outlist)
    return(outlist)
  } else {
    if (add.info) {
      info <- matrix(ncol = J, nrow = 2)
      info[, 1] <- c(expl, sig.print)
      store.mat <- rbind(store.mat, info)
    }
    if (print) print(store.mat)
    return(as.data.frame(store.mat))
  }
}

#' @rdname summary.fmlogit
#' @export summary.fmlogit.wtp
#' @exportS3Method

summary.fmlogit.wtp <- function(object, varlist = NULL, sepline = FALSE,
                                digits = 3, sigcode = c(0.05, 0.01, 0.001),
                                print = FALSE, ...) {
  if (!inherits(object, "fmlogit.wtp")) stop("'object' must be an fmlogit.wtp object.")
  has.se <- identical(colnames(object$wtp), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  xnames <- if (has.se) rownames(object$wtp) else colnames(object$wtp)
  if (length(varlist) == 0) varlist <- xnames
  var.col.no <- match(varlist, xnames)
  if (anyNA(var.col.no)) warning("Some requested covariates were not found and have been omitted.")
  var.col.no <- var.col.no[!is.na(var.col.no)]
  n.vars <- length(var.col.no)
  if (n.vars == 0L) stop("No matching variables in 'varlist'.")
  varlist <- xnames[var.col.no]
  if (!has.se) {
    store.mat <- object$wtp[, var.col.no, drop = FALSE]
    if (print) print(store.mat)
    return(store.mat)
  }
  sig.print <- paste("Significance code: ",
                     paste0(sigcode, " ", asterisk(sigcode*0.99999999, sigcode), collapse = "  "))
  tab <- object$wtp[var.col.no, , drop = FALSE]
  est <- matrix(tab[, 1L], ncol = 1L, dimnames = list(varlist, NULL))
  se.mat <- matrix(tab[, 2L], ncol = 1L)
  p.values <- matrix(tab[, 4L], ncol = 1L)
  store.mat <- formatSummary(est, se.mat, p.values, sepline, digits, sigcode)
  if (!sepline) {
    store.mat <- as.data.frame(store.mat)
    colnames(store.mat) <- NULL
  } else {
    store.mat <- stats::setNames(as.vector(store.mat), rownames(store.mat))
  }
  if (print) {
    print(store.mat)
    print(sig.print)
  }
  return(store.mat)
}

# Format variables in rows and outcomes in columns, keeping p-values unrounded.
formatSummary <- function(est, se.mat, p.values, sepline, digits, sigcode) {
  est <- signif(est, digits = digits)
  if (is.null(se.mat)) return(est)

  stars <- matrix(asterisk(p.values, sigcode), nrow = nrow(est), ncol = ncol(est))
  se.text <- matrix(paste0("(", signif(se.mat, digits = digits), ")"),
                    nrow = nrow(est), ncol = ncol(est))
  if (!sepline) {
    return(matrix(paste0(est, se.text, stars), nrow = nrow(est), ncol = ncol(est),
                  dimnames = dimnames(est)))
  }

  n.vars <- nrow(est)
  store.mat <- matrix("", nrow = 2L*n.vars, ncol = ncol(est),
                      dimnames = list(rep("", 2L*n.vars), colnames(est)))
  for (k in seq_len(n.vars)) {
    store.mat[2L*k-1L, ] <- paste0(est[k, ], stars[k, ])
    store.mat[2L*k, ] <- se.text[k, ]
    rownames(store.mat)[2L*k-1L] <- rownames(est)[k]
  }
  store.mat
}
