#' Plot partial effects and user-valued aggregates
#'
#' @name plot.fmlogit
#' @aliases plot.fmlogit.margins
#' @description
#' Plot observation-level marginal or discrete effects, or their user-valued
#' aggregates from \code{wtp()}, using base R graphics.
#'
#' @param x An \code{"fmlogit.margins"} object computed with
#'   \code{marg.type = "aveacr"}, or an \code{"fmlogit"} fit for aggregate plots.
#' @param varlist Covariate names to plot, in the desired order. Defaults to
#'   all available effects, excluding the intercept for a fitted model.
#' @param x_data,y Optional design and observed-share matrices, such as
#'   \code{fit$x} and \code{fit$y}, aligned with the retained observations.
#'   Used by the margins method for named axes or groups.
#' @param against For margins, a numeric vector of horizontal coordinates.
#'   For a fitted model, a non-intercept design-column name: vary that column
#'   from its observed minimum to maximum, fixing others at their means.
#'   The default in both methods is to plot the observed effects by row number.
#' @param against.x,against.y Column name in \code{x_data} or \code{y} to use as the
#'   horizontal axis. Supply only one of the three axis arguments.
#' @param group Optional grouping vector for the margins method.
#' @param group.x Column name in \code{x_data} supplying groups.
#' @param group.algebra Optional function applied to the grouping values,
#'   for example \code{function(z) z > 0}.
#' @param mfrow Optional two-element panel layout, as in \code{par(mfrow = ...)}.
#'   By default, use a layout close to square. Graphics settings are restored.
#' @param type Base graphics plot type; defaults to points for margins and
#'   lines for aggregates. Lines connect increasing horizontal coordinates within groups.
#' @param plot.show Logical; draw the plots? Plotting data are returned either way.
#' @param xlab Optional horizontal-axis label.
#' @param wtp.vec Numeric outcome values, in the order of the fitted outcome columns.
#' @param t Number of evaluation points when varying \code{against}.
#' @param effect Either \code{"discrete"} (the default) or \code{"marginal"}.
#' @param se Compute pointwise aggregate SEs and draw plus/minus one SE?
#'   Defaults to \code{FALSE}.
#' @param R Number of Krinsky-Robb draws per evaluation point when \code{se = TRUE}.
#' @param ... Base graphics options, including \code{col}, \code{pch},
#'   \code{lty}, \code{main}, \code{ylab}, and \code{ylim}.
#'
#' @details
#' The margins method plots the stored observation-level effects. Their
#' standard errors are unavailable: the stored SEs describe average effects,
#' not the uncertainty of each observation-level effect.
#'
#' The fitted-model method aggregates effects using \code{wtp.vec}. Optional
#' pointwise SEs use \code{effects.fmlogit()} and the existing independence
#' approximation in \code{wtp()}; they omit covariance between outcome effects.
#' They condition on the evaluation design and can be slow to simulate.
#' Use \code{set.seed()} for reproducible SEs. The bars represent one SE,
#' rather than a particular confidence level.
#' Effects refer to processed design columns; varying one column does not
#' update related factor indicators, powers, or interactions. See
#' \code{effects.fmlogit()} for the reference-profile and averaging definitions.
#'
#' @return Invisibly, a list with \code{x}, the horizontal coordinates;
#'   \code{effects}, a matrix with one row per coordinate and one column per
#'   panel; \code{se}, a matching matrix or \code{NULL}; and \code{group}, a
#'   factor or \code{NULL}. Margins columns are named \code{variable:outcome};
#'   Aggregate columns are covariate names. Rows retain their input order even when
#'   lines are drawn in coordinate order.
#'
#' @examples
#' x <- data.frame(income = seq(-1, 1, length.out = 40), urban = rep(0:1, 20))
#' base <- 0.3 + 0.05 * x$income - 0.02 * x$urban
#' y <- cbind(base = base, other = 1 - base)
#' fit <- fmlogit(y, x, print.level = 0)
#' margins <- effects.fmlogit(fit, marg.type = "aveacr")
#' plot(margins, x_data = fit$x, against.x = "income", group.x = "urban")
#' curves <- plot(fit, wtp.vec = c(100, 200), varlist = "income",
#'                against = "income", t = 20)
#' head(curves$effects)
#' observed <- plot(margins, plot.show = FALSE)
#' head(observed$effects)
#'
#' @rdname plot.fmlogit
#' @export plot.fmlogit.margins
#' @exportS3Method
plot.fmlogit.margins <- function(x, varlist = NULL, x_data = NULL, y = NULL,
                                 against = NULL, against.x = NULL, against.y = NULL,
                                 group.x = NULL, group.algebra = NULL, mfrow = NULL,
                                 type = "p", plot.show = TRUE, group = NULL,
                                 xlab = NULL, ...) {
  object <- x
  if (!inherits(object, "fmlogit.margins") || is.null(object$marg.list)) {
    stop("'object' must contain margins computed with 'marg.type = \"aveacr\"'.")
  }
  N <- nrow(object$marg.list[[1L]])
  J <- nrow(object$effects)
  x <- x_data
  xnames <- colnames(object$effects)
  ynames <- rownames(object$effects)
  if (is.null(varlist)) varlist <- xnames
  varlist <- plot_variables(varlist, xnames)

  if (sum(c(!is.null(against), !is.null(against.x), !is.null(against.y))) > 1L) {
    stop("Supply only one of 'against', 'against.x', and 'against.y'.")
  }
  axis_label <- "Observation"
  if (!is.null(against.x)) {
    against <- plot_column(x, against.x, N)
    axis_label <- against.x
  } else if (!is.null(against.y)) {
    against <- plot_column(y, against.y, N)
    axis_label <- against.y
  } else if (!is.null(against)) {
    axis_label <- "against"
  }
  if (is.null(against)) against <- seq_len(N)
  if (!is.numeric(against) || length(against) != N || any(!is.finite(against))) {
    stop("'against' must contain one finite numeric coordinate per observation.")
  }

  if (!is.null(group.x)) {
    if (!is.null(group)) stop("Supply only one of 'group' and 'group.x'.")
    group <- plot_column(x, group.x, N)
  }
  if (!is.null(group.algebra)) {
    if (!is.function(group.algebra) || is.null(group)) {
      stop("'group.algebra' must be a function applied to 'group' or 'group.x'.")
    }
    group <- group.algebra(group)
  }
  if (!is.null(group)) {
    if (length(group) != N || anyNA(group)) stop("'group' must have one nonmissing value per observation.")
    group <- factor(group)
  }

  effects_mat <- do.call(cbind, object$marg.list[varlist])
  colnames(effects_mat) <- paste(rep(varlist, each = J), rep(ynames, times = length(varlist)), sep = ":")
  plot_data <- list(x = against, effects = effects_mat, se = NULL, group = group)
  if (plot.show) {
    if (is.null(xlab)) xlab <- axis_label
    draw_effects(against, effects_mat, group = group, mfrow = mfrow, type = type,
                 xlab = xlab, ...)
  }
  invisible(plot_data)
}

#' @rdname plot.fmlogit
#' @export plot.fmlogit
#' @exportS3Method
plot.fmlogit <- function(x, wtp.vec, varlist = NULL, against = NULL,
                         mfrow = NULL, t = 500L, effect = c("discrete", "marginal"),
                         type = "l", plot.show = TRUE, se = FALSE, R = 1000L,
                         xlab = NULL, ...) {
  object <- x
  if (!inherits(object, "fmlogit")) stop("'object' must be an fmlogit fit.")
  effect <- match.arg(effect)
  x <- object$x
  y <- object$y
  N <- nrow(y)
  J <- ncol(y)
  K <- ncol(x)
  xnames <- colnames(x)
  varlist <- plot_variables(varlist, xnames)
  if (!is.numeric(wtp.vec) || !is.null(dim(wtp.vec)) ||
      length(wtp.vec) != J || any(!is.finite(wtp.vec))) {
    stop("'wtp.vec' must be a finite numeric vector with one value per outcome.")
  }
  if (!is.logical(se) || length(se) != 1L || is.na(se)) stop("'se' must be TRUE or FALSE.")

  x_eval <- x
  against_values <- seq_len(N)
  if (!is.null(against)) {
    k <- match(against, xnames)
    if (length(k) != 1L || is.na(k) || k == 1L) stop("'against' must name a non-intercept design column.")
    if (!is.numeric(t) || length(t) != 1L || !is.finite(t) || t < 1 || t != floor(t)) {
      stop("'t' must be a positive integer.")
    }
    against_values <- seq(min(x[, k]), max(x[, k]), length.out = t)
    x_eval <- matrix(rep(colMeans(x), each = t), nrow = t, ncol = K,
                     dimnames = list(NULL, xnames))
    x_eval[, k] <- against_values
  }

  se_mat <- NULL
  if (is.null(against) && !se) {
    # Compute observation-level effects in a single call.
    margins <- effects.fmlogit(object, effect = effect, marg.type = "aveacr", varlist = varlist)
    effects_mat <- wtp(margins, wtp.vec, indv.obs = TRUE)$wtp.obs
  } else {
    effects_mat <- matrix(NA_real_, nrow(x_eval), length(varlist),
                    dimnames = list(rownames(x_eval), varlist))
    if (se) se_mat <- effects_mat
    for (i in seq_len(nrow(x_eval))) {
      margins <- effects.fmlogit(object, effect = effect, varlist = varlist,
                                at = x_eval[i, -1L], se = se, R = R)
      result <- wtp(margins, wtp.vec)
      effects_mat[i, ] <- if (se) result$wtp[, "estimate"] else result$wtp[1L, ]
      if (se) se_mat[i, ] <- result$wtp[, "std"]
    }
  }
  plot_data <- list(x = against_values, effects = effects_mat, se = se_mat, group = NULL)
  if (plot.show) {
    if (is.null(xlab)) xlab <- if (is.null(against)) "Observation" else against
    draw_effects(against_values, effects_mat, se = se_mat, mfrow = mfrow, type = type,
                 xlab = xlab, ...)
  }
  invisible(plot_data)
}

# Shared selection and plotting keep both methods on the same data/graphics path.
plot_variables <- function(varlist, xnames) {
  if (is.null(varlist)) varlist <- setdiff(xnames, "(Intercept)")
  if (!is.character(varlist) || !length(varlist) || anyNA(varlist) ||
      anyDuplicated(varlist) || anyNA(match(varlist, xnames))) {
    stop("'varlist' must contain distinct, recognized covariate names.")
  }
  varlist
}

plot_column <- function(data, column, N) {
  if (is.null(data) || nrow(data) != N || length(column) != 1L ||
      !column %in% colnames(data)) stop("Supply aligned data containing column '", column, "'.")
  data[, column]
}

draw_effects <- function(against, effects_mat, se = NULL, group = NULL, mfrow = NULL,
                         type = "p", xlab = "Observation", ylab = "Effect",
                         main = colnames(effects_mat), ylim = NULL, col = NULL,
                         pch = 1, lty = 1, cex = 1, ...) {
  n_panels <- ncol(effects_mat)
  if (is.null(mfrow)) {
    side <- ceiling(sqrt(n_panels))
    mfrow <- c(ceiling(n_panels / side), side)
  }
  if (!is.numeric(mfrow) || length(mfrow) != 2L || any(!is.finite(mfrow)) ||
      any(mfrow < 1 | mfrow != floor(mfrow))) stop("'mfrow' must contain two positive integers.")
  if (is.null(type)) type <- "l"
  type <- match.arg(type, c("p", "l", "b", "o", "h", "s", "S", "n"))
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = mfrow)

  show_legend <- !is.null(group)
  if (is.null(group)) group <- factor(rep("All", length(against)))
  n_groups <- nlevels(group)
  if (is.null(col)) col <- seq_len(n_groups)
  col <- rep(col, length.out = n_groups)
  pch <- rep(pch, length.out = n_groups)
  lty <- rep(lty, length.out = n_groups)
  main <- rep(main, length.out = n_panels)

  for (panel in seq_len(n_panels)) {
    limits <- if (is.null(ylim)) {
      range(0, effects_mat[, panel], if (!is.null(se)) c(effects_mat[, panel] - se[, panel],
                                                 effects_mat[, panel] + se[, panel]), finite = TRUE)
    } else ylim
    graphics::plot(against, effects_mat[, panel], type = "n", xlab = xlab, ylab = ylab,
                   main = main[panel], ylim = limits, ...)
    graphics::abline(h = 0, col = "grey70", lty = 3)
    for (level in seq_len(n_groups)) {
      i <- which(group == levels(group)[level])
      i <- i[order(against[i])]
      if (!is.null(se)) {
        graphics::segments(against[i], effects_mat[i, panel] - se[i, panel],
                           against[i], effects_mat[i, panel] + se[i, panel], col = col[level])
      }
      graphics::points(against[i], effects_mat[i, panel], type = type,
                       col = col[level], pch = pch[level], lty = lty[level], cex = cex)
    }
    if (show_legend) {
      graphics::legend("topright", legend = levels(group), col = col,
                       pch = if (type %in% c("l", "s", "S")) NA else pch,
                       lty = if (type == "p") 0 else lty, bty = "n")
    }
  }
}
