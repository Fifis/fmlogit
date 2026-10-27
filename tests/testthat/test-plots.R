test_that("Margins plots preserve effects, grouping and graphics settings", {
  income <- matrix(c(-0.2, -0.1, 0.1, 0.2, 0.2, 0.1, -0.1, -0.2), 4L,
                   dimnames = list(NULL, c("base", "other")))
  individual <- list(income = income, indicator = 2 * income)
  object <- structure(list(effects = vapply(individual, colMeans, numeric(2L)),
                           marg.list = individual,
                           se = matrix(0.5, 2L, 2L), R = 100L),
                      class = "fmlogit.margins")
  x <- cbind(income = c(2, 0, 3, 1), indicator = c(0, 1, 0, 1))

  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  mfrow <- graphics::par("mfrow")
  result <- plot(x = object, varlist = c("indicator", "income"), x_data = x,
                 against.x = "income", group.x = "indicator",
                 group.algebra = function(z) z > 0, mfrow = c(2, 2))

  expected <- cbind(individual$indicator, individual$income)
  colnames(expected) <- c("indicator:base", "indicator:other",
                         "income:base", "income:other")
  expect_equal(result$x, x[, "income"])
  expect_equal(result$effects, expected)
  expect_equal(result$group, factor(x[, "indicator"] > 0))
  # Standard errors of average margins do not describe individual effects.
  expect_null(result$se)
  expect_equal(graphics::par("mfrow"), mfrow)
})

test_that("WTP plots agree with individual effects and inference at grid points", {
  set.seed(19)
  N <- 80L
  x <- cbind(income = seq(-2, 2, length.out = N),
             indicator = rep(c(0, 1), N / 2L))
  p_other <- stats::plogis(0.2 + 0.5 * x[, "income"] - 0.3 * x[, "indicator"])
  share <- stats::rbinom(N, 20L, p_other) / 20
  y <- cbind(base = 1 - share, other = share)
  object <- fmlogit(y, x, method = "BFGS", print.level = 0L)
  varlist <- c("indicator", "income")
  wtp.vec <- c(0, 100)

  result <- plot(object, wtp.vec, varlist, effect = "marginal", plot.show = FALSE)
  margins <- effects.fmlogit(object, effect = "marginal", marg.type = "aveacr",
                             varlist = varlist)
  expected <- wtp(margins, wtp.vec, indv.obs = TRUE)$wtp.obs
  expect_equal(result$x, seq_len(N))
  expect_equal(result$effects, expected)
  expect_null(result$se)

  set.seed(7)
  result <- plot(object, wtp.vec, varlist, against = "income", t = 3L,
                 effect = "marginal", se = TRUE, R = 30L, plot.show = FALSE)
  at <- colMeans(object$x[, -1L, drop = FALSE])
  at["income"] <- min(x[, "income"])
  set.seed(7)
  margins <- effects.fmlogit(object, effect = "marginal", varlist = varlist,
                             at = at, se = TRUE, R = 30L)
  expected <- wtp(margins, wtp.vec)$wtp
  expect_equal(result$x, seq(min(x[, "income"]), max(x[, "income"]), length.out = 3L))
  expect_equal(unname(result$effects[1L, ]), unname(expected[, "estimate"]))
  expect_equal(unname(result$se[1L, ]), unname(expected[, "std"]))
})
