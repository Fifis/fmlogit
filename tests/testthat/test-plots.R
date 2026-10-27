
test_that("Margins plots preserve effects, grouping and graphics settings", {
  s <- simulatedShares(N = 100L)
  fit <- fmlogit(s$y.multinom, s$x, print.level = 0L)
  varlist <- c("x2", "x1")
  object <- effects.fmlogit(fit, effect = "marginal", marg.type = "aveacr", varlist = varlist)

  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({
    grDevices::dev.off()
    unlink(file)
  }, add = TRUE)
  mfrow <- graphics::par("mfrow")

  result <- plot(object, varlist = varlist, x.data = fit$x, against.x = "x1", group.x = "x2",
    group.algebra = function(z) z > 0, mfrow = c(2, 3))

  expected <- do.call(cbind, object$marg.list[varlist])
  colnames(expected) <- unlist(lapply(varlist, function(v) {
    paste(v, colnames(fit$y), sep = ":")
  }), use.names = FALSE)

  expect_equal(result$x, fit$x[, "x1"])
  expect_equal(result$effects, expected)
  expect_equal(result$group, factor(fit$x[, "x2"] > 0))

  # Standard errors of average margins do not describe individual effects.
  expect_null(result$se)
})



test_that("WTP plots agree with individual effects and inference at grid points", {
  s <- simulatedShares(N = 100L)
  object <- fmlogit(s$y.multinom, s$x, print.level = 0L)
  N <- nrow(object$x)
  varlist <- c("x2", "x1")
  wtp.vec <- c(0, 50, 100)

  result <- plot(object, wtp.vec = wtp.vec, varlist = varlist,
                 effect = "marginal", plot.show = FALSE)

  margins <- effects.fmlogit(object, effect = "marginal", marg.type = "aveacr", varlist = varlist)
  expected <- wtp(margins, wtp.vec, indv.obs = TRUE)$wtp.obs

  expect_equal(result$x, seq_len(N))
  expect_equal(result$effects, expected)
  expect_null(result$se)

  # Compare the first evaluation point and its standard errors.
  set.seed(1)
  result <- plot(object, wtp.vec = wtp.vec, varlist = varlist, against = "x1", t = 3L,
    effect = "marginal", se = TRUE, R = 30L, plot.show = FALSE)

  weights <- object$weights
  at <- apply(object$x[, -1L, drop = FALSE], 2L, stats::weighted.mean, w = weights)
  support <- object$x[weights > 0, "x1"]
  at["x1"] <- min(support)

  set.seed(1)
  margins <- effects.fmlogit(object, effect = "marginal", varlist = varlist,
    at = at, se = TRUE, R = 30L)

  expected <- wtp(margins, wtp.vec)$wtp

  expect_equal(result$x, seq(min(support), max(support), length.out = 3L))
  expect_equal(unname(result$effects[1L, ]), unname(expected[, "Estimate"]))
  expect_equal(unname(result$se[1L, ]), unname(expected[, "Std. Error"]))
})

