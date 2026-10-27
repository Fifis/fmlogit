summary_model <- function() {
  x <- data.frame(income = seq(-1, 1, length.out = 40), urban = rep(0:1, 20))
  y <- cbind(base = 0.3 + 0.05 * x$income - 0.02 * x$urban)
  y <- cbind(y, other = 1 - y[, "base"])
  fmlogit(y, x, print.level = 0L)
}

test_that("Summaries preserve baseline conventions and requested covariate order", {
  fit <- summary_model()
  table <- summary(fit, add.info = FALSE)$estimates
  expect_identical(rownames(table), colnames(fit$x))
  expect_identical(colnames(table), colnames(fit$y)[-1L])

  varlist <- c("urban", "income")
  table <- summary(fit, varlist = varlist, add.info = FALSE)$estimates
  expect_identical(rownames(table), varlist)

  margins <- effects(fit)
  table <- summary(margins, varlist = varlist, add.info = FALSE)$estimates
  expect_identical(rownames(table), varlist)
  expect_identical(colnames(table), colnames(fit$y))

  aggregate <- wtp(margins, c(0, 100))
  expect_identical(colnames(summary(aggregate, varlist = varlist)), varlist)
})

test_that("Rounding p-values cannot give a nonsignificant estimate a star", {
  fit <- summary_model()
  set.seed(123)
  margins <- effects(fit, varlist = "income", se = TRUE, R = 10L)
  aggregate <- wtp(margins, c(0, 100))

  # A p-value just above 5% must remain nonsignificant after table formatting.
  fit$estimates[[1L]][, 4L] <- 0.05001
  margins$ztable$income[, "Pr(>|z|)"] <- 0.05001
  aggregate$wtp[, "Pr(>|z|)"] <- 0.05001
  sigcode <- c(0.05, 0.01, 0.001)
  tables <- list(
    summary(fit, add.info = FALSE, sigcode = sigcode)$estimates,
    summary(margins, add.info = FALSE, sigcode = sigcode)$estimates,
    summary(aggregate, sigcode = sigcode)
  )
  for (table in tables) expect_false(any(grepl("*", as.matrix(table), fixed = TRUE)))
})
