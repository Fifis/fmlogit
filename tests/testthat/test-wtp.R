test_that("Aggregate uncertainty retains outcome covariance and adding-up identities", {
  effects <- matrix(c(-0.03, 0.01, 0.02), nrow = 3L,
                    dimnames = list(c("base", "second", "third"), "income"))
  # Cov(second, third) = 0.03; the baseline effect offsets their sum.
  cov.mat <- matrix(c(0.19, -0.07, -0.12,
                      -0.07, 0.04, 0.03,
                      -0.12, 0.03, 0.09), nrow = 3L)
  margins <- structure(list(effects = effects, cov.list = list(income = cov.mat),
                            R = 100L), class = "fmlogit.margins")
  values <- c(0, 10, 20)
  result <- wtp(margins, values)$wtp
  expect_equal(unname(result[, "Estimate"]), sum(values * effects))
  expected.se <- sqrt(10^2*0.04 + 20^2*0.09 + 2*10*20*0.03)
  expect_equal(unname(result[, "Std. Error"]), expected.se)
  expect_equal(wtp(margins, values + 100)$wtp, result)

  constant <- wtp(margins, rep(7, 3L))$wtp
  expect_equal(unname(constant[, c("Estimate", "Std. Error")]), c(0, 0))
})

test_that("Binary aggregate SEs agree with the single independent share effect", {
  set.seed(12)
  x <- data.frame(income = seq(-1, 1, length.out = 80L))
  share <- rbinom(nrow(x), 20L, plogis(0.2 + 0.5*x$income))/20
  fit <- fmlogit(cbind(base = 1 - share, other = share), x, print.level = 0L)
  margins <- effects(fit, varlist = "income", marg.type = "aveacr", se = TRUE, R = 50L)
  result <- wtp(margins, c(5, 25))$wtp

  # With two shares, a = (v_2 - v_1) m_2, so its SE scales by |v_2 - v_1|.
  expect_equal(unname(result[, "Estimate"]), 20*unname(margins$effects[2L, ]))
  expect_equal(unname(result[, "Std. Error"]), 20*unname(margins$se[2L, ]))

  margins$cov.list <- NULL
  expect_error(wtp(margins, c(5, 25)), "recompute margins")
})
