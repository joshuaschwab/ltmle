context("summary.ltmleEffectMeasures: RR and OR on the original outcome scale")

sim <- function(Y.fun, n = 1000, seed = 1) {
  set.seed(seed)
  W <- rnorm(n); A <- rbinom(n, 1, plogis(0.3 * W))
  data.frame(W, A, Y = Y.fun(A, W, n))
}
fit <- function(d, ...) suppressMessages(ltmle(d, Anodes = "A", Ynodes = "Y", abar = list(1, 0), estimate.time = FALSE, ...))
est <- function(s, name) s$effect.measures[[name]]$estimate

test_that("binary outcome: RR and OR match the internal-scale formulas", {
  r <- fit(sim(function(A, W, n) rbinom(n, 1, plogis(-0.5 + 0.7 * A + 0.5 * W))))
  s <- summary(r); b <- r$beta
  expect_equal(est(s, "RR"), plogis(b[1] + b[2]) / plogis(b[1]), check.attributes = FALSE)
  expect_equal(est(s, "OR"), exp(b[2]), check.attributes = FALSE)
})

test_that("binary outcome with Yrange = c(0, 1) runs and matches Yrange = NULL", {
  d <- sim(function(A, W, n) rbinom(n, 1, plogis(-0.5 + 0.7 * A + 0.5 * W)))
  s1 <- summary(fit(d)); s2 <- summary(fit(d, Yrange = c(0, 1)))
  expect_equal(s1$effect.measures, s2$effect.measures)
})

test_that("rescaled nonnegative continuous outcome: RR is the ratio of original-scale means, no OR", {
  s <- summary(fit(sim(function(A, W, n) pmax(20 + 4 * A + 3 * W + rnorm(n, 0, 4), 0))))
  expect_equal(est(s, "RR"), est(s, "treatment") / est(s, "control"))
  expect_null(s$effect.measures$OR)
})

test_that("rescaled outcome with min(Yrange) > 0: RR is not scaled by the range", {
  s <- summary(fit(sim(function(A, W, n) 50 + 10 * A + 5 * W + rnorm(n, 0, 5)), Yrange = c(10, 100)))
  expect_equal(est(s, "RR"), est(s, "treatment") / est(s, "control"))
})

test_that("outcome in [0, 1] (not rescaled): RR and OR from the means", {
  s <- summary(fit(sim(function(A, W, n) pmin(pmax(0.1 + 0.05 * A + 0.03 * W + rnorm(n, 0, 0.05), 0), 1))))
  m1 <- est(s, "treatment"); m0 <- est(s, "control")
  expect_equal(est(s, "RR"), m1 / m0)
  expect_equal(est(s, "OR"), (m1 / (1 - m1)) / (m0 / (1 - m0)))
})

test_that("outcome that can be negative: no RR or OR", {
  s <- summary(fit(sim(function(A, W, n) rnorm(n, A + W))))
  expect_null(s$effect.measures$RR)
  expect_null(s$effect.measures$OR)
})
