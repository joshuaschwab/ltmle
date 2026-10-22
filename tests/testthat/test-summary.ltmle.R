context("Test summary.ltmle") 

test_that("which effect measures are returned depends on the range of Y", {
  n <- 100
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  # non-binary Y in [0, 1]: RR and OR are both returned
  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=list(1, 0), variance.method="ic")
  s <- summary(r1)
  expect_equal(names(s$effect.measures), c("treatment", "control", "ATE", "RR", "OR"))

  # binary Y: RR and OR are both returned
  t1 <- ltmle(transform(data, Y=round(Y)), Anodes="A", Ynodes="Y", abar=list(1, 0), survivalOutcome=FALSE, estimate.time = FALSE, variance.method = "ic")
  u <- summary(t1)
  expect_equal(names(u$effect.measures), c("treatment", "control", "ATE", "RR", "OR"))

  # nonnegative Y outside [0, 1]: RR only
  v1 <- ltmle(transform(data, Y=Y*10), Anodes="A", Ynodes="Y", abar=list(1, 0), estimate.time = FALSE, variance.method = "ic")
  v <- summary(v1)
  expect_equal(names(v$effect.measures), c("treatment", "control", "ATE", "RR"))

  # Y can be negative: neither RR nor OR
  w1 <- ltmle(transform(data, Y=Y-0.5), Anodes="A", Ynodes="Y", abar=list(1, 0), estimate.time = FALSE, variance.method = "ic")
  w <- summary(w1)
  expect_equal(names(w$effect.measures), c("treatment", "control", "ATE"))
})

test_that("CI is truncated if outcome is binary", {
  n <- 1000
  W <- rnorm(n)
  A <- rexpit(W * 10)
  Y <- rbinom(n, 1, 0.99)
  data <- data.frame(W, A, Y)
  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, variance.method="iptw", estimate.time = FALSE)
  s1 <- summary(r1)
  expect_true(max(s1$treatment$CI) < 1.0001)
})

test_that("CI is not truncated if outcome is not binary", {
  n <- 100
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rnorm(n, mean = 5))
  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, variance.method="ic", estimate.time = FALSE)
  s1 <- summary(r1)
  expect_true(max(s1$treatment$CI) > 1)
})