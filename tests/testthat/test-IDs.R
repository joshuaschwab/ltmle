context("id")

n <- 500
W <- rnorm(n)
A <- rexpit(W)
Y <- rexpit(W + A)
data <- data.frame(W, A, Y)

test_that("invalid id throws error", {
  expect_error(ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = FALSE, id=0:n), "id must be a vector with length") 
})

test_that("unique ids same as NULL", {
  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = FALSE, variance.method = "ic")
  r2 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = FALSE, variance.method = "ic", id=10000 + sample(1:n, n, replace = FALSE))
  r1$call <- r2$call <- NULL
  expect_equal(summary(r1), summary(r2))
})

test_that("multiple subjects per household increases variance", {
  n <- 500
  num.households <- 20
  u <- rep(rnorm(num.households), length.out=n)
  id <- rep(1:num.households, length.out=n)
  W <- rnorm(n) + u
  A <- rexpit(W)
  Y <- rexpit(W + A + u)
  data <- data.frame(W, A, Y)
  
  
  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = FALSE, variance.method = "ic")
  r2 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = FALSE, variance.method = "ic", id=id)
  for (est in c("tmle", "iptw")) {
    expect_true(summary(r1, estimator=est)$treatment$std.dev < summary(r2, estimator=est)$treatment$std.dev)
    expect_equal(r1$estimates[est], r2$estimates[est])
  }
  r3 <- ltmleMSM(data, Anodes="A", Ynodes="Y", regimes=list(function (x) 1, function (x) 0), summary.measures=NULL, working.msm = "Y~1", variance.method = "ic", estimate.time = FALSE)
  r4 <- ltmleMSM(data, Anodes="A", Ynodes="Y", regimes=list(function (x) 1, function (x) 0), summary.measures=NULL, working.msm = "Y~1", variance.method = "ic", id=id, estimate.time = FALSE)
  
  for (est in c("tmle", "iptw")) {
    expect_true(summary(r3, estimator=est)$cmat[, 2] < summary(r4, estimator=est)$cmat[, 2])
    expect_equal(summary(r3, estimator=est)$cmat[, 1], summary(r4, estimator=est)$cmat[, 1])
  }
})
