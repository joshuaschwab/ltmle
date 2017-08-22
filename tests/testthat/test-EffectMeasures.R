context("effect measures")

expect_equals <- function(...) expect_equal(..., tolerance=1e-8, scale=1, check.attributes=FALSE)

test_that("effect measures with continuous Y is same as separate calls", {
  n <- 200
  W <- rnorm(n)
  A <- rbinom(n, 1, plogis(W))
  Y <- 10 + W + A + rnorm(n)
  data <- data.frame(W, A, Y)
  variance.method <- "ic"
  
  ltmle.fit0 <- ltmle(data, Anodes="A", Ynodes="Y", abar=0, variance.method = variance.method)
  ltmle.fit1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, variance.method = variance.method)
  ltmle.fit <- ltmle(data, Anodes="A", Ynodes="Y", abar=list(1, 0), variance.method = variance.method)
  
  s <- summary(ltmle.fit)
  
  s1 <- sqrt(var((ltmle.fit1$IC$tmle - ltmle.fit0$IC$tmle)) /n) 
  s2 <- s$effect.measures$ATE$std.dev
  
  expect_equals(s1, s2)
  expect_equals(s$effect.measures$treatment$estimate, ltmle.fit1$estimates["tmle"])
  expect_equals(s$effect.measures$control$estimate, ltmle.fit0$estimates["tmle"])
  expect_equals(s$effect.measures$ATE$estimate, ltmle.fit1$estimates["tmle"] - ltmle.fit0$estimates["tmle"])
})