context("Test summary.ltmle") 

test_that("OR and RR are not calculated when Y is not binary", {
	n <- 10
  set.seed(50)
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=list(1, 0), IC.variance.only=TRUE)

  s <- summary(r1)

  expect_that(names(s$effect.measures), equals(c("treatment", "control", "ATE")))

  t1 <- ltmle(transform(data, Y=round(Y)), Anodes="A", Ynodes="Y", abar=list(1, 0), survivalOutcome=FALSE)
  u <- summary(t1)

  expect_that(names(u$effect.measures), equals(c("treatment", "control", "ATE", "RR", "OR")))

})