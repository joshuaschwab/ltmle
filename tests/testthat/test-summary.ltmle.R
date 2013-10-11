context("Test summary.ltmle") 

test_that("OR and RR are not calculated when Y is not binary", {
	n <- 10
  set.seed(50)
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1)
  r0 <- ltmle(data, Anodes="A", Ynodes="Y", abar=0)  

  s <- summary(r1, r0)

  expect_that(names(s$effect.measures), equals(c("ATE")))

  t1 <- ltmle(transform(data, Y=round(Y)), Anodes="A", Ynodes="Y", abar=1, survivalOutcome=FALSE)
  t0 <- ltmle(transform(data, Y=round(Y)), Anodes="A", Ynodes="Y", abar=0, survivalOutcome=FALSE)  

  u <- summary(t1, t0)

  expect_that(names(u$effect.measures), equals(c("ATE", "RR", "OR")))

})