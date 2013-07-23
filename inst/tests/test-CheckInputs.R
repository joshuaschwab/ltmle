context("Test CheckInputs")

test_that("Lnodes, Anodes, Cnodes, Ynodes include all columns after the first one", {
  n <- 100 
  set.seed(2345)
  dat <- data.frame(W=rbinom(n, 1, .5),
                    A1=rbinom(n, 1, .5),
                    L1=rbinom(n, 1, .5),
                    A2=rbinom(n, 1, .5),
                    Y=rbinom(n, 1, .5)
                    )

  expect_that(ltmle(dat, Anodes=c("A1", "A2"),
     Lnodes = NULL, 
     Ynodes = "Y", 
     abar = c(1,1)),
     throws_error("All nodes after the A-, C-, L-, or Ynodes must be in A-, C-, L-, or Ynodes"))
  })
