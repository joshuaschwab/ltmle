context("Test CheckInputs")

test_that("Lnodes, Anodes, Cnodes, Ynodes include all columns after the first one", {
  n <- 10 
  set.seed(2345)
  dat <- data.frame(W=rbinom(n, 1, .5),
                    A1=rbinom(n, 1, .5),
                    L2=rbinom(n, 1, .5),
                    A2=rbinom(n, 1, .5),
                    Y=rbinom(n, 1, .5)
                    )

  expect_that(ltmle(dat, Anodes=c("A1", "A2"),
     Lnodes = NULL, 
     Ynodes = "Y", 
     abar = c(1,1)),
     throws_error("All nodes after the A-, C-, L-, or Ynodes must be in A-, C-, L-, or Ynodes"))
  })

test_that("In formulas RHS variables are parents of the current node", {

n <- 10 
  set.seed(2345)
  dat <- data.frame(W=rbinom(n, 1, .5),
                    A1=rbinom(n, 1, .5),
                    L2=rbinom(n, 1, .5),
                    A2=rbinom(n, 1, .5),
                    Y=rbinom(n, 1, .5)
                    )  

  bad.gform <- c(A1="A1~W", A2="A2~Y+L2+A1+W")
  bad.Qform <- c(L2="Qkplus1~A1+W", Y="Qkplus1~L2+A2+Y")

expect_that(ltmle(dat, Anodes=c("A1", "A2"),
     Lnodes = "L2", 
     Ynodes = "Y", 
     abar = c(1,1)),
     gform = bad.gform,  
     throws_error("something something"))

expect_that(ltmle(dat, Anodes=c("A1", "A2"),
     Lnodes = "L2", 
     Ynodes = "Y", 
     abar = c(1,1)),
     Qform = bad.Qform,  
     throws_error("something something"))
  })