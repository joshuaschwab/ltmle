context("Test CheckInputs")

test_that("Lnodes, Anodes, Cnodes, Ynodes include all columns after the first one", {
  n <- 30 
  dat <- data.frame(W=rbinom(n, 1, .5),
                    A1=rbinom(n, 1, .5),
                    L2=rbinom(n, 1, .5),
                    A2=rbinom(n, 1, .5),
                    Y=rbinom(n, 1, .5)
                    )

  expect_that(ltmle(dat, Anodes=c("A1", "A2"),
     Lnodes = NULL, 
     Ynodes = "Y", 
     abar = c(1,1),
     survivalOutcome = FALSE,
     estimate.time = FALSE),
     throws_error("All nodes after the first A/C node must be in A-, C-, L-, or Ynodes"))
  })

test_that("In formulas RHS variables are parents of the current node", {
  n <- 30 
  dat <- data.frame(W=rbinom(n, 1, .5),
                    A1=rbinom(n, 1, .5),
                    L2=rbinom(n, 1, .5),
                    A2=rbinom(n, 1, .5),
                    Y=rbinom(n, 1, .5)
                    )  

  expect_that(ltmle(dat, Anodes=c("A1", "A2"),
                    Lnodes = "L2", 
                    Ynodes = "Y", 
                    abar = c(1,1),
                    gform = c("A1~W", "A2~Y+L2+A1+W"),
                    survivalOutcome = FALSE,
                    estimate.time = FALSE),  
       throws_error("Some nodes in gform\\[2\\] are not parents of A2"))

  expect_that(ltmle(dat, Anodes=c("A1", "A2"),
                    Lnodes = "L2", 
                    Ynodes = "Y", 
                    abar = c(1,1),
                    Qform = c(L2="Q.kplus1~A1+W", Y="Q.kplus1~L2+A2+Y"),
                    survivalOutcome = FALSE,
                    estimate.time = FALSE),  
       throws_error("Some nodes in Qform\\[2\\] are not parents of Y"))

  #No baseline covars before first A node
  dat2 <- dat[, -1]
  expect_that(ltmle(dat2, Anodes=c("A1", "A2"),
                    Lnodes = "L2", 
                    Ynodes = "Y", 
                    abar = c(1,1),
                    gform = c("A1~A1", "A2~1"),
                    survivalOutcome = FALSE,
                    estimate.time = FALSE),
       throws_error("Some nodes in gform\\[1\\] are not parents of A1")) 

  expect_that(ltmle(dat2, Anodes=c("A1", "A2"),
                    Lnodes = "L2", 
                    Ynodes = "Y", 
                    abar = c(1,1),
                    gform = c("A1~1", "A2~1"),
                    survivalOutcome = FALSE,
                    estimate.time = FALSE),  
       is_a("ltmle"))            
})

test_that("Non-binary outcome can be in [0, 1]", {
  n <- 30
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  expect_that(ltmle(data, Anodes="A", Ynodes="Y", abar=1, variance.method="ic", estimate.time = FALSE), is_a("ltmle"))
})

test_that("survivalOutcome=TRUE requires binary outcomes.", {
  n <- 30
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  expect_that(ltmle(data, Anodes="A", Ynodes="Y", abar=1, survivalOutcome=TRUE, variance.method="ic", estimate.time = FALSE), 
    throws_error("When survivalOutcome is TRUE, all Ynodes should be 0, 1, or NA"))
})

test_that("binaryOutcome flag set correctly", {
  n <- 30
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = rbinom(n, 2, .5)/2)

  expect_false(ltmle(data, Anodes="A", Ynodes="Y", abar=1, survivalOutcome=FALSE, variance.method="ic",  estimate.time = FALSE)$binaryOutcome)
  expect_true(ltmle(transform(data, Y=round(Y)), Anodes="A", Ynodes="Y", abar=1, survivalOutcome=FALSE)$binaryOutcome) 
})


test_that("Y outside of [0, 1] is scaled", {
  n <- 30
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y = runif(n)*10)

  expect_that(l <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, variance.method="ic",  estimate.time = FALSE), 
    shows_message("Some Ynodes are not in \\[0, 1\\], and Yrange was NULL, so all Y nodes are being\\ntransformed to \\(Y-min.of.all.Ys\\)/range.of.all.Ys"))  
  expect_true(l$transformOutcome)
  expect_that(ltmle(data, Anodes="A", Ynodes="Y", abar=1, Yrange=c(0, 3), variance.method="ic", estimate.time = FALSE), 
    shows_message("Some Ynodes are not in \\[Yrange\\[1\\], Yrange\\[2\\]\\], Y values are truncated"))
  expect_false(ltmle(transform(data, Y=Y/30), Anodes="A", Ynodes="Y", abar=1, variance.method="ic", estimate.time = FALSE)$transformOutcome)
})

test_that("survivalOutcome required if outcome is binary", {
  n <- 30
  data <- data.frame(W = rnorm(n),
                     A = rbinom(n, 1, .5), 
                     Y1 = rbinom(n, 1, .5),
                     Y2 = rbinom(n, 1, .5))

  expect_that(ltmle(data, Anodes="A", Ynodes=c("Y1", "Y2"), abar=1, estimate.time = FALSE),
    throws_error("All Ynodes are 0, 1, or NA; the outcome is treated as binary. The 'survivalOutcome' argument must be specified."))
  expect_error(ltmle(data, Anodes="A", Ynodes=c("Y1", "Y2"), abar=1, estimate.time = FALSE, Yrange = c(0, 0.5), survivalOutcome = TRUE), "All Ynodes are 0, 1, or NA, but Yrange is something other")
})

test_that("cvControl requires correct names and makes a difference", {
  set.seed(1)
  n <- 50
  W <- rnorm(n)
  A <- rbinom(n, 1, plogis(W))
  Y <- rbinom(n, 1, plogis(W + A))
  data <- data.frame(W, A, Y)
  expect_error(ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = FALSE, SL.library = "default", SL.cvControl = 5), "SL.cvControl must be a list")
  expect_error(ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = FALSE, SL.library = "default", SL.cvControl = list(validRows = 1)), "The valid names for SL.cvControl are V, stratifyCV, shuffle. validRows is not currently supported.")
  r1 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = FALSE, SL.library = c("SL.step", "SL.glm", "SL.mean", "SL.lm"))
  r2 <- ltmle(data, Anodes="A", Ynodes="Y", abar=1, estimate.time = FALSE, SL.library = c("SL.step", "SL.glm", "SL.mean", "SL.lm"), SL.cvControl = list(V = 5, shuffle = F, stratifyCV = T))
  expect_gt(abs(r1$estimates["tmle"] - r2$estimates["tmle"]), 0.001)
})
