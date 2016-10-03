context("abar/regimes tests")

n <- 50
W <- rnorm(n)
C.binary <- rexpit(W)
C <- BinaryToCensoring(is.uncensored = C.binary)
A <- rexpit(W)
Y <- rexpit(W+A)
data <- data.frame(W, C, A, Y)


test_that("logical abar is the same as numeric",{
  r1 <- ltmle(data, Anodes="A", Cnodes="C", Ynodes="Y", abar=as.matrix(W > 0), estimate.time = F)
  r2 <- ltmle(data, Anodes="A", Cnodes="C", Ynodes="Y", abar=as.matrix(as.numeric(W > 0)), estimate.time = F)
  r1$call <- r2$call <- NULL
  expect_equal(r1, r2)
})

test_that("list of rules works", {
  rule1 <- function (row) row["W"] > 0
  rule2 <- function (row) row["W"] > 1
  r1 <- ltmle(data, Anodes="A", Cnodes="C", Ynodes="Y", rule=list(rule1, rule2), estimate.time = F)
  r2 <- ltmle(data, Anodes="A", Cnodes="C", Ynodes="Y", rule=list(rule2, rule1), estimate.time = F)
  expect_equal(summary(r1)$effect.measures$ATE$estimate, -1 * summary(r2)$effect.measures$ATE$estimate, tolerance = 1e-4)
})

test_that("abar can be NULL if no Anodes, C nodes can be binary or factor", {
  r1 <- ltmle(data, Anodes=NULL, Lnodes="A", Cnodes="C", Ynodes="Y", abar=NULL, estimate.time = F)
  data2 <- data.frame(W, C=C.binary, A, Y)
  r2 <- ltmle(data2, Anodes=NULL, Lnodes="A", Cnodes="C", Ynodes="Y", abar=NULL, estimate.time = F)
  r1$call <- r2$call <- NULL
  expect_equal(r1, r2)
  r3 <- ltmleMSM(data, Anodes=NULL, Lnodes="A", Cnodes="C", Ynodes="Y", regimes=NULL, estimate.time = F, summary.measures = NULL, working.msm = "Y~1")
})

test_that("IPTW is NA if no intervention_match", {
  data.all1 <- data.frame(W, A=1, Y=rexpit(W))
  expect_warning(r <- ltmle(data.all1, Anodes="A", Ynodes="Y", abar=0, estimate.time = F), "no rows uncensored and matching regimes/abar - IPTW returns NA")
  expect_true(is.na(r$estimates["iptw"]))
})

test_that("abar is same as equivalent rule", {
  A2 <- rexpit(W)
  data2 <- data.frame(W, C, A, A2, Y)
  r1 <- ltmle(data2, Anodes=c("A", "A2"), Cnodes="C", Ynodes="Y", abar=cbind(W < 0, W > 1), estimate.time = F)
  r2 <- ltmle(data2, Anodes=c("A", "A2"), Cnodes="C", Ynodes="Y", rule=function (row) c(row["W"] < 0, row["W"] > 1), estimate.time = F)
  r1$call <- r2$call <- NULL
  expect_equal(r1, r2)
})

test_that("regimes is same as equivalent rule", {
  A2 <- rexpit(W)
  data2 <- data.frame(W, C, A, A2, Y)
  rule1 <- function (row) c(row["W"] < 0, row["W"] > 1)
  rule2 <- function (row) c(0, 1)
  rule3 <- function (row) c(1, 1)
  
  regimes <- array(1, dim=c(n, 2, 3))
  regimes[, , 1] <- cbind(W < 0, W > 1)
  regimes[, , 2] <- cbind(rep(0, n), rep(1, n))
    
  r1 <- ltmleMSM(data2, Anodes=c("A", "A2"), Cnodes="C", Ynodes="Y", regimes=regimes, estimate.time = F, summary.measures = NULL, working.msm = "Y~1")
  r2 <- ltmleMSM(data2, Anodes=c("A", "A2"), Cnodes="C", Ynodes="Y", regimes=list(rule1, rule2, rule3), estimate.time = F, summary.measures = NULL, working.msm = "Y~1")
  r1$call <- r2$call <- NULL
  expect_equal(r1, r2)
  
  r3 <- ltmle(data, Anodes="A", Cnodes="C", Ynodes="Y", abar=as.matrix(W > 0), estimate.time = F)
  r4 <- ltmle(data, Anodes="A", Cnodes="C", Ynodes="Y", rule=function (row) row$W > 0, estimate.time = F)
  r3$call <- r4$call <- NULL
  expect_equal(r3, r4)
  
  rule1 <- function (row) row["W"] < 0  
  rule2 <- function (row) 0
  rule3 <- function (row) 1
  regimes <- array(1, dim=c(n, 1, 3))
  regimes[, 1, 1] <- W < 0
  regimes[, 1, 2] <- 0
  r1 <- ltmleMSM(data, Anodes="A", Cnodes="C", Ynodes="Y", regimes=regimes, estimate.time = F, summary.measures = NULL, working.msm = "Y~1")
  r2 <- ltmleMSM(data, Anodes="A", Cnodes="C", Ynodes="Y", regimes=list(rule1, rule2, rule3), estimate.time = F, summary.measures = NULL, working.msm = "Y~1")
  r1$call <- r2$call <- NULL
  expect_equal(r1, r2)
})

test_that("NA in regimes after censoring is OK", {
  abar <- rep(NA, n)
  abar[C=="uncensored"] <- 1
 
  data2 <- data.frame(W, C, L=rnorm(n), A, Y)
  ltmle(data2, Anodes="A", Cnodes="C", Ynodes="Y", Lnodes="L", estimate.time = F, abar=AsMatrix(abar)) 
})