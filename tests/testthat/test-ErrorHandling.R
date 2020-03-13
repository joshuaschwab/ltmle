context("error handling")

test_that("invalid inputs throw errors 1", {
  n <- 50
  W <- rnorm(n)
  A <- rexpit(W)
  Y <- rexpit(W + A)
  data <- data.frame(W, A, Y)
  expect_error(ltmle(data, Anodes = "A", Ynodes = "Y", estimate.time = FALSE, abar = function (row) 1), "abar should be a vector or matrix, not a function")
  expect_error(ltmle(data, Anodes = "A", Ynodes = "Y", estimate.time = FALSE, abar = 1, variance.method = TRUE), "variance.method must be one of")
  expect_error(ltmle(data, Anodes = "A", Ynodes = "Y", estimate.time = FALSE, abar = 1, gform = "W ~ A"), "The LHS of gform")
  expect_error(ltmle(data, Anodes = "AA", Ynodes = "Y", estimate.time = FALSE, abar = 1), "named node")
  expect_error(ltmleMSM(data, Anodes = "A", Ynodes = "Y", estimate.time = FALSE, working.msm = "Y~1", summary.measures = NULL, regimes = NULL), "regimes must not be NULL")
  expect_error(ltmleMSM(data, Anodes = "A", Ynodes = "Y", estimate.time = FALSE, working.msm = "Y~1", summary.measures = NULL, regimes = matrix(1, n, 1)), "regimes must be an array with 3 dimensions")
  expect_error(expect_message(ltmle(data, Anodes = "A", Ynodes = "Y", estimate.time = T, abar = 1, deterministic.g.function = function() 1), "Timing estimate unavailable"))
})

test_that("invalid inputs throw errors 2", {
  data(sampleDataForLtmleMSM)
  Anodes <- grep("^A", names(sampleDataForLtmleMSM$data))
  Lnodes <- c("CD4_1", "CD4_2")
  Ynodes <- grep("^Y", names(sampleDataForLtmleMSM$data))
  s <- sampleDataForLtmleMSM$summary.measures
  s[, 1, ] <- s[, 2, ]  
  expect_warning(r <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=T, regimes=sampleDataForLtmleMSM$regimes, summary.measures=s, final.Ynodes=Ynodes, working.msm="Y ~ time + switch.time", estimate.time=FALSE, variance.method="ic"), "rcond")
  expect_warning(print(summary(r)), "Unable to compute standard errors")
  expect_error(ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=T, regimes=sampleDataForLtmleMSM$regimes, summary.measures=NULL, working.msm="Y~1", msm.weights = array(1, dim=c(1, 1, 1))), "msm.weights must be")
})

test_that("invalid inputs throw errors 3", {
  n <- 100
  W <- rnorm(n)
  A <- rexpit(W)
  Y <- rexpit(W + A)
  data <- data.frame(W, A1=rep(1, n), A2=A, Y)
  expect_error(ltmle(data, Anodes = c("A1", "A2"), Ynodes = "Y", estimate.time = FALSE, abar = c(0, 1), deterministic.g.function = function (data, current.node, nodes) if (current.node < 3) NULL else list(is.deterministic=data$A1==0, prob1=0)), "This error occured while calling deterministic.g.function on data where Anodes are set to abar")
  
  data <- data.frame(W, C=BinaryToCensoring(is.censored = rexpit(W)), A, Y)
  abar <- matrix(1, n)
  abar[1:floor(n/2)] <- NA
  expect_error(ltmle(data, Anodes="A", Cnodes="C", Ynodes="Y", abar=abar, estimate.time = FALSE), "NA in regimes/abar not allowed")
  
  abar <- matrix(1, n)
  abar[min(which(data$C=="censored"))] <- NA
  expect_warning(expect_error(ltmle(data, Anodes="A", Cnodes="C", Ynodes="Y", abar=abar, estimate.time = FALSE), "NA in tempC"), "NA in regimes/abar before the first L/Y node will probably cause an error")

  expect_error(ltmle(data.frame(W, C=rep("uncensored", n), A, Y, stringsAsFactors = F), abar=1, Anodes="A", Cnodes="C", Ynodes="Y"), "in data, all Cnodes should be factors with two levels")
  
  expect_error(ltmle(data.frame(W, C=rep("???", n), A, Y, stringsAsFactors = T), abar=1, Anodes="A", Cnodes="C", Ynodes="Y"), "all levels of data")
})

test_that("rarely used lines run", {
  n <- 20
  data <- data.frame(Y=rexpit(rnorm(n)), X1=rep(1, n), X2=rep(2, n))
  ltmle.glm("Y ~ X1 + X2", family=binomial(), data = data, weights = NULL)
  
  expect_warning(MakePosDef(matrix(-c(1, 2, 2, 4), 2, 2)), "Covariance matrix from EstimateVariance not positive definite")
  
  is.deterministic <- matrix(T, n, 1)
  deterministic.Q <- matrix(0, n, 1)
  h.g.ratio <- array(1, dim=c(n, 1, 1))
  intervention.match <- matrix(T, n, 1)
  regimes.with.positive.weight <- 1
  off <- rep(0, n)
  X <- matrix(1, n, 1)
  Qstar.kplus1 <- matrix(runif(n), n, 1)
  Qstar <- matrix(runif(n), n, 1)
  
  uncensored <- rep(T, n)
  expect_error(FixScoreEquation(Qstar.kplus1, h.g.ratio, uncensored, intervention.match, is.deterministic, deterministic.Q, off, X, regimes.with.positive.weight), "minimizer failed")
})