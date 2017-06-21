context("Test variance estimation")

test_that("TMLE based variance estimate > IC based variance estimate under positivity",{
  skip_on_cran()
  niter <- 40
  std.dev.ratio.tmle <- std.dev.ratio.iptw <- numeric(niter)
  for (i in 1:niter) {
    n <- 1000
    W <- rnorm(n)
    A <- rbinom(n, 1, plogis(5*W))
    Y <- rbinom(n, 1, 0.5)
    expect_warning({
      r.var.tmle <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE)
      r.var.ic <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, variance.method="ic")
      r.var.iptw <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, variance.method="iptw")}
      , "Variance estimate is based on")
    std.dev.ratio.tmle[i] <- summary(r.var.tmle)$treatment$std.dev / summary(r.var.ic)$treatment$std.dev
    std.dev.ratio.iptw[i] <- summary(r.var.iptw)$treatment$std.dev / summary(r.var.ic)$treatment$std.dev
  }
  expect_true(mean(std.dev.ratio.tmle) > 1.2)
  expect_true(mean(std.dev.ratio.tmle) > mean(std.dev.ratio.iptw))
  expect_true(min(std.dev.ratio.iptw) >= 1)
  expect_true(min(std.dev.ratio.tmle) >= 1)
})

test_that("All variance estimation methods run", {
  n <- 300
  W1 <- rnorm(n)
  W2 <- rnorm(n)
  A1 <- rexpit(W1 + W2)
  L1 <- rnorm(n) + A1
  C1 <- rexpit(W1 + W2)
  Y1 <- rexpit(W1 + W2 + A1)
  A2 <- rexpit(W1 + W2 + A1)
  L2 <- rnorm(n) + A1 + L1
  C2 <- rexpit(W1 + W2 + L1 + L2)
  Y2 <- rexpit(W1 + W2 - L1 - L2 + A1 - 2*A2)
  Y2[Y1 == 1] <- 1
  data <- data.frame(W1, W2, A1, L1, C1=BinaryToCensoring(is.uncensored = C1), Y1, A2, L2, C2=BinaryToCensoring(is.uncensored = C2), Y2)
  Anodes <- c("A1", "A2")
  Lnodes <- c("L1", "L2")
  Cnodes <- c("C1", "C2")
  Ynodes <- c("Y1", "Y2")
  r1 <- ltmle(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, abar=c(0,1), estimate.time = F, survivalOutcome = T)
  r2 <- ltmle(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, abar=c(0,1), estimate.time = F, survivalOutcome = T, variance.method = "iptw")
  expect_warning(r3 <- ltmle(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, abar=c(0,1), estimate.time = F, survivalOutcome = T, variance.method = "ic"), "Variance estimate is based on")
  regimes <- list(function(row) c(1,1), function(row) c(0,0))
  summary.measures <- array(1:0, dim=c(2, 1, 1))
  colnames(summary.measures) <- "A"
  r4 <- ltmleMSM(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, estimate.time = F, survivalOutcome = T, regimes=regimes, summary.measures=summary.measures, working.msm="Y~W1+W2+A", final.Ynodes = "Y2") 
  r5 <- ltmleMSM(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, estimate.time = F, survivalOutcome = T, regimes=regimes, summary.measures=summary.measures, working.msm="Y~W1+W2+A", final.Ynodes = "Y2", variance.method = "iptw") 
  expect_warning(r6 <- ltmleMSM(data=data, Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes, Ynodes=Ynodes, estimate.time = F, survivalOutcome = T, regimes=regimes, summary.measures=summary.measures, working.msm="Y~W1+W2+A", final.Ynodes = "Y2", variance.method = "ic"), "Variance estimate is based on") 
  
})

test_that("warning message prints when variance estimate ratio > 100", {
  n <- 1000
  W <- rnorm(n)
  A <- rexpit(8 * W)
  Y <- rexpit(W + A)
  data <- data.frame(W, A, Y)
  expect_warning(print(summary(ltmle(data=data, Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, gbounds = c(0, 1)))), "When this ratio is greater than 100, both variance estimates are less likely to be accurate.")
})