context("Test print.ltmle, print.summary.ltmle, print.ltmleEffectMeasures, print.summary.ltmleEffectMeasures")

test_that("Print methods work", {
  n <- 30
  W <- rnorm(n)
  A <- rexpit(W)
  Y <- rexpit(W + A)
  r <- ltmle(data=data.frame(W, A, Y), Anodes="A", Ynodes="Y", estimate.time = F, abar=1)
  r.eff.meas <- ltmle(data=data.frame(W, A, Y), Anodes="A", Ynodes="Y", estimate.time = F, abar=list(1, 0))
  r.gcomp <- ltmle(data=data.frame(W, A, Y), Anodes="A", Ynodes="Y", estimate.time = F, abar=1, gcomp=T)
  print(r)
  print(summary(r))
  print(r.eff.meas)
  print(summary(r.eff.meas))
  print(r, estimator="iptw")
  print(summary(r), estimator="iptw")
  print(r.eff.meas, estimator="iptw")
  print(summary(r.eff.meas), estimator="iptw")
  
  print(r.gcomp)
  print(summary(r.gcomp))
})

