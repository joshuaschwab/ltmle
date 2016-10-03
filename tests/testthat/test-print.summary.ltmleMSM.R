context("Test print.summary.ltmleMSM") #also print.ltmleMSM

test_that("transformOutcome NOTE is printed by print.summary.ltmleMSM", {
  
  data(sampleDataForLtmleMSM)
  Anodes <- grep("^A", names(sampleDataForLtmleMSM$data))
  Lnodes <- c("CD4_1", "CD4_2")
  Ynodes <- grep("^Y", names(sampleDataForLtmleMSM$data))
  
  data <- sampleDataForLtmleMSM$data
  data[,Ynodes] <- matrix(runif(length(Ynodes)*nrow(data)), nrow(data))*10
  data <- as.data.frame(apply(data, 2, function(x) {
    x[is.na(x)] <- max(x, na.rm=TRUE)
    x
  }))
  
  expect_warning(
    result <- ltmleMSM(data, Anodes=Anodes, Lnodes=Lnodes,
                       Ynodes=Ynodes, survivalOutcome=FALSE,
                       regimes=sampleDataForLtmleMSM$regimes, 
                       summary.measures=sampleDataForLtmleMSM$summary.measures,
                       final.Ynodes=Ynodes, 
                       working.msm="Y ~ time + I(pmax(time - switch.time, 0))", 
                       estimate.time=FALSE, variance.method="ic")
    , "Variance estimate is based on")
  
  expect_that(print(result), prints_text("NOTE: The MSM is modeling the transformed outcome"))
  expect_that(print(summary(result)), prints_text("NOTE: The MSM is modeling the transformed outcome"))    
  
  expect_warning(
    result2 <- ltmleMSM(data, Anodes=Anodes, Lnodes=Lnodes,
                        Ynodes=Ynodes, survivalOutcome=FALSE,
                        regimes=sampleDataForLtmleMSM$regimes, 
                        summary.measures=sampleDataForLtmleMSM$summary.measures,
                        final.Ynodes=Ynodes, 
                        working.msm="Y ~ time + I(pmax(time - switch.time, 0))", 
                        estimate.time=FALSE, variance.method="ic", gcomp=TRUE)
    , "Variance estimate is based on")
  expect_that(print(result2), prints_text("NOTE: The MSM is modeling the transformed outcome"))
  expect_that(print(summary(result2)), prints_text("Warning: inference for gcomp is not accurate!"))
  
  expect_warning(
    result3 <- ltmle(data, Anodes=Anodes, Lnodes=Lnodes,
                     Ynodes=Ynodes, survivalOutcome=FALSE,
                     abar=list(c(1,1,1), c(0,0,0)), 
                     estimate.time=FALSE, gcomp=TRUE, variance.method = "ic")
    , "Variance estimate is based on")
  expect_that(print(summary(result3)), prints_text("Warning: inference for gcomp is not accurate!"))
})