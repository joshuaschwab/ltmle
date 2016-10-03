context("Test MaintainControl and MaintainTreatment")

test_that("MaintainControl and MaintainTreatment work", {
  n <- 100
  W <- rnorm(n)
  data <- data.frame(W, A1 = rexpit(W), A2 = rexpit(W), Y = rexpit(W))
  nodes <- list(A = 2:3, AC = 2:3)
  expect_equal(MaintainControl(data, current.node = 3, nodes), list(is.deterministic=data$A1==0, prob1=0))
  expect_equal(MaintainControl(data, current.node = 2, nodes), NULL)
  expect_equal(MaintainTreatment(data, current.node = 3, nodes), list(is.deterministic=data$A1==1, prob1=1))
  expect_equal(MaintainTreatment(data, current.node = 2, nodes), NULL)
  
})