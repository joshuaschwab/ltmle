context("Tests with random data")

library(tmle)

expect_equals <- function(...) expect_equal(..., tolerance=0.0001, scale=1, check.attributes=FALSE)

test_that("treatment specific mean point treatment matches Susan Gruber tmle package", {
  set.seed(NULL)
  niter <- 10
  for (i in 1:niter) {
    n <- 1000
    W <- matrix(rnorm(n), ncol=1)
    A <- rbinom(n, 1, plogis(W))
    Y <- rbinom(n, 1, plogis(A + W))
    Qform <- c("Y ~ A + W")
    gform <- "A ~ W"
    lgbound <- 0.01
    gbounds <- c(lgbound, 1)
    r1 <- tmle(Y, A, W, Qform = Qform, gform = gform, family = "binomial", Qbounds=c(0,1), alpha=0.9999, gbound=lgbound)
    
    data <- data.frame(W, A, Y)
    r2.1 <- ltmle(data, Anodes="A", Ynodes="Y", Qform=c(Y="Q.kplus1 ~ A + W"), gform=gform, abar=1, gbounds=gbounds, survivalOutcome=TRUE, estimate.time=FALSE)
    r2.0 <- ltmle(data, Anodes="A", Ynodes="Y", Qform=c(Y="Q.kplus1 ~ A + W"), gform=gform, abar=0, gbounds=gbounds, survivalOutcome=TRUE, estimate.time=FALSE)
    s <- summary(r2.1, r2.0)
    
    expect_equals(r1$estimates$ATE$psi, s$effect.measures$ATE$estimate)
    expect_equals(sqrt(r1$estimates$ATE$var.psi), s$effect.measures$ATE$std.dev)
    expect_equals(r1$estimates$ATE$pvalue, s$effect.measures$ATE$pvalue)
    
    expect_equals(r1$estimates$RR$psi, s$effect.measures$RR$estimate)
    expect_equals(r1$estimates$RR$pvalue, s$effect.measures$RR$pvalue)
    
    expect_equals(r1$estimates$OR$psi, s$effect.measures$OR$estimate)
    expect_equals(r1$estimates$OR$pvalue, s$effect.measures$OR$pvalue)
  }
  
})

test_that("simple longitudinal data matches code from Susan Gruber paper", {
  set.seed(NULL)
  niter <- 10
  for (i in 1:niter) {
    n <- 1000
    W1 <- rbinom(n,1, .5)    
    W2 <- rbinom(n, 1, .5) 
    W3 <- rnorm(n, 4, 1)
    A0 <- rbinom(n, 1, .5)   
    logitA1 <- .1 +.5*W1 + W2 - .1*(W3) + A0  
    A1 <- rbinom(n, 1, plogis(logitA1))
    L1 <- 3 + A0 - 2*W1*W2 - .5*W3  + rnorm(n, 1)
    logitA2 <- (-1.2  - .2*W2 +.1* W3 + .4*L1)
    A2 <- rbinom(n, 1, plogis(logitA2))
    A2[A0==0] <- 0   # no switching from ctl to treatment
    logitA3 <-  1.8 + .1*W2 - .05*(W3) - .4*L1 + 1.5*A2 
    A3 <- rbinom(n, 1, plogis(logitA3))
    Y <- rbinom(n, 1, plogis(3 - .3*A0 - .5*A2 + .1*W2 - .5*L1 + rnorm(n)))
    
    data <- data.frame(W1, W2, W3, A0, A1, L1, A2, A3, Y)
    data$A1[data$A0==0] <- 0
    data$A2[data$A1==0] <- 0
    data$A3[data$A2==0] <- 0
    data$Y[data$A3==0] <- 0
    data$L1[data$A1==0] <- 0
    
    Anodes <- grep("^A", names(data))
    Qform <- c(L1="Q.kplus1 ~ W1 + W2 + W3 + A0", Y="Q.kplus1 ~ A0 + A2 + W2 + L1")
    gform <- c("A0 ~ 1", "A1 ~ W1 + W2 + W3", "A2 ~ W2 + W3 + L1", "A3 ~ W2 + W3 + L1 + A2")
    lgbound <- 0.01
    
    r1 <- SuppressGivenWarnings(ltmle.sg(data, Inodes=Anodes, Lnodes=c(6, 9), Ynodes=9, Qform=Qform, gform=gform, gbd=lgbound),
                                "prediction from a rank-deficient fit may be misleading")
    r2 <- ltmle(data, Anodes=Anodes, Lnodes=6, Ynodes=9, Qform=Qform, gform=gform, abar=c(1, 1, 1, 1), stratify=TRUE, gbounds=c(lgbound, 1), survivalOutcome=TRUE, estimate.time=FALSE)
    
    expect_equals(r1["iptw"], r2$estimates["iptw"])
    expect_equals(r1["tmle"], r2$estimates["tmle"])
    expect_equals(sqrt(r1["var.tmle"]), summary(r2)$treatment$std.dev)
  }
})