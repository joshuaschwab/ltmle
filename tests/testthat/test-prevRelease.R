context("Test that results match previous release")

expect_equals <- function(...) expect_equal(..., tolerance=0.001, scale=1, check.attributes=FALSE)

if (T) {
test_that("tests from 'create tests to compare versions.R'", {
  load("~/Dropbox (UC Berkeley Biostat)/Josh-Berkeley/ltmle development/PreviousReleaseTests.RData")
  for (j in seq_along(btests)) {
    args <- btests[[j]]$args
    withr::with_envvar(c("LTMLE.REGRESSION.TESTING.SEED" = 12345), result <- do.call(btests[[j]]$fun, args))
    current <- btests[[j]]$compareFun(result)
    
    if (is.null(btests[[j]]$result$binaryOutcome)) btests[[j]]$result$binaryOutcome <- FALSE
    prev <- btests[[j]]$compareFun(btests[[j]]$result)
    expect_equals(current, prev, info=paste(btests[[j]]$info, "j = ", j, "current = ", paste(current, collapse=" "), "prev = ", paste(prev, collapse=" ")))
  }
})

test_that("estimate time tests work", {
  #slow tests - run to get 100% coverage once everything else is working
  cat("- est time - ", date(), "\n", file = "~/Dropbox (UC Berkeley Biostat)/Josh-Berkeley/ltmle development/log.txt", append = F)
  EstTime <- function(n) {
    cat("- in function - ", date(), "\n", file = "~/Dropbox (UC Berkeley Biostat)/Josh-Berkeley/ltmle development/log.txt", append=T)
    end.time <- 3
    alive <- uncensored <- rep(T, n)
    prevA <- rexpit(rnorm(n) + 2)
    W <- rnorm(n)
    data <- data.frame(W)
    for (t in 1:end.time) {
      L <- A <- C <- Y <- rep(NA, n)
      L[alive & uncensored] <- -1 + rnorm(sum(alive & uncensored)) + prevA[alive & uncensored]
      A[prevA == 0 & alive & uncensored] <- 0
      A[prevA == 1 & alive & uncensored] <- rexpit(2 + W[prevA == 1 & alive & uncensored] + L[prevA == 1 & alive & uncensored])
      prevA <- A
      C[alive & uncensored] <- rexpit(-3 + -2*L[alive & uncensored] - A[alive & uncensored])
      uncensored[alive] <- uncensored[alive] & !C[alive]
      Y[alive & uncensored] <- rexpit(-3 + W[alive & uncensored] + L[alive & uncensored] + A[alive & uncensored])
      Y[!alive] <- 1
      alive[uncensored] <- alive[uncensored] & Y[uncensored]==0
      data <- data.frame(data, L, A, BinaryToCensoring(is.censored = C), Y)
      names(data)[ncol(data)] <- paste0("Y", t)
      names(data)[ncol(data) - 1] <- paste0("C", t)
      names(data)[ncol(data) - 2] <- paste0("A", t)
      names(data)[ncol(data) - 3] <- paste0("L", t)
    }
    regimes <- array(1, dim=c(n, end.time, 2))
    for (ii in 1:2) {
      for (t in 2:end.time) {
        regimes[, t, ii] <- regimes[, t - 1, ii] * rexpit(rnorm(n) + ii) 
      }
    }
    summary.measures <- array(c(1, 2, 1, 2), dim=c(2, 1, 2))
    colnames(summary.measures) <- "R"
    r <- ltmleMSM(data=data, Anodes=paste0("A", 1:end.time), Cnodes=paste0("C", 1:end.time), Lnodes=paste0("L", 2:end.time), Ynodes=paste0("Y", 1:end.time), final.Ynodes=c("Y1", paste0("Y", end.time)), regimes=regimes, summary.measures=summary.measures, working.msm="Y~R", estimate.time=T, survivalOutcome = T)
    return(r)
  }
  if (T) {
    expect_equal(class(EstTime(4000)), "ltmleMSM")
    expect_equal(class(EstTime(10000)), "ltmleMSM")
  }
})

test_that("invalid inputs throw errors SL", {
  #some of these won't work in check("ltmle-dev") because they depend on libraries for SL
  set.seed(2)
  n <- 100
  W <- rnorm(n)
  A <- rexpit(W)
  Y <- rexpit(W + A)
  data <- data.frame(W, W2=rnorm(n), A1=A, Y1=Y, A2=A, Y2=as.numeric(rexpit(W) | Y))
  print(head(data))
  expect_warning(ltmle(data, Anodes = c("A1", "A2"), Ynodes = c("Y1", "Y2"), estimate.time = FALSE, abar = c(1, 1), SL.library = list(g=NULL, Q=c("SL.mean", "SL.knn")), survivalOutcome = T), "SuperLearner returned predicted.values > 1 or < 0")
  expect_error(ltmle(data, Anodes = c("A1", "A2"), Ynodes = c("Y1", "Y2"), estimate.time = FALSE, abar = c(1, 1), SL.library = list(g=NULL, Q=c("SL.xxxx")), survivalOutcome = T), "Error occured during call to SuperLearner")
  expect_error(ltmle(data, Anodes = c("A1", "A2"), Ynodes = c("Y1", "Y2"), estimate.time = FALSE, abar = c(1, 1), SL.library = list(g=NULL, Q=c("SL.glmnet")), survivalOutcome = T), "Note that some SuperLeaner libraries crash when called with continuous dependent variables")
  expect_error(ltmle(data, Anodes = c("A1", "A2"), Ynodes = c("Y1", "Y2"), estimate.time = FALSE, abar = c(1, 1), SL.library = list(g=NULL, Q=c("SL.mean", "SL.cforest")), survivalOutcome = T), "SuperLearner returned all NAs")
  expect_error(ltmle(data, Anodes = c("A1", "A2"), Ynodes = c("Y1", "Y2"), estimate.time = FALSE, abar = c(1, 1), survivalOutcome = T, deterministic.Q.function = function(data, ...) list(is.deterministic=rep(T, nrow(data)), Q.value=0.5)), "inconsistent deterministic Q at node")
})
}