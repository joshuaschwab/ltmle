context("Test that results match previous release")

expect_equals <- function(...) expect_equal(..., tolerance=0.0001, scale=1, check.attributes=FALSE)

IsSurvival <- function(data) {
  finalY <- data[, ncol(data)]
  return(all(finalY %in% c(0, 1, NA)))
}

test_that("tests from 'create tests to compare versions.R'", {
  data(PreviousReleaseTests)
  for (j in seq_along(btests)) {
    additional.args <- list(survivalOutcome=IsSurvival(btests[[j]]$args$data))

    args <- c(btests[[j]]$args, additional.args)
    if (! is.null(args$regimens)) {
      args$regimes <- args$regimens #regimens were previously referred to as regimes
      args$regimens <- NULL
    }
    set.seed(1) #keep superlearner synced
    result <- btests[[j]]$compareFun(do.call(btests[[j]]$fun, args))
    
    current <- result
    prev <- btests[[j]]$result
    expect_equals(current, prev, info=paste(btests[[j]]$info, "j = ", j, "current = ", paste(current, collapse=" "), "prev = ", paste(prev, collapse=" ")))
  }
})
