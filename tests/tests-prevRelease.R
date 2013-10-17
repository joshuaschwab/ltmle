expect_equals <- function(...) expect_equal(..., tolerance=0.0001, scale=1, check.attributes=FALSE)

IsSurvival <- function(data) {
  finalY <- data[, ncol(data)]
  return(all(finalY %in% c(0, 1, NA)))
}

test_that("results match previous release", {
  data(PreviousReleaseTests)
  for (j in seq_along(btests)) {
    additional.args <- list(survivalOutcome=IsSurvival(btests[[j]]$args$data))

    args <- c(btests[[j]]$args, additional.args)
    set.seed(1) #keep superlearner synced
    result <- do.call(btests[[j]]$fun, args)
    
    current <- btests[[j]]$compareFun(result)
    prev <- btests[[j]]$compareFun(btests[[j]]$result)
    expect_equals(current, prev, info=paste(btests[[j]]$info, "j = ", j, "current = ", paste(current, collapse=" "), "prev = ", paste(prev, collapse=" ")))
  }
})
