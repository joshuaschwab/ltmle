context("Test that results match previous release")

expect_equals <- function(...) expect_equal(..., tolerance=0.001, scale=1, check.attributes=FALSE)

test_that("tests from 'create tests to compare versions.R'", {
  skip_on_cran() #these are slow
  data(PreviousReleaseTests)
  for (j in seq_along(btests)) {
    cat("in test_that, btests[[k]], k=", j, "\n")
    args <- btests[[j]]$args

    result <- do.call(btests[[j]]$fun, args)
    current <- btests[[j]]$compareFun(result)
    if (is.null(btests[[j]]$result$binaryOutcome)) btests[[j]]$result$binaryOutcome <- FALSE
    prev <- btests[[j]]$compareFun(btests[[j]]$result)
    expect_equals(current, prev, info=paste(btests[[j]]$info, "j = ", j, "current = ", paste(current, collapse=" "), "prev = ", paste(prev, collapse=" ")))
  }
})
