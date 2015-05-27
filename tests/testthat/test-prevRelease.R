context("Test that results match previous release")

expect_equals <- function(...) expect_equal(..., tolerance=0.0001, scale=1, check.attributes=FALSE)

IsSurvival <- function(data) {
  finalY <- data[, ncol(data)]
  return(all(finalY %in% c(0, 1, NA)))
}

VersionHasNewVariance <- function(ver) {
  if (ver %in% c("0.9", "0.9.3")) return(FALSE)
  if (ver %in% c("0.9.4")) return(TRUE) #this hasn't actually been used yet
  stop("unexpected ver")
}


test_that("tests from 'create tests to compare versions.R'", {
  skip_on_cran() #these are slow
  data(PreviousReleaseTests)
  for (j in seq_along(btests)) {
    additional.args <- list(survivalOutcome=IsSurvival(btests[[j]]$args$data), IC.variance.only=!VersionHasNewVariance(btests[[j]]$ver))

    args <- c(btests[[j]]$args, additional.args)
    args$regimes <- args$regimens #regimens were previously referred to as regimes
    args$regimens <- NULL
    args$mhte.iptw <- NULL #mhte.iptw was removed
    prev.seed <- .Random.seed
    set.seed(1) #keep superlearner (and rnorm in FixScoreEquation) synced 
    result <- do.call(btests[[j]]$fun, args)
    .Random.seed <<- prev.seed
    current <- btests[[j]]$compareFun(result)
    
    if (is.null(btests[[j]]$result$binaryOutcome)) btests[[j]]$result$binaryOutcome <- FALSE
    prev <- btests[[j]]$compareFun(btests[[j]]$result)
    expect_equals(current, prev, info=paste(btests[[j]]$info, "j = ", j, "current = ", paste(current, collapse=" "), "prev = ", paste(prev, collapse=" ")))
  }
})
