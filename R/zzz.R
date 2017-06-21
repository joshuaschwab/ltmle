#...
#search for fixme, browser
#run: create tests to compare versions.R (check that it works with niter=1, then run again with niter=10)
#move test-prevRelease to ltmle-dev/tests/testthat
#test("ltmle-dev")
#run: check coverage 
#source('~/Dropbox (UC Berkeley Biostat)/Josh-Berkeley/ltmle development/check coverage.R'
#move test-prevRelease.R back to /ltmle development
#run: CheckBeforeSubmit
#update cran-comments.rd
#update in github
#build_win
#submit
#once accepted, tag in github as accepted on CRAN

release_questions <- function() {
  c('Have you updated the R version in cran-comments.rd', 'Have you run check coverage.R (in /ltmle development)', 'Have you moved test-prevRelease.R back to /ltmle development?', 'Have you run "create tests to compare versions.R" (in /ltmle development)?', 'Have you run CheckBeforeSubmit()?  (in utils/GeneralUtilities.R)') # nocov
}
