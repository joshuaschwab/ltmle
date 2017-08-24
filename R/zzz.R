#...
#search for fixme, browser
#run: create tests to compare versions.R (check that it works with niter=1, then run again with niter=10, then back to niter=1)
#move test-prevRelease to ltmle-dev/tests/testthat
#test("ltmle-dev")
#run: check coverage (done, but not sure about above)
#source('~/Dropbox (UC Berkeley Biostat)/Josh-Berkeley/ltmle development/check coverage.R')
#move test-prevRelease.R back to /ltmle development
#document("ltmle-dev")
#run: CheckBeforeSubmit:
#  source('~/Dropbox (UC Berkeley Biostat)/Josh-Berkeley/utils/GeneralUtilities.R')
#  CheckBeforeSubmit("ltmle-dev")
#revdep_check("ltmle-dev")
#update cran-comments.rd
#update in github
#build_win("ltmle-dev", version = "R-release")
#build_win("ltmle-dev", version = "R-devel")
#release("ltmle-dev)
#once accepted, tag in github as accepted on CRAN

release_questions <- function() {
  c('Have you updated the R version in cran-comments.rd', 'Have you run check coverage.R (in /ltmle development)', 'Have you moved test-prevRelease.R back to /ltmle development?', 'Have you run "create tests to compare versions.R" (in /ltmle development)?', 'Have you run CheckBeforeSubmit()?  (in utils/GeneralUtilities.R)') # nocov
}
