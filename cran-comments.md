## Test environments
* local OS X install, R 3.4.2
* win-builder (devel) [win-builder release seems to be down]

## R CMD check results
There were no ERRORs or WARNINGs. There is one NOTE:  
  New maintainer:
  Joshua Schwab <jschwab77@berkeley.edu>    [this is my new correct email]
Old maintainer(s):
  Joshua Schwab <joshuaschwab@yahoo.com>    [I have sent an email from this address confirming the change]

Found the following (possibly) invalid URLs: [see below]
  URL: http://doi.org/10.18637/jss.v081.i01
    From: man/ltmle-package.Rd
    Status: 404
    Message: Not Found

Found the following (possibly) invalid DOIs: [see below]
  DOI: 10.18637/jss.v081.i01
    From: inst/CITATION
    Status: Not Found
    Message: 404

## Downstream dependencies
I have run devtools::revdep_check on ltmle. There is 1 downstream dependency and no failures were reported.

The DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN. (This is the message JSS asked me to send.)
