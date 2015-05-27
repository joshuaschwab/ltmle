## Test environments
* local OS X install, R 3.2.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Joshua Schwab <joshuaschwab@yahoo.com>’
(this is my correct email)  
  
## Downstream dependencies
I have run devtools::revdep_check on ltmle. There is 1 downstream dependency and no failures were reported.

## Apologies for frequent submission
I apologize for submitting this version less than one month after the previous version. Several bugs were reported by users that I failed to catch and now need to fix.

## Resubmission
I have updated to current R (3.2.0).

I believe I have fixed the problems that caused errors on solaris-sparc.  
