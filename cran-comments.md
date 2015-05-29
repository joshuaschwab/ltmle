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
I have added comments for \donttest sections.

I am unable to reproduce the error reported in tests/test-all.R on either my default local OS X install or an installation configured using --disable-long-double. The error reported is:
 Error in solve.default(C.old, t(IC)) :
    system is computationally singular: reciprocal condition number = 0
but the code checks rcond(C.old) first. Nonetheless, I have added a try() around the line. If you can suggest a way I can reproduce this error without submitting to CRAN, I would be happy to do so. 
Thanks,
Josh

