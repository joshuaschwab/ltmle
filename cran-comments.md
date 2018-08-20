## Test environments
* local OS X install, R 3.5.0
* win-builder (devel and release) 

## R CMD check results
There were no ERRORs or WARNINGs. There is one NOTE - the previous version was archived on CRAN because I did not include a package used in /tests in Suggests. This is now fixed. My apologies for not fixing this earlier.  
   
## Downstream dependencies
I have run devtools::revdep_check on ltmle. There is 1 downstream dependency and no failures were reported.

