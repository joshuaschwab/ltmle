ltmle was removed from CRAN because it depended on speedglm which has been removed.
I have removed the dependency on speedglm. 

I submitted a new version to CRAN on March 29 but it received a note that checktime is 13 minutes on r-devel-windows-x86_64. I have now removed the vignettes and replaced them with articles so the check time is faster.

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔


