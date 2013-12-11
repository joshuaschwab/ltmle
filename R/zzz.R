.onAttach <- function(...) {
  packageStartupMessage(
"The influence curve based variance estimator implemented in this package provides 
asymptotically conservative standard error estimates when the treatment mechanism is 
estimated with MLE using a correctly specific parametric model. However, we have found 
in simulations that in finite samples this variance estimator may be substantially
anti-conservative in settings where the outcome is rare or there are positivity 
violations. An improved variance estimator addressing these challenges is currently 
under development and will be released in a subsequent version of the package.\n")
}