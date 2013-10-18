# Helper function to generate censoring as factors
#' @export
BinaryToCensoring <- function(is.censored, is.uncensored) {
  if (! xor(missing(is.censored), missing(is.uncensored))) stop("exactly one of is.censored and is.uncensored must be passed")
  calling.name <- names(sys.call(0))[2]
  if (length(calling.name) == 0 || ! calling.name %in% c("is.censored", "is.uncensored")) stop("the argument to BinaryToCensoring must be completely named - see ?BinaryToCensoring")
  if (missing(is.uncensored)) {
    is.uncensored <- ! is.censored
  }
  if (! all(is.uncensored %in% c(0, 1, NA))) stop("the argument to BinaryToCensoring should be binary (0, 1, or NA) or logical")
  y <- character(length(is.uncensored))
  y[is.uncensored == 0] <- "censored"
  y[is.uncensored == 1] <- "uncensored"
  y[is.na(is.uncensored)] <- NA
  return(factor(y))
}