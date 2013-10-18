\name{BinaryToCensoring}
\alias{BinaryToCensoring}

\title{
BinaryToCensoring
}
\description{
Helper function for creating censoring columns as factors.
}
\usage{
BinaryToCensoring(is.censored, is.uncensored)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{is.censored}{binary vector: 0=uncensored, 1=censored}
  \item{is.uncensored}{binary vector: 0=censored, 1=uncensored}
}
\details{
 Exactly one of \code{is.censored} and \code{is.uncensored} must be specified as a \emph{named} argument.
 All elements of the input vector must be 0, 1, or \code{NA}
}
\value{
 an object of class "\code{factor}" with levels "censored" and "uncensored"
}

\author{
Joshua Schwab \email{joshuaschwab@yahoo.com}
}

\seealso{
\code{\link{factor}}
}
\examples{
 BinaryToCensoring(is.censored=c(0, 1, 1, 0, NA))
 BinaryToCensoring(is.uncensored=c(1, 0, 0, 1, NA))   #the same
 
 \dontrun{
 BinaryToCensoring(c(0, 1))   #error because the input must be named
 }
}

