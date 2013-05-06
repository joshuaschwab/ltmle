\name{summary.ltmle}
\alias{summary.ltmle}
\alias{print.ltmle}
\alias{print.summary.ltmle}
\alias{summary.ltmleMSM}
\alias{print.ltmleMSM}
\alias{print.summary.ltmleMSM}
\title{
Summarizing results from Longitudinal Targeted Maximum Likelihood Estimation (ltmle)
}
\description{
These functions are methods for class \code{ltmle} or \code{summary.ltmle} objects.
}
\usage{
\method{summary}{ltmle}(object, control.object = NULL, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...)
\method{summary}{ltmleMSM}(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...)
\method{print}{ltmle}(x, ...)
\method{print}{ltmleMSM}(x, ...)
}
\arguments{
  \item{object}{
an object of class "\code{ltmle}" or "\code{ltmleMSM}", usually a result of a call to \code{\link{ltmle}} or \code{\link{ltmleMSM}}. 
}
  \item{x}{
an object of class "\code{summary.ltmle}" or "\code{summary.ltmleMSM}", usually a result of a call to \code{\link{summary.ltmle}} or \code{\link{summary.ltmleMSM}}.
}
  \item{control.object}{
an object of class "\code{ltmle}", usually a result of a call to \code{\link{ltmle}}. If calling with \code{control.object}, \code{object} should be the 'treatment' result, while \code{control.object} should be the 'control' result
}
  \item{estimator}{
character; one of "tmle", "iptw", "gcomp", or "naive" ["naive" is available for ltmle but not ltmleMSM]. The estimator for which to get summary measures. "tmle" is valid iff the original ltmle/ltmleMSM call used gcomp=FALSE. "gcomp" is valid iff the original ltmle/ltmleMSM call used gcomp=TRUE
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
If \code{summary.ltmle} is called with \code{control.object = NULL} (the default), the parameter value of the estimator, the estimated variance, a 95 percent confidence interval, and a p-value are returned. If \code{control.object != NULL}, the additive treatment effect for \code{object} (treatment) vs \code{control.object} (control), relative risk, and odds ratio are also returned, along with the variance, confidence interval, and p-value for each.
}
\value{
\code{summary.ltmle} returns an object of class "\code{summary.ltmle}", a list with components
\item{treatment}{a list with components summarizing the estimate of \code{object}
  \itemize{
    \item \code{estimate} - the parameter estimate of \eqn{E[Y_a]}
    \item \code{std.dev} - estimated standard deviation of parameter
    \item \code{p.value} - two-sided p-value
    \item \code{CI} - vector of length 2 with 95 percent confidence interval
  }
}
\item{control}{a list with components summarizing the estimate of \code{control.object}, same components as \code{treatment}, above}
\item{effect.measures}{a list with components
  \itemize{
  \item\code{ATE} - list with components \code{est}: estimate of average treatment effect, \code{IC}: influence curve values for average treatment effect
  \item \code{RR} - list with components \code{est}: estimate of relative risk, \code{IC}: influence curve values for relative risk
  \item \code{OR} - list with components \code{est}: estimate of odds ratio, \code{IC}: influence curve values for odds ratio
  }
}
\item{treatment.call}{the matched call to \code{ltmle} for the treatment object}
\item{control.call}{the matched call to \code{ltmle} for the control object}
\item{estimator}{the \code{estimator} input argument}

\code{summary.ltmleMSM} returns an object of class "\code{summary.ltmleMSM}", a matrix with rows for each MSM parameter and columns for the point estimate, standard error, 2.5percent confidence interval, 97.5percent confidence interval, and p-value.
}

\seealso{
\code{\link{ltmle}}, \code{\link{summary}}
}
\examples{
## For examples see example(ltmle)
}

