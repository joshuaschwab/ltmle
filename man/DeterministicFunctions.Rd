\name{deterministic.g.function_template}
\alias{MaintainTreatment}
\alias{MaintainControl}
\alias{deterministic.g.function_template}
\alias{deterministic.Q.function_template}


\title{
Deterministic g/Q functions - examples and templates
}
\description{
 \code{deterministic.g.function_template} and \code{deterministic.Q.function_template} are templates for the
 \code{deterministic.g.function} and \code{deterministic.Q.function} arguments to \code{\link{ltmle}} or \code{\link{ltmleMSM}}.
 
 \code{MaintainTreatment} and \code{MaintainControl} are two commonly used \code{deterministic.g.function}s.
}
\usage{
MaintainTreatment(data, current.node, nodes)
MaintainControl(data, current.node, nodes)
deterministic.g.function_template(data, current.node, nodes)
deterministic.Q.function_template(data, current.node, nodes, called.from.estimate.g)
}

\arguments{
  \item{data}{the 'data' data.frame passed to \code{ltmle} or \code{ltmleMSM}}
  \item{current.node}{the column index of data corresponding to the A or C node (for g) or L or Y node (for Q)}
  \item{nodes}{list of column indicies, components: 
    \itemize{
      \item \code{A} Anodes (treatment)
      \item \code{C} Cnodes (censoring)
      \item \code{L} Lnodes (time-varying covariates)
      \item \code{Y} Ynodes (events)
      \item \code{AC} Anodes and Cnodes combined and sorted
      \item \code{LY} Lnodes and Ynodes combined, sorted, "blocks" removed - see \code{\link{ltmle}}
    }
  }
  \item{called.from.estimate.g}{TRUE or FALSE - your function will be called with \code{called.from.estimate.g=TRUE} during 
     estimation of g and \code{called.from.estimate.g=FALSE} during estimation of Q.}
}
\details{
   The intended use of the templates is for the user to copy and paste the function arguments and body and then
   fill in the required sections. They will not run as-is.
   
   MaintainTreatment and MaintainControl may be passed as-is for the \code{deterministic.g.function} argument to \code{\link{ltmle}} or \code{\link{ltmleMSM}}

   Note that censoring nodes in \code{data} may be passed as binaries but they are converted to the
   preferred format of factors with levels "censored" and "uncensored" before deterministic functions are called.
   Also note that nodes may be passed to ltmle as either the names of nodes or numerical column indicies, but they
   are all converted to numerical indicies before deterministic functions are called.
   If the \code{survivalFunction} argument to \code{ltmle} or \code{ltmleMSM} is \code{TRUE}, the package automatically
   assumes that once Y jumps to 1, all future Y nodes stay 1 and treatment does not change. It is not necessary to specify this 
   in deterministic functions.
}
\value{
 A deterministic.g.function should return a list with components:
   \item{is.deterministic }{vector of logicals, length=nrow(data)}
   \item{prob1 }{the probability that data[is.deterministic, current.node] == 1, 
                  vector of length 1 or length(which(is.deterministic))}
 A deterministic.Q.function should return a list with components:
   \item{is.deterministic }{vector of logicals, length=nrow(data)}
   \item{Q.value}{the iterated expectation of the final Y, 
                  vector of length 1 or length(which(is.deterministic))}
                  
NOTE: The \code{Q.value} component is not used or required when \code{called.from.estimate.g} is \code{TRUE}
   
}
\author{
Joshua Schwab \email{joshuaschwab@yahoo.com}
}


\seealso{
\code{\link{ltmle}}, \code{\link{ltmleMSM}}
}

\examples{
 #Show a template for a deterministic.g.function
 deterministic.g.function_template
 
 #Show a template for a deterministic.Q.function
 deterministic.Q.function_template
 
 #Use MaintainControl
 set.seed(1)
 rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
 n <- 100
 W <- rnorm(n)
 A1 <- rexpit(W)
 A2 <- as.numeric(rexpit(W) | A1)  #treatment at time 1 implies treatment at time 2
 Y <- rexpit(W + A1 + A2 + rnorm(n))
 data <- data.frame(W, A1, A2, Y)

 result <- ltmle(data, Anodes=c("A1", "A2"), Ynodes="Y", abar=c(1, 1), 
   SL.library=c("SL.glm", "SL.step", "SL.mean"), 
   deterministic.g.function=MaintainTreatment)
}
