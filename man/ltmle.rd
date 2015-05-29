\name{ltmle}
\alias{ltmle}
\alias{ltmleMSM}
\title{
Longitudinal Targeted Maximum Likelihood Estimation
}
\description{
\code{ltmle} is Targeted Maximum Likelihood Estimation (TMLE) of treatment/censoring specific mean outcome for point-treatment and longitudinal data. \code{ltmleMSM} adds Marginal Structural Models. Both always provide Inverse Probability of Treatment/Censoring Weighted estimate (IPTW) as well. Maximum likelihood based G-computation estimate (G-comp) can be obtained instead of TMLE. \code{ltmle} can be used to calculate additive treatment effect, risk ratio, and odds ratio.
}
\usage{
ltmle(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=NULL, Qform=NULL, 
 gform=NULL, abar, rule=NULL, gbounds=c(0.01, 1), Yrange=NULL, 
 deterministic.g.function=NULL, stratify=FALSE, SL.library=NULL, 
 estimate.time=TRUE, gcomp=FALSE, iptw.only=FALSE, 
 deterministic.Q.function=NULL, IC.variance.only=FALSE, observation.weights=NULL)
ltmleMSM(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=NULL, Qform=NULL,
 gform=NULL, gbounds=c(0.01, 1), Yrange=NULL, deterministic.g.function=NULL, 
 SL.library=NULL, regimes, working.msm, summary.measures, 
 final.Ynodes=NULL, stratify=FALSE, msm.weights="empirical", 
 estimate.time=TRUE, gcomp=FALSE,  
 iptw.only=FALSE, deterministic.Q.function=NULL, memoize=TRUE, 
 IC.variance.only=FALSE, observation.weights=NULL)
}
\arguments{
  \item{data}{data frame following the time-ordering of the nodes. See 'Details'.}
  \item{Anodes}{column names or indicies in \code{data} of treatment nodes}
  \item{Cnodes}{column names or indicies in \code{data} of censoring nodes}
  \item{Lnodes}{column names or indicies in \code{data} of time-dependent covariate nodes}
  \item{Ynodes}{column names or indicies in \code{data} of outcome nodes}
  \item{survivalOutcome}{If \code{TRUE}, then Y nodes are indicators of an event, and if Y at some time point is 1, then all following should be 1. Required to be \code{TRUE} or \code{FALSE} if outcomes are binary and there are multiple Ynodes.}
  \item{Qform}{character vector of regression formulas for \eqn{Q}. See 'Details'.}
  \item{gform}{character vector of regression formulas for \eqn{g} or a matrix/array of prob(A=1). See 'Details'.}
  \item{abar}{binary vector (numAnodes x 1) or matrix (n x numAnodes) of counterfactual treatment or a list of length 2. See 'Details'.}
  \item{rule}{a function to be applied to each row (a named vector) of \code{data} that returns a numeric vector of length numAnodes or a list of length 2. See 'Details'.}
  \item{gbounds}{lower and upper bounds on estimated cumulative probabilities for g-factors. Vector of length 2, order unimportant.}
  \item{Yrange}{NULL or a numerical vector where the min and max of \code{Yrange} specify the range of all Y nodes. See 'Details'.}
  \item{deterministic.g.function}{optional information on A and C nodes that are given deterministically. See 'Details'. Default \code{NULL} indicates no deterministic links.}
  \item{stratify}{if \code{TRUE} stratify on following \code{abar} when estimating Q and g. If \code{FALSE}, pool over \code{abar}.}
  \item{SL.library}{optional character vector of libraries to pass to \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. \code{NULL} indicates \link{glm} should be called instead of \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. '\code{default}' indicates a standard set of libraries. May be separately specified for \eqn{Q} and \eqn{g}. See 'Details'.}
  \item{estimate.time}{if \code{TRUE}, run an initial estimate using only 50 observations and use this to print a very rough estimate of the total time to completion. No action if there are fewer than 50 observations.}
  \item{gcomp}{if \code{TRUE}, run the maximum likelihood based G-computation estimate \emph{instead} of TMLE}
  \item{regimes}{binary array: n x numAnodes x numRegimes of counterfactual treatment or a list of 'rule' functions}
  \item{working.msm}{character formula for the working marginal structural model}
  \item{summary.measures}{array: num.regimes x num.summary.measures x num.final.Ynodes - measures summarizing the regimes that will be used on the right hand side of \code{working.msm} (baseline covariates may also be used in the right hand side of \code{working.msm} and do not need to be included in \code{summary.measures})}
  \item{final.Ynodes}{vector subset of Ynodes - used in MSM to pool over a set of outcome nodes}
  \item{msm.weights}{projection weights for the working MSM. If "empirical", weight by empirical proportions of rows matching each regime for each final.Ynode, with duplicate regimes given zero weight. If \code{NULL}, no weights. Or an array of user-supplied weights with dimensions c(n, num.regimes, num.final.Ynodes) or c(num.regimes, num.final.Ynodes).}
  \item{iptw.only}{by default (\code{iptw.only = FALSE}), both TMLE and IPTW are run in \code{ltmle} and \code{ltmleMSM}. If \code{iptw.only = TRUE}, only IPTW is run, which is faster.}
  \item{deterministic.Q.function}{optional information on Q given deterministically. See 'Details'. Default \code{NULL} indicates no deterministic links.}
  \item{memoize}{If \code{TRUE}, glm regressions will be memoized. It is recommended to leave this as \code{TRUE} (the default), especially if there are multiple \code{final.Ynodes}, because the code is not written as efficiently as it should be and will end up repeating the same glm call. Will be fixed in a future release.}
  \item{observation.weights}{observation (sampling) weights. Vector of length n. If \code{NULL}, assumed to be all 1.}
  \item{IC.variance.only}{If \code{FALSE}, compute both the robust variance estimate using TMLE and the influence curve based variance estimate (use the larger of the two). If \code{TRUE}, only compute the influence curve based variance estimate, which is faster, but may be substantially anti-conservative if there are positivity violations or rare outcomes.}
}
\details{
  The estimates returned by \code{ltmle} are of a treatment specific mean, \eqn{E[Y_{\bar{a}}]}, the mean of the final treatment node, where all treatment nodes, \eqn{A}, are set to \eqn{\bar{a}} (\code{abar}) and all censoring nodes \eqn{C} are set to 1 (uncensored). The estimates returned by \code{ltmleMSM} are similar but are the parameters in a working marginal structural model.
    
  \code{data} should be a data frame where the order of the columns corresponds to the time-ordering of the model. 
  \itemize{
      \item in censoring columns (Cnodes): factor with two levels: "censored" and "uncensored". The helper function \code{CensoringToBinary} can be used to create these factors.
      \item in treatment columns (Anodes): 1 = treated, 0 = untreated (must be binary)
      \item in event columns (Ynodes): If \code{survivalOutcome} is \code{TRUE}, then Y nodes are treated as indicators of a one-time event. See details for \code{survivalOutocme}. If \code{survivalOutcome} is \code{FALSE}, Y nodes are treated as binary if all values are 0 or 1, and are treated as continuous otherwise. If Y nodes are continuous, they may be automatically scaled. See details for \code{Yrange}.
      \item time-dependent covariate columns (Lnodes): can be any numeric data
      \item  Data in \code{Cnodes}, \code{Anodes}, \code{Lnodes} and \code{Ynodes} are not used after (to the right of) censoring (or an event when \code{survivalOutcome==TRUE}) and may be coded as \code{NA} or any other value.
      \item Columns in \code{data} that are before (to the left of) the first of \code{Cnodes} or \code{Anodes} are treated as baseline variables, even if they are specified as \code{Lnodes}. 
      \item After the first of \code{Cnodes}, \code{Anodes}, \code{Ynodes}, or \code{Lnodes}, every column must be in one of \code{Cnodes}, \code{Anodes}, \code{Ynodes}, or \code{Lnodes}. 
      }
  
  If \code{survivalOutcome} is \code{TRUE},  all Y values are indicators of an event (e.g. death) at or before the current time, where 1 = event and 0 = no event. The events in Ynodes must be of the form where once Y jumps to 1, Y remains 1 at subsequent nodes. 
  
  For continuous outcomes, (\code{survivalOutcome==FALSE} and some Y nodes are not 0 or 1,) Y values are truncated at the minimum and maximum of \code{Yrange} if specified, and then transformed and scaled to be in [0,1]. That is, transformed to \code{(Y-min(Yrange))/(max(Yrange)-min(Yrange))}. If \code{Yrange} is \code{NULL}, it is set to the range of all Y nodes. In that case, Y nodes are only scaled if any values fall outside of [0,1]. For intervention specific means (\code{ltmle}), parameter estimates are transformed back based \code{Yrange}. 
 
  
  \code{Qform} should be \code{NULL}, in which case all parent nodes of each L and Y node will be used as regressors, or a named character vector that can be coerced to class "\code{\link{formula}}". The length of \code{Qform} must be equal to \code{length(Lnodes) + length(Ynodes)}** and the names and order of the formulas must be the same as the names and order of the L and Y nodes in \code{data}. The left hand side of each formula should be "\code{Q.kplus1}". If \code{SL.library} is \code{NULL}, \code{glm} will be called using the elements of \code{Qform}. If \code{SL.library} is specified, \code{\link[SuperLearner:SuperLearner]{SuperLearner}} will be called and all variables appearing on the right hand side of a formula will be passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. The actual functional form of the formula is unimportant if \code{\link[SuperLearner:SuperLearner]{SuperLearner}} is called.
  
  ** If there is a "block" of L and Y nodes not separated by A or C nodes, only one regression is required at the first L/Y node in a block. You can pass regression formulas for the other L/Y nodes, but they will be ignored (with a message). See example 5.
  
  \code{gform} should be \code{NULL}, in which case all parent nodes of each L and Y node will be used as regressors, or a character vector that can be coerced to class "\code{\link{formula}}", or a matrix/array of Prob(A=1). If \code{gform} is a character vector, the length of \code{gform} must be equal to \code{length(Anodes) + length(Cnodes)} and the order of the formulas must be the same as the order the A and C nodes appear in \code{data}. The left hand side of each formula should be the name of the Anode or Cnode. If \code{SL.library} is \code{NULL}, \code{glm} will be called using the elements of \code{gform}. If \code{SL.library} is specified, \code{\link[SuperLearner:SuperLearner]{SuperLearner}} will be called and all variables appearing on the right hand side of a formula will be passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. The actual functional form of the formula is unimportant if \code{\link[SuperLearner:SuperLearner]{SuperLearner}} is called. 
  
  In \code{ltmle}, \code{gform} can also be a n x numACnodes matrix where entry (i, j) is the probability that the ith observation of the jth A/C node is 1 (if an Anode) or uncensored (if a Cnode), conditional on following abar up to that node. In \code{ltmleMSM}, \code{gform} can similarly be a n x numACnodes x numRegimes array, where entry (i, j, k) is the probability that the ith observation of the jth A/C node is 1 (if an Anode) or uncensored (if a Cnode), conditional on following regime k up to that node. If \code{gform} is a matrix/array, \code{deterministic.g.function} will not be used and should be \code{NULL}.
  
  \code{abar} specifies the counterfactual values of the Anodes, using the order they appear in \code{data} and should have the same length (if abar is a vector) or number of columns (if abar is a matrix) as \code{Anodes}. 

  \code{rule} can be used to specify a dynamic treatment rule. \code{rule} is a function applied to each row of \code{data} which returns the a numeric vector of the same length as \code{Anodes}.
  
  \code{abar} and \code{rule} cannot both be specified. If one of them if a list of length 2, additive treatment effect, risk ratio, and odds ratio can be computed using \code{\link{summary.ltmleEffectMeasures}}. 

  \code{regimes} can be a binary array: n x numAnodes x numRegimes of counterfactual treatment or a list of 'rule' functions as described above for the \code{rule} parameter for the \code{ltmle} function
  
  \code{deterministic.g.function} can be a function used to specify model knowledge about value of Anodes and/or Cnodes that are set deterministically. For example, it may be the case that once a patient starts treatment, they always stay on treatment. For details on the form of the function and examples, see \code{\link{deterministic.g.function_template}}

  \code{deterministic.Q.function} can be a function used to specify model knowledge about the final event state. For example, it may be the case that a patient can complete the study at some intermediate time point, in which case the probability of death is 0 (assuming they have not died already). For details on the form of the function and examples, see \code{\link{deterministic.Q.function_template}}

\code{SL.library} may be a character vector of libraries (or \code{NULL} or '\code{default}'), in which case these libraries are     used to estimate both \eqn{Q} and \eqn{g} OR a list with two components, \code{Q} and \code{g}, where each is a character vector of libraries (or \code{NULL} or '\code{default}'). 
\code{NULL} indicates \link{glm} should be called instead of \code{\link[SuperLearner:SuperLearner]{SuperLearner}}
If \code{SL.library} is the string '\code{default}', \code{SL.library} is set to \code{list("SL.glm", "SL.stepAIC", "SL.bayesglm", c("SL.glm", "screen.corP"), c("SL.step", "screen.corP"), c("SL.step.forward", "screen.corP"), c("SL.stepAIC", "screen.corP"), c("SL.step.interaction", "screen.corP"), c("SL.bayesglm", "screen.corP")}. 
Note that the default set of libraries consists of main terms models. It may be advisable to include squared terms, interaction terms, etc in \code{data} or include libraries that consider non-linear terms.

The print method for \code{ltmle} objects only prints the tmle estimates. 
}
\value{
\code{ltmle} returns an object of class "\code{ltmle}" (unless \code{abar} or \code{rule} is a list, in which case it returns an object of class \code{ltmleSummaryMeasures}, which has the same components as \code{ltmleMSM}.)
The function \code{\link{summary}} (i.e. \code{\link{summary.ltmle}}) can be used to obtain or print a summary of the results.
An object of class "\code{ltmle}" is a list containing the following components:
\item{estimates}{a named vector of length 4 with elements, each an estimate of \eqn{E[Y_{bar{a}}]}:
  \itemize{
  \item \code{tmle} - Targeted Maximum Likelihood Estimate [NULL if \code{gcomp} is \code{TRUE}]
  \item \code{iptw} - Inverse Probability of Treatment/Censoring Weighted estimate 
  \item \code{gcomp} - maximum likelihood based G-computation estimate [NULL if \code{gcomp} is \code{FALSE}]
  }
}
\item{IC}{a list with the following components of Influence Curve values}
  \itemize{
  \item \code{tmle} - vector of influence curve values for Targeted Maximum Likelihood Estimate [NULL if \code{gcomp} is \code{TRUE}]
  \item \code{iptw} - vector of influence curve values for Inverse Probability of Treatment/Censoring Weighted estimate 
  \item \code{gcomp} - vector of influence curve values for Targeted Maximum Likelihood Estimate without updating [NULL if \code{gcomp} is \code{FALSE}]
  }
\item{cum.g}{cumulative g, after bounding: for ltmle, n x numACnodes, for ltmleMSM, n x numACnodes x num.regimes}
\item{cum.g.unbounded}{cumulative g, before bounding: for ltmle, n x numACnodes, for ltmleMSM, n x numACnodes x num.regimes}
\item{call}{the matched call}
\item{gcomp}{the \code{gcomp} input}
\item{formulas}{a \code{list} with elements \code{Qform} and \code{gform}}
\item{fit}{a list with the following components}
  \itemize{
  \item \code{g} - list of length numACnodes - \code{glm} or \code{SuperLearner} return objects from fitting g regressions
  \item \code{Q} - list of length numLYnodes - \code{glm} or \code{SuperLearner} return objects from fitting Q regressions
  \item \code{Qstar} - list of length numLYnodes - \code{glm} (or numerical optimization if \code{glm} fails to solve the score equation) return objects from updating the Q fit
  }

\code{ltmleMSM} returns an object of class "\code{ltmleMSM}"
The function \code{\link{summary}} (i.e. \code{\link{summary.ltmleMSM}}) can be used to obtain or print a summary of the results.
An object of class "\code{ltmleMSM}" is a list containing the following components:
\item{beta}{parameter estimates for working.msm using TMLE (GCOMP if \code{gcomp} input is \code{TRUE})}
\item{beta.iptw}{parameter estimates for working.msm using IPTW}
\item{IC}{matrix, n x numBetas - influence curve values for TMLE (without updating if \code{gcomp} input is \code{TRUE})}
\item{IC.iptw}{matrix, n x numBetas - influence curve values for IPTW}
\item{msm}{object of class glm - the result of fitting the working.msm}
\item{cum.g}{array, n x numACnodes x numRegimes - cumulative g, after bounding}
\item{cum.g.unbounded}{array, n x numACnodes x numRegimes - cumulative g, before bounding}
\item{call}{the matched call}
\item{gcomp}{the \code{gcomp} input}
\item{formulas}{a \code{list} with elements \code{Qform} and \code{gform}}
\item{fit}{a list with the following components}
  \itemize{
  \item \code{g} - list of length numRegimes of list of length numACnodes - \code{glm} or \code{SuperLearner} return objects from fitting g regressions
  \item \code{Q} - list of length numLYnodes - \code{glm} or \code{SuperLearner} return objects from fitting Q regressions
  \item \code{Qstar} - list of length numLYnodes - \code{glm} (or numerical optimization if \code{glm} fails to solve the score equation) return objects from updating the Q fit
  }
}

\references{
Lendle, Samuel, Schwab, Joshua, Petersen, Maya and and van der Laan, Mark J "ltmle: An R Package Implementing Targeted Minimum Loss-based Estimation for Longitudinal Data", Forthcoming  
	
Petersen, Maya, Schwab, Joshua and van der Laan, Mark J, "Targeted Maximum Likelihood Estimation of Marginal Structural Working Models for Dynamic Treatments Time-Dependent Outcomes", Forthcoming

van der Laan, Mark J. and Gruber, Susan, "Targeted Minimum Loss Based Estimation of an Intervention Specific Mean Outcome" (August 2011). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 290.
\url{http://biostats.bepress.com/ucbbiostat/paper290}

van der Laan, Mark J. and Rose, Sherri, "Targeted Learning: Causal Inference for Observational and Experimental Data" New York: Springer, 2011.
}

\author{
Joshua Schwab \email{joshuaschwab@yahoo.com}, Samuel Lendle, Maya Petersen, and Mark van der Laan
}

\seealso{
\code{\link{summary.ltmle}}, \code{\link{summary.ltmleMSM}}, \code{\link[SuperLearner:SuperLearner]{SuperLearner}}, \code{\link{deterministic.g.function_template}}, \code{\link{deterministic.Q.function_template}}
}
\examples{
rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
                 
# Example 1: Single time point example. True value of E[Y_1] (expected value of
#   Y setting A to 1) is approximately 0.5939.
set.seed(2)
n <- 1000
W1 <- rnorm(n)
W2 <- rbinom(n, size=1, prob=0.3)   
W3 <- rnorm(n)
A <- rexpit(-1 + 2 * W1^2)
Y <- rexpit(-0.5 + 2 * W1^2 + 0.5 * W2 - 0.5 * A + 0.2 * W3 * A 
       - 1.1 * W3 + 0.2 * rnorm(n))

data <- data.frame(W1, W2, W3, A, Y)


\donttest{ #This takes about 4 seconds to run
library(SuperLearner)

#SuperLearner semiparametric estimation using all parents as regressors 
result1 <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", abar=1, 
  SL.library=c("SL.glm", "SL.step", "SL.mean"))
summary(result1)
summary(result1, estimator="iptw")

#SuperLearner semiparametric estimation using (incorrectly) specified regressors
#note: The functional form for Qform and gform is unimportant if 
# using SuperLearner - see 'Details'
result1a <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
 Qform=c(Y="Q.kplus1 ~ W1 + W3 + A"), gform="A ~ W1", abar=1, 
 SL.library=c("SL.glm", "SL.step", "SL.mean"))
summary(result1a)
}

#glm using correctly specified Qform and gform
result.abar1 <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
 Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ I(W1^2)", 
 abar=1, SL.library=NULL)

\donttest{ #This takes about 18 seconds to run
#Get summary measures (additive treatment effect, odds ratio, relative risk) 
#  for abar=1 vs abar=0
result.compare <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
                      Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ I(W1^2)", 
                      abar=list(1, 0), SL.library=NULL)
summary(result.compare)


# Example 2: Longitudinal example. Includes informative censoring and treatment. 
# Time ordering of data is W, C1, L1, A1, Y1, C2, L2, A2, Y2
# True value of E[Y_(1,1,1,1)] (expected value of Y setting C1, A1, C2, A2 all to 1)
#  is approximately 0.413.
# A1 is known to always be 1 if L1 < -2, and is 1 with probability 0.1 if L1 > 2 
# A2 is known to always be 1 if A1 is 1 
# We incorporate this knowledge using deterministic.g.function

# Generate data:
ua <- rep(TRUE, n)   #ua = uncensored and alive
L1 <- A1 <- Y1 <- C2.binary <- L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
W <- rnorm(n)
C1 <- BinaryToCensoring(is.uncensored=rexpit(2 + W))
ua <- ua & C1 == "uncensored"
L1[ua] <- rnorm(n)[ua] + W[ua]
A1[ua] <- rexpit(L1[ua])
A1[ua & L1 < -2] <- 1
A1[ua & L1 >  2] <- rbinom(n, size=1, prob=0.1)[ua & L1 >  2]
Y1[ua] <- rexpit((W + L1 - A1)[ua])
ua <- ua & !Y1
C2.binary[ua] <- rexpit((1 + 0.7 * L1 - A1)[ua])
C2 <- BinaryToCensoring(is.uncensored=C2.binary)
ua <- ua & C2 == "uncensored"
L2[ua] <- (0.5 * L1 - 0.9 * A1 + rnorm(n))[ua]
A2[ua] <- rexpit((0.5 * L1 + 0.8 * L2)[ua]) | A1[ua]
Y2[ua] <- rexpit((0.7 * L1 + L2 - 0.8 * A1 - A2)[ua])
Y2[Y1 == 1] <- 1  # if a patient dies at time 1, record death at time 2 as well
data <- data.frame(W, C1, L1, A1, Y1, C2, L2, A2, Y2)

deterministic.g.function <- function(data, current.node, nodes) {
  if (names(data)[current.node] == "A1") {
    det <- (data$L1 < -2 | data$L1 > 2) & !is.na(data$L1)
    prob1 <- ((data$L1 < -2) * 1 + (data$L1 > 2) * 0.1)[det]
  } else if (names(data)[current.node] == "A2") {
    det <- data$A1 == 1 & !is.na(data$A1)
    prob1 <- 1
  } else if (names(data[current.node]) \%in\% c("C1", "C2")){
    return(NULL)  #this returns the default of no deterministic links 
    #note that it is not necessary to specify that prior censoring indicates future censoring
  } else {
    stop("unexpected current.node")
  }
  return(list(is.deterministic=det, prob1=prob1))  
}

result2 <- ltmle(data, Anodes=c("A1","A2"), Cnodes=c("C1", "C2"), 
                Lnodes=c("L1", "L2"), Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
                deterministic.g.function=deterministic.g.function, survivalOutcome=TRUE)
summary(result2) 
 
# Example 3: Dynamic treatment, observation weights
# W -> A1 -> L -> A2 -> Y
# Treatment regime of interest is: Always treat at time 1 (A1 = 1), 
#   treat at at time 2 (A2 = 1), iff L > 0
# Weight by pmax(W + 2, 0)

n <- 1000
W <- rnorm(n)
A1 <- rexpit(W)
L <- 0.3 * W + 0.2 * A1 + rnorm(n)
A2 <- rexpit(W + A1 + L)
Y <- rexpit(W - 0.6 * A1 + L - 0.8 * A2)
data <- data.frame(W, A1, L, A2, Y)

abar <- matrix(nrow=n, ncol=2)
abar[, 1] <- 1
abar[, 2] <- L > 0

result3 <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", 
  survivalOutcome=TRUE, abar=abar, observation.weights = pmax(W + 2, 0))
summary(result3)

# Example 3.1: The regime can also be specified as a rule function

rule <- function(row) c(1, row["L"] > 0)

result.rule <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", 
  survivalOutcome=TRUE, rule=rule)
# This should be the same as the above result
summary(result.rule)

# Example 4: Deterministic Q function
# W -> A1 -> Y1 -> L2 -> A2 -> Y2
n <- 200
L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
W <- rnorm(n)
A1 <- rexpit(W)
Y1 <- rexpit(W - A1)
alive <- Y1 == 0
L2[alive] <- (0.5 * W - 0.9 * A1 + rnorm(n))[alive]
completed.study <- alive & L2 > 0

#Specify that Q is deterministically 0 when L2 is in the history of the 
# current Q regression and L2 > 0
#Note 1: det.Q.fun doesn't condition on called.from.estimate.g so g will also be set 
#        deterministically after L2 > 0 
#Note 2: It is not necessary to specify that Q is deterministically 1 if Y1 is 1; this is automatic
det.Q.fun.4a <- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "L2")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data$L2 > 0 & !is.na(data$L2)
  return(list(is.deterministic=is.deterministic, Q.value=0))
}

#patients don't change treatment after leaving study; leave their A2 as NA
A2[alive & !completed.study] <- rexpit((0.5 * W + 0.8 * L2)[alive & !completed.study])

Y2[alive & !completed.study] <- rexpit((L2 - 0.8 * A1 - A2)[alive & !completed.study])
Y2[alive & completed.study] <- 0
Y2[!alive] <- 1  # if a patient dies at time 1, record death at time 2 as well
data <- data.frame(W, A1, Y1, L2, A2, Y2)

result4a <- ltmle(data, Anodes=c("A1","A2"), Lnodes="L2", Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
  SL.library=NULL, estimate.time=FALSE, deterministic.Q.function=det.Q.fun.4a, survivalOutcome=TRUE,
  IC.variance.only=TRUE)
  #IC.variance.only=FALSE is not currently compatible with deterministic.Q.function
#note: You will get the same result if you pass Lnodes=NULL (see next example)
summary(result4a)

#In this variation, suppose that treatment can still change after a patient leaves the study

det.Q.fun.4b <- function(data, current.node, nodes, called.from.estimate.g) {
  #there is no deterministic information when calculating g - treatment may still change
  if (called.from.estimate.g) return(NULL)  
  
  L2.index <- which(names(data) == "L2")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data$L2 > 0 & !is.na(data$L2)
  return(list(is.deterministic=is.deterministic, Q.value=0))
}

A2[alive] <- rexpit((0.5 * W + 0.8 * L2)[alive])  #patients can change treatment after leaving study
Y2[alive & !completed.study] <- rexpit((L2 - 0.8 * A1 - A2)[alive & !completed.study])
Y2[alive & completed.study] <- 0
Y2[!alive] <- 1  # if a patient dies at time 1, record death at time 2 as well
data <- data.frame(W, A1, Y1, L2, A2, Y2)

result4b <- ltmle(data, Anodes=c("A1","A2"), Lnodes="L2", Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
 SL.library=NULL, estimate.time=FALSE, deterministic.Q.function=det.Q.fun.4b, survivalOutcome=TRUE,
 IC.variance.only=TRUE) 
 #IC.variance.only=FALSE is not currently compatible with deterministic.Q.function
summary(result4b)

# Example 5: Multiple time-dependent covariates and treatments at each time point, 
#            continuous Y values
# age -> gender -> A1 -> L1a -> L1b -> Y1 -> A2 -> L2a -> L2b -> Y2
n <- 100
age <- rbinom(n, 1, 0.5)
gender <- rbinom(n, 1, 0.5)
A1 <- rexpit(age + gender)
L1a <- 2*age - 3*gender + 2*A1 + rnorm(n)
L1b <- rexpit(age + 1.5*gender - A1)
Y1 <- plogis(age - gender + L1a + 0.7*L1b + A1 + rnorm(n))
A2 <- rexpit(age + gender + A1 - L1a - L1b)
L2a <- 2*age - 3*gender + 2*A1 + A2 + rnorm(n)
L2b <- rexpit(age + 1.5*gender - A1 - A2)
Y2 <- plogis(age - gender + L1a + L1b + A1 + 1.8*A2 + rnorm(n))
data <- data.frame(age, gender, A1, L1a, L1b, Y1, A2, L2a, L2b, Y2)

#also show some different ways of specifying the nodes:
result5 <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), Ynodes=
 grep("^Y", names(data)), abar=c(1, 0), SL.library=NULL, estimate.time=FALSE, 
 survivalOutcome=FALSE, IC.variance.only=TRUE) 
 #IC.variance.only=FALSE is not currently compatible with non-binary outcomes
summary(result5)

result5a <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
 Ynodes=grep("^Y", names(data)), abar=c(1, 0), SL.library=NULL, estimate.time=FALSE, 
 survivalOutcome=FALSE, gform=c("A1 ~ gender", "A2 ~ age"), 
 IC.variance.only=TRUE) 
 #IC.variance.only=FALSE is not currently compatible with non-binary outcomes
summary(result5a)

#Usually you would specify a Qform for all of the Lnodes and Ynodes but in this case 
# L1a, L1b, Y1 is a "block" of L/Y nodes not separated by Anodes or Cnodes (the same is true for 
# L2a, L2b, Y2). Only one regression is required at the first L/Y node in a block. You can pass 
# regression formulas for the other L/Y nodes, but they'll be ignored.
result5b <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
 Ynodes=grep("^Y", names(data)), abar=c(1, 0), estimate.time=FALSE, survivalOutcome=FALSE, 
 gform=c("A1 ~ gender", "A2 ~ age"), Qform=c(L1a="Q.kplus1 ~ 1", L2a="Q.kplus1 ~ 1"), 
 IC.variance.only=TRUE) 
 #IC.variance.only=FALSE is not currently compatible with non-binary outcomes
summary(result5b)


#Gives the same result but prints a message saying some regression formulas will be dropped:
result5c <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
 Ynodes=grep("^Y", names(data)), abar=c(1, 0), estimate.time=FALSE, survivalOutcome=FALSE, 
 gform=c("A1 ~ gender", "A2 ~ age"), Qform=c(L1a="Q.kplus1 ~ 1", L1b="Q.klus1~A1", 
 Y1="Q.kplus1~L1a", L2a="Q.kplus1 ~ 1", L2b="Q.klus1~A1", Y2="Q.kplus1~A2 + gender"), 
 IC.variance.only=TRUE) 
 #IC.variance.only=FALSE is not currently compatible with non-binary outcomes

summary(result5c)


#If there were a Anode or Cnode between L1b and Y1, Y1 would also need a Q regression formula


# Example 6: MSM
# Given data over 3 time points where A switches to 1 once and then stays 1. We want to know
# how death varies as a function of gender, time and an indicator of whether a patient's 
# intended regime was to switch before time.
# Note that working.msm includes time and switch.time, which are columns of 
# summary.measures; working.msm also includes male, which is ok because it is a baseline
# covariate (it comes before any A/C/L/Y nodes).
data(sampleDataForLtmleMSM)
Anodes <- grep("^A", names(sampleDataForLtmleMSM$data))
Lnodes <- c("CD4_1", "CD4_2")
Ynodes <- grep("^Y", names(sampleDataForLtmleMSM$data))
msm.weights <- matrix(1:12, nrow=4, ncol=3) #just an example (can also use a 200x3x4 array), 
                                            #or NULL (for no weights), or "empirical" (the default)

result6 <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
                   survivalOutcome=TRUE,
                   regimes=sampleDataForLtmleMSM$regimes, 
                   summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, 
                   working.msm="Y ~ male + time + I(pmax(time - switch.time, 0))", 
                   msm.weights=msm.weights, estimate.time=FALSE)
print(summary(result6))


# Example 6.1: regimes can also be specified as a list of rule functions

regimesList <- list(function(row) c(1,1,1),
                     function(row) c(0,1,1),
                     function(row) c(0,0,1),
                     function(row) c(0,0,0))
result.regList <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
                   survivalOutcome=TRUE, regimes=regimesList, 
                   summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, 
                   working.msm="Y ~ male + time + I(pmax(time - switch.time, 0))", 
                   msm.weights=msm.weights, estimate.time=FALSE)
# This should be the same as the above result
print(summary(result.regList))         


# Example 7: variance estimation
# A simple point treatment problem W, A, Y. But there is a positivity problem - 
# for small values of W, Prob(A = 1) is very small.
# The true parameter value, E[Y_1] is approximately 0.697
# The true TMLE standard deviation is approximately 0.064, 
# the true IPTW standard deviation is approximately 0.058.
set.seed(2)
n <- 1000
W <- rnorm(n)
A <- rexpit(8 * W)
Y <- rexpit(W + A)
r1 <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE)
r2 <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, 
 IC.variance.only=TRUE)
print(summary(r1))
print(summary(r2))
print(summary(r1, "iptw"))
print(summary(r2, "iptw")) #the same - IC.variance.only only affects TMLE
}
}