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
ltmle(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=FALSE, Qform=NULL, gform=NULL, abar, rule=NULL,
 gbounds=c(0.01, 1), deterministic.acnode.map=NULL, stratify=FALSE, SL.library=NULL, 
 estimate.time=nrow(data) > 50, gcomp=FALSE, mhte.iptw=FALSE, iptw.only=FALSE, deterministic.Q.map=NULL)
ltmleMSM(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=FALSE, Qform=NULL, gform=NULL, 
 gbounds=c(0.01, 1), deterministic.acnode.map=NULL, SL.library=NULL, regimens, working.msm, 
 summary.measures, summary.baseline.covariates=NULL, final.Ynodes=NULL, pooledMSM=TRUE, 
 stratify=FALSE, weight.msm=TRUE, estimate.time=nrow(data) > 50, gcomp=FALSE, mhte.iptw=FALSE, 
 iptw.only=FALSE, deterministic.Q.map=NULL, memoize=TRUE)
}
\arguments{
  \item{data}{data frame following the time-ordering of the nodes. See 'Details'.}
  \item{Anodes}{column names or indicies in \code{data} of treatment nodes}
  \item{Cnodes}{column names or indicies in \code{data} of censoring nodes}
  \item{Lnodes}{column names or indicies in \code{data} of time-dependent covariate nodes}
  \item{Ynodes}{column names or indicies in \code{data} of outcome nodes}
  \item{survivalOutcome}{If \code{TRUE}, then Y nodes are indicators of an event, and if Y at some time point is 1, then all following should be 1.}
  \item{Qform}{character vector of regression formulas for \eqn{Q}. See 'Details'.}
  \item{gform}{character vector of regression formulas for \eqn{g}. See 'Details'.}
  \item{abar}{binary vector (numAnodes x 1) or matrix (n x numAnodes) of counterfactual treatment}
  \item{rule}{a function to be applied to each row (a named vector) of \code{data} that returns a numeric vector of length numAnodes}
  \item{gbounds}{lower and upper bounds on estimated probabilities for g-factors. Vector of length 2, order unimportant.}
  \item{deterministic.acnode.map}{optional information on A and C nodes that are given deterministically. See 'Details'. Default \code{NULL} indicates no deterministic links.}
  \item{stratify}{if \code{TRUE} stratify on following \code{abar} when estimating Q and g. If \code{FALSE}, pool over \code{abar}.}
  \item{SL.library}{optional character vector of libraries to pass to \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. \code{NULL} indicates \link{glm} should be called instead of \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. '\code{default}' indicates a standard set of libraries. May be separately specified for \eqn{Q} and \eqn{g}. See 'Details'.}
  \item{estimate.time}{if \code{TRUE}, run an initial estimate using only 50 observations and use this to print a very rough estimate of the total time to completion. No action if there are fewer than 50 observations.}
  \item{gcomp}{if \code{TRUE}, run the maximum likelihood based G-computation estimate \emph{instead} of TMLE}
  \item{regimens}{binary array: n x numAnodes x numRegimens of counterfactual treatment or a list of 'rule' functions}
  \item{working.msm}{character formula for the working marginal structural model}
  \item{summary.measures}{array: num.regimens x num.summary.measures x num.final.Ynodes - measures summarizing the regimens that will be used on the right hand side of working.msm}
  \item{summary.baseline.covariates}{NOT FULLY IMPLEMENTED YET - leave as NULL (default)}
  \item{final.Ynodes}{vector subset of Ynodes - used in MSM to pool over a set of outcome nodes}
  \item{pooledMSM}{if \code{TRUE}, the TMLE targeted step will pool across regimens}
  \item{weight.msm}{if \code{TRUE}, the working.msm will be weighted by the empirical probability of each regimen [in the future more flexible weighting may be possible]} 
  \item{mhte.iptw}{if \code{TRUE}, IPTW is calculated using the modified Horvitz-Thompson estimator (normalizes by sum of the inverse weights)}
  \item{iptw.only}{by default (\code{iptw.only = FALSE}), both TMLE and IPTW are run in \code{ltmle} and \code{ltmleMSM(pooledMSM=FALSE)}. If \code{iptw.only = TRUE}, only IPTW is run, which is faster. This parameter is not used in \code{ltmleMSM(pooledMSM=FALSE)} since IPTW is not run.}
  \item{deterministic.Q.map}{optional information on Q given deterministically. See 'Details'. Default \code{NULL} indicates no deterministic links.}
  \item{memoize}{If \code{TRUE}, glm regressions will be memoized. It is recommended to leave this as \code{TRUE} (the default), especially if there are multiple \code{final.Ynodes}, because the code is not written as efficiently as it should be and will end up repeating the same glm call. Will be fixed in a future release.}
}
\details{
  The estimates returned by \code{ltmle} are of a treatment specific mean, \eqn{E[Y_{\bar{a}}]}, the mean of the final treatment node, where all treatment nodes, \eqn{A}, are set to \eqn{\bar{a}} (\code{abar}) and all censoring nodes \eqn{C} are set to 1 (uncensored). The estimates returned by \code{ltmleMSM} are similar but are the parameters in a working marginal structural model.
  
  By calling \code{ltmle} twice, using two different values of \code{abar}, additive treatment effect, risk ratio, and odds ratio can be computed using \code{\link{summary.ltmle}}. 
  
  \code{data} should be a data frame where the order of the columns corresponds to the time-ordering of the model. 
  \itemize{
      \item in censoring columns (Cnodes): factor with two levels: "censored" and "uncensored". 
      \item in treatment columns (Anodes): 1 = treated, 0 = untreated (must be binary)
      \item in event columns (Ynodes): binary where 1 = event (e.g. death) 0 = no event OR continuous in [0, 1). Note that Y may not be exactly equal to 1 in the case of continuous Y. Genral continous Y may be scaled to [0, 0.999999] and then unscaled.
      \item time-dependent covariate columns (Lnodes): can be any numeric data
  }
  
  There may be columns in \code{data} that are not in any of \code{Cnodes}, \code{Anodes}, \code{Ynodes}, and \code{Lnodes}.
  
  The events in Ynodes must be of the form where once Y jumps to 1, Y remains 1 at subsequent nodes. Future releases of this package may relax this assumption. 
  
  Data in Cnodes, Anodes, and Lnodes is not used after an event (e.g. death) or censoring and may be coded as NA or any other value.
  
  A node should only be classified as an Lnode if it comes after a Cnode or Anode. If it does not, it is considered a a baseline covariate, not a time-dependent covariate. 
  
  \code{Qform} should be \code{NULL}, in which case all parent nodes of each L and Y node will be used as regressors, or a named character vector that can be coerced to class "\code{\link{formula}}". The length of \code{Qform} must be equal to \code{length(Lnodes) + length(Ynodes)}** and the names and order of the formulas must be the same as the names and order of the L and Y nodes in \code{data}. The left hand side of each formula should be "\code{Q.kplus1}". If \code{SL.library} is \code{NULL}, \code{glm} will be called using the elements of \code{Qform}. If \code{SL.library} is specified, \code{\link[SuperLearner:SuperLearner]{SuperLearner}} will be called and all variables appearing on the right hand side of a formula will be passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. The actual functional form of the formula is unimportant if \code{\link[SuperLearner:SuperLearner]{SuperLearner}} is called.
  
  ** If there is a "block" of L and Y nodes not separated by A or C nodes, only one regression is required at the first L/Y node in a block. You can pass regression formulas for the other L/Y nodes, but they will be ignored (with a message). See example 5.
  
  \code{gform} should be \code{NULL}, in which case all parent nodes of each L and Y node will be used as regressors, or a character vector that can be coerced to class "\code{\link{formula}}". The length of \code{gform} must be equal to \code{length(Anodes) + length(Cnodes)} and the order of the formulas must be the same as the order the A and C nodes appear in \code{data}. The left hand side of each formula should be the name of the Anode or Cnode. If \code{SL.library} is \code{NULL}, \code{glm} will be called using the elements of \code{gform}. If \code{SL.library} is specified, \code{\link[SuperLearner:SuperLearner]{SuperLearner}} will be called and all variables appearing on the right hand side of a formula will be passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. The actual functional form of the formula is unimportant if \code{\link[SuperLearner:SuperLearner]{SuperLearner}} is called.
  
  \code{abar} specifies the counterfactual values of the Anodes, using the order they appear in \code{data} and should have the same length (if abar is a vector) or number of columns (if abar is a matrix) as \code{Anodes}.

  \code{rule} can be used to specify a dynamic treatment rule. \code{rule} is a function applied to each row of \code{data} which returns the a numeric vector of the same length as \code{Anodes}.

  \code{regimens} can be a binary array: n x numAnodes x numRegimens of counterfactual treatment or a list of 'rule' functions as described above for the \code{rule} parameter for the \code{ltmle} function
  
  \code{deterministic.acnode.map} can be used to specify model knowledge about value of Anodes and/or Cnodes that are set deterministically. For example, it may be the case that once a patient starts treatment, they always stay on treatment. \code{deterministic.acnode.map} is a list [1 ... number of deterministic nodes] of lists with the following components:
    \itemize{
      \item \code{node} - the name or index of the node \emph{that we wish to set the probability of treatment/censoring for}
      \item \code{is.deterministic} - logical vector, length equal to \code{nrow(data)}, \code{TRUE} if a patient's status for this node is deterministic
      \item \code{prob1} - probability that node is 1 (treated/uncensored), length equal to 1 or \code{length(which(is.deterministic))}
}

  \code{deterministic.Q.map} can be used to specify model knowledge about the final event state. For example, it may be the case that a patient can complete the study at some intermediate time point, in which case the probability of death is 0 (assuming they have not died already). \code{deterministic.Q.map} is a list [1 ... number of deterministic nodes] of lists with the following components:
    \itemize{
      \item \code{node} - the name or index of the node. \emph{If this node is in the history of the L or Y node at the current Q regression, Q will be set} 
      \item \code{is.deterministic} - logical vector, length equal to \code{nrow(data)}, \code{TRUE} if a patient's status for this node is deterministic
      \item \code{Q.value} - value to set Q, length equal to 1 or \code{length(which(is.deterministic))}
      \item \code{implies.deterministic.g} - logical length 1 - if \code{TRUE}, after a deterministic Q event, subsequent values of treatment and censoring nodes are known to be fixed (e.g. in the case of completing a study) 
}

  \emph{Important! "node" means different things in \code{deterministic.acnode.map} and \code{deterministic.Q.map}.} In \code{deterministic.acnode.map}, \code{node="A1"} means that for patients where \code{is.deterministic = TRUE}, the probability (g) that A1 == 1 is set according to \code{prob1}. In \code{deterministic.Q.map}, \code{node="A1"} means that in Q regressions at all L or Y nodes after A1, for patients where \code{is.deterministic = TRUE}, Q (iterated expectation of final Y) is set according to \code{Q.value}.

\code{SL.library} may be a character vector of libraries (or \code{NULL} or '\code{default}'), in which case these libraries are     used to estimate both \eqn{Q} and \eqn{g} OR a list with two components, \code{Q} and \code{g}, where each is a character vector of libraries (or \code{NULL} or '\code{default}'). 
\code{NULL} indicates \link{glm} should be called instead of \code{\link[SuperLearner:SuperLearner]{SuperLearner}}
If \code{SL.library} is the string '\code{default}', \code{SL.library} is set to \code{list("SL.glm", "SL.glmnet", "SL.stepAIC", "SL.bayesglm", c("SL.glm", "screen.corP"), c("SL.glmnet", "screen.corP"), c("SL.step", "screen.corP"), c("SL.step.forward", "screen.corP"), c("SL.stepAIC", "screen.corP"), c("SL.step.interaction", "screen.corP"), c("SL.bayesglm", "screen.corP"))}. Note that the default set of libraries consists of main terms models. It may be advisable to include squared terms, interaction terms, etc in \code{data} or include libraries that consider non-linear terms.

The print method for \code{ltmle} objects only prints the tmle estimates. 
}
\value{
\code{ltmle} returns an object of class "\code{ltmle}"
The function \code{\link{summary}} (i.e. \code{\link{summary.ltmle}}) can be used to obtain or print a summary of the results.
An object of class "\code{ltmle}" is a list containing the following components:
\item{estimates}{a named vector of length 4 with elements, each an estimate of \eqn{E[Y_{bar{a}}]}:
  \itemize{
  \item \code{tmle} - Targeted Maximum Likelihood Estimate [NULL if \code{gcomp} is \code{TRUE}]
  \item \code{iptw} - Inverse Probability of Treatment/Censoring Weighted estimate 
  \item \code{gcomp} - maximum likelihood based G-computation estimate [NULL if \code{gcomp} is \code{FALSE}]
  \item \code{naive} - naive estimate \eqn{E[Y|a=\bar{a}]}
  }
}
\item{IC}{a list with the following components of Influence Curve values}
  \itemize{
  \item \code{tmle} - vector of influence curve values for Targeted Maximum Likelihood Estimate [NULL if \code{gcomp} is \code{TRUE}]
  \item \code{iptw} - vector of influence curve values for Inverse Probability of Treatment/Censoring Weighted estimate 
  \item \code{gcomp} - vector of influence curve values for Targeted Maximum Likelihood Estimate without updating [NULL if \code{gcomp} is \code{FALSE}]
  }
\item{call}{the matched call}
\item{gcomp}{the \code{gcomp} input}
\item{formulas}{a \code{list} with elements \code{Qform} and \code{gform}}

\code{ltmleMSM} returns an object of class "\code{ltmleMSM}"
The function \code{\link{summary}} (i.e. \code{\link{summary.ltmleMSM}}) can be used to obtain or print a summary of the results.
An object of class "\code{ltmleMSM}" is a list containing the following components:
\item{beta}{parameter estimates for working.msm using TMLE (GCOMP if \code{gcomp} input is \code{TRUE})}
\item{beta.iptw}{parameter estimates for working.msm using IPTW (\code{NULL} if \code{pooledMSM} is \code{TRUE})}
\item{IC}{matrix, n x numBetas - influence curve values for TMLE (without updating if \code{gcomp} input is \code{TRUE})}
\item{IC.iptw}{matrix, n x numBetas - influence curve values for IPTW (\code{NULL} if \code{pooledMSM} is \code{TRUE})}
\item{msm}{object of class glm - the result of fitting the working.msm}
\item{cum.g}{matrix, n x numACnodes - cumulative g, after bounding}
\item{call}{the matched call}
\item{gcomp}{the \code{gcomp} input}
\item{formulas}{a \code{list} with elements \code{Qform} and \code{gform}}
}

\references{
van der Laan, Mark J. and Gruber, Susan, "Targeted Minimum Loss Based Estimation of an Intervention Specific Mean Outcome" (August 2011). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 290.
\url{http://biostats.bepress.com/ucbbiostat/paper290}

van der Laan, Mark J. and Rose, Sherri, "Targeted Learning: Causal Inference for Observational and Experimental Data" New York: Springer, 2011.

Maya Petersen, Joshua Schwab, Susan Gruber, Nello Blaser, Michael Schomaker, and Mark van der Laan. "Targeted Maximum Likelihood Estimation for Longitudinal Marginal Structural Working Models". Forthcoming.
}
\author{
Joshua Schwab \email{joshuaschwab@yahoo.com}, Maya Petersen, and Mark van der Laan
}

\seealso{
\code{\link{summary.ltmle}}, \code{\link{summary.ltmleMSM}}, \code{\link[SuperLearner:SuperLearner]{SuperLearner}}
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


\dontrun{
# The SuperLearner examples are a little slow.
library(SuperLearner)

#SuperLearner semiparametric estimation using all parents as regressors 
result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", abar=1, 
  SL.library=c("SL.glm", "SL.step", "SL.mean"))
summary(result)
summary(result, estimator="iptw")

#SuperLearner semiparametric estimation using (incorrectly) specified regressors
#note: The functional form for Qform and gform is unimportant if 
# using SuperLearner - see 'Details'
result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
 Qform=c(Y="Q.kplus1 ~ W1 + W3 + A"), gform="A ~ W1", abar=1, 
 SL.library=c("SL.glm", "SL.step", "SL.mean"))
summary(result)
}

#glm using correctly specified Qform and gform
result.abar1 <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
 Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ I(W1^2)", 
 abar=1, SL.library=NULL)

#Get summary measures (additive treatment effect, odds ratio, relative risk) 
#  for abar=1 vs abar=0
result.abar0 <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
 Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ I(W1^2)", 
 abar=0, SL.library=NULL)
summary(result.abar1, result.abar0)


# Example 2: Longitudinal example. Includes informative censoring and treatment. 
# Time ordering of data is W, C1, L1, A1, Y1, C2, L2, A2, Y2
# True value of E[Y_(1,1,1,1)] (expected value of Y setting C1, A1, C2, A2 all to 1)
#  is approximately 0.413.
# A1 is known to always be 1 if L1 < -2, and is 1 with probability 0.1 if L1 > 2 
# A2 is known to always be 1 if A1 is 1 
# We incorporate this knowledge using deterministic.acnode.map

# Helper function to generate censoring as factors
BinaryToCensoring <- function(x) {
  if (! all(x \%in\% c(0, 1, NA))) stop("x should be binary")
  y <- character(length(x))
  y[x == 0] <- "censored"
  y[x == 1] <- "uncensored"
  y[is.na(x)] <- NA
  return(factor(y))
}

# Generate data:
ua <- rep(TRUE, n)   #ua = uncensored and alive
L1 <- A1 <- Y1 <- C2.binary <- L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
W <- rnorm(n)
C1 <- BinaryToCensoring(rexpit(2 + W))
ua <- ua & C1 == "uncensored"
L1[ua] <- rnorm(n)[ua] + W[ua]
A1[ua] <- rexpit(L1[ua])
A1[ua & L1 < -2] <- 1
A1[ua & L1 >  2] <- rbinom(n, size=1, prob=0.1)[ua & L1 >  2]
Y1[ua] <- rexpit((W + L1 - A1)[ua])
ua <- ua & !Y1
C2.binary[ua] <- rexpit((1 + 0.7 * L1 - A1)[ua])
C2 <- BinaryToCensoring(C2.binary)
ua <- ua & C2 == "uncensored"
L2[ua] <- (0.5 * L1 - 0.9 * A1 + rnorm(n))[ua]
A2[ua] <- rexpit((0.5 * L1 + 0.8 * L2)[ua]) | A1[ua]
Y2[ua] <- rexpit((0.7 * L1 + L2 - 0.8 * A1 - A2)[ua])
Y2[Y1 == 1] <- 1  # if a patient dies at time 1, record death at time 2 as well
data <- data.frame(W, C1, L1, A1, Y1, C2, L2, A2, Y2)

# Analyze data:
deterministic.acnode.map <- function(data, current.ac.node, Anodes, Cnodes) {
  if (names(data[current.ac.node])=="A1") {
    det <- (data$L1 < -2 | data$L1 > 2) & !is.na(L1)
    prob1 <- ((data$L1 < -2) * 1 + (data$L1 > 2) * 0.1)[det]
  } else if (names(data[current.ac.node])=="A2") {
    det <- data$A1 == 1 & !is.na(data$A1)
    prob1 <- 1
  } else if (names(data[current.ac.node]) \%in\% c("C1", "C2")){
    return(NULL)
  } else {
    stop("unexpected current.ac.node")
  }
  return(list(is.deterministic=det, prob1=prob1))  
}

result <- ltmle(data, Anodes=c("A1","A2"), Cnodes=c("C1", "C2"), 
                Lnodes=c("L1", "L2"), Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, abar=c(1, 1), 
                deterministic.acnode.map=deterministic.acnode.map, SL.library=NULL)
summary(result) 
 
# Example 3: Dynamic treatment
# W -> A1 -> L -> A2 -> Y
# Treatment regimen of interest is: Always treat at time 1 (A1 = 1), treat at at time 2 (A2 = 1), iff L > 0
# True value of E[Y_d] is approximately 0.346

n <- 10000
W <- rnorm(n)
A1 <- rexpit(W)
L <- 0.3 * W + 0.2 * A1 + rnorm(n)
A2 <- rexpit(W + A1 + L)
Y <- rexpit(W - 0.6 * A1 + L - 0.8 * A2)
data <- data.frame(W, A1, L, A2, Y)

abar <- matrix(nrow=n, ncol=2)
abar[, 1] <- 1
abar[, 2] <- L > 0

result <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", survivalOutcome=TRUE, abar=abar, SL.library=NULL)
summary(result)

# Example 3.1: The regimen can also be specified as a rule function

rule <- function(row) c(1, row["L"] > 0)

result.rule <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", survivalOutcome=TRUE, rule=rule, SL.library=NULL)
# This should be the same as the above result
summary(result.rule)

# Example 4: Deterministic Q map
# W -> A1 -> Y1 -> L2 -> A2 -> Y2
n <- 1000
L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
deterministic.Q.map <- NULL
W <- rnorm(n)
A1 <- rexpit(W)
Y1 <- rexpit(W - A1)
alive <- Y1 == 0
L2[alive] <- (0.5 * W - 0.9 * A1 + rnorm(n))[alive]
completed.study <- alive & L2 > 0

#specify that Q is deterministically 0 when L2 is in the history of the current Q regression and Y1==0 and L2 > 0
#it is not necessary to specify that Q is deterministically 1 if Y1 is 1 - this is automatic
deterministic.Q.map[[1]] <- list(node="L2", #node can also be the column index of the node in data
                             is.deterministic=completed.study,
                             Q.value=0, #Q.value can also be a vector of length which(is.deterministic)
                             implies.deterministic.g=TRUE)
A2[alive & !completed.study] <- rexpit((0.5 * W + 0.8 * L2)[alive & !completed.study])
A2[alive & completed.study] <- A1[alive & completed.study] #this doesn't actually change the results because the package assumes that patients follow abar after a deterministic event
Y2[alive & !completed.study] <- rexpit((L2 - 0.8 * A1 - A2)[alive & !completed.study])
Y2[alive & completed.study] <- 0
Y2[!alive] <- 1  # if a patient dies at time 1, record death at time 2 as well
data <- data.frame(W, A1, Y1, L2, A2, Y2)

#result <- ltmle(data, Anodes=c("A1","A2"), Lnodes="L2", Ynodes=c("Y1", "Y2"), survivalOutcome=TRUE, abar=c(1, 1), SL.library=NULL, estimate.time=FALSE, deterministic.Q.map=deterministic.Q.map) #needs fixing!
#note: You will get the same result if you pass Lnodes=NULL (see next example)
summary(result)


# Example 5: Multiple time-dependent covariates and treatments at each time point, continuous Y values
# age -> gender -> A1 -> L1a -> L1b -> Y1 -> A2 -> L2a -> L2b -> Y2
n <- 1000
Y1 <- A2 <- L2a <- L2b <- C2 <- Y2 <- as.numeric(rep(NA, n))
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
result <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), Ynodes=grep("^Y", names(data)), abar=c(1, 0), SL.library=NULL, estimate.time=FALSE)
summary(result)

result <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), Ynodes=grep("^Y", names(data)), abar=c(1, 0), SL.library=NULL, estimate.time=FALSE, gform=c("A1 ~ gender", "A2 ~ age"))
summary(result)

#Usually you would specify a Qform for all of the Lnodes and Ynodes but in this case L1a, L1b, Y1 is a "block" of L/Y nodes not separated by Anodes or Cnodes (the same is true for L2a, L2b, Y2). Only one regression is required at the first L/Y node in a block. You can pass regression formulas for the other L/Y nodes, but they'll be ignored.
result1 <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), Ynodes=grep("^Y", names(data)), abar=c(1, 0), SL.library=NULL, estimate.time=FALSE, gform=c("A1 ~ gender", "A2 ~ age"), Qform=c(L1a="Q.kplus1 ~ 1", L2a="Q.kplus1 ~ 1"))
summary(result1)


\dontrun{
#Gives the same result but prints a message saying some regression formulas will be dropped:
result2 <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), Ynodes=grep("^Y", names(data)), abar=c(1, 0), SL.library=NULL, estimate.time=FALSE, gform=c("A1 ~ gender", "A2 ~ age"), Qform=c(L1a="Q.kplus1 ~ 1", L1b="Q.klus1~A1", Y1="Q.kplus1~L1a", L2a="Q.kplus1 ~ 1", L2b="Q.klus1~A1", Y2="Q.kplus1~A2 + gender"))
summary(result2)
}


#If there were a Anode or Cnode between L1b and Y1, Y1 would also need a Q regression formula



# Example 6: MSM
# Given data over 3 time points where A switches to 1 once and then stays 1. We want to know
# how death varies as a function of time and an indicator of whether a patient's intended
# regimen was to switch before time.
data(sampleDataForLtmleMSM)
Anodes <- grep("^A", names(sampleDataForLtmleMSM$data))
Lnodes <- c("CD4_2", "CD4_3")
Ynodes <- grep("^Y", names(sampleDataForLtmleMSM$data))

result <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE,
                   regimens=sampleDataForLtmleMSM$regimens, 
                   summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, 
                   working.msm="Y ~ time + I(switch.time <= time)", estimate.time=FALSE)
print(summary(result))


# Example 6.1: regimens can also be specified as a list of rule functions

regimensList <- list(function(row) c(1,1,1),
                     function(row) c(0,1,1),
                     function(row) c(0,0,1),
                     function(row) c(0,0,0))
result.regList <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=TRUE,
                   regimens=regimensList, 
                   summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, 
                   working.msm="Y ~ time + I(switch.time <= time)", estimate.time=FALSE)
# This should be the same as the above result
print(summary(result.regList))         

}
