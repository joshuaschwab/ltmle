#' Longitudinal Targeted Maximum Likelihood Estimation
#' 
#' \code{ltmle} is Targeted Maximum Likelihood Estimation (TMLE) of
#' treatment/censoring specific mean outcome for point-treatment and
#' longitudinal data. \code{ltmleMSM} adds Marginal Structural Models. Both
#' always provide Inverse Probability of Treatment/Censoring Weighted estimate
#' (IPTW) as well. Maximum likelihood based G-computation estimate (G-comp) can
#' be obtained instead of TMLE. \code{ltmle} can be used to calculate additive
#' treatment effect, risk ratio, and odds ratio.
#' 
#' The estimates returned by \code{ltmle} are of a treatment specific mean,
#' \eqn{E[Y_{\bar{a}}]}, the mean of the final treatment node, where all
#' treatment nodes, \eqn{A}, are set to \eqn{\bar{a}} (\code{abar}) and all
#' censoring nodes \eqn{C} are set to 1 (uncensored). The estimates returned by
#' \code{ltmleMSM} are similar but are the parameters in a working marginal
#' structural model.
#' 
#' \code{data} should be a data frame where the order of the columns
#' corresponds to the time-ordering of the model.  \itemize{ \item in censoring
#' columns (Cnodes): factor with two levels: "censored" and "uncensored". The
#' helper function \code{BinaryToCensoring} can be used to create these
#' factors.  \item in treatment columns (Anodes): 1 = treated, 0 = untreated
#' (must be binary) \item in event columns (Ynodes): If \code{survivalOutcome}
#' is \code{TRUE}, then Y nodes are treated as indicators of a one-time event.
#' See details for \code{survivalOutocme}. If \code{survivalOutcome} is
#' \code{FALSE}, Y nodes are treated as binary if all values are 0 or 1, and
#' are treated as continuous otherwise. If Y nodes are continuous, they may be
#' automatically scaled. See details for \code{Yrange}.  \item time-dependent
#' covariate columns (Lnodes): can be any numeric data \item Data in
#' \code{Cnodes}, \code{Anodes}, \code{Lnodes} and \code{Ynodes} are not used
#' after (to the right of) censoring (or an event when
#' \code{survivalOutcome==TRUE}) and may be coded as \code{NA} or any other
#' value.  \item Columns in \code{data} that are before (to the left of) the
#' first of \code{Cnodes} or \code{Anodes} are treated as baseline variables,
#' even if they are specified as \code{Lnodes}.  \item After the first of
#' \code{Cnodes}, \code{Anodes}, \code{Ynodes}, or \code{Lnodes}, every column
#' must be in one of \code{Cnodes}, \code{Anodes}, \code{Ynodes}, or
#' \code{Lnodes}.  }
#' 
#' If \code{survivalOutcome} is \code{TRUE}, all Y values are indicators of an
#' event (e.g. death) at or before the current time, where 1 = event and 0 = no
#' event. The events in Ynodes must be of the form where once Y jumps to 1, Y
#' remains 1 at subsequent nodes.
#' 
#' For continuous outcomes, (\code{survivalOutcome==FALSE} and some Y nodes are
#' not 0 or 1,) Y values are truncated at the minimum and maximum of
#' \code{Yrange} if specified, and then transformed and scaled to be in [0,1].
#' That is, transformed to \code{(Y-min(Yrange))/(max(Yrange)-min(Yrange))}. If
#' \code{Yrange} is \code{NULL}, it is set to the range of all Y nodes. In that
#' case, Y nodes are only scaled if any values fall outside of [0,1]. For
#' intervention specific means (\code{ltmle}), parameter estimates are
#' transformed back based \code{Yrange}.
#' 
#' \code{Qform} should be \code{NULL}, in which case all parent nodes of each L
#' and Y node will be used as regressors, or a named character vector that can
#' be coerced to class "\code{\link{formula}}". The length of \code{Qform} must
#' be equal to \code{length(Lnodes) + length(Ynodes)}** and the names and order
#' of the formulas must be the same as the names and order of the L and Y nodes
#' in \code{data}. The left hand side of each formula should be
#' "\code{Q.kplus1}". If \code{SL.library} is \code{NULL}, \code{glm} will be
#' called using the elements of \code{Qform}. If \code{SL.library} is
#' specified, \code{\link[SuperLearner:SuperLearner]{SuperLearner}} will be
#' called after a design matrix is created using \code{Qform.} 
#' 
#' ** If there is a "block" of L and Y nodes not separated by A or C nodes,
#' only one regression is required at the first L/Y node in a block. You can
#' pass regression formulas for the other L/Y nodes, but they will be ignored
#' (with a message). See example 5.
#' 
#' \code{gform} should be \code{NULL}, in which case all parent nodes of each L
#' and Y node will be used as regressors, or a character vector that can be
#' coerced to class "\code{\link{formula}}", or a matrix/array of Prob(A=1). If
#' \code{gform} is a character vector, the length of \code{gform} must be equal
#' to \code{length(Anodes) + length(Cnodes)} and the order of the formulas must
#' be the same as the order the A and C nodes appear in \code{data}. The left
#' hand side of each formula should be the name of the Anode or Cnode. If
#' \code{SL.library} is \code{NULL}, \code{glm} will be called using the
#' elements of \code{gform}. If \code{SL.library} is specified,
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}} will be called after a 
#' design matrix is created using \code{gform}. 
#' 
#' In \code{ltmle}, \code{gform} can also be a n x numACnodes matrix where
#' entry (i, j) is the probability that the ith observation of the jth A/C node
#' is 1 (if an Anode) or uncensored (if a Cnode), conditional on following abar
#' up to that node. In \code{ltmleMSM}, \code{gform} can similarly be a n x
#' numACnodes x numRegimes array, where entry (i, j, k) is the probability that
#' the ith observation of the jth A/C node is 1 (if an Anode) or uncensored (if
#' a Cnode), conditional on following regime k up to that node. If \code{gform}
#' is a matrix/array, \code{deterministic.g.function} will not be used and
#' should be \code{NULL}.
#' 
#' \code{abar} specifies the counterfactual values of the Anodes, using the
#' order they appear in \code{data} and should have the same length (if abar is
#' a vector) or number of columns (if abar is a matrix) as \code{Anodes}.
#' 
#' \code{rule} can be used to specify a dynamic treatment rule. \code{rule} is
#' a function applied to each row of \code{data} which returns the a numeric
#' vector of the same length as \code{Anodes}.
#' 
#' \code{abar} and \code{rule} cannot both be specified. If one of them if a
#' list of length 2, additive treatment effect, risk ratio, and odds ratio can
#' be computed using \code{\link{summary.ltmleEffectMeasures}}.
#' 
#' \code{regimes} can be a binary array: n x numAnodes x numRegimes of
#' counterfactual treatment or a list of 'rule' functions as described above
#' for the \code{rule} parameter for the \code{ltmle} function
#' 
#' \code{deterministic.g.function} can be a function used to specify model
#' knowledge about value of Anodes and/or Cnodes that are set
#' deterministically. For example, it may be the case that once a patient
#' starts treatment, they always stay on treatment. For details on the form of
#' the function and examples, see
#' \code{\link{deterministic.g.function_template}}
#' 
#' \code{deterministic.Q.function} can be a function used to specify model
#' knowledge about the final event state. For example, it may be the case that
#' a patient can complete the study at some intermediate time point, in which
#' case the probability of death is 0 (assuming they have not died already).
#' For details on the form of the function and examples, see
#' \code{\link{deterministic.Q.function_template}}
#' 
#' \code{SL.library} may be a character vector of libraries (or '\code{glm}' or
#' '\code{default}'), in which case these libraries are used to estimate both
#' \eqn{Q} and \eqn{g} OR a list with two components, \code{Q} and \code{g},
#' where each is a character vector of libraries (or '\code{glm}' or
#' '\code{default}').  '\code{glm}' indicates \link{glm} should be called
#' instead of \code{\link[SuperLearner:SuperLearner]{SuperLearner}} If
#' \code{SL.library} is the string '\code{default}', \code{SL.library} is set
#' to \code{list("SL.glm", "SL.stepAIC", "SL.bayesglm", c("SL.glm",
#' "screen.corP"), c("SL.step", "screen.corP"), c("SL.step.forward",
#' "screen.corP"), c("SL.stepAIC", "screen.corP"), c("SL.step.interaction",
#' "screen.corP"), c("SL.bayesglm", "screen.corP")}.  Note that the default set
#' of libraries consists of main terms models. It may be advisable to include
#' squared terms, interaction terms, etc in \code{gform} and \code{Qform} or 
#' include libraries that consider non-linear terms.
#' 
#' If \code{attr(SL.library, "return.fit") == TRUE}, then \code{fit$g} and 
#' \code{fit$Q} will return full \code{SuperLearner} or \code{speedglm} objects. 
#' If not, only a summary matrix will be returned to save memory.
#' 
#' The print method for \code{ltmle} objects only prints the tmle estimates.
#' 
#' @aliases ltmle ltmleMSM
#' @param data data frame following the time-ordering of the nodes. See
#' 'Details'.
#' @param Anodes column names or indicies in \code{data} of treatment nodes
#' @param Cnodes column names or indicies in \code{data} of censoring nodes
#' @param Lnodes column names or indicies in \code{data} of time-dependent
#' covariate nodes
#' @param Ynodes column names or indicies in \code{data} of outcome nodes
#' @param survivalOutcome If \code{TRUE}, then Y nodes are indicators of an
#' event, and if Y at some time point is 1, then all following should be 1.
#' Required to be \code{TRUE} or \code{FALSE} if outcomes are binary and there
#' are multiple Ynodes.
#' @param Qform character vector of regression formulas for \eqn{Q}. See
#' 'Details'.
#' @param gform character vector of regression formulas for \eqn{g} or a
#' matrix/array of prob(A=1). See 'Details'.
#' @param abar binary vector (numAnodes x 1) or matrix (n x numAnodes) of
#' counterfactual treatment or a list of length 2. See 'Details'.
#' @param rule a function to be applied to each row (a named vector) of
#' \code{data} that returns a numeric vector of length numAnodes. See 'Details'.
#' @param gbounds lower and upper bounds on estimated cumulative probabilities
#' for g-factors. Vector of length 2, order unimportant.
#' @param Yrange NULL or a numerical vector where the min and max of
#' \code{Yrange} specify the range of all Y nodes. See 'Details'.
#' @param deterministic.g.function optional information on A and C nodes that
#' are given deterministically. See 'Details'. Default \code{NULL} indicates no
#' deterministic links.
#' @param stratify if \code{TRUE} stratify on following \code{abar} when
#' estimating Q and g. If \code{FALSE}, pool over \code{abar}.
#' @param SL.library optional character vector of libraries to pass to
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. \code{NULL} indicates
#' \link{glm} should be called instead of
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}}. '\code{default}'
#' indicates a standard set of libraries. May be separately specified for
#' \eqn{Q} and \eqn{g}. See 'Details'.
#' @param estimate.time if \code{TRUE}, run an initial estimate using only 50
#' observations and use this to print a very rough estimate of the total time
#' to completion. No action if there are fewer than 50 observations.
#' @param gcomp if \code{TRUE}, run the maximum likelihood based G-computation
#' estimate \emph{instead} of TMLE
#' @param regimes binary array: n x numAnodes x numRegimes of counterfactual
#' treatment or a list of 'rule' functions
#' @param working.msm character formula for the working marginal structural
#' model
#' @param summary.measures array: num.regimes x num.summary.measures x
#' num.final.Ynodes - measures summarizing the regimes that will be used on the
#' right hand side of \code{working.msm} (baseline covariates may also be used
#' in the right hand side of \code{working.msm} and do not need to be included
#' in \code{summary.measures})
#' @param final.Ynodes vector subset of Ynodes - used in MSM to pool over a set
#' of outcome nodes
#' @param msm.weights projection weights for the working MSM. If "empirical",
#' weight by empirical proportions of rows matching each regime for each
#' final.Ynode, with duplicate regimes given zero weight. If \code{NULL}, no
#' weights. Or an array of user-supplied weights with dimensions c(n,
#' num.regimes, num.final.Ynodes) or c(num.regimes, num.final.Ynodes).
#' @param iptw.only by default (\code{iptw.only = FALSE}), both TMLE and IPTW
#' are run in \code{ltmle} and \code{ltmleMSM}. If \code{iptw.only = TRUE},
#' only IPTW is run, which is faster.
#' @param deterministic.Q.function optional information on Q given
#' deterministically. See 'Details'. Default \code{NULL} indicates no
#' deterministic links.
#' @param observation.weights observation (sampling) weights. Vector of length
#' n. If \code{NULL}, assumed to be all 1.
#' @param variance.method Method for estimating variance of TMLE. 
#' One of "ic", "tmle", "iptw". If "tmle", compute both the robust variance
#' estimate using TMLE and the influence curve based variance estimate (use the
#' larger of the two). If "iptw", compute both the robust variance
#' estimate using IPTW and the influence curve based variance estimate (use the
#' larger of the two). If "ic", only compute the influence curve based
#' variance estimate. "ic" is fastest, but may be substantially
#' anti-conservative if there are positivity violations or rare outcomes. "tmle" is
#' slowest but most robust if there are positivity violations or rare outcomes. 
#' "iptw" is a compromise between speed and robustness.
#' variance.method="tmle" or "iptw" are not yet available with non-binary outcomes, 
#' gcomp=TRUE, stratify=TRUE, or deterministic.Q.function.
#' @param id Household or subject identifiers. Vector of length n or \code{NULL}. 
#' Integer, factor, or character recommended, but any type that can be coerced 
#' to factor will work. \code{NULL} means all distinct ids.
#' 
#' @return \code{ltmle} returns an object of class "\code{ltmle}" (unless
#' \code{abar} or \code{rule} is a list, in which case it returns an object of
#' class \code{ltmleSummaryMeasures}, which has the same components as
#' \code{ltmleMSM}.) The function \code{\link{summary}} (i.e.
#' \code{\link{summary.ltmle}}) can be used to obtain or print a summary of the
#' results. An object of class "\code{ltmle}" is a list containing the
#' following components: \item{estimates}{a named vector of length 4 with
#' elements, each an estimate of \eqn{E[Y_{bar{a}}]}: \itemize{ \item
#' \code{tmle} - Targeted Maximum Likelihood Estimate [NULL if \code{gcomp} is
#' \code{TRUE}] \item \code{iptw} - Inverse Probability of Treatment/Censoring
#' Weighted estimate \item \code{gcomp} - maximum likelihood based
#' G-computation estimate [NULL if \code{gcomp} is \code{FALSE}] } }
#' \item{IC}{a list with the following components of Influence Curve values}
#' \itemize{ \item \code{tmle} - vector of influence curve values for Targeted
#' Maximum Likelihood Estimate [NULL if \code{gcomp} is \code{TRUE}] \item
#' \code{iptw} - vector of influence curve values for Inverse Probability of
#' Treatment/Censoring Weighted estimate \item \code{gcomp} - vector of
#' influence curve values for Targeted Maximum Likelihood Estimate without
#' updating [NULL if \code{gcomp} is \code{FALSE}] } \item{cum.g}{cumulative g,
#' after bounding: for ltmle, n x numACnodes, for ltmleMSM, n x numACnodes x
#' num.regimes} \item{cum.g.unbounded}{cumulative g, before bounding: for
#' ltmle, n x numACnodes, for ltmleMSM, n x numACnodes x num.regimes} 
#' \item{cum.g.used}{binary - TRUE if an entry of cum.g was used in the updating
#' step (note: even if cum.g.used is FALSE, a small value of cum.g.unbounded may
#' still indicate a positivity problem): for ltmle, n x numACnodes, 
#' for ltmleMSM, n x numACnodes x num.regimes} 
#' \item{call}{the matched call} \item{gcomp}{the \code{gcomp} input}
#' \item{formulas}{a \code{list} with elements \code{Qform} and \code{gform}}
#' \item{fit}{a list with the following components} \itemize{ \item \code{g} -
#' list of length numACnodes - \code{glm} or \code{SuperLearner} (see Details) 
#' return objects from fitting g regressions 
#' \item \code{Q} - list of length numLYnodes - \code{glm} or \code{SuperLearner} 
#' (see Details) return objects from fitting Q regressions
#' \item \code{Qstar} - list of length numLYnodes - \code{glm} (or numerical
#' optimization if \code{glm} fails to solve the score equation) return objects
#' from updating the Q fit }
#' 
#' \code{ltmleMSM} returns an object of class "\code{ltmleMSM}" The function
#' \code{\link{summary}} (i.e. \code{\link{summary.ltmleMSM}}) can be used to
#' obtain or print a summary of the results. An object of class
#' "\code{ltmleMSM}" is a list containing the following components:
#' \item{beta}{parameter estimates for working.msm using TMLE (GCOMP if
#' \code{gcomp} input is \code{TRUE})} \item{beta.iptw}{parameter estimates for
#' working.msm using IPTW} \item{IC}{matrix, n x numBetas - influence curve
#' values for TMLE (without updating if \code{gcomp} input is \code{TRUE})}
#' \item{IC.iptw}{matrix, n x numBetas - influence curve values for IPTW}
#' \item{msm}{object of class glm - the result of fitting the working.msm}
#' \item{cum.g}{array, n x numACnodes x numRegimes - cumulative g, after
#' bounding} \item{cum.g.unbounded}{array, n x numACnodes x numRegimes -
#' cumulative g, before bounding} \item{call}{the matched call}
#' \item{gcomp}{the \code{gcomp} input} \item{formulas}{a \code{list} with
#' elements \code{Qform} and \code{gform}} 
#' \item{fit}{a list with the following components} 
#' \itemize{ \item \code{g} - list of length numRegimes of list of length 
#' numACnodes - \code{glm} or \code{SuperLearner} (see Details) return objects from
#' fitting g regressions \item \code{Q} - list of length numLYnodes -
#' \code{glm} or \code{SuperLearner} (see Details) return objects from fitting Q 
#' regressions
#' \item \code{Qstar} - list of length numLYnodes - \code{glm} (or numerical
#' optimization if \code{glm} fails to solve the score equation) return objects
#' from updating the Q fit }
#' @author Joshua Schwab \email{jschwab77@berkeley.edu}, Samuel Lendle, Maya
#' Petersen, and Mark van der Laan
#' @seealso \code{\link{summary.ltmle}}, \code{\link{summary.ltmleMSM}},
#' \code{\link[SuperLearner:SuperLearner]{SuperLearner}},
#' \code{\link{deterministic.g.function_template}},
#' \code{\link{deterministic.Q.function_template}}
#' @examples
#' 
#' rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
#'                  
#' # Example 1: Single time point example. True value of E[Y_1] (expected value of
#' #   Y setting A to 1) is approximately 0.5939.
#' set.seed(2)
#' n <- 1000
#' W1 <- rnorm(n)
#' W2 <- rbinom(n, size=1, prob=0.3)   
#' W3 <- rnorm(n)
#' A <- rexpit(-1 + 2 * W1^2)
#' Y <- rexpit(-0.5 + 2 * W1^2 + 0.5 * W2 - 0.5 * A + 0.2 * W3 * A 
#'        - 1.1 * W3 + 0.2 * rnorm(n))
#' 
#' data <- data.frame(W1, W2, W3, A, Y)
#' 
#' 
#' \donttest{ #This takes about 4 seconds to run
#' library(SuperLearner)
#' 
#' #SuperLearner semiparametric estimation using all parents as regressors 
#' result1 <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", abar=1, 
#'   SL.library=c("SL.glm", "SL.step", "SL.mean"))
#' summary(result1)
#' summary(result1, estimator="iptw")
#' 
#' #SuperLearner semiparametric estimation using correctly specified regressors
#' result1a <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
#'  Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ I(W1^2)", abar=1, 
#'  SL.library=c("SL.glm", "SL.step", "SL.mean"))
#' summary(result1a)
#' }
#' 
#' #glm using correctly specified Qform and gform
#' result.abar1 <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
#'  Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ I(W1^2)", 
#'  abar=1, SL.library=NULL)
#' 
#' \donttest{ #This takes about 18 seconds to run
#' #Get summary measures (additive treatment effect, odds ratio, relative risk) 
#' #  for abar=1 vs abar=0
#' result.compare <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
#'                       Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ I(W1^2)", 
#'                       abar=list(1, 0), SL.library=NULL)
#' summary(result.compare)
#' 
#' 
#' # Example 2: Longitudinal example. Includes informative censoring and treatment. 
#' # Time ordering of data is W, C1, L1, A1, Y1, C2, L2, A2, Y2
#' # True value of E[Y_(1,1,1,1)] (expected value of Y setting C1, A1, C2, A2 all to 1)
#' #  is approximately 0.413.
#' # A1 is known to always be 1 if L1 < -2, and is 1 with probability 0.1 if L1 > 2 
#' # A2 is known to always be 1 if A1 is 1 
#' # We incorporate this knowledge using deterministic.g.function
#' 
#' # Generate data:
#' set.seed(2)
#' ua <- rep(TRUE, n)   #ua = uncensored and alive
#' L1 <- A1 <- Y1 <- C2.binary <- L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
#' W <- rnorm(n)
#' C1 <- BinaryToCensoring(is.uncensored=rexpit(2 + W))
#' ua <- ua & C1 == "uncensored"
#' L1[ua] <- rnorm(n)[ua] + W[ua]
#' A1[ua] <- rexpit(L1[ua])
#' A1[ua & L1 < -2] <- 1
#' A1[ua & L1 >  2] <- rbinom(n, size=1, prob=0.1)[ua & L1 >  2]
#' Y1[ua] <- rexpit((W + L1 - A1)[ua])
#' ua <- ua & !Y1
#' C2.binary[ua] <- rexpit((1 + 0.7 * L1 - A1)[ua])
#' C2 <- BinaryToCensoring(is.uncensored=C2.binary)
#' ua <- ua & C2 == "uncensored"
#' L2[ua] <- (0.5 * L1 - 0.9 * A1 + rnorm(n))[ua]
#' A2[ua] <- rexpit((0.5 * L1 + 0.8 * L2)[ua]) | A1[ua]
#' Y2[ua] <- rexpit((0.7 * L1 + L2 - 0.8 * A1 - A2)[ua])
#' Y2[Y1 == 1] <- 1  # if a patient dies at time 1, record death at time 2 as well
#' data <- data.frame(W, C1, L1, A1, Y1, C2, L2, A2, Y2)
#' 
#' deterministic.g.function <- function(data, current.node, nodes) {
#'   if (names(data)[current.node] == "A1") {
#'     det <- (data$L1 < -2 | data$L1 > 2) & !is.na(data$L1)
#'     prob1 <- ((data$L1 < -2) * 1 + (data$L1 > 2) * 0.1)[det]
#'   } else if (names(data)[current.node] == "A2") {
#'     det <- data$A1 == 1 & !is.na(data$A1)
#'     prob1 <- 1
#'   } else if (names(data[current.node]) %in% c("C1", "C2")){
#'     return(NULL)  #this returns the default of no deterministic links 
#'     #note that it is not necessary to specify that prior censoring indicates future censoring
#'   } else {
#'     stop("unexpected current.node")
#'   }
#'   return(list(is.deterministic=det, prob1=prob1))  
#' }
#' 
#' result2 <- ltmle(data, Anodes=c("A1","A2"), Cnodes=c("C1", "C2"), 
#'                 Lnodes=c("L1", "L2"), Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
#'                 deterministic.g.function=deterministic.g.function, survivalOutcome=TRUE)
#' summary(result2) 
#'  
#' # Example 3: Dynamic treatment, observation weights
#' # W -> A1 -> L -> A2 -> Y
#' # Treatment regime of interest is: Always treat at time 1 (A1 = 1), 
#' #   treat at at time 2 (A2 = 1), iff L > 0
#' # Weight by pmax(W + 2, 0)
#' 
#' set.seed(2)
#' n <- 1000
#' W <- rnorm(n)
#' A1 <- rexpit(W)
#' L <- 0.3 * W + 0.2 * A1 + rnorm(n)
#' A2 <- rexpit(W + A1 + L)
#' Y <- rexpit(W - 0.6 * A1 + L - 0.8 * A2)
#' data <- data.frame(W, A1, L, A2, Y)
#' 
#' abar <- matrix(nrow=n, ncol=2)
#' abar[, 1] <- 1
#' abar[, 2] <- L > 0
#' 
#' result3 <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", 
#'   survivalOutcome=TRUE, abar=abar, observation.weights = pmax(W + 2, 0))
#' summary(result3)
#' 
#' # Example 3.1: The regime can also be specified as a rule function
#' 
#' rule <- function(row) c(1, row["L"] > 0)
#' 
#' result.rule <- ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", 
#'   survivalOutcome=TRUE, rule=rule, observation.weights = pmax(W + 2, 0))
#' # This should be the same as the above result
#' summary(result.rule)
#' 
#' # Example 4: Deterministic Q function
#' # W -> A1 -> Y1 -> L2 -> A2 -> Y2
#' set.seed(2)
#' n <- 200
#' L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
#' W <- rnorm(n)
#' A1 <- rexpit(W)
#' Y1 <- rexpit(W - A1)
#' alive <- Y1 == 0
#' L2[alive] <- (0.5 * W - 0.9 * A1 + rnorm(n))[alive]
#' completed.study <- alive & L2 > 0
#' 
#' #Specify that Q is deterministically 0 when L2 is in the history of the 
#' # current Q regression and L2 > 0
#' #Note 1: det.Q.fun doesn't condition on called.from.estimate.g so g will also be set 
#' #        deterministically after L2 > 0 
#' #Note 2: It is not necessary to specify that Q is deterministically 1 if Y1 is 1; this is automatic
#' det.Q.fun.4a <- function(data, current.node, nodes, called.from.estimate.g) {
#'   L2.index <- which(names(data) == "L2")
#'   stopifnot(length(L2.index) == 1)
#'   L2.in.history <- L2.index < current.node
#'   if (! L2.in.history) return(NULL)
#'   
#'   is.deterministic <- data$L2 > 0 & !is.na(data$L2)
#'   return(list(is.deterministic=is.deterministic, Q.value=0))
#' }
#' 
#' #patients don't change treatment after leaving study; leave their A2 as NA
#' A2[alive & !completed.study] <- rexpit((0.5 * W + 0.8 * L2)[alive & !completed.study])
#' 
#' Y2[alive & !completed.study] <- rexpit((L2 - 0.8 * A1 - A2)[alive & !completed.study])
#' Y2[alive & completed.study] <- 0
#' Y2[!alive] <- 1  # if a patient dies at time 1, record death at time 2 as well
#' data <- data.frame(W, A1, Y1, L2, A2, Y2)
#' 
#' result4a <- ltmle(data, Anodes=c("A1","A2"), Lnodes="L2", Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
#'   SL.library=NULL, estimate.time=FALSE, deterministic.Q.function=det.Q.fun.4a, survivalOutcome=TRUE)
#' #note: You will get the same result if you pass Lnodes=NULL (see next example)
#' summary(result4a)
#' 
#' #In this variation, suppose that treatment can still change after a patient leaves the study
#' 
#' det.Q.fun.4b <- function(data, current.node, nodes, called.from.estimate.g) {
#'   #there is no deterministic information when calculating g - treatment may still change
#'   if (called.from.estimate.g) return(NULL)  
#'   
#'   L2.index <- which(names(data) == "L2")
#'   stopifnot(length(L2.index) == 1)
#'   L2.in.history <- L2.index < current.node
#'   if (! L2.in.history) return(NULL)
#'   
#'   is.deterministic <- data$L2 > 0 & !is.na(data$L2)
#'   return(list(is.deterministic=is.deterministic, Q.value=0))
#' }
#' 
#' A2[alive] <- rexpit((0.5 * W + 0.8 * L2)[alive])  #patients can change treatment after leaving study
#' Y2[alive & !completed.study] <- rexpit((L2 - 0.8 * A1 - A2)[alive & !completed.study])
#' Y2[alive & completed.study] <- 0
#' Y2[!alive] <- 1  # if a patient dies at time 1, record death at time 2 as well
#' data <- data.frame(W, A1, Y1, L2, A2, Y2)
#' 
#' result4b <- ltmle(data, Anodes=c("A1","A2"), Lnodes="L2", Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
#'  SL.library=NULL, estimate.time=FALSE, deterministic.Q.function=det.Q.fun.4b, survivalOutcome=TRUE)
#' summary(result4b)
#' 
#' # Example 5: Multiple time-dependent covariates and treatments at each time point, 
#' #            continuous Y values
#' # age -> gender -> A1 -> L1a -> L1b -> Y1 -> A2 -> L2a -> L2b -> Y2
#' set.seed(2)
#' n <- 100
#' age <- rbinom(n, 1, 0.5)
#' gender <- rbinom(n, 1, 0.5)
#' A1 <- rexpit(age + gender)
#' L1a <- 2*age - 3*gender + 2*A1 + rnorm(n)
#' L1b <- rexpit(age + 1.5*gender - A1)
#' Y1 <- plogis(age - gender + L1a + 0.7*L1b + A1 + rnorm(n))
#' A2 <- rexpit(age + gender + A1 - L1a - L1b)
#' L2a <- 2*age - 3*gender + 2*A1 + A2 + rnorm(n)
#' L2b <- rexpit(age + 1.5*gender - A1 - A2)
#' Y2 <- plogis(age - gender + L1a + L1b + A1 + 1.8*A2 + rnorm(n))
#' data <- data.frame(age, gender, A1, L1a, L1b, Y1, A2, L2a, L2b, Y2)
#' 
#' #Note that gform is not correctly specified in these examples.
#' 
#' #Also show some different ways of specifying the nodes:
#' 
#' result5a <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
#'  Ynodes=grep("^Y", names(data)), abar=c(1, 0), SL.library=NULL, estimate.time=FALSE, 
#'  survivalOutcome=FALSE, gform=c("A1 ~ gender", "A2 ~ age")) 
#' summary(result5a)
#' 
#' #Usually you would specify a Qform for all of the Lnodes and Ynodes but in this case 
#' # L1a, L1b, Y1 is a "block" of L/Y nodes not separated by Anodes or Cnodes (the same is true for 
#' # L2a, L2b, Y2). Only one regression is required at the first L/Y node in a block. You can pass 
#' # regression formulas for the other L/Y nodes, but they'll be ignored.
#' result5b <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
#'  Ynodes=grep("^Y", names(data)), abar=c(1, 0), estimate.time=FALSE, survivalOutcome=FALSE, 
#'  gform=c("A1 ~ gender", "A2 ~ age"), Qform=c(L1a="Q.kplus1 ~ 1", L2a="Q.kplus1 ~ 1"))
#' summary(result5b)
#' 
#' 
#' #Gives the same result but prints a message saying some regression formulas will be dropped:
#' result5c <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
#'  Ynodes=grep("^Y", names(data)), abar=c(1, 0), estimate.time=FALSE, survivalOutcome=FALSE, 
#'  gform=c("A1 ~ gender", "A2 ~ age"), Qform=c(L1a="Q.kplus1 ~ 1", L1b="Q.klus1~A1", 
#'  Y1="Q.kplus1~L1a", L2a="Q.kplus1 ~ 1", L2b="Q.klus1~A1", Y2="Q.kplus1~A2 + gender"))
#' 
#' summary(result5c)
#' 
#' 
#' #If there were a Anode or Cnode between L1b and Y1, Y1 would also need a Q regression formula
#' 
#' 
#' # Example 6: MSM
#' # Given data over 3 time points where A switches to 1 once and then stays 1. We want to know
#' # how death varies as a function of gender, time and an indicator of whether a patient's 
#' # intended regime was to switch before time.
#' # Note that working.msm includes time and switch.time, which are columns of 
#' # summary.measures; working.msm also includes male, which is ok because it is a baseline
#' # covariate (it comes before any A/C/L/Y nodes).
#' data(sampleDataForLtmleMSM)
#' Anodes <- grep("^A", names(sampleDataForLtmleMSM$data))
#' Lnodes <- c("CD4_1", "CD4_2")
#' Ynodes <- grep("^Y", names(sampleDataForLtmleMSM$data))
#' msm.weights <- matrix(1:12, nrow=4, ncol=3) #just an example (can also use a 200x3x4 array), 
#'                                             #or NULL (for no weights), or "empirical" (the default)
#' 
#' result6 <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
#'                    survivalOutcome=TRUE,
#'                    regimes=sampleDataForLtmleMSM$regimes, 
#'                    summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, 
#'                    working.msm="Y ~ male + time + I(pmax(time - switch.time, 0))", 
#'                    msm.weights=msm.weights, estimate.time=FALSE)
#' print(summary(result6))
#' 
#' 
#' # Example 6.1: regimes can also be specified as a list of rule functions
#' 
#' regimesList <- list(function(row) c(1,1,1),
#'                      function(row) c(0,1,1),
#'                      function(row) c(0,0,1),
#'                      function(row) c(0,0,0))
#' result.regList <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
#'                    survivalOutcome=TRUE, regimes=regimesList, 
#'                    summary.measures=sampleDataForLtmleMSM$summary.measures, final.Ynodes=Ynodes, 
#'                    working.msm="Y ~ male + time + I(pmax(time - switch.time, 0))", 
#'                    msm.weights=msm.weights, estimate.time=FALSE)
#' # This should be the same as the above result
#' print(summary(result.regList))         
#' 
#' 
#' # Example 7: variance estimation
#' # A simple point treatment problem W, A, Y. But there is a positivity problem - 
#' # for small values of W, Prob(A = 1) is very small.
#' # The true parameter value, E[Y_1] is approximately 0.697
#' # The true TMLE standard deviation is approximately 0.064, 
#' # the true IPTW standard deviation is approximately 0.058.
#' set.seed(2)
#' n <- 1000
#' W <- rnorm(n)
#' A <- rexpit(8 * W)
#' Y <- rexpit(W + A)
#' r1 <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE)
#' r2 <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, 
#'  variance.method="ic")
#' r3 <- ltmle(data.frame(W, A, Y), Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, 
#'  variance.method="iptw")
#' print(summary(r1))
#' print(summary(r2))
#' print(summary(r3))
#' print(summary(r1, estimator="iptw"))
#' print(summary(r2, estimator="iptw")) #the same - variance.method only affects TMLE
#' print(summary(r3, estimator="iptw")) #the same - variance.method only affects TMLE
#' }
#' 
#' @export ltmle
ltmle <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=NULL, Qform=NULL, gform=NULL, abar, rule=NULL, gbounds=c(0.01, 1), Yrange=NULL, deterministic.g.function=NULL, stratify=FALSE, SL.library="glm", estimate.time=TRUE, gcomp=FALSE, iptw.only=FALSE, deterministic.Q.function=NULL, variance.method="tmle", observation.weights=NULL, id=NULL) {
  msm.inputs <- GetMSMInputsForLtmle(data, abar, rule, gform)
  inputs <- CreateInputs(data=data, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=survivalOutcome, Qform=Qform, gform=msm.inputs$gform, Yrange=Yrange, gbounds=gbounds, deterministic.g.function=deterministic.g.function, SL.library=SL.library, regimes=msm.inputs$regimes, working.msm=msm.inputs$working.msm, summary.measures=msm.inputs$summary.measures, final.Ynodes=msm.inputs$final.Ynodes, stratify=stratify, msm.weights=msm.inputs$msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, variance.method=variance.method, observation.weights=observation.weights, id=id) 
  result <- LtmleFromInputs(inputs)
  result$call <- match.call()
  return(result)
}
# General code flow:
#  ltmle -> CreateInputs -> LtmleFromInputs -> LtmleMSMFromInputs -> ...
#  ltmleMSM -> CreateInputs -> LtmleMSMFromInputs -> ...

RegimesFromAbar <- function(data, abar, rule) {
  if (!is.null(rule)) {
    if (!(missing(abar) || is.null(abar))) stop("'abar' should not be specified when using a 'rule' function")
    rule.output.length <- length(as.numeric(rule(data[1, ])))
    if (all(sapply(data, is.numeric))) {
      abar <- apply(data, 1, rule)
      if (rule.output.length == 1) {
        abar <- AsMatrix(abar)
      } else {
        abar <- t(abar)
      }
    } else {
      #apply doesn't work well if there are factors (apply calls as.matrix and everything becomes a character)
      if (nrow(data) > 10000) warning("Using a rule input may be significantly slower than using abar/regimes if your data contains censoring nodes or is otherwise not all numeric.")
      abar <- matrix(nrow=nrow(data), ncol=rule.output.length)
      for (i in 1:nrow(data)) {
        abar[i, ] <- as.numeric(rule(data[i, ]))
      }
    }
  }
  if (is.vector(abar)) {
    abar <- matrix(rep(abar, each=nrow(data)), nrow=nrow(data))
  } else if (is.null(abar)) {
    abar <- matrix(numeric(0), nrow=nrow(data), ncol=0)
  } else if (is.function(abar)) {
    stop("abar should be a vector or matrix, not a function. Use the 'rule' paramter to pass a function instead.")
  }
  regimes <- abar
  dim(regimes) <- c(nrow(regimes), ncol(regimes), 1)
  return(regimes)
}

# ltmle is a special case of ltmleMSM - get the arguments used by ltmleMSM for the special case
GetMSMInputsForLtmle <- function(data, abar, rule, gform) {
  if ((!missing(abar) && is.list(abar)) || is.list(rule)) {
    if (is.list(rule)) {
      if (length(rule) != 2) stop("If rule is a list, it must be of length 2")
      regimes1 <- RegimesFromAbar(data, abar, rule[[1]])
      regimes0 <- RegimesFromAbar(data, abar, rule[[2]])
    } else {
      if (length(abar) != 2) stop("If abar is a list, it must be of length 2")
      regimes1 <- RegimesFromAbar(data, abar[[1]], rule)
      regimes0 <- RegimesFromAbar(data, abar[[2]], rule)
    }
    if (ncol(regimes1) != ncol(regimes0)) stop("If abar or rule is a list, both elements must give a matrix with the same number of columns")
    if (nrow(regimes1) != nrow(regimes0)) stop("If abar or rule is a list, both elements must give a matrix with the same number of rows")
    regimes <- c(regimes1, regimes0)
    dim(regimes) <- c(nrow(regimes1), ncol(regimes1), 2)
    summary.measures <- array(1:0, dim=c(2, 1, 1))
    colnames(summary.measures) <- "A"
    working.msm <- "Y ~ A"
    msm.weights <- matrix(1, nrow=2, ncol=1)
  } else {
    regimes <- RegimesFromAbar(data, abar, rule)
    working.msm <- "Y ~ 1"
    msm.weights <- matrix(1, nrow=1, ncol=1)
    summary.measures <- array(dim=c(1, 0, 1))
  }
  msm.inputs <- list(regimes=regimes, working.msm=working.msm, summary.measures=summary.measures, gform=gform, final.Ynodes=NULL, msm.weights=msm.weights)
  return(msm.inputs)
}

# run ltmle from the ltmleInputs object
LtmleFromInputs <- function(inputs) {
  msm.result <- LtmleMSMFromInputs(inputs)
  num.regimes <- dim(inputs$regimes)[3]
  stopifnot(num.regimes %in% 1:2)
  if (num.regimes == 2) {
    class(msm.result) <- "ltmleEffectMeasures"
    return(msm.result)
  }
  names(msm.result$beta.iptw) <- names(msm.result$beta) <- NULL
  iptw <- plogis(msm.result$beta.iptw)
  iptw.list <- list(iptw.estimate=iptw, iptw.IC=iptw*(1-iptw)*msm.result$IC.iptw[, 1])
  
  r <- list()
  if (inputs$iptw.only) {
    tmle <- NA
    tmle.IC <- rep(NA, nrow(inputs$data))
  } else {
    tmle <- plogis(msm.result$beta)
    tmle.IC <- msm.result$IC[, 1] #only one regime
  }
  r$estimates <- c(tmle=tmle, iptw=iptw.list$iptw.estimate)
  r$IC <- list(tmle=tmle.IC * tmle * (1 - tmle), iptw=iptw.list$iptw.IC)
  if (!is.null(msm.result$variance.estimate)) {
    stopifnot(length(msm.result$variance.estimate) == 1)
    r$variance.estimate <- msm.result$variance.estimate[1] * (tmle * (1 - tmle))^2 
  }
  
  if (inputs$gcomp) {
    names(r$estimates)[1] <- names(r$IC)[1] <- "gcomp"
  }
  
  r$cum.g <- AsMatrix(msm.result$cum.g[, , 1]) #only one regime
  r$cum.g.unbounded <- AsMatrix(msm.result$cum.g.unbounded[, , 1]) #only one regime
  r$cum.g.used <- AsMatrix(msm.result$cum.g.used[, , 1]) #only one regime
  r$gcomp <- inputs$gcomp
  r$fit <- msm.result$fit
  r$fit$g <- r$fit$g[[1]]  #only one regime
  r$fit$Q <- r$fit$Q[[1]]  #only one regime 
  r$Qstar <- msm.result$Qstar[, 1, 1] #1 regime, 1 final.Ynode
  
  r$formulas <- msm.result$formulas
  r$binaryOutcome <- msm.result$binaryOutcome
  r$transformOutcome <- msm.result$transformOutcome==TRUE #Want to store transformOutcome flag without attributes
  
  if (msm.result$transformOutcome) {
    Yrange <- attr(msm.result$transformOutcome, "Yrange")
    #back transform estimate and IC
    r$estimates <- r$estimates*diff(Yrange) + min(Yrange)  
    r$IC <- lapply(r$IC, function (IC) IC * diff(Yrange))
    r$variance.estimate <- r$variance.estimate * (diff(Yrange))^2 
  }
  class(r) <- "ltmle"
  return(r)
}

#' @describeIn ltmle Longitudinal Targeted Maximum Likelihood Estimation for a Marginal Structural Model
#' @export 
ltmleMSM <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=NULL, Qform=NULL, gform=NULL, gbounds=c(0.01, 1), Yrange=NULL, deterministic.g.function=NULL, SL.library="glm", regimes, working.msm, summary.measures, final.Ynodes=NULL, stratify=FALSE, msm.weights="empirical", estimate.time=TRUE, gcomp=FALSE, iptw.only=FALSE, deterministic.Q.function=NULL, variance.method="tmle", observation.weights=NULL, id=NULL) {
  inputs <- CreateInputs(data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, estimate.time, gcomp, iptw.only, deterministic.Q.function, variance.method, observation.weights, id)
  result <- LtmleMSMFromInputs(inputs)
  result$call <- match.call()
  class(result) <- "ltmleMSM"
  return(result) 
}

# run ltmleMSM from ltmleInputs object
LtmleMSMFromInputs <- function(inputs) {  
  if (inputs$estimate.time) EstimateTime(inputs)
  result <- MainCalcs(inputs)
  result$gcomp <- inputs$gcomp
  result$formulas <- list(Qform=inputs$Qform, gform=inputs$gform)
  result$binaryOutcome <- inputs$binaryOutcome
  result$transformOutcome <- inputs$transformOutcome
  result$survivalOutcome <- inputs$survivalOutcome
  return(result)
}

# create the ltmleInputs object used by many other functions - fills in defaults and does error checking
CreateInputs <- function(data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, estimate.time, gcomp, iptw.only, deterministic.Q.function, variance.method, observation.weights, id) {
  if (is.list(regimes)) {
    if (!all(sapply(regimes, is.function))) stop("If 'regimes' is a list, then all elements should be functions.")
    regimes <- simplify2array(lapply(regimes, function (rule) drop3(RegimesFromAbar(data, rule = rule))), higher=TRUE)
  }
  if (!(is.null(regimes) || length(dim(regimes)) == 3)) {
    stop("regimes must be an array with 3 dimensions (unless Anodes is NULL, in which case regimes can be NULL)")
  }
  if (is.null(regimes) || dim(regimes)[3]==0) {
    if (length(Anodes) != 0) {
      stop("regimes must not be NULL (or have dim(regimes)[3]==0) unless Anodes is also NULL")
    }
    regimes <- array(numeric(0), dim=c(nrow(data), 0, 1))
  }
  num.regimes <- dim(regimes)[3]
  if (is.logical(regimes)) {
    regimes <- regimes * 1
    message("abar or regimes was passed as logical and was converted to numeric")
  }
  all.nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes)
  Qform <- CreateLYNodes(data, all.nodes, check.Qform=TRUE, Qform=Qform)$Qform
  data <- ConvertCensoringNodes(data, Cnodes, has.deterministic.functions=!is.null(deterministic.g.function) && is.null(deterministic.Q.function))
  if (is.null(final.Ynodes)) {
    final.Ynodes <- max(all.nodes$Y)
  } else {
    final.Ynodes <- NodeToIndex(data, final.Ynodes)
  }
  
  #Using get to avoid the "no visible binding for global variable" note in R CMD check
  if (identical(SL.library, 'default')) SL.library <- get("Default.SL.Library")
  SL.library.Q <- GetLibrary(SL.library, "Q")
  SL.library.g <- GetLibrary(SL.library, "g")
  
  if (is.null(summary.measures)) {
    summary.measures <- matrix(nrow=num.regimes, ncol=0)
  }
  if (length(dim(summary.measures)) == 2) {
    num.final.Ynodes <- length(final.Ynodes)
    summary.measures <- array(repmat(summary.measures, m=1, n=num.final.Ynodes), dim=c(nrow(summary.measures), ncol(summary.measures), num.final.Ynodes), dimnames=list(rownames(summary.measures), colnames(summary.measures), NULL))
  }
  if (is.null(observation.weights)) observation.weights <- rep(1, nrow(data))
  
  if (is.matrix(gform)) {
    if (num.regimes > 1 && variance.method != "ic") stop("If there is more than one regime (using ltmle with list abar or ltmleMSM) and numeric gform and variance.method != 'ic', then gform must be an array, not a matrix.")
    gform <- array(gform, dim=c(nrow(gform), ncol(gform), num.regimes))
  }
  
  #error checking (also get value for survivalOutcome if NULL)
  check.results <- CheckInputs(data, all.nodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, deterministic.Q.function, observation.weights, gcomp, variance.method, id) 
  survivalOutcome <- check.results$survivalOutcome
  
  if (!isTRUE(attr(data, "called.from.estimate.variance", exact=TRUE)) && !isTRUE(attr(data, "skip.clean.data", exact=TRUE))) { 
    data <- CleanData(data, all.nodes, deterministic.Q.function, survivalOutcome)
  }
  transform.list <- TransformOutcomes(data, all.nodes, Yrange)
  data <- transform.list$data
  transformOutcome <- transform.list$transformOutcome
  binaryOutcome <- check.results$binaryOutcome
  
  if (length(Qform) == 0) Qform <- GetDefaultForm(data, all.nodes, is.Qform=TRUE, stratify, survivalOutcome, showMessage=TRUE)
  if (length(gform) == 0) gform <- GetDefaultForm(data, all.nodes, is.Qform=FALSE, stratify, survivalOutcome, showMessage=TRUE)
  
  # Several functions in the pooled version are only written to accept main terms MSM
  # Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S1 + S3 + S4" where 
  # S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
  main.terms <- ConvertToMainTerms(data, working.msm, summary.measures, all.nodes)
  intervention.match <- CalcInterventionMatchArray(data, regimes, all.nodes$A)
  
  inputs <- list(data=data, all.nodes=all.nodes, survivalOutcome=survivalOutcome, Qform=Qform, gform=gform, gbounds=gbounds, Yrange=Yrange, deterministic.g.function=deterministic.g.function, SL.library.Q=SL.library.Q, SL.library.g=SL.library.g, regimes=regimes, working.msm=main.terms$msm, combined.summary.measures=main.terms$summary.measures, final.Ynodes=final.Ynodes, stratify=stratify, msm.weights=msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, binaryOutcome=binaryOutcome, transformOutcome=transformOutcome, variance.method=variance.method, observation.weights=observation.weights, baseline.column.names=main.terms$baseline.column.names, beta.names=main.terms$beta.names, uncensored=check.results$uncensored, intervention.match=intervention.match, id=id)
  class(inputs) <- "ltmleInputs"
  if (length(all.nodes$AC) == 0) inputs$variance.method <- "ic" #no warning for this case
  if (inputs$variance.method != "ic" && !is.null(VarianceAvailableWarning(inputs))) inputs$variance.method <- "ic"
  return(inputs)
}

# Prevent misspelled argument after $ from returning NULL
`$.ltmleInputs` <- function(x, name) {
  if (! (name %in% names(x))) stop(paste(name, "is not an element of x"))
  return(x[[name]])
}

# Loop over final Ynodes, run main calculations
MainCalcs <- function(inputs) {
  num.final.Ynodes <- length(inputs$final.Ynodes)
  num.betas <- dim(inputs$combined.summary.measures)[2]
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes       note: num.measures is summary measures and baseline covariates, converted to main terms
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  Qstar <- array(dim=c(n, num.regimes, num.final.Ynodes))
  all.msm.weights <- GetMsmWeights(inputs) #n x num.regimes x num.final.Ynodes
  new.var.y <- array(dim=c(num.betas, num.betas, num.final.Ynodes))
  IC <- matrix(0, n, num.betas)
  #store IC for each final Ynode, compare var(IC) to sum(var(IC.ynode))
  IC.y <- array(dim=c(n, num.betas, num.final.Ynodes))
  
  g.list <- EstimateG(inputs)
  iptw <- CalcIPTW(inputs, g.list$cum.g, all.msm.weights)
  fit <- list(g=g.list$fit) 
  if (inputs$iptw.only) {
    beta <- rep(NA, length(iptw$beta))
    fitted.msm <- NULL
    variance.estimate <- NULL
    fixed.tmle <- list(cum.g.used=array(NA, dim=dim(g.list$cum.g)))
  } else {
    for (j in 1:num.final.Ynodes) {
      fixed.tmle <- FixedTimeTMLE(inputs, nodes = SubsetNodes(inputs$all.nodes, final.Ynode=inputs$final.Ynodes[j]), msm.weights = drop3(all.msm.weights[, , j, drop=FALSE]), combined.summary.measures = dropn(inputs$combined.summary.measures[, , , j, drop=FALSE], n=4), g.list = g.list)
      IC <- IC + fixed.tmle$IC
      IC.y[, , j] <- fixed.tmle$IC
      Qstar[, , j] <- fixed.tmle$Qstar # n x num.regimes
      new.var.y[, , j] <- fixed.tmle$est.var 
    }
    fit <- c(fit, fixed.tmle$fit)
    if (isTRUE(attr(inputs$data, "called.from.estimate.variance", exact=TRUE))) { 
      return(list(IC=matrix(NA, 1, 1), msm=NULL, beta=qlogis(mean(Qstar)), cum.g=g.list$cum.g, cum.g.unbounded=g.list$cum.g.unbounded, fit=fit, variance.estimate=NULL, beta.iptw=iptw$beta, IC.iptw=iptw$IC, Qstar=Qstar, cum.g.used=fixed.tmle$cum.g.used))
    }
    fitted.msm <- FitPooledMSM(inputs$working.msm, Qstar, inputs$combined.summary.measures, all.msm.weights * inputs$observation.weights) 
    IC <- FinalizeIC(IC, inputs$combined.summary.measures, Qstar, fitted.msm$m.beta, all.msm.weights, inputs$observation.weights, inputs$id) #n x num.betas
    C.old <- NormalizeIC(IC, inputs$combined.summary.measures, fitted.msm$m.beta, all.msm.weights, inputs$observation.weights, g.ratio = NULL) #C without using g.ratio 
    g.ratio <- CalcGUnboundedToBoundedRatio(g.list, inputs$all.nodes, inputs$final.Ynodes)
    CheckForVarianceWarning(inputs, g.ratio)
    if (inputs$variance.method == "ic") {   
      variance.estimate <- NULL
    } else {
      new.var <- matrix(NA, num.betas, num.betas)
      for (i in 1:num.betas) {
        for (j in 1:num.betas) {
          if (num.final.Ynodes > 1) {
            cov.IC <- cov(IC.y[, i, ], IC.y[, j, ])
            diag(cov.IC) <- new.var.y[i, j, ]
            new.var[i, j] <- sum(cov.IC)
          } else {
            new.var[i, j] <- new.var.y[i, j, 1]
          }
        }
      }
      
      C <- NormalizeIC(IC, inputs$combined.summary.measures, fitted.msm$m.beta, all.msm.weights, inputs$observation.weights, g.ratio)
      variance.estimate <- safe.solve(C) %*% new.var %*% t(safe.solve(C))
    }
    IC <- t(safe.solve(C.old, t(IC))) #IC %*% solve(C) 
    beta <- coef(fitted.msm$m)
    names(beta) <- inputs$beta.names
  }
  
  return(list(IC=IC, msm=fitted.msm$m, beta=beta, cum.g=g.list$cum.g, cum.g.unbounded=g.list$cum.g.unbounded, fit=fit, variance.estimate=variance.estimate, beta.iptw=iptw$beta, IC.iptw=iptw$IC, Qstar=Qstar, cum.g.used=fixed.tmle$cum.g.used)) #note: only returns cum.g and fit and cum.g.used for the last final.Ynode
}

VarianceAvailableWarning <- function(inputs) {
  if (!inputs$binaryOutcome) return("Robust variance estimate is not currently available with non binary outcomes")
  if (!is.null(inputs$deterministic.Q.function)) return("Robust variance estimate is not currently available with deterministic.Q.function")
  if (inputs$gcomp) return("Robust variance estimate is not currently available with gcomp")
  if (inputs$stratify) return("Robust variance estimate is not currently available with stratify=TRUE")
  if (!is.null(inputs$id)) return("Robust variance estimate is not currently available with non-null id")
  return(NULL)
}

CheckForVarianceWarning <- function(inputs, g.ratio) {
  if (inputs$variance.method == "ic") {
    positivity <- mean(g.ratio < 1, na.rm=TRUE) > 0.01
    rare.events <- inputs$binaryOutcome && (colMeans(inputs$data[, inputs$final.Ynodes, drop=FALSE], na.rm=TRUE) < 0.03) 
    if (positivity || rare.events) {
      variance.available.warning <- VarianceAvailableWarning(inputs)
      warning.msg <- "Variance estimate is based on influence curve only, which may be significantly anticonservative because your data appears to contain"
      if (positivity) warning.msg <- paste(warning.msg, "positivity violations")
      if (positivity && rare.events) warning.msg <- paste(warning.msg, "and")
      if (rare.events) warning.msg <- paste(warning.msg, "rare events")
      if (is.null(variance.available.warning)) {
        warning.msg <- paste0(warning.msg, ". It is recommended to use variance.method='tmle' or variance.method='iptw' to obtain a more robust variance estimate (but run time may be significantly longer). See variance.method details in ?ltmle")
      } else {
        warning.msg <- paste0(warning.msg, ". ", variance.available.warning, " but this will be addressed in a future release.")
      }
      warning(warning.msg)
    }
  }
  invisible(NULL)
}

CalcIPTW <- function(inputs, cum.g, msm.weights) {
  if (isTRUE(attr(inputs$data, "called.from.estimate.variance", exact=TRUE))) { 
    return(list(beta=NA, IC=matrix(NA, 1, 1)))
  }
  nodes <- inputs$all.nodes
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  num.final.Ynodes <- length(inputs$final.Ynodes)
  Y.vec <- X.mat <- weight.vec <- NULL
  save.xy <- list()
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- inputs$final.Ynodes[j]
    intervention.match <- InterventionMatch(inputs$intervention.match, nodes$A, cur.node=final.Ynode)
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, final.Ynode)
    for (i in 1:num.regimes) {
      index <- uncensored & intervention.match[, i]
      col.index <- which.max(nodes$AC[nodes$AC < final.Ynode]) 
      Y <- inputs$data[index, final.Ynode]
      if (length(col.index > 0)) {
        g <- cum.g[index, col.index, i] 
      } else {
        g <- 1
      }
      X <- inputs$combined.summary.measures[index, , i, j]
      if (is.vector(X)) { #if only one summary.measure or sum(index)==1, X is dropped to vector
        dim(X) <- c(sum(index), ncol(inputs$combined.summary.measures))
      }
      weight <- msm.weights[index, i, j] * inputs$observation.weights[index] / g
      weight[msm.weights[index, i, j] == 0 | inputs$observation.weights[index] == 0] <- 0 #avoid problems where weight and g are both 0
      
      save.xy[[length(save.xy) + 1]] <- list(X=X, Y=Y, weight=weight, index=index)
      Y.vec <- c(Y.vec, Y)
      X.mat <- rbind(X.mat, X)
      weight.vec <- c(weight.vec, weight) 
    }
  }
  colnames(X.mat) <- colnames(inputs$combined.summary.measures)
  
  if (nrow(X.mat) == 0) {
    #this happens if there are no rows uncensored and intervention.match
    warning("no rows uncensored and matching regimes/abar - IPTW returns NA")
    num.beta <- ncol(inputs$combined.summary.measures)
    return(list(beta=rep(NA, num.beta), IC=matrix(nrow=n, ncol=num.beta)))
  }
  m.glm <- ltmle.glm(formula(inputs$working.msm), family=quasibinomial(), data=data.frame(Y=Y.vec, X.mat, weight.vec), weights=as.vector(scale(weight.vec, center=FALSE))) #note: scale weights because there were rare problems where large weights caused convergence problems
  beta <- coef(m.glm)
  IC <- matrix(0, nrow=n, ncol=length(beta))  #n x num.betas
  m.beta <- array(dim=c(n, num.regimes, num.final.Ynodes)) 
  cnt <- 0
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- inputs$final.Ynodes[j]
    for (i in 1:num.regimes) {
      newdata <- data.frame(inputs$combined.summary.measures[, , i, j])
      colnames(newdata) <- colnames(inputs$combined.summary.measures) #needed if only one summary measure
      SuppressGivenWarnings(m.beta[, i, j] <- predict(m.glm, newdata=newdata, type="response"), "prediction from a rank-deficient fit may be misleading")
      
      cnt <- cnt + 1
      XY.list <- save.xy[[cnt]]
      IC[XY.list$index, ] <- IC[XY.list$index, ] + XY.list$weight * XY.list$X * (XY.list$Y - m.beta[XY.list$index, i, j]) #recycles weight, Y, m.beta
    }
  }
  
  C <- NormalizeIC(IC, inputs$combined.summary.measures, m.beta, msm.weights, observation.weights=inputs$observation.weights, g.ratio=NULL)
  normalized.IC <- t(safe.solve(C, t(IC)))  
  household.IC <- HouseholdIC(normalized.IC, inputs$id)
  names(beta) <- inputs$beta.names
  return(list(beta=beta, IC=household.IC))
}

# ltmleMSM for a single final.Ynode
FixedTimeTMLE <- function(inputs, nodes, msm.weights, combined.summary.measures, g.list) {
  #combined.summary.measures: n x num.measures x num.regimes   (num.measures=num.summary.measures + num.baseline.covariates)
  data <- inputs$data
  
  num.regimes <- dim(inputs$regimes)[3]
  n <- nrow(data)
  num.betas <- ncol(combined.summary.measures)
  tmle <- rep(NA, num.regimes)
  IC <- matrix(0, nrow=n, ncol=num.betas)
  cum.g.used <- array(FALSE, dim=dim(g.list$cum.g)) #n x num.AC.nodes x num.regimes
  
  est.var <- matrix(0, num.betas, num.betas)  
  regimes.with.positive.weight <- which(apply(msm.weights > 0, 2, any))
  if (length(regimes.with.positive.weight) == 0) stop("All regimes have weight 0 (one possible reason is that msm.weights='emipirical' and no data rows match any of the regimes and are uncensored)")
  fit.Qstar <- fit.Q <- vector("list", length(nodes$LY))
  names(fit.Qstar) <- names(fit.Q) <- names(data)[nodes$LY]
  Qstar.kplus1 <- matrix(data[, max(nodes$Y)], nrow=n, ncol=num.regimes)
  mean.summary.measures <- apply(abs(combined.summary.measures), 2, mean)
  if (length(nodes$LY) > 0) {
    for (LYnode.index in length(nodes$LY):1) {
      cur.node <- nodes$LY[LYnode.index]
      deterministic.list.origdata <- IsDeterministic(data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=FALSE, inputs$survivalOutcome)
      uncensored <- IsUncensored(inputs$uncensored, nodes$C, cur.node)
      intervention.match <- InterventionMatch(inputs$intervention.match, nodes$A, cur.node)
      if (inputs$stratify) {
        subs <- uncensored & intervention.match & !deterministic.list.origdata$is.deterministic #n x num.regimes
      } else {
        subs <- uncensored & !deterministic.list.origdata$is.deterministic #vector
      }
      Q.est <- Estimate(inputs, form = inputs$Qform[LYnode.index], Qstar.kplus1=if (LYnode.index == length(nodes$LY)) Qstar.kplus1[, 1] else Qstar.kplus1, family=quasibinomial(), subs=subs, type="link", nodes=nodes, called.from.estimate.g=FALSE, calc.meanL=FALSE, cur.node=cur.node, regimes.meanL=NULL, regimes.with.positive.weight=regimes.with.positive.weight) #if this is the last node, only pass the first column as a vector
      logitQ <- Q.est$predicted.values
      fit.Q[[LYnode.index]] <- Q.est$fit 
      ACnode.index  <- which.max(nodes$AC[nodes$AC < cur.node])
      SuppressGivenWarnings(update.list <- UpdateQ(Qstar.kplus1, logitQ, combined.summary.measures, g.list$cum.g[, ACnode.index, ], inputs$working.msm, uncensored, intervention.match, deterministic.list.origdata$is.deterministic, msm.weights, inputs$gcomp, inputs$observation.weights), GetWarningsToSuppress(update.step = TRUE))
      if (length(ACnode.index) > 0) cum.g.used[, ACnode.index, ] <- cum.g.used[, ACnode.index, ] | update.list$cum.g.used 
      Qstar <- update.list$Qstar
      Qstar[Q.est$is.deterministic] <- Q.est$deterministic.Q[Q.est$is.deterministic] #matrix indexing
      curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
      curIC.relative.error <- abs(colSums(curIC)) 
      curIC.relative.error[mean.summary.measures > 0] <- curIC.relative.error[mean.summary.measures > 0] / mean.summary.measures[mean.summary.measures > 0]
      if (any(curIC.relative.error > 0.001) && !inputs$gcomp) {
        SetSeedIfRegressionTesting()
        fix.score.list <- FixScoreEquation(Qstar.kplus1, update.list$h.g.ratio, uncensored, intervention.match, Q.est$is.deterministic, Q.est$deterministic.Q, update.list$off, update.list$X, regimes.with.positive.weight)
        Qstar <- fix.score.list$Qstar
        curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
        update.list$fit <- fix.score.list$fit      
      }
      est.var <- est.var + EstimateVariance(inputs, nodes, combined.summary.measures, regimes.with.positive.weight, uncensored, alive=!deterministic.list.origdata$is.deterministic, Qstar, Qstar.kplus1, cur.node, msm.weights, LYnode.index, ACnode.index, g.list$cum.g, g.list$prob.A.is.1, g.list$cum.g.meanL, g.list$cum.g.unbounded, g.list$cum.g.meanL.unbounded, inputs$observation.weights, is.last.LYnode=(LYnode.index==length(nodes$LY)), intervention.match) 
      IC <- IC + curIC 
      Qstar.kplus1 <- Qstar
      fit.Qstar[[LYnode.index]] <- update.list$fit
    }
  } else {
    Qstar <- Qstar.kplus1 #if there are no LY nodes (can happen if no A/C nodes before a final.Ynode)
  }
  #tmle <- colMeans(Qstar)
  return(list(IC=IC, Qstar=Qstar, est.var=est.var, fit=list(Q=ReorderFits(fit.Q), Qstar=fit.Qstar), cum.g.used=cum.g.used)) 
}

EstimateVariance <- function(inputs, nodes, combined.summary.measures, regimes.with.positive.weight, uncensored, alive, Qstar, Qstar.kplus1, cur.node, msm.weights, LYnode.index, ACnode.index, cum.g, prob.A.is.1, cum.g.meanL, cum.g.unbounded, cum.g.meanL.unbounded, observation.weights, is.last.LYnode, intervention.match) { 
  if (inputs$variance.method == "ic") return(NA)
  est.var.iptw <- inputs$variance.method == "iptw"
  TmleOfVariance <- function(Z, Z.meanL) {
    if (all(is.na(Z))) stop("all Z are NA in EstimateVariance")
    #length(Z.meanL) == 0 --> point treatment case, return mean(Z); all Z can be zero when summary.measure is 0
    if (length(Z.meanL) == 0 || all(Z==0 | is.na(Z))) {
      Qstar <- Scale(Z, 0, 1)
      return(list(EZd1 = mean(Z, na.rm=T), Qstar = Qstar))
    }
    if (est.var.iptw) {
      index <- uncensored & intervention.match[, d1]
      g <- cum.g[index, ACnode.index, d1]
      Y <- Scale(Z, 0, 1)[index]
      iptw.estimate <- sum(Y / g) / sum(1 / g) 
      return(list(EZd1=iptw.estimate * diff(range(Z, na.rm=T)) + min(Z, na.rm=T), Qstar=rep(NA, length(Z))))
    }
    
    sparsity.data <- inputs$data[, 1:cur.node]
    sparsity.data[, cur.node] <- Scale(Z, 0, 1)
    temp.nodes <- lapply(nodes, function (x) x[x <= cur.node])
    if (cur.node %in% temp.nodes$L) {
      #if the current node is an L node, make current node a Y node (last node has to be a Y node)
      temp.nodes$L <- setdiff(temp.nodes$L, cur.node)
      temp.nodes$Y <- c(temp.nodes$Y, cur.node)
    }
    stratify <- FALSE
    Qform <- paste(GetDefaultForm(sparsity.data[, 1:cur.node], nodes=temp.nodes, is.Qform=TRUE, stratify=stratify, survivalOutcome=FALSE, showMessage=FALSE), paste0("+ sparityAdj_Z.meanL_", 1:length(temp.nodes$LY)))
    Qform[length(Qform)] <- "IDENTITY"
    
    Z.meanL <- apply(AsMatrix(Z.meanL), 2, LogitScale)
    sparsity.data <- cbind(Z.meanL, sparsity.data)
    names(sparsity.data)[sseq(1, ncol(Z.meanL))] <- paste0("sparityAdj_Z.meanL_", sseq(1, ncol(Z.meanL)))
    temp.nodes <- lapply(temp.nodes, function (x) x + ncol(Z.meanL))
    
    names(Qform) <- names(sparsity.data)[temp.nodes$LY]
    attr(sparsity.data, "called.from.estimate.variance") <- TRUE
    var.tmle <- ltmle(sparsity.data, Anodes=temp.nodes$A, Cnodes=temp.nodes$C, Lnodes=temp.nodes$L, Ynodes=temp.nodes$Y, survivalOutcome=FALSE, Qform=Qform, gform=drop3(prob.A.is.1[, 1:ACnode.index, d1, drop=FALSE]), abar= GetABar(inputs$regimes, d1, temp.nodes$A), gbounds=inputs$gbounds, stratify=stratify, estimate.time=FALSE, deterministic.Q.function=det.q.function, variance.method="ic", observation.weights=observation.weights) 
    
    EZd1 <- var.tmle$estimates["tmle"] * diff(range(Z, na.rm=T)) + min(Z, na.rm=T)
    return(list(EZd1 = EZd1, Qstar=var.tmle$Qstar))
  }
  
  EqualRegimesIndex <- function(dd1, dd2) {
    #index of each observation where regime d1 matches regime d2 up to cur.node
    if (!any(nodes$A <= cur.node)) return(rep(TRUE, n))
    eq <- rowAlls(AsMatrix(inputs$regimes[, which(nodes$A <= cur.node), dd1]) == AsMatrix(inputs$regimes[, which(nodes$A <= cur.node), dd2]))
    eq[is.na(eq)] <- FALSE #regimes can be NA after censoring or death - after death Z is zero, after censoring Z is not used so it doesn't really matter whether eq is TRUE or FALSE, but FALSE can be faster (TRUE implies not static treatment)
    return(eq) 
  }
  
  IsStaticTreatment <- function() {
    #static = for all observations, regime d1 matches regime d2 only when d1=d2
    for (dd1 in regimes.with.positive.weight) {  
      for (dd2 in regimes.with.positive.weight[regimes.with.positive.weight > dd1]) {
        if (any(EqualRegimesIndex(dd1, dd2))) return(FALSE)
      }
    }
    return(TRUE)
  }
  
  num.regimes <- dim(inputs$regimes)[3]
  num.betas <- ncol(combined.summary.measures)
  n <- nrow(inputs$data)
  
  #used in ltmle call below
  if (inputs$survivalOutcome) {
    det.q.function <- function(data, current.node, nodes, called.from.estimate.g) {
      if (!any(nodes$Y < current.node)) return(NULL)
      prev.Y <- data[, nodes$Y[nodes$Y < current.node], drop=F]
      prev.Y[is.na(prev.Y)] <- 0
      is.deterministic <- rowAnys(prev.Y == 1)
      Q.value <- data[is.deterministic, max(nodes$Y)] #this is 0 before scaling but may be nonzero after scaling
      return(list(is.deterministic=is.deterministic, Q.value=Q.value))   
    }
  } else {
    det.q.function <- NULL
  }
  static.treatment <- IsStaticTreatment()
  variance.estimate <- matrix(0, num.betas, num.betas)
  Sigma <- array(dim=c(n, num.regimes, num.regimes))
  if (!is.last.LYnode) Q.data <- inputs$data[alive, 1:cur.node, drop=F]
  #Sigma is not exactly symmetric due to SetA on d1, but could probably save time with approximate symmetry
  for (d1 in regimes.with.positive.weight) {
    if (static.treatment) {
      d2.regimes <- d1 #only need diagonal elements
    } else {
      d2.regimes <- regimes.with.positive.weight #need all elements
    }
    for (d2 in d2.regimes) { 
      if (is.last.LYnode) {
        Sigma[, d1, d2] <- Qstar[, d1] * (1 - Qstar[, d1]) 
      } else {
        if (any(alive)) { #if none alive, skip this to avoid warnings on range
          resid.sq <- (Qstar.kplus1[alive, d1] - Qstar[alive, d1]) * (Qstar.kplus1[alive, d2] - Qstar[alive, d2]) 
          resid.sq.range <- range(resid.sq, na.rm=T)
          if (diff(resid.sq.range) > 0.0001) {
            Q.data[, cur.node] <- (resid.sq - resid.sq.range[1]) / diff(resid.sq.range)
            names(Q.data)[cur.node]  <- "Q.kplus1" #to match with Qform
            m <- ltmle.glm(formula = formula(inputs$Qform[LYnode.index]), family = quasibinomial(), data = Q.data, weights=NULL) 
            Q.newdata <- SetA(data = Q.data, regimes = inputs$regimes[alive, , d1, drop=F], Anodes = nodes$A, cur.node = cur.node)
            SuppressGivenWarnings(Q.resid.sq.pred <- predict(m, newdata = Q.newdata, type = "response"), "prediction from a rank-deficient fit may be misleading")
            Sigma[alive, d1, d2] <- Q.resid.sq.pred * diff(resid.sq.range) + resid.sq.range[1]
          } else {
            resid.sq.value <- min(resid.sq, na.rm = T) #all values are the same, just get one non-NA
            Sigma[alive, d1, d2] <- resid.sq.value
          }
        }
        Sigma[!alive, d1, d2] <- 0
      }
    }
  }
  
  if (est.var.iptw) Z.without.sum.meas.meanL <- Z.meanL <- NA 
  no.V <- length(inputs$baseline.column.names) == 0
  if ((!est.var.iptw && static.treatment) || (est.var.iptw && static.treatment && no.V)) { #est.var.iptw doesn't work with the V code because we don't have var.tmle$Qstar
    for (d1 in regimes.with.positive.weight) {   
      Z.without.sum.meas <- Sigma[, d1, d1] / cum.g[, ACnode.index, d1] * cum.g.unbounded[, ACnode.index, d1] / cum.g[, ACnode.index, d1] * msm.weights[, d1]^2 * observation.weights^2
      if (!est.var.iptw) Z.without.sum.meas.meanL <- 1 / cum.g.meanL[, ACnode.index, d1, ] * cum.g.meanL.unbounded[, ACnode.index, d1, ] / cum.g.meanL[, ACnode.index, d1, ] * msm.weights[, d1]^2 * observation.weights^2
      var.tmle <- TmleOfVariance(Z.without.sum.meas, Z.without.sum.meas.meanL)
      if (no.V) {
        variance.estimate <- variance.estimate + (combined.summary.measures[1, , d1] %*% t(combined.summary.measures[1, , d1])) * var.tmle$EZd1
      } else {
        #has V (so combined.summary.measures varies)
        baseline.msm <- paste("Qstar ~", paste(inputs$baseline.column.names, collapse=" + "), "+", paste0("I(", inputs$baseline.column.names, "^2)", collapse=" + "))
        V.data <- data.frame(Qstar=var.tmle$Qstar, inputs$data[, inputs$baseline.column.names, drop=FALSE])
        m <- ltmle.glm(formula(baseline.msm), family = quasibinomial(), data=V.data, weights = NULL)
        SuppressGivenWarnings(pred.Qstar <- predict(m, type = "response", newdata = V.data) * diff(range(Z.without.sum.meas, na.rm=T)) + min(Z.without.sum.meas, na.rm=T), "prediction from a rank-deficient fit may be misleading")  #n x 1
        variance.estimate.sum <- crossprod(combined.summary.measures[, , d1], combined.summary.measures[, , d1] * pred.Qstar) #equivalent to:  for (i in 1:n) variance.estimate.sum <- variance.estimate.sum + (combined.summary.measures[i, , d1] %*% t(combined.summary.measures[i, , d1])) * pred.Qstar[i]
        variance.estimate <- variance.estimate + variance.estimate.sum / n
      }
    }
  } else {
    for (beta.index2 in 1:num.betas) {
      for (d1 in regimes.with.positive.weight) {          
        Z.base <- rep(0, n)  #Z without h1(d1, V, beta.index1)
        if (!est.var.iptw) Z.base.meanL <- matrix(0, n, dim(cum.g.meanL)[4])
        for (d2 in regimes.with.positive.weight) {
          equal.regimes.index <- EqualRegimesIndex(d1, d2) #index of each observation where regime d1 matches regime d2 
          h1 <- combined.summary.measures[, beta.index2, d2] * msm.weights[, d2]
          Z.base[equal.regimes.index] <- Z.base[equal.regimes.index] + h1[equal.regimes.index] * Sigma[equal.regimes.index, d1, d2] / cum.g[equal.regimes.index, ACnode.index, d1] * observation.weights[equal.regimes.index] #this is equivalent to using cum.g.unbounded in the denominator and multiplying by phi=cum.g.unbounded/cum.g.bounded  
          if (!est.var.iptw) Z.base.meanL[equal.regimes.index, ] <- Z.base.meanL[equal.regimes.index, ] + h1[equal.regimes.index] * 1 / cum.g.meanL[equal.regimes.index, ACnode.index, d1, ] * observation.weights[equal.regimes.index] #recycles
        }        
        for (beta.index1 in 1:num.betas) {  
          if (beta.index1 >= beta.index2) {
            Z <- combined.summary.measures[, beta.index1, d1] * msm.weights[, d1] * cum.g.unbounded[, ACnode.index, d1] / cum.g[, ACnode.index, d1] * observation.weights * Z.base 
            if (!est.var.iptw) Z.meanL <- combined.summary.measures[, beta.index1, d1] * msm.weights[, d1] * cum.g.meanL.unbounded[, ACnode.index, d1, ] / cum.g.meanL[, ACnode.index, d1, ] * observation.weights * Z.base.meanL
            var.tmle <- TmleOfVariance(Z, Z.meanL)
            variance.estimate[beta.index1, beta.index2] <- variance.estimate[beta.index1, beta.index2] + var.tmle$EZd1
          } else {
            variance.estimate[beta.index1, beta.index2] <- variance.estimate[beta.index2, beta.index1] #use symmetry
          }
        }
      } 
    }
  }
  if (max(abs(variance.estimate - t(variance.estimate))) > 1e-5) stop("not symmetric")  
  if (any(eigen(variance.estimate, only.values = TRUE)$values < -1e-8)) {
    #this happens very rarely
    variance.estimate <- MakePosDef(variance.estimate)
  }
  return(variance.estimate)
}

MakePosDef <- function(variance.estimate) {
  orig.variance.estimate <- variance.estimate
  try.result <- try({
    near.pd <- Matrix::nearPD(variance.estimate) #may cause an error if variance.estimate is negative definite
    variance.estimate <- as.matrix(near.pd$mat)
  }, silent = TRUE)
  if (inherits(try.result, "try-error") || !near.pd$converged || any(abs(orig.variance.estimate - variance.estimate) > 0.001 & (abs(orig.variance.estimate - variance.estimate) / orig.variance.estimate) > 0.1)) {
    warning("Covariance matrix from EstimateVariance not positive definite, unable to compute standard errors. You may want to try variance.method='ic'.") 
    variance.estimate <- matrix(nrow=nrow(variance.estimate), ncol=ncol(variance.estimate))
  }
  return(variance.estimate)
}

CalcGUnboundedToBoundedRatio <- function(g.list, nodes, final.Ynodes) {
  CalcForFinalYNode <- function(num.AC.nodes) {
    if (num.AC.nodes == 0) return(1)
    if (! anyNA(g.list$cum.g)) return(AsMatrix(g.list$cum.g.unbounded[, num.AC.nodes, ] / g.list$cum.g[, num.AC.nodes, ]))
    #cum.g is NA after censoring - for censored observations use cum.g.meanL  
    #If censored at node j, set all nodes > j to meanl. 
    #[,,k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L (na.rm=T)
    g.ratio1 <- matrix(NA, n, num.regimes)
    for (i in 1:num.regimes) {
      g.ratio.temp <- cbind(g.list$cum.g.meanL.unbounded[, num.AC.nodes, i, ] / g.list$cum.g.meanL[, num.AC.nodes, i, ], g.list$cum.g.unbounded[, num.AC.nodes, i] / g.list$cum.g[, num.AC.nodes, i])
      index <- max.col(!is.na(g.ratio.temp), "last")
      g.ratio1[, i] <- g.ratio.temp[sub2ind(1:n, col = index, num.rows = n)]
    }
    return(g.ratio1)
  }
  #calc for each final.ynode - num.AC.nodes varies
  n <- dim(g.list$cum.g)[1]
  num.regimes <- dim(g.list$cum.g)[3]
  num.final.Ynodes <- length(final.Ynodes)
  g.ratio <- array(dim=c(n, num.regimes, num.final.Ynodes))
  for (j in 1:num.final.Ynodes) {
    num.AC.nodes <- sum(nodes$AC < final.Ynodes[j])
    g.ratio[, , j] <- CalcForFinalYNode(num.AC.nodes)
  }
  return(g.ratio)
}

# remove any nodes after final.Ynode
SubsetNodes <- function(nodes, final.Ynode) {
  return(lapply(nodes, function (x) x[x <= final.Ynode]))
}

# Fit the MSM
FitPooledMSM <- function(working.msm, Qstar, combined.summary.measures, msm.weights) {
  #Qstar: n x num.regimes x num.final.Ynodes
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
  #msm.weights: n x num.regimes x num.final.Ynodes
  
  n <- dim(Qstar)[1]
  num.regimes <- dim(Qstar)[2]
  num.final.Ynodes <- dim(Qstar)[3]
  num.summary.measures <- dim(combined.summary.measures)[2]
  
  X <- apply(combined.summary.measures, 2, rbind) 
  Y <- as.vector(Qstar)
  weight.vec <- as.vector(msm.weights)
  data.pooled <- data.frame(Y, X)
  positive.weight <- weight.vec > 0 #speedglm crashes if Y is NA even if weight is 0
  
  m <- ltmle.glm(formula(working.msm), data=data.pooled[positive.weight, ], family=quasibinomial(), weights=weight.vec[positive.weight])
  SuppressGivenWarnings(m.beta <- predict(m, newdata=data.pooled, type="response"), "prediction from a rank-deficient fit may be misleading")
  dim(m.beta) <- dim(Qstar)
  return(list(m=m, m.beta=m.beta))
}

#convert from recordIC (num.records x num.measures) to householdIC (num.households x num.measures)
HouseholdIC <- function(recordIC, id) {
  if (is.null(id)) return(recordIC)
  householdIC <- as.matrix(aggregate(recordIC, list(id=id), sum)[, -1, drop=FALSE])
  num.records <- nrow(recordIC)
  num.households <- nrow(householdIC)
  householdIC <- householdIC * num.households / num.records
  return(householdIC)
}

#final step in calculating TMLE influence curve
FinalizeIC <- function(IC, combined.summary.measures, Qstar, m.beta, msm.weights, observation.weights, id) {
  #mBeta, Qstar: n x num.regimes x num.final.Ynodes
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
  
  #summary.measures: num.regimes x num.summary.measures x num.final.Ynodes
  #msm.weights: n x num.regimes x num.final.Ynodes
  
  num.betas <- ncol(IC)
  n <- nrow(Qstar)
  num.regimes <- ncol(Qstar)
  num.final.Ynodes <- dim(Qstar)[3]
  
  stopifnot(num.betas == ncol(combined.summary.measures))
  
  finalIC <- matrix(0, nrow=n, ncol=num.betas)
  for (j in 1:num.final.Ynodes) {
    for (i in 1:num.regimes) {
      if (any(msm.weights[, i, j] > 0)) {
        m1 <- matrix(Qstar[, i, j] - m.beta[, i, j], ncol=1)   #n x 1
        for (k in 1:num.betas) {
          m2 <- combined.summary.measures[, k, i, j] # n x 1
          finalIC[, k] <- finalIC[, k] + msm.weights[, i, j] * observation.weights * (m1 * m2) 
        }
      }
    }  
  }
  IC <- IC + finalIC
  # return(IC)
  return(HouseholdIC(IC, id))
}

# Normalize the influence curve matrix
NormalizeIC <- function(IC, combined.summary.measures, m.beta, msm.weights, observation.weights, g.ratio) {    
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
  #g.ratio = g.unbounded / g.bounded : n x num.regimes x num.final.Ynodes
  n <- dim(combined.summary.measures)[1]
  num.betas <- dim(combined.summary.measures)[2]
  num.regimes <- dim(combined.summary.measures)[3]
  num.final.Ynodes <- dim(combined.summary.measures)[4]
  
  if (is.null(g.ratio)) {
    g.ratio <- array(1, dim=c(n, num.regimes, num.final.Ynodes)) #if variance.method="ic", g.ratio should be NULL 
  }
  C <- array(0, dim=c(num.betas, num.betas))
  for (j in 1:num.final.Ynodes) { 
    for (i in 1:num.regimes) {
      tempC <- crossprod(combined.summary.measures[, , i, j] * g.ratio[, i, j], combined.summary.measures[, , i, j] * g.ratio[, i, j] * msm.weights[, i, j] * m.beta[, i, j] * (1 - m.beta[, i, j]) * observation.weights)
      if (anyNA(tempC)) stop("NA in tempC")
      C <- C + tempC
    }
  }
  C <- C / n
  
  # For easier reading, this is equivalent (but much slower):
  # C2 <- array(0, dim=c(num.betas, num.betas, n))
  # for (j in 1:num.final.Ynodes) {
  #   for (i in 1:num.regimes) {
  #     positive.msm.weights <- which(msm.weights[, i, j] > 0)
  #     tempC <- array(0, dim=c(num.betas, num.betas, n))
  #     for (k in positive.msm.weights) {
  #       m.beta.temp <- m.beta[k, i, j]  
  #       h <- matrix(combined.summary.measures[k, , i, j], ncol=1) * msm.weights[k, i, j] * g.ratio[k, i, j]
  #       tempC[, , k] <- h %*% t(h) * m.beta.temp * (1 - m.beta.temp) * observation.weights[k] / msm.weights[k, i, j]
  #     }
  #     if (anyNA(tempC)) stop("NA in tempC")
  #     C2 <- C2 + tempC
  #   }
  # }
  # C2 <- apply(C2, c(1, 2), mean)
  
  if (rcond(C) < 1e-12) {
    C <- matrix(NA, nrow=num.betas, ncol=num.betas)
    warning("rcond(C) near 0, standard errors not available")
  } 
  return(C)
}

# Get a single regime from the regimes array
GetABar <- function(regimes, regime.index, Anodes) {
  abar <- AsMatrix(regimes[, seq_along(Anodes), regime.index]) #if there's only 1 Anode, make sure abar comes back as a matrix
  return(abar)
}

# Targeting step - update the initial fit of Q using clever covariates
UpdateQ <- function(Qstar.kplus1, logitQ, combined.summary.measures, cum.g, working.msm, uncensored, intervention.match, is.deterministic, msm.weights, gcomp, observation.weights) { 
  #logitQ, Qstar.kplus1: n x num.regimes
  #cum.g: n x num.regimes (already indexed for this node)
  #uncensored: n x 1
  #intervention.match: n x num.regimes
  #is.deterministic: n x 1
  #summary.measures: num.regimes x num.summary.measures
  #baseline.covariates: names/indicies: num.baseline.covariates x 1
  #msm.weights: n x num.regimes
  #stacked.summary.measures: (n*num.regimes) x num.measures
  #combined.summary.measures: n x num.measures x num.regimes   (num.measures=num.summary.measures + num.baseline.covariates)
  #h.g.ratio: n x num.regimes x num.measures
  #observation.weights: n x 1
  n <- nrow(logitQ)
  num.regimes <- ncol(logitQ)
  off <- as.vector(logitQ)
  Y <- as.vector(Qstar.kplus1)
  
  if (length(cum.g) == 0) cum.g <- 1 #if there are no A/C nodes 
  stacked.summary.measures <- apply(combined.summary.measures, 2, rbind)
  
  subs.vec <- uncensored & !is.deterministic & as.vector(intervention.match) #recycles uncensored and is.deterministic
  weight.vec <- numeric(n * num.regimes)
  weight.vec[subs.vec] <- (observation.weights * as.vector(msm.weights) / as.vector(cum.g))[subs.vec] #recycles observation.weights (subsetting avoids problems with NA in cum.g) 
  if (anyNA(weight.vec)) stop("NA in weight.vec")
  
  f <- as.formula(paste(working.msm, "+ offset(off)"))
  data.temp <- data.frame(Y, stacked.summary.measures, off)
  if (gcomp) {
    Qstar <- plogis(logitQ)
    m <- "no Qstar fit because gcomp=TRUE (so no updating step)"
  } else {
    if (any(weight.vec > 0)) {
      m <- ltmle.glm(f, data=data.temp[weight.vec > 0, ], family=quasibinomial(), weights=as.vector(scale(weight.vec[weight.vec > 0], center=FALSE))) #this should include the indicators; note: scale converts to matrix, as.vector needed to convert back 
      Qstar <- matrix(predict(m, newdata=data.temp, type="response"), nrow=nrow(logitQ))
    } else {
      Qstar <- plogis(logitQ)
      m <- "no Qstar fit because no subjects alive, uncensored, following intervention"
    }
    
  }
  indicator <- matrix(uncensored * observation.weights, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) * matrix(intervention.match, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) #I(A=rule and uncensored) * observation.weights
  h.g.ratio <- stacked.summary.measures / matrix(cum.g, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) * indicator # I() * h * observation.weights / g
  dim(h.g.ratio) <- c(n, num.regimes, ncol(h.g.ratio))
  for (i in 1:num.regimes) {
    h.g.ratio[, i, ] <- h.g.ratio[, i, ] * msm.weights[, i] #recycles msm.weights
    weight.zero.index <- msm.weights[, i] == 0
    h.g.ratio[weight.zero.index, i, ] <- 0  #cum.g is 0 so X is NA so h.g.ratio is NA when weight is 0
  }
  cum.g.used <- weight.vec > 0 & msm.weights > 0
  dim(cum.g.used) <- c(n, num.regimes)
  return(list(Qstar=Qstar, h.g.ratio=h.g.ratio, X=stacked.summary.measures, off=off, fit=m, cum.g.used=cum.g.used)) 
}

# Sometimes GLM doesn't converge and the updating step of TMLE doesn't solve the score equation (sum of TMLE influence curve not equal to zero). This function attempts to solve the score equation directly
FixScoreEquation <- function(Qstar.kplus1, h.g.ratio, uncensored, intervention.match, is.deterministic, deterministic.Q, off, X, regimes.with.positive.weight) {
  CalcScore <- function(e) {
    Qstar <- QstarFromE(e)
    ICtemp <- CalcIC(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
    return(sum(colSums(ICtemp) ^ 2)) #each column has to add to zero
  }
  
  QstarFromE <- function(e) {
    Qstar <- plogis(off + X %*% e) #X: n x (num.summary.measures + num.baseline.covariates) (which should be num.beta);  e: num.beta x 1 
    dim(Qstar) <- dim(Qstar.kplus1)
    Qstar[is.deterministic] <- deterministic.Q[is.deterministic] #matrix indexing
    return(Qstar)
  }
  
  FindMin <- function(minimizer) {
    num.tries <- 20 
    init.e <- numeric(num.betas) #first try an initial estimate of epsilon=0
    for (i in 1:num.tries) {
      m <- nlminb(start=init.e, objective=CalcScore, control=list(abs.tol=max.objective, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
      e <- m$par
      obj.val <- m$objective
      if (obj.val < max.objective) {
        m$ltmle.msg <- "updating step using glm failed to solve score equation; solved using nlminb"
        return(list(e=e, solved=TRUE, m=m))
      }
      init.e <- rnorm(num.betas) #if the first try didn't work, try a random initial estimate of epsilon 
    }
    return(list(e=numeric(num.betas), solved=FALSE, m="score equation not solved!")) #nocov - return Q (not updated)
  }
  max.objective <- 0.0001 ^ 2
  num.betas <- ncol(X)
  for (offset.lbound in c(1e-8, 0.0001, 0.001, 0.01)) {
    off <- Bound(off, qlogis(c(offset.lbound, 1-offset.lbound)))
    l <- FindMin("nlminb")
    if (l$solved) break
  }
  if (! l$solved) stop("minimizer failed")
  Qstar <- QstarFromE(l$e)
  return(list(Qstar=Qstar, fit=l$m))
}

# Estimate how long it will take to run ltmleMSM
EstimateTime <- function(inputs) {
  sample.size <- 50
  if (nrow(inputs$data) < sample.size) {
    message(paste("Timing estimate unavailable when n <", sample.size))
    return(NULL)
  }
  sample.index <- sample(nrow(inputs$data), size=sample.size)
  small.inputs <- inputs
  small.inputs$data <- small.inputs$data[sample.index, ]
  small.inputs$regimes <- small.inputs$regimes[sample.index, , , drop=F]
  small.inputs$observation.weights <- small.inputs$observation.weights[sample.index]
  small.inputs$uncensored <- small.inputs$uncensored[sample.index, , drop=F]
  small.inputs$intervention.match <- small.inputs$intervention.match[sample.index, , , drop=F]
  small.inputs$combined.summary.measures <- small.inputs$combined.summary.measures[sample.index, , , , drop=F]
  if (!is.null(small.inputs$id)) small.inputs$id <- small.inputs$id[sample.index]
  if (is.numeric(inputs$gform)) small.inputs$gform <- small.inputs$gform[sample.index, , , drop=F]
  if (length(dim(inputs$msm.weights)) == 3) small.inputs$msm.weights <- small.inputs$msm.weights[sample.index, , , drop=F]
  start.time <- Sys.time()
  try.result <- suppressWarnings(try(MainCalcs(small.inputs), silent=TRUE))
  if (inherits(try.result, "try-error")) {
    message("Timing estimate unavailable")
  } else {
    elapsed.time <- Sys.time() - start.time 
    est.time1 <- round(sqrt(as.double(elapsed.time, units="mins") * nrow(inputs$data) / sample.size), digits=0)
    est.time2 <- round(as.double(elapsed.time, units="mins") * nrow(inputs$data) / sample.size, digits=0)
    if (est.time2 == 0) {
      est.time.str <- "< 1 minute"
    } else if (est.time2 == 1) {
      est.time.str <- "1 minute" 
    } else {
      est.time.str <- paste(est.time1, "to", est.time2, "minutes") 
    }
    message("Estimate of time to completion: ", est.time.str)
  }
  return(NULL)
}

#' Get standard error, p-value, and confidence interval for one ltmle object 
#' Summarizing results from Longitudinal Targeted Maximum Likelihood Estimation
#' (ltmle)
#' 
#' These functions are methods for class \code{ltmle} or \code{summary.ltmle}
#' objects.
#' 
#' \code{summary.ltmle} returns the parameter value of the estimator, the
#' estimated variance, a 95 percent confidence interval, and a p-value.
#' 
#' \code{summary.ltmleEffectMeasures} returns the additive treatment effect for
#' each of the two objects in the \code{abar} list passed to \code{ltmle}.
#' Relative risk, and odds ratio are also returned, along with the variance,
#' confidence interval, and p-value for each.
#' 
#' \code{summary.ltmleMSM} returns a matrix of MSM parameter estimates.
#' 
#' @aliases summary.ltmle print.ltmle print.summary.ltmle summary.ltmleMSM
#' print.ltmleMSM print.summary.ltmleMSM summary.ltmleEffectMeasures
#' print.ltmleEffectMeasures print.summary.ltmleEffectMeasures
#' @param object an object of class "\code{ltmle}" or "\code{ltmleMSM}" or
#' "\code{ltmleEffectMeasures}", usually a result of a call to
#' \code{\link{ltmle}} or \code{\link{ltmleMSM}}.
#' @param x an object of class "\code{summary.ltmle}" or
#' "\code{summary.ltmleMSM}" or "\code{ltmleEffectMeasures}", usually a result
#' of a call to \code{\link{summary.ltmle}} or \code{\link{summary.ltmleMSM}}.
#' @param estimator character; one of "tmle", "iptw", "gcomp". The estimator
#' for which to get effect measures. "tmle" is valid iff the original
#' ltmle/ltmleMSM call used gcomp=FALSE. "gcomp" is valid iff the original
#' ltmle/ltmleMSM call used gcomp=TRUE
#' @param digits the number of significant digits to use when printing.
#' @param signif.stars logical. If \code{TRUE}, significance stars are printed
#' for each coefficient.
#' @param \dots further arguments passed to or from other methods.
#' @return \code{summary.ltmle} returns an object of class
#' "\code{summary.ltmle}", a list with components \item{treatment}{a list with
#' components summarizing the estimate of \code{object} \itemize{ \item
#' \code{estimate} - the parameter estimate of \eqn{E[Y_d]} \item
#' \code{std.dev} - estimated standard deviation of parameter \item
#' \code{p.value} - two-sided p-value \item \code{CI} - vector of length 2 with
#' 95 percent confidence interval } }
#' 
#' \item{call}{the matched call to \code{ltmle} for \code{object}}
#' \item{estimator}{the \code{estimator} input argument}
#' \item{variance.estimate.ratio}{ratio of the TMLE based variance estimate to
#' the influence curve based variance estimate}
#' 
#' \code{summary.ltmleEffectMeasures} returns an object of class
#' "\code{summary.ltmleEffectMeasures}", a list with same components as
#' \code{summary.ltmle} above, but also includes: \item{effect.measures}{a list
#' with components, each with the same components as \code{treatment} in
#' \code{summary.ltmle} above \itemize{ \item \code{treatment} - corresponds to
#' the first in the list \code{abar} (or \code{rule}) passed to \code{ltmle}
#' \item \code{control} - corresponds to the second in the list \code{abar} (or
#' \code{rule}) passed to \code{ltmle} \item \code{ATE} - average treatment
#' effect \item \code{RR} - relative risk \item \code{OR} - odds ratio } }
#' 
#' \code{summary.ltmleMSM} returns an object of class
#' "\code{summary.ltmleMSM}", a matrix with rows for each MSM parameter and
#' columns for the point estimate, standard error, 2.5percent confidence
#' interval, 97.5percent confidence interval, and p-value.
#' @seealso \code{\link{ltmle}}, \code{\link{summary}}
#' @examples
#' 
#' rexpit <- function(x) rbinom(n = length(x), size = 1, prob = plogis(x))
#' 
#' # Compare the expected outcomes under two counterfactual plans: Treatment plan:
#' # set A1 to 1 if W > 0, set A2 to 1 if W > 1.5, always set A3 to 1 Control plan:
#' # always set A1, A2, and A3 to 0
#' W <- rnorm(1000)
#' A1 <- rexpit(W)
#' A2 <- rexpit(W + 2 * A1)
#' A3 <- rexpit(2 * A1 - A2)
#' Y <- rexpit(W - A1 + 0.5 * A2 + 2 * A3)
#' data <- data.frame(W, A1, A2, A3, Y)
#' treatment <- cbind(W > 0, W > 1.5, 1)
#' control <- matrix(0, nrow = 1000, ncol = 3)
#' result <- ltmle(data, Anodes = c("A1", "A2", "A3"), Ynodes = "Y", abar = list(treatment, 
#'     control))
#' print(summary(result))
#' 
#' ## For examples of summary.ltmle and summary.ltmleMSM, see example(ltmle)
#' 
#' @export 
summary.ltmle <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  if ("control.object" %in% names(list(...))) stop("The control.object parameter has been deprecated. To obtain additive treatment effect, risk ratio, and relative risk, call ltmle with abar=list(treatment, control). See ?ltmle and ?summary.ltmleEffectMeasures.")
  if (! estimator[1] %in% c("tmle", "iptw", "gcomp")) stop("estimator should be one of: tmle, iptw, gcomp. If you are trying to use control.object, the control.object parameter has been deprecated. To obtain additive treatment effect, risk ratio, and relative risk, call ltmle with abar=list(treatment, control). See ?ltmle and ?summary.ltmleEffectMeasures.")
  if (estimator == "tmle" && object$gcomp) stop("estimator 'tmle' is not available because ltmleMSM was called with gcomp=TRUE")
  if (estimator == "gcomp" && !object$gcomp) stop("estimator 'gcomp' is not available because ltmleMSM was called with gcomp=FALSE")
  
  IC.variance <- var(object$IC[[estimator]])
  if (estimator=="tmle" && !is.null(object$variance.estimate)) {
    v <- max(IC.variance, object$variance.estimate) 
  } else {
    v <- IC.variance
  }
  variance.estimate.ratio=v/IC.variance
  
  if (object$binaryOutcome) {
    CIBounds <- c(0, 1)
  } else {
    CIBounds <- c(-Inf, Inf)  #could truncate at Yrange, but it's not clear that's right
  }
  treatment <- GetSummary(list(long.name=NULL, est=object$estimates[estimator], gradient=1, log.std.err=FALSE, CIBounds=CIBounds), v, n=length(object$IC[[estimator]]))
  ans <- list(treatment=treatment, call=object$call, estimator=estimator, variance.estimate.ratio=variance.estimate.ratio)
  class(ans) <- "summary.ltmle"
  return(ans)
}

# Get standard errors, p-values, confidence intervals for an ltmleEffectMeasures object: treatment EYd, control EYd, additive effect, relative risk, odds ratio
#' @rdname summary.ltmle
#' @export 
summary.ltmleEffectMeasures <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  info <- GetSummaryLtmleMSMInfo(object, estimator)
  beta <- info$estimate
  IC <- info$IC
  y0 <- plogis(beta[1])
  y1 <- plogis(beta[1] + beta[2])
  names(y0) <- names(y1) <- NULL
  eff.list <- list(
    treatment=list(long.name="Treatment Estimate", est=y1, gradient=c(y1*(1-y1), y1*(1-y1)), log.std.err=FALSE, CIBounds=0:1),  
    control=list(long.name="Control Estimate", est=y0, gradient=c(y0*(1-y0), 0), log.std.err=FALSE, CIBounds=0:1), 
    ATE=list(long.name="Additive Treatment Effect", est=y1 - y0, gradient=c(y1*(1-y1) - y0*(1-y0), y1*(1-y1)), log.std.err=FALSE, CIBounds=c(-1, 1)), 
    RR=list(long.name="Relative Risk", est=y1 / y0, gradient=c(y0-y1, 1-y1), log.std.err=TRUE, CIBounds=c(0, Inf)), 
    OR=list(long.name="Odds Ratio", est=exp(beta[2]), gradient=c(0, 1), log.std.err=TRUE, CIBounds=c(0, Inf)))
  if (!object$binaryOutcome) {
    eff.list$RR <- eff.list$OR <- NULL #not valid if non-binary outcome
  }
  n <- nrow(IC)
  
  measures.IC <- lapply(eff.list, GetSummary, var(IC), n)
  if (is.null(object$variance.estimate)) {
    measures.variance.estimate <- NULL #if variance.method="ic"
  } else {
    measures.variance.estimate <- lapply(eff.list, GetSummary, object$variance.estimate, n)  
  }
  measures.max <- measures.IC
  for (i in seq_along(measures.variance.estimate)) {
    std.dev.diff <- measures.variance.estimate[[i]]$std.dev - measures.IC[[i]]$std.dev
    if (!is.na(std.dev.diff) && (std.dev.diff > 0)) { #can be NA if all Y_d are near 0 or 1
      measures.max[[i]] <- measures.variance.estimate[[i]]
    }
  }
  if (object$transformOutcome) {
    #transform back to original scale
    Yrange <- attr(object$transformOutcome, "Yrange")
    measures.max <- lapply(measures.max, function (x) {
      x$estimate <- x$estimate * diff(Yrange)
      x$std.dev <- x$std.dev * diff(Yrange)
      x$CI <- x$CI * diff(Yrange)
      if (x$long.name %in% c("Treatment Estimate", "Control Estimate")) {
        x$estimate <- x$estimate + min(Yrange)
        x$CI <- x$CI + min(Yrange)
      } else {
        stopifnot(x$long.name == "Additive Treatment Effect")
      }
      return(x)
    })
  }
  ans <- list(call=object$call, effect.measures=measures.max, variance.estimate.ratio=info$variance.estimate.ratio, estimator=estimator)
  class(ans) <- "summary.ltmleEffectMeasures"
  return(ans) 
}

# Do some error checking and get basic info about the ltmleMSM object
GetSummaryLtmleMSMInfo <- function(object, estimator) {
  if (! estimator %in% c("tmle", "iptw", "gcomp")) stop("estimator should be one of: tmle, iptw, gcomp")
  if (estimator == "tmle") {
    if (object$gcomp) stop("estimator 'tmle' is not available because ltmleMSM was called with gcomp=TRUE")
    estimate <- object$beta
    IC <- object$IC
  } else if (estimator == "iptw") {
    estimate <- object$beta.iptw
    IC <- object$IC.iptw
  } else if (estimator == "gcomp") {
    if (!object$gcomp) stop("estimator 'gcomp' is not available because ltmleMSM was called with gcomp=FALSE")
    estimate <- object$beta
    IC <- object$IC
  }
  IC.variance <- apply(IC, 2, var)
  if (is.null(object$variance.estimate)) { 
    v <- IC.variance 
  } else {
    v <- pmax(diag(object$variance.estimate), IC.variance)
  }
  variance.estimate.ratio <- v / IC.variance
  return(list(estimate=estimate, IC=IC, variance.estimate.ratio=variance.estimate.ratio, v=v))
}

# Get summary measures for MSM parameters (standard errors, p-values, confidence intervals)
#' @rdname summary.ltmle
#' @export 
summary.ltmleMSM <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  info <- GetSummaryLtmleMSMInfo(object, estimator)
  estimate <- info$estimate
  v <- info$v
  n <- nrow(info$IC)
  std.dev <- sqrt(v/n)
  pval <- 2 * pnorm(-abs(estimate / std.dev))
  CI <- GetCI(estimate, std.dev)  
  cmat <- cbind(estimate, std.dev, CI, pval)
  dimnames(cmat) <- list(names(estimate), c("Estimate", "Std. Error", "CI 2.5%", "CI 97.5%", "p-value"))
  ans <- list(cmat=cmat, estimator=estimator, transformOutcome=object$transformOutcome, variance.estimate.ratio=info$variance.estimate.ratio) 
  class(ans) <- "summary.ltmleMSM"
  return(ans)
}

# Print method for summary.ltmleMSM
#' @rdname summary.ltmle
#' @export 
print.summary.ltmleMSM <- function(x, digits = max(3, getOption("digits") - 3), signif.stars = getOption("show.signif.stars"), ...) {
  cat("Estimator: ", x$estimator, "\n")
  if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
  printCoefmat(x$cmat, digits = digits, signif.stars = signif.stars, 
               na.print = "NA", has.Pvalue=TRUE, ...)
  if (x$transformOutcome) {
    Yrange <- attr(x$transformOutcome, "Yrange")
    cat("NOTE: The MSM is modeling the transformed outcome ( Y -", min(Yrange),
        ")/(", max(Yrange),"-", min(Yrange),")")
  }
  CheckVarianceEstimateRatio(x)
  invisible(x)
}

# Print method for summary.ltmle
#' @rdname summary.ltmle
#' @export 
print.summary.ltmle <- function(x, ...) {
  cat("Estimator: ", x$estimator, "\n")
  if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
  PrintCall(x$call)
  PrintSummary(x$treatment)
  CheckVarianceEstimateRatio(x)
  invisible(x)
}

# Print method for ltmleEffectMeasures
#' @rdname summary.ltmle
#' @export 
print.ltmleEffectMeasures <- function(x, ...) {
  PrintCall(x$call)
  cat("Use summary(...) to get estimates, standard errors, p-values, and confidence intervals for treatment EYd, control EYd, additive effect, relative risk, and odds ratio.\n")
  invisible(x)
}

# Print method for summary.ltmleEffectMeasures
#' @rdname summary.ltmle
#' @export 
print.summary.ltmleEffectMeasures <- function(x, ...) {
  cat("Estimator: ", x$estimator, "\n")
  if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
  PrintCall(x$call)
  lapply(x$effect.measures, PrintSummary)
  CheckVarianceEstimateRatio(x)
  invisible(x)
}

# Print a warning message if the TMLE based variance estimate is much greater than the IC based variance estimate 
CheckVarianceEstimateRatio <- function(summary.obj) {
  if (anyNA(summary.obj$variance.estimate.ratio)) {
    warning("Unable to compute standard errors.") 
    return(NULL)
  }
  if (any(summary.obj$variance.estimate.ratio > 100)) {
    warning.msg <- paste0("max(TMLE based variance estimate / IC based variance estimate) = ", floor(max(summary.obj$variance.estimate.ratio)), ".\nWhen this ratio is greater than 100, both variance estimates are less likely to be accurate.")
    warning(warning.msg)
  }
}

# Print method for ltmleMSM
#' @rdname summary.ltmle
#' @export 
print.ltmleMSM <- function(x, ...) {
  PrintCall(x$call)
  if (x$gcomp) {
    cat("GCOMP Beta Estimates: \n")
  } else {
    cat("TMLE Beta Estimates: \n")
  }
  print(x$beta)
  if (x$transformOutcome) {
    Yrange <- attr(x$transformOutcome, "Yrange")
    cat("NOTE: The MSM is modeling the transformed outcome ( Y -", min(Yrange),
        ")/(", max(Yrange),"-", min(Yrange),")")
  }
  invisible(x)
}

# Print method for ltmle
#' @rdname summary.ltmle
#' @export 
print.ltmle <- function(x, ...) {
  PrintCall(x$call)
  if (x$gcomp) {
    cat("GCOMP Estimate: ", x$estimates["gcomp"], "\n")
  } else {
    cat("TMLE Estimate: ", x$estimates["tmle"], "\n")
  }  
  invisible(x)
}

# Print a call
PrintCall <- function(cl) {
  cat("Call:\n", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}

# Print estimate, standard error, p-value, confidence interval
PrintSummary <- function(x) {
  if (!is.null(x$long.name)) cat(x$long.name, ":\n", sep="")
  cat("   Parameter Estimate: ", signif(x$estimate, 5), "\n")
  if (x$log.std.err) {
    if (x$long.name == "Relative Risk") {
      param.abbrev <- "RR"
    } else if (x$long.name == "Odds Ratio") {
      param.abbrev <- "OR"
    } else {
      stop("unexpected x$long.name") # nocov (should never occur - ignore in code coverage checks)
    }
    cat("  Est Std Err log(", param.abbrev, "):  ", sep="")
  } else {
    cat("    Estimated Std Err:  ")
  }
  cat(signif(x$std.dev, 5), "\n")
  cat("              p-value: ", ifelse(x$pvalue <= 2*10^-16, "<2e-16",signif(x$pvalue, 5)), "\n")
  cat("    95% Conf Interval:",paste("(", signif(x$CI[1], 5), ", ", signif(x$CI[2], 5), ")", sep=""),"\n\n")
  invisible(x)
}

#Calculate estimate, standard deviation, p-value, confidence interval
GetSummary <- function(eff.list, cov.mat, n) {
  estimate <- eff.list$est
  v <- t(eff.list$gradient) %*% cov.mat %*% eff.list$gradient
  stopifnot(length(v) == 1)
  std.dev <- sqrt(v[1, 1] / n)
  
  if (eff.list$log.std.err) {
    pvalue <- 2 * pnorm(-abs(log(estimate) / std.dev))
    CI <- exp(GetCI(log(estimate), std.dev))
  } else {
    pvalue <- 2 * pnorm(-abs(estimate / std.dev))
    CI <- GetCI(estimate, std.dev)
  }
  CI <- Bound(CI, eff.list$CIBounds) 
  return(list(long.name=eff.list$long.name, estimate=estimate, std.dev=std.dev, pvalue=pvalue, CI=CI, log.std.err=eff.list$log.std.err))
}

# Calculate 95% confidence interval
GetCI <- function(estimate, std.dev) {
  x <- qnorm(0.975) * std.dev
  CI <- cbind("2.5%"=estimate - x, "97.5%"=estimate + x)
  return(CI)
}

# Parametric or SuperLeaner estimation of each g-factor
EstimateG <- function(inputs) {
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  nodes <- inputs$all.nodes
  
  g <- cum.g <- cum.g.unbounded <- prob.A.is.1 <- array(NaN, dim=c(n, length(nodes$AC), num.regimes))
  if (inputs$variance.method == "ic") {
    cum.g.meanL <- cum.g.meanL.unbounded <- NULL
  } else {
    g.meanL <- cum.g.meanL <- cum.g.meanL.unbounded <- array(NaN, dim=c(n, length(nodes$AC), num.regimes, length(nodes$LY)-1)) 
  }
  fit <- vector("list", length(nodes$AC))
  names(fit) <- names(inputs$data)[nodes$AC]
  
  if (inputs$variance.method != "ic" && anyNA(inputs$regimes)) {
    regimes.meanL <- inputs$regimes
    for (i in seq_along(nodes$A)) {
      for (regime.index in 1:num.regimes) {
        regimes.meanL[is.na(regimes.meanL[, i, regime.index]), i, regime.index] <- Mode(inputs$regimes[, i, regime.index], na.rm = TRUE)
      }
    }
  } else {
    regimes.meanL <- NULL
  }
  
  for (i in seq_along(nodes$AC)) {
    cur.node <- nodes$AC[i]
    uncensored <- IsUncensored(inputs$uncensored, nodes$C, cur.node)
    deterministic.origdata <- IsDeterministic(inputs$data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=TRUE, inputs$survivalOutcome)$is.deterministic #deterministic due to death or Q.function
    if (is.numeric(inputs$gform)) {
      if (!is.null(inputs$deterministic.g.function)) stop("deterministic.g.function is not compatible with numeric gform")
      prob.A.is.1[, i, ] <- inputs$gform[, i, ]
      g.est <- list(is.deterministic = deterministic.origdata) 
      fit[[i]] <- "no fit due to numeric gform"
    } else {
      form <- inputs$gform[i]
      deterministic.g.list.origdata <- IsDeterministicG(inputs$data, cur.node, inputs$deterministic.g.function, nodes, using.newdata=F) #deterministic due to deterministic.g.function - using original data
      deterministic.g.origdata <- deterministic.g.list.origdata$is.deterministic
      if (inputs$stratify) {
        intervention.match <- InterventionMatch(inputs$intervention.match, nodes$A, cur.node=nodes$AC[i])
        subs <- uncensored & intervention.match & !deterministic.origdata & !deterministic.g.origdata
      } else {
        subs <- uncensored & !deterministic.origdata & !deterministic.g.origdata
      }
      g.est <- Estimate(inputs, form=form, Qstar.kplus1=NULL, subs=subs, family=quasibinomial(), type="response", nodes=nodes, called.from.estimate.g=TRUE, calc.meanL=inputs$variance.method != "ic", cur.node=cur.node, regimes.meanL=regimes.meanL, regimes.with.positive.weight=1:num.regimes) #assume all regimes have positive weight for some final.Ynode 
      prob.A.is.1[, i, ] <- g.est$predicted.values
      fit[[i]] <- g.est$fit
    }
    #prob.A.is.1 is prob(a=1), gmat is prob(a=abar)
    #cur.abar can be NA after censoring/death if treatment is dynamic
    if (cur.node %in% nodes$A) {
      cur.abar <- AsMatrix(inputs$regimes[, nodes$A == cur.node, ])
      if (is.null(regimes.meanL)) {
        cur.abar.meanL <- cur.abar
      } else {
        cur.abar.meanL <- AsMatrix(regimes.meanL[, nodes$A == cur.node, ])
      }
    } else {
      cur.abar <- cur.abar.meanL <- matrix(1, nrow(inputs$data), num.regimes)  #if this is a cnode, abar is always 1 (uncensored)
    }
    g[, i, ] <- CalcG(AsMatrix(prob.A.is.1[, i, ]), cur.abar, g.est$is.deterministic)
    if (inputs$variance.method != "ic") {
      if (is.numeric(inputs$gform)) {
        if (anyNA(g[, i, ])) stop("Error - NA in numeric gform. There may not be NA values in gform (including after censoring if variance.method is 'tmle' or 'iptw'.")
        g.meanL[, i, , ] <- g[, i, ] #recyles
      } else {
        for (j in sseq(1, dim(g.meanL)[4])) {
          g.meanL[, i, , j] <- CalcG(AsMatrix(g.est$prob.A.is.1.meanL[, , j]), cur.abar.meanL, g.est$is.deterministic)
        }
      }
    }
    if (anyNA(g[uncensored, i, ])) stop("Error - NA in g. g should only be NA after censoring. If you passed numeric gform, make sure there are no NA values except after censoring. Otherwise something has gone wrong.")
  }
  
  
  for (regime.index in 1:num.regimes) {
    cum.g.list <- CalcCumG(AsMatrix(g[, , regime.index]), inputs$gbounds)
    cum.g[, , regime.index] <- cum.g.list$bounded
    cum.g.unbounded[, , regime.index] <- cum.g.list$unbounded
    if (inputs$variance.method != "ic") {
      for (j in sseq(1, dim(g.meanL)[4])) {
        cum.g.list <- CalcCumG(AsMatrix(g.meanL[, , regime.index, j]), inputs$gbounds)
        cum.g.meanL[, , regime.index, j] <- cum.g.list$bounded
        cum.g.meanL.unbounded[, , regime.index, j] <- cum.g.list$unbounded
      }
    }
  }
  return(list(cum.g=cum.g, cum.g.unbounded=cum.g.unbounded, cum.g.meanL=cum.g.meanL, fit=ReorderFits(fit), prob.A.is.1=prob.A.is.1, cum.g.meanL.unbounded=cum.g.meanL.unbounded))
}

CalcG <- function(prob.A.is.1, cur.abar, deterministic.newdata) {
  g <- matrix(NA_real_, nrow(prob.A.is.1), ncol(prob.A.is.1))
  g[!is.na(cur.abar) & cur.abar == 1] <- prob.A.is.1[!is.na(cur.abar) & cur.abar == 1] #matrix indexing
  g[!is.na(cur.abar) & cur.abar == 0] <- 1 - prob.A.is.1[!is.na(cur.abar) & cur.abar == 0] #matrix indexing    
  g[deterministic.newdata] <- 1  #a=abar deterministically after death or other deterministic Q (matrix indexing)
  return(g)
}

# Truncate values within supplied bounds
Bound <- function(x, bounds) {
  stopifnot(length(bounds) == 2 && !anyNA(bounds))
  x[x < min(bounds)] <- min(bounds)
  x[x > max(bounds)] <- max(bounds)
  return(x)
}

# Change from list (num.nodes) of list (num.regimes) to list (num.regimes) of list (num.nodes)
ReorderFits <- function(l1) {
  if (length(l1) == 0) l1 <- list("no fit due to no A/C nodes")
  num.regimes <- length(l1[[1]])
  num.nodes <- length(l1)
  
  l2 <- vector("list", num.regimes)
  for (i in 1:num.regimes) {
    l2[[i]] <- vector("list", num.nodes)
    names(l2[[i]]) <- names(l1)
    for (j in 1:num.nodes) {
      l2[[i]][[j]] <- l1[[j]][[i]]
    }
  }
  return(l2)
}
# Convert named nodes to indicies of nodes
NodeToIndex <- function(data, node) {
 # if (! is.data.frame(data)) stop("data must be a data frame")
  if (is.numeric(node) || is.null(node)) return(node)
  if (! is.character(node)) stop("nodes must be numeric, character, or NULL")
  index <- match(node, names(data))
  if (anyNA(index)) {
    stop(paste("\nnamed node(s) not found:", node[is.na(index)]))
  }
  return(index)
}

 # check if glm should be run instead of SuperLearner
is.glm <- function(SL.library) {
  is.equal(SL.library, "glm", check.attributes = FALSE)
}

# Run GLM or SuperLearner
Estimate <- function(inputs, form, subs, family, type, nodes, Qstar.kplus1, cur.node, calc.meanL, called.from.estimate.g, regimes.meanL, regimes.with.positive.weight) {
  FitAndPredict <- function() {
    if (length(Y.subset) < 2) stop("Estimation failed because there are fewer than 2 observations to fit")
    Y.subset.range <- range(Y.subset)
    if (anyNA(Y.subset.range[1])) stop("Internal error - NA in Y during Estimate")
    if (Y.subset.range[1] < -0.0001 || Y.subset.range[2] > 1.0001) stop("Internal error - Y negative or greater than 1 in Estimate")
    if (Y.subset.range[2] - Y.subset.range[1] < 0.0001) { 
      # all equal Y values causes errors in some SL libraries (and is a waste of time anyway)
      Y.value <- Y.subset.range[2]
      m <- list("no estimation occured because all Y values are the same", Y.value=Y.value)
      predicted.values <- ValuesByType(rep(Y.value, nrow(newdata)))
      class(m) <- "no.Y.variation"
    } else {
      if (use.glm) {
        #estimate using GLM
        # cat(" fit: ", form, "\n")
        SuppressGivenWarnings({
          m <- ltmle.glm.fit(y=Y.subset, x=X.subset, family=family, weights=observation.weights.subset, offset=offst, intercept=intercept)
          m$terms <- tf
          predicted.values <- predict(m, newdata=newdata, type=type)
        }, GetWarningsToSuppress())
      } else {
        #estimate using SuperLearner
        newX.list <- GetNewX(newdata)
        SetSeedIfRegressionTesting()
        try.result <- try({
          SuppressGivenWarnings(m <- SuperLearner::SuperLearner(Y=Y.subset, X=X.subset, SL.library=SL.library, verbose=FALSE, family=family, newX=newX.list$newX, obsWeights=observation.weights.subset, id=id.subset, env = environment(SuperLearner::SuperLearner)), c("non-integer #successes in a binomial glm!", "prediction from a rank-deficient fit may be misleading")) 
        })
        if (!inherits(try.result, "try-error") && all(is.na(m$SL.predict))) { #there's a bug in SuperLearner - if a library returns NAs, it gets coef 0 but the final prediction is still all NA; predict(..., onlySL = TRUE) gets around this
          m$SL.predict <- predict(m, newX.list$newX, X.subset, Y.subset, onlySL = TRUE)$pred
        }
        predicted.values <- ProcessSLPrediction(m$SL.predict, newX.list$new.subs, try.result)
      }
    }
    return(list(m = m, predicted.values = predicted.values))
  }
  GetSLStopMsg <- function(Y) {
    ifelse(all(Y %in% c(0, 1, NA)), "", "\n Note that some SuperLeaner libraries crash when called with continuous dependent variables, as in the case of initial Q regressions with continuous Y or subsequent Q regressions even if Y is binary.") 
  }
  
  ProcessSLPrediction <- function(pred, new.subs, try.result) {
    if (inherits(try.result, "try-error")) {
      stop(paste("\n\nError occured during call to SuperLearner:\n", form, GetSLStopMsg(Y.subset), "\n The error reported is:\n", try.result)) 
    }
    if (all(is.na(pred))) {
      stop(paste("\n\n Unexpected error: SuperLearner returned all NAs during regression:\n", form, GetSLStopMsg(Y.subset))) #nocovr
    }
    predicted.values <- rep(NA, nrow(newdata))
    predicted.values[new.subs] <- pred
    if (max(predicted.values, na.rm=T) > 1 || min(predicted.values, na.rm=T) < 0) {
      msg <- paste("SuperLearner returned predicted.values > 1 or < 0: [min, max] = [", min(predicted.values, na.rm=T), ",", max(predicted.values, na.rm=T), "]. Bounding to [0,1]") 
      warning(msg)
      predicted.values <- Bound(predicted.values, bounds=c(0, 1))
    }
    return(ValuesByType(predicted.values))
  }
  
  PredictOnly <- function(newdata1) {
    if (class(m)[1] == "no.Y.variation") return(rep(m$Y.value, nrow(newdata1)))
    if (use.glm) {
      SuppressGivenWarnings(pred <- predict(m, newdata1, type), "prediction from a rank-deficient fit may be misleading")
    } else {
      newX.list <- GetNewX(newdata1)
      pred <- ProcessSLPrediction(predict(m, newX.list$newX, X.subset, Y.subset, onlySL = TRUE)$pred, newX.list$new.subs, try.result=NULL)
    }
    return(pred)
  }
  
  ValuesByType <- function(x) {
    if (type == "link") {
      stopifnot(family$family %in% c("binomial", "quasibinomial"))
      qlogis(Bound(x, bounds=c(0.0001, 0.9999)))
    } else {
      x
    }
  }
  
  GetNewX <- function(newdata1) {
    new.mod.frame <- model.frame(f, data = newdata1, drop.unused.levels = TRUE, na.action = na.pass)
    newX.temp <- model.matrix(terms(f), new.mod.frame)
    if (!use.glm) {
      colnames(newX.temp) <- paste0("Xx.", 1:ncol(newX.temp)) #change to temp colnames to avoid problems in some SL libraries; SL.gam has problems with names like (Intercept) 
    }
    new.subs <- !rowAnyMissings(newX.temp) #remove NA values from newdata - these will output to NA anyway and cause errors in SuperLearner
    newX <- as.data.frame(newX.temp[new.subs, , drop=FALSE])
    if (ncol(X) == 1) { 
      #SuperLearner crashes if there are screening algorithms and only one column - add a constant
      X.subset <<- cbind(X.subset, ltmle.added.constant=1) 
      newX <- cbind(newX, ltmle.added.constant=1)
    }
    return(list(newX=newX, new.subs=new.subs))
  }
  
  PredictProbAMeanL <- function() {
    #Predict the probability that A=1 if L and Y nodes are set to their mean (or median) values
    
    #probAis1.meanL is n x num.LYnodes - 1
    #probAis1.meanL[, k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L
    
    #somewhat inefficient - for W A.1 L.2 A.2 L.3 A.3 Y, does P(A.1=1) setting L.3 to mean and then L.2 and L.3 to mean, but none of these can be used in P(A.1=1) because they're after A.1
    
    #A is already set to abar in data
    probAis1.meanL <- matrix(NaN, nrow(inputs$data), length(nodes$LY) - 1)
    if (ncol(probAis1.meanL) == 0) return(probAis1.meanL)
    all.LY.nodes <- sort(union(nodes$L, nodes$Y)) #not the same as nodes$LY, which removes blocks
    LYindex <- length(nodes$LY)
    for (i in length(all.LY.nodes):1) { 
      regression.node <- all.LY.nodes[i]
      L <- data[single.subs, regression.node]
      if (is.numeric(L) && !IsBinary(L)) {
        meanL <- mean(L, na.rm = TRUE)
      } else {
        meanL <- Mode(L, na.rm = TRUE) #for factors and binaries
      }
      newdata.meanL[, regression.node] <- meanL
      if (regression.node %in% nodes$LY[1:length(nodes$LY)-1]) {
        LYindex <- LYindex - 1
        probAis1.meanL[, LYindex] <- PredictOnly(newdata = newdata.meanL)
      }  
    }
    if (anyNA(probAis1.meanL[, 1])) stop("NA in probAis1.meanL[, 1]")
    return(probAis1.meanL)
  }
  stopifnot(type %in% c("link", "response"))
  num.regimes <- dim(inputs$regimes)[3]
  if (form == "IDENTITY") {
    stopifnot(is.vector(Qstar.kplus1) == 1)
    predicted.values <- ValuesByType(matrix(Qstar.kplus1, nrow = nrow(inputs$data), ncol = num.regimes))
    fit <- as.list(rep("no fit because form == IDENTITY", num.regimes))
    deterministic.list.olddata <- IsDeterministic(inputs$data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g, inputs$survivalOutcome)
    is.deterministic <- matrix(deterministic.list.olddata$is.deterministic, nrow=nrow(inputs$data), ncol=num.regimes)
    deterministic.Q <- matrix(NA, nrow(inputs$data), num.regimes)
    deterministic.Q[is.deterministic, ] <- deterministic.list.olddata$Q
    return(list(predicted.values=predicted.values, fit=fit, is.deterministic=is.deterministic, deterministic.Q=deterministic.Q, prob.A.is.1.meanL=NULL))
  }
  data <- inputs$data
  if (cur.node %in% nodes$C) {
    data[, cur.node] <- ConvertCensoringNodeToBinary(data[, cur.node]) #convert factors to binaries for compatability with glm and some SL librarie
  }
  f <- as.formula(form)
  SL.library <- if (called.from.estimate.g) inputs$SL.library.g else inputs$SL.library.Q
  use.glm <- (is.glm(SL.library) || length(RhsVars(f)) == 0)  #in a formula like "Y ~ 1", call glm
  
  first.regime <- min(regimes.with.positive.weight)
  if (is.null(Qstar.kplus1)) {
    data.with.Qstar <- data
  } else {
    if (is.matrix(Qstar.kplus1)) {
      data.with.Qstar <- cbind(data, Q.kplus1=Qstar.kplus1[, first.regime])
    } else {
      data.with.Qstar <- cbind(data, Q.kplus1=Qstar.kplus1)
    }
  }
  mod.frame <- model.frame(f, data = data.with.Qstar, drop.unused.levels = TRUE, na.action = na.pass)
  Y <- mod.frame[[1]]
  tf <- terms(f)
  X <- model.matrix(tf, mod.frame)
  offst <- model.offset(mod.frame)
  intercept <- attributes(tf)$intercept
  
  if (!use.glm) {
    if (is.equal(family, quasibinomial())) family <- binomial()
    if (!is.null(offst)) stop("offset in formula not supported with SuperLearner")
    colnames(X) <- paste0("Xx.", 1:ncol(X)) #change to temp colnames to avoid problems in some SL libraries; SL.gam has problems with names like (Intercept) 
    X <- as.data.frame(X)
  }
  
  fit <- vector("list", num.regimes)
  predicted.values <- deterministic.Q <- matrix(NA, nrow(data), num.regimes)
  is.deterministic <- matrix(FALSE, nrow(data), num.regimes)
  Qstar.index <- subs.index <- 1
  fit.and.predict <- NULL
  multiple.subs <- is.matrix(subs)
  multiple.Qstar <- is.matrix(Qstar.kplus1)
  if (calc.meanL) {
    prob.A.is.1.meanL <- array(NaN, dim=c(nrow(inputs$data), num.regimes, length(nodes$LY)-1))
    Anode.index <- which(nodes$A < cur.node)
  } else {
    prob.A.is.1.meanL <- NULL
  }
  for (regime.index in regimes.with.positive.weight) {
    newdata <- SetA(data = data.with.Qstar, regimes = inputs$regimes, Anodes = nodes$A, regime.index = regime.index, cur.node = cur.node)
    if (calc.meanL) {
      if (!is.null(regimes.meanL)) {
        newdata.meanL <- SetA(data = data.with.Qstar, regimes = regimes.meanL, Anodes = nodes$A, regime.index = regime.index, cur.node = cur.node)
      } else {
        newdata.meanL <- newdata
      }
    }
    
    deterministic.list.newdata <- IsDeterministic(newdata, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g, inputs$survivalOutcome)
    if (called.from.estimate.g && !is.null(inputs$deterministic.g.function)) {
      newdata.with.current <- newdata
      stopifnot(cur.node %in% nodes$AC)
      if (cur.node %in% nodes$A) {
        newdata.with.current[, cur.node] <- inputs$regimes[, which(nodes$A == cur.node), regime.index] #set current node to regime for consistency checking in IsDeterministicG
      } else {
        newdata.with.current <- newdata
      }
      deterministic.g.list.newdata <- IsDeterministicG(newdata.with.current, cur.node, inputs$deterministic.g.function, nodes, using.newdata=T) #deterministic g - using data modified so A = abar
    } else {
      deterministic.g.list.newdata <- list(is.deterministic = rep(FALSE, nrow(data)), prob1 = NULL)
    }
    if (regime.index > first.regime && multiple.Qstar) {
      Y <- Qstar.kplus1[, Qstar.index]
    }
    if (regime.index == first.regime || multiple.subs) {
      single.subs <- if (multiple.subs) subs[, subs.index] else subs
      X.subset <- X[single.subs, , drop=FALSE]
      id.subset <- inputs$id[single.subs]
      if (any(single.subs)) X.subset[, colAlls(X.subset == 0)] <- 1 #if there is a column of all zeros, speedglm may crash - replace with column of 1s [any(single.subs) is needed because X.subset==0 is dropped to NULL and causes an error if !any(single.subs)]
      observation.weights.subset <- inputs$observation.weights[single.subs]
      offst.subset <- offst[single.subs]
    }
    if (regime.index == first.regime || multiple.subs || multiple.Qstar) {
      Y.subset <- Y[single.subs]
      if (anyNA(Y.subset)) stop("NA in Estimate")
    }
    if (!all(deterministic.list.newdata$is.deterministic | deterministic.g.list.newdata$is.deterministic)) {
      if (is.null(fit.and.predict) || multiple.Qstar || multiple.subs) {
        fit.and.predict <- FitAndPredict()
        m <- fit.and.predict$m
        predicted.values[, regime.index] <- fit.and.predict$predicted.values
      } else {
        #just predict
        predicted.values[, regime.index] <- PredictOnly(newdata)
      }
      if (calc.meanL) prob.A.is.1.meanL[, regime.index, ] <- PredictProbAMeanL()
    } else {
      m <- "all rows are deterministic, no estimation took place" 
    }
    predicted.values[deterministic.g.list.newdata$is.deterministic, regime.index] <- deterministic.g.list.newdata$prob1 
    if (calc.meanL) prob.A.is.1.meanL[deterministic.g.list.newdata$is.deterministic, regime.index, ] <- deterministic.g.list.newdata$prob1 
    is.deterministic[, regime.index] <- deterministic.list.newdata$is.deterministic
    if (!called.from.estimate.g) deterministic.Q[deterministic.list.newdata$is.deterministic, regime.index] <- deterministic.list.newdata$Q
    if (isTRUE(attr(SL.library, "return.fit", exact = TRUE))) {
      fit[[regime.index]] <- m  
    } else {
      if (use.glm) {
        if (class(m)[1] %in% c("speedglm", "glm")) {
          fit[[regime.index]] <- summary(m)$coefficients #only return matrix to save space (unless return.fit attr)
        } else {
          stopifnot(class(m)[1] %in% c("no.Y.variation", "character"))
          fit[[regime.index]] <- m #there was no fit because all determinsitic or all Y the same
        }
      } else {
        capture.output(print.m <- print(m))
        fit[[regime.index]] <- print.m #only return print matrix to save space (unless return.fit attr)
      }
    } 
    if (multiple.subs) subs.index <- subs.index + 1
    if (multiple.Qstar) Qstar.index <- Qstar.index + 1
  }
  return(list(predicted.values=predicted.values, fit=fit, is.deterministic=is.deterministic, deterministic.Q=deterministic.Q, prob.A.is.1.meanL=prob.A.is.1.meanL))
}

# Calculate bounded cumulative G
CalcCumG <- function(g, gbounds) {
  cum.g <- rowCumprods(g)
  return(list(unbounded=cum.g, bounded=Bound(cum.g, gbounds)))
}

# Calc logical matrix - n x numCnodes = is uncensored up to and including Cnode[i]
CalcUncensoredMatrix <- function(data, Cnodes) {
  uncensored <- matrix(nrow=nrow(data), ncol=length(Cnodes))
  cum.uncensored <- rep(TRUE, nrow(data))
  for (Cnode.index in seq_along(Cnodes)) {
    cum.uncensored <- cum.uncensored & (data[, Cnodes[Cnode.index]] %in% c("uncensored", NA))
    uncensored[, Cnode.index] <- cum.uncensored
  }
  return(uncensored)
}

# Determine which patients are uncensored
#return vector of [numDataRows x 1] I(C=uncensored) from Cnodes[1] to the Cnode just before cur.node
# note: if calling from outside ltmle:::, cur.node needs to be the node index, not a string!
IsUncensored <- function(uncensored.matrix, Cnodes, cur.node) {
  index <- which.max(Cnodes[Cnodes < cur.node])
  if (length(index) == 0) return(rep(TRUE, nrow(uncensored.matrix)))
  return(uncensored.matrix[, index])
}

# Calc logical array - n x num.Anodes x num.regimes = follows regime[j] up to and including Anode[i]
CalcInterventionMatchArray <- function(data, regimes, Anodes) {
  num.regimes <- dim(regimes)[3]
  intervention.match <- array(dim=c(nrow(data), length(Anodes), num.regimes))
  cum.intervention.match <- matrix(TRUE, nrow(data), num.regimes)
  for (Anode.index in seq_along(Anodes)) {
    cum.intervention.match <- cum.intervention.match & ((data[, Anodes[Anode.index]] == regimes[, Anode.index, ]) %in% c(TRUE, NA)) #recycles regimes
    intervention.match[, Anode.index, ] <- cum.intervention.match
  }
  return(intervention.match)
}

# Determine which patients are following specified treatment regime
#return matrix of [numObservations x numRegimes] I(A==abar) from Anodes[1] to the Anode just before cur.node
# note: if calling from outside ltmle:::, cur.node needs to be the node index, not a string!
InterventionMatch <- function(intervention.match.array, Anodes, cur.node) {
  index <- which.max(Anodes[Anodes < cur.node])
  if (length(index) == 0) return(matrix(TRUE, nrow(intervention.match.array), dim(intervention.match.array)[3]))
  return(AsMatrix(intervention.match.array[, index, ]))
}

# Determine which patients have died or have Q set deterministically by user function before cur.node
# return list:
#    is.deterministic: vector of [numObservations x 1] - true if patient is already dead before cur.node or set by deterministic.Q.function
#    Q.value: vector of [which(is.deterministic) x 1] - value of Q
IsDeterministic <- function(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g, survivalOutcome) {
  #set Q.value to 1 if previous y node is 1
  if (survivalOutcome && any(nodes$Y < cur.node)) {
    last.Ynode <- max(nodes$Y[nodes$Y < cur.node])
    is.deterministic <- data[, last.Ynode] %in% TRUE
  } else {
    is.deterministic <- rep(FALSE, nrow(data))
  }
  
  #get Q values from deterministic.Q.function
  default <- list(is.deterministic=is.deterministic, Q.value=1)
  if (is.null(deterministic.Q.function)) return(default)
  #put this in a try-catch?
  det.list <- deterministic.Q.function(data=data, current.node=cur.node, nodes=nodes, called.from.estimate.g=called.from.estimate.g)
  if (is.null(det.list)) return(default)
  if (called.from.estimate.g) {
    #it's ok if Q.value isn't returned if called.from.estimate.g
    if (!is.list(det.list) || !("is.deterministic" %in% names(det.list)) || !(length(det.list) %in% 1:2)) stop("deterministic.Q.function should return a list with names: is.deterministic, Q.value")
  } else {
    if (!is.list(det.list) || !setequal(names(det.list), c("is.deterministic", "Q.value")) || length(det.list) != 2) stop("deterministic.Q.function should return a list with names: is.deterministic, Q.value")
  }
  
  if (! length(det.list$Q.value) %in% c(1, length(which(det.list$is.deterministic)))) stop("the length of the 'Q.value' element of deterministic.Q.function's return argument should be either 1 or length(which(det.list$is.deterministic))")
  
  #check that these observations where Q.value is 1 due to death (previous y is 1) aren't set to anything conflicting by deterministic.Q.function
  Q.value.from.function <- rep(NA, nrow(data))
  Q.value.from.function[det.list$is.deterministic] <- det.list$Q.value
  set.by.function.and.death <- is.deterministic & det.list$is.deterministic
  if (any(Q.value.from.function[set.by.function.and.death] != 1)) {
    stop(paste("inconsistent deterministic Q at node:", names(data)[cur.node])) 
  }
  finalY <- data[, max(nodes$Y)]
  inconsistent.rows <- (det.list$Q.value %in% c(0,1)) & (det.list$Q.value != finalY[det.list$is.deterministic]) & !is.na(finalY[det.list$is.deterministic])
  if (any(inconsistent.rows)) stop(paste("At node:",names(data)[cur.node], "deterministic.Q.function is inconsistent with data - Q.value is either 0 or 1 but this does not match the final Y node value\nCheck data rows:", paste(head(rownames(data)[det.list$is.deterministic][inconsistent.rows]), collapse=" ")))
  
  #return combined values
  Q.value <- rep(NA, nrow(data))
  Q.value[is.deterministic] <- 1
  Q.value[det.list$is.deterministic] <- det.list$Q.value
  is.deterministic <- is.deterministic | det.list$is.deterministic
  Q.value <- Q.value[is.deterministic]
  if (anyNA(is.deterministic) || anyNA(Q.value)) stop("NA in is.deterministic or Q.value")
  return(list(is.deterministic=is.deterministic, Q.value=Q.value))
}

# Determine which patients have an Anode value which is deterministic 
# For example, deterministic.g.function may be used to specify that once a patient starts treatment, they stay on treatment and this should be taken into consideration during estimation of G
IsDeterministicG <- function(data, cur.node, deterministic.g.function, nodes, using.newdata) {
  stopifnot(cur.node %in% nodes$AC)
  default <- list(is.deterministic=rep(FALSE, nrow(data)), prob1=NULL)
  if (is.null(deterministic.g.function)) return(default)
  #put this in a try-catch?
  det.list <- deterministic.g.function(data=data, current.node=cur.node, nodes=nodes)
  if (is.null(det.list)) return(default)
  if (!is.list(det.list) || !setequal(names(det.list), c("is.deterministic", "prob1")) || length(det.list) != 2) stop("deterministic.g.function should return a list with names: is.deterministic, prob1")
  if (! length(det.list$prob1) %in% c(1, length(which(det.list$is.deterministic)))) stop("the length of the 'prob1' element of deterministic.g.function's return argument should be either 1 or length(which(det.list$is.deterministic))")
  
  inconsistent.rows <- (det.list$prob1 %in% c(0,1)) & (det.list$prob1 != data[det.list$is.deterministic, cur.node]) & !is.na(data[det.list$is.deterministic, cur.node])
  if (any(inconsistent.rows)) {
    err.msg <- paste("At node:",names(data)[cur.node], "deterministic.g.function is inconsistent with data - prob1 is either 0 or 1 but this does not match the node value.\nCheck data rows:", paste(head(rownames(data)[det.list$is.deterministic][inconsistent.rows]), collapse=" ")) 
    if (using.newdata) {
      err.msg <- paste(err.msg, "\n This error occured while calling deterministic.g.function on data where Anodes are set to abar.")
      cat("deterministic.g.function is inconsistent with data.\nAfter setting Anodes to abar, the data looks like this:\n")
      print(head(data[det.list$is.deterministic[inconsistent.rows], ]))
    }
    stop(err.msg)
  }
  return(det.list)
}

# Calculate the TMLE influence curve for one node
CalcIC <- function(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight) {
  n <- nrow(Qstar)
  num.regimes <- ncol(Qstar)
  num.betas <- dim(h.g.ratio)[3] #h.g.ratio: n x num.regimes x num.betas
  
  IC <- matrix(0, nrow=n, ncol=num.betas)
  for (i in regimes.with.positive.weight) {
    index <- uncensored & intervention.match[, i]
    if (any(h.g.ratio[index, i, ] != 0)) {
      IC[index, ] <- IC[index, ] + (Qstar.kplus1[index, i] - Qstar[index, i]) * h.g.ratio[index, i, ]
    }
  }
  return(IC)
}

#Set the Anodes of data to regime[, , regime.index] up to cur.node
SetA <- function(data, regimes, Anodes, regime.index, cur.node) {
  Anode.index <- which(Anodes < cur.node)
  for (a in Anode.index) {
    data[, Anodes[a]] <- regimes[, a, regime.index]
  }
  return(data)
}

# Return the left hand side variable of formula f as a character
LhsVars <- function(f) {
  f <- as.formula(f)
  return(as.character(f[[2]]))
}

# Return the right hand side variables of formula f as a character vector
RhsVars <- function(f) {
  f <- as.formula(f)
  return(all.vars(f[[3]]))
}

# Error checking for inputs
CheckInputs <- function(data, nodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, deterministic.Q.function, observation.weights, gcomp, variance.method, id) {
  stopifnot(length(dim(regimes)) == 3)
  num.regimes <- dim(regimes)[3]
  if (!is.glm(GetLibrary(SL.library, "Q")) || !is.glm(GetLibrary(SL.library, "g"))) {
    if (!requireNamespace("SuperLearner")) stop("SuperLearner package is required if SL.library is not NULL or 'glm'")
  } 
  #each set of nodes should be sorted - otherwise causes confusion with gform, Qform, abar
  if (is.unsorted(nodes$A, strictly=TRUE)) stop("Anodes must be in increasing order")
  if (is.unsorted(nodes$C, strictly=TRUE)) stop("Cnodes must be in increasing order")
  if (is.unsorted(nodes$L, strictly=TRUE)) stop("Lnodes must be in increasing order")
  if (is.unsorted(nodes$Y, strictly=TRUE)) stop("Ynodes must be in increasing order")
  if (is.unsorted(final.Ynodes, strictly=TRUE)) stop("final.Ynodes must be in increasing order")
  
  if (length(nodes$L) > 0) {
    if (max(nodes$L) > max(nodes$Y)) stop("Lnodes are not allowed after the final Y node")
  }
  
  all.nodes <- c(nodes$A, nodes$C, nodes$L, nodes$Y)
  if (length(all.nodes) > length(unique(all.nodes))) stop("A node cannot be listed in more than one of Anodes, Cnodes, Lnodes, Ynodes")
  if (is.null(nodes$Y)) stop("Ynodes cannot be null")
  
  if (length(nodes$AC) > 0 && !all(sseq(min(nodes$AC), ncol(data)) %in% all.nodes)) {
    stop("All nodes after the first A/C node must be in A-, C-, L-, or Ynodes")
  }
  
  for (reserved.name in c("observation.weights", "Q.kplus1", "Qstar")) {
    #these might cause conflicts
    if (reserved.name %in% names(data)) stop(paste(reserved.name, "is reserved and may not be used as a column name of data"))
  }
  
  if (length(variance.method) != 1 || !(variance.method %in% c("ic", "tmle", "iptw"))) {
    stop("variance.method must be one of 'ic', 'tmle', 'iptw'") 
  }
  
  #If gform is NULL, it will be set by GetDefaultForm; no need to check here
  if (length(gform) > 0) {
    if (is.character(gform)) {
      if (length(gform) != length(nodes$AC)) stop("length(gform) != length(c(Anodes, Cnodes))")
      for (i in 1:length(gform)) {
        if (LhsVars(gform[i]) != names(data)[nodes$AC[i]]) {
          stop("The LHS of gform[", i, "] should be the name of the ", i, "th A or C node")
        }
        parents <- if(nodes$AC[i] > 1) {
          names(data)[1:(nodes$AC[i]-1)]
        } else {
          NULL
        }
        if (!all(RhsVars(gform[i]) %in% parents)) {
          stop("Some nodes in gform[", i, "] are not parents of ", LhsVars(gform[i]))
        }
        if (any(RhsVars(gform[i]) %in% names(data)[nodes$C])) stop("Cnodes should not be used as RHS variables in gform (regressions are only run on uncensored observations so including a Cnode has no effect and slows down regressions)")
      }
    } else {
      if (! is.numeric(gform)) stop("gform should be a character vector or numeric")
      if (nrow(gform) != nrow(data)) stop("if gform is numeric, it should have the same number of rows as data")
      if (ncol(gform) != length(nodes$AC)) stop("if gform is numeric, it should have the same number of columns as length(c(Anodes, Cnodes))")
      if (length(dim(gform)) != 3 || dim(gform)[3] != num.regimes) stop("if gform is numeric, dim[3] should be num.regimes (gform can also be a matrix if variance.method == 'ic')")
      if (!is.null(deterministic.g.function)) stop("if gform is numeric, deterministic.g.function must be NULL")
      if (max(gform, na.rm=T) > 1 || min(gform, na.rm=T) < 0) stop("if gform is numeric, all values should be probabilities")
      if (!is.null(deterministic.Q.function) && !isTRUE(attr(data, "called.from.estimate.variance", exact=TRUE))) warning("If gform is numeric and deterministic.Q.function is not NULL, deterministic.Q.function will only affect g based on the observed values of the Anodes, not the counterfactual values. If your deterministic.Q.function does depends on the values of the Anodes, it is recommended to not use numeric gform.")
    }
  }
  
  #If Qform is NULL, it will be set by GetDefaultForm; no need to check here
  if (length(Qform) > 0) {
    if (! is.character(Qform)) stop("Qform should be a character vector")
    if (length(Qform) != length(nodes$LY)) stop("length of Qform is not equal to number of L/Y nodes")
    for (i in 1:length(Qform)) {
      if (length(names(Qform[i])) == 0) stop("Each element of Qform must be named. The name must match the name of the corresponding L/Y node in data.")
      if (names(Qform[i]) != names(data)[nodes$LY[i]]) stop("The name of each element of Q must match the name of the corresponding L/Y node in data.")
      if (Qform[i] != "IDENTITY") {
        #IDENTITY is only meant to be used in a call by ltmle:::EstimateVariance
        if (LhsVars(Qform[i]) != "Q.kplus1") stop("LHS of each Qform should be Q.kplus1")
        parents <- names(data)[1:(nodes$LY[i]-1)]
        if (!all(RhsVars(Qform[i]) %in% parents)) {
          stop("Some nodes in Qform[", i, "] are not parents of ", names(Qform[i]))
        }
        if (any(RhsVars(Qform[i]) %in% names(data)[nodes$C])) stop("Cnodes should not be used as RHS variables in Qform (regressions are only run on uncensored observations so including a Cnode has no effect and slows down regressions)")
      }
    }
  }
  
  if (length(gbounds) != 2) stop("gbounds should have length 2")
  if (! (is.null(deterministic.g.function) || is.function(deterministic.g.function))) stop("deterministic.g.function should be a function or NULL")
  
  if (! all(unlist(data[, nodes$A]) %in% c(0, 1, NA))) stop("in data, all Anodes should be binary")
  #note: Cnodes are checked in ConvertCensoringNodes
  if (! all(sapply(data[, c(nodes$A, nodes$Y), drop=F], is.numeric))) stop("in data, all Anodes and Ynodes should be numeric (not, for instance, logical)")
  all.Y <- unlist(data[, nodes$Y])
  
  binaryOutcome <- all(all.Y %in% c(0, 1, NA))
  
  if (binaryOutcome) {
    if (is.null(survivalOutcome)) {
      if (length(nodes$Y) == 1) {
        survivalOutcome <- FALSE #doesn't matter 
      } else {
        stop("All Ynodes are 0, 1, or NA; the outcome is treated as binary. The 'survivalOutcome' argument must be specified if there are multiple Ynodes.")
      }
    }    
    if (!is.null(Yrange) && !is.equal(Yrange, c(0L, 1L))) {
      stop("All Ynodes are 0, 1, or NA, but Yrange is something other than NULL or c(0, 1)") 
    }
  } else {
    if (is.null(survivalOutcome)) survivalOutcome <- FALSE
    if (survivalOutcome) stop("When survivalOutcome is TRUE, all Ynodes should be 0, 1, or NA")
  }
  
  uncensored.array <- CalcUncensoredMatrix(data, nodes$C)  
  for (Ynode in nodes$Y) {
    uncensored <- IsUncensored(uncensored.array, nodes$C, cur.node=Ynode)
    deterministic <- IsDeterministic(data, cur.node=Ynode, deterministic.Q.function=NULL, nodes, called.from.estimate.g=FALSE, survivalOutcome)$is.deterministic #pass deterministic.Q.function=NULL so we're only picking up deaths (if surivalOutcome=FALSE, deterministic will all be FALSE)
    if (anyNA(data[deterministic, Ynode]) || ! all(data[deterministic, Ynode] == 1)) stop("For survival outcomes, once a Ynode jumps to 1 (e.g. death), all subsequent Ynode values should be 1.")    
    if (anyNA(data[uncensored, Ynode])) stop("Ynodes may not be NA except after censoring")
  }
  
  if (! is.equal(dim(regimes)[1:2], c(nrow(data), length(nodes$A)))) stop("Problem with abar or regimes:\n   In ltmleMSM, regimes should have dimensions n x num.Anodes x num.regimes\n   In ltmle, abar should be a matrix with dimensions n x num.Anodes or a vector with length num.Anodes")
  stopifnot(num.regimes == nrow(summary.measures))
  if (!all(regimes %in% c(0, 1, NA))) stop("all regimes should be binary")
  for (Anode.index in seq_along(nodes$A)) {
    first.LYnode <- min(nodes$LY)
    cur.node <- nodes$A[Anode.index]
    uncensored <- IsUncensored(uncensored.array, Cnodes = nodes$C, cur.node = cur.node)
    deterministic <- IsDeterministic(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g=TRUE, survivalOutcome)$is.deterministic
    if (anyNA(regimes[uncensored & !deterministic, Anode.index, ])) {
      stop("NA in regimes/abar not allowed (except after censoring/death)") 
    }
    if (cur.node < first.LYnode && anyNA(regimes[!deterministic, Anode.index, ])) {
      warning("NA in regimes/abar before the first L/Y node will probably cause an error") #may not cause an error if this A is not part of any Qform and variance.method="ic"  
    }
  }
  
  num.final.Ynodes <- length(final.Ynodes) 
  if ((length(dim(summary.measures)) != 3) || ! is.equal(dim(summary.measures)[c(1, 3)], c(num.regimes, num.final.Ynodes))) stop("summary.measures should be an array with dimensions num.regimes x num.summary.measures x num.final.Ynodes")
  if (class(working.msm) != "character") stop("class(working.msm) must be 'character'")
  if (LhsVars(working.msm) != "Y") stop("the left hand side variable of working.msm should always be 'Y' [this may change in future releases]")
  if (!is.vector(observation.weights) || length(observation.weights) != nrow(data) || anyNA(observation.weights) || any(observation.weights < 0) || max(observation.weights) == 0) stop("observation.weights must be NULL or a vector of length nrow(data) with no NAs, no negative values, and at least one positive value")
  
  if (!(is.null(msm.weights) || is.equal(msm.weights, "empirical") || is.equal(dim(msm.weights), c(nrow(data), num.regimes, num.final.Ynodes)) || is.equal(dim(msm.weights), c(num.regimes, num.final.Ynodes)))) {
    stop("msm.weights must be NULL, 'empirical', or an array with dim(msm.weights) = c(n, num.regimes, num.final.Ynodes) or c(num.regimes, num.final.Ynodes)") 
  }
  
  if (!(is.null(id) || (is.vector(id) && length(id) == nrow(data)))) stop("id must be a vector with length nrow(data) or be NULL")
  return(list(survivalOutcome=survivalOutcome, binaryOutcome=binaryOutcome, uncensored=uncensored.array))
}

TransformOutcomes <- function(data, nodes, Yrange) {
  all.Y <- unlist(data[, nodes$Y])
  transformOutcome <- FALSE
  if (!is.null(Yrange)) {
    #if Yrange was specified
    rng <- range(all.Y, na.rm=TRUE)
    if (min(rng) < min(Yrange) || max(rng) > max(Yrange)) {
      #Truncate if Y vals are outside Yrange
      message("Some Ynodes are not in [Yrange[1], Yrange[2]], Y values are truncated")
      data[,nodes$Y][data[,nodes$Y] < min(Yrange)]<- min(Yrange)
      data[,nodes$Y][data[,nodes$Y] > max(Yrange)] <- max(Yrange)       
    } 
    #Then transform
    transformOutcome <- TRUE
  } else {
    #if Yrange was not specified, get it
    Yrange <- range(all.Y, na.rm=TRUE)
    if (min(Yrange) < 0 || max(Yrange) > 1) {
      #And see if we need to transform
      transformOutcome <- TRUE
      message("Some Ynodes are not in [0, 1], and Yrange was NULL, so all Y nodes are being\ntransformed to (Y-min.of.all.Ys)/range.of.all.Ys") 
    }
  }
  if (transformOutcome) {
    attr(transformOutcome, 'Yrange') <- Yrange 
    data[,nodes$Y] <- (data[, nodes$Y]-min(Yrange))/diff(Yrange) 
  }
  return(list(data=data, transformOutcome=transformOutcome))
}

# Set all nodes (except Y) to NA after death or censoring; Set Y nodes to 1 after death
CleanData <- function(data, nodes, deterministic.Q.function, survivalOutcome, showMessage=TRUE) {
  #make sure binaries have already been converted before calling this function
  is.nan.df <- function (x) {
    y <- if (length(x)) {
      do.call("cbind", lapply(x, "is.nan"))
    } else {
      matrix(FALSE, length(row.names(x)), 0)
    }
  }
  is.na.strict <- function (x) is.na(x) & !is.nan.df(x)  #only for data.frames
  changed <- FALSE
  ua <- rep(TRUE, nrow(data))  #uncensored and alive
  if (ncol(data) == 1) return(data)
  deterministic.Q.function.depends.on.called.from.estimate.g <- !is.null(deterministic.Q.function) && length(grep("called.from.estimate.g", as.character(body(deterministic.Q.function)))) > 0
  for (i in 1:(ncol(data)-1)) {
    if (anyNA(data[ua, i])) stop("NA values are not permitted in data except after censoring or a survival event")
    if (i %in% c(nodes$L, nodes$Y, nodes$AC)) { #don't use nodes$LY - this can leave out some Y nodes because of blocks
      is.deterministic <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=TRUE, survivalOutcome=survivalOutcome)$is.deterministic #check determinisitic including node i 
      
      if (deterministic.Q.function.depends.on.called.from.estimate.g) {
        is.deterministic.Q <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=FALSE, survivalOutcome=survivalOutcome)$is.deterministic 
        if (any(is.deterministic[ua] & !is.deterministic.Q[ua])) stop("Any row set deterministic by deterministic.Q.function(..., called.from.estimate.g=TRUE) must imply that the row is also set deterministic by deterministic.Q.function(..., called.from.estimate.g=FALSE)") #det.Q.fun(T) should imply det.Q.fun(F)
      }
      
      ua[ua] <- !is.deterministic[ua]
      if (anyNA(ua)) stop("internal ltmle error - ua should not be NA in CleanData")
      if (! all(is.na.strict(data[is.deterministic, setdiff((i+1):ncol(data), nodes$Y), drop=FALSE]))) {
        data[is.deterministic, setdiff((i+1):ncol(data), nodes$Y)] <- NA #if deterministic, set all nodes except Y to NA
        changed <- TRUE
      }
      
      if (i %in% nodes$C) {
        censored <- data[, i] == "censored" & ua
        if (! all(is.na.strict(data[censored, (i+1):ncol(data), drop=FALSE]))) {
          data[censored, (i+1):ncol(data)] <- NA  #if censored, set all nodes (including Y) to NA
          changed <- TRUE
        }
        ua[ua] <- !censored[ua] 
        if (anyNA(ua)) stop("internal ltmle error - ua should not be NA in CleanData")
      } 
      if (changed && showMessage) {
        message("Note: for internal purposes, all nodes after a censoring event are set to NA and \n all nodes (except Ynodes) are set to NA after Y=1 if survivalFunction is TRUE (or if specified by deterministic.Q.function).\n Your data did not conform and has been adjusted. This may be relevant if you are \n writing your own deterministic function(s) or debugging ltmle.")
      }
    }
  }
  return(data)
}

# Get the default Q or g formula - each formula consists of all parent nodes except censoring and event nodes [also except A nodes if stratifying]
GetDefaultForm <- function(data, nodes, is.Qform, stratify, survivalOutcome, showMessage) {
  if (is.Qform) {
    lhs <- rep("Q.kplus1", length(nodes$LY))
    node.set <- nodes$LY
  } else {
    lhs <- names(data)[nodes$AC]
    node.set <- nodes$AC
  }
  if (stratify) {
    stratify.nodes <- c(nodes$C, nodes$A)
  } else {
    stratify.nodes <- c(nodes$C)
  }
  if (survivalOutcome) {
    stratify.nodes <- c(stratify.nodes, nodes$Y)
  }
  form <- NULL
  for (i in seq_along(node.set)) {
    cur.node <- node.set[i]
    if (cur.node == 1) {
      form[i] <- paste(lhs[i], "~ 1")  #no parent nodes
    } else {
      parent.node.names <- names(data)[setdiff(1:(cur.node - 1), stratify.nodes)]
      if (length(parent.node.names) == 0) {
        form[i] <- paste(lhs[i], "~ 1")
      } else {
        form[i] <- paste(lhs[i], "~", paste(parent.node.names, collapse=" + "))     
      }
    }
    names(form)[i] <- names(data)[cur.node]
  }
  
  if (showMessage) {
    #Prints formulas with automatic wrapping thanks to print.formula
    message(ifelse(is.Qform, "Qform", "gform"),
            " not specified, using defaults:")
    lapply(seq_along(form), function(i, names) {
      message("formula for ", names[i], ":")
      #Using print on a formula because it nicely wraps
      message(capture.output(print(as.formula(form[i]), showEnv=FALSE)))
    }, names=names(form))
    message("")
  }
  return(form)
}

# Organize nodes
CreateNodes <- function(data, Anodes, Cnodes, Lnodes, Ynodes) {  
  Anodes <- NodeToIndex(data, Anodes)
  Cnodes <- NodeToIndex(data, Cnodes)
  Lnodes <- NodeToIndex(data, Lnodes)
  Ynodes <- NodeToIndex(data, Ynodes)
  nodes <- SuppressGivenWarnings(list(A=Anodes, C=Cnodes, L=Lnodes, Y=Ynodes, AC=sort(c(Anodes, Cnodes))), "is.na() applied to non-(list or vector) of type 'NULL'") #suppress warnings if no A/C nodes
  nodes$baseline <- sseq(1, min(c(nodes$A, nodes$L, nodes$C, nodes$Y)) - 1)
  nodes$LY <- CreateLYNodes(data, nodes, check.Qform=FALSE)
  return(nodes)
}

# Get the LY nodes but don't include "blocks" of L/Y nodes uninterrupted by A/C nodes
CreateLYNodes <- function(data, nodes, check.Qform, Qform) {
  LYnodes <- sort(c(nodes$L, nodes$Y))
  SuppressGivenWarnings(nodes.to.remove <- LYnodes[LYnodes < min(nodes$AC)], "no non-missing arguments to min; returning Inf") #no warning if no AC nodes
  
  #if there are no A/C nodes between two or more LY nodes, only the first LY node in the block is considered an LY node
  if (length(LYnodes) > 1) {
    for (i in 1:(length(LYnodes) - 1)) {
      cur.node <- LYnodes[i]
      next.node <- LYnodes[i + 1]
      if (! any(cur.node:next.node %in% nodes$AC)) {
        nodes.to.remove <- c(nodes.to.remove, next.node)
      }
    }
  }
  new.LYnodes <- setdiff(LYnodes, nodes.to.remove)
  if (check.Qform) {
    removed.Qform.index <- NULL
    for (i in nodes.to.remove) {
      index <- which(names(Qform) == names(data)[i])
      if (length(index) > 0) {
        removed.Qform.index <- c(removed.Qform.index, index)
      }
    }
    if (! is.null(removed.Qform.index)) {
      message("L/Y nodes (after removing blocks)  : ", names(data)[new.LYnodes], "\n")
      message("Qform names                        : ", names(Qform), "\n")
      message(paste("The following nodes are not being considered as L/Y nodes because they are part of a block of L/Y nodes. They are being dropped from Qform:\n"), paste(names(Qform)[removed.Qform.index], "\n", collapse=" "))
      Qform <- Qform[-removed.Qform.index]
    }
    return(list(LYnodes=new.LYnodes, Qform=Qform))
  }
  return(new.LYnodes)
}

# SL.library can be a character vector of library or a list with two separate vectors, one for Q and one for g
GetLibrary <- function(SL.library, estimate.type) {
  NullToGlm <- function(libr) if (is.null(libr)) "glm" else libr
  if (is.null(names(SL.library))) return(NullToGlm(SL.library))
  if (! identical(sort(names(SL.library)), sort(c("Q", "g")))) stop("If SL.library has names, it must have two names: Q and g")
  if (! estimate.type %in% c("Q", "g")) stop("bad estimate.type") 
  if (length(setdiff(names(attributes(SL.library)), c("names", "return.fit"))) > 0) stop("If SL.library has attributes, the only valid attributes are name and return.fit")
  lib <- SL.library[[estimate.type]]
  attr(lib, "return.fit") <- attr(SL.library, "return.fit", exact = TRUE)
  return(NullToGlm(SL.library[[estimate.type]]))
}

GetMsmWeights <- function(inputs) {
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  stopifnot(num.regimes >= 1)
  num.final.Ynodes <- length(inputs$final.Ynodes)
  if (is.equal(inputs$msm.weights, "empirical")) {
    #default is probability of following abar given alive, uncensored; conditioning on past treatment/no censoring, but not L, W; duplicates get weight 0
    msm.weights <- matrix(nrow=num.regimes, ncol=num.final.Ynodes)
    
    for (j in 1:num.final.Ynodes) {
      final.Ynode <- inputs$final.Ynodes[j]
      regimes.subset <- inputs$regimes[, inputs$all.nodes$A < final.Ynode, , drop = FALSE]
      if (ncol(regimes.subset) > 0) {
        is.duplicate <- duplicated(regimes.subset, MARGIN=3)
      } else {
        is.duplicate <- c(FALSE, rep(TRUE, num.regimes - 1))  #in case there are C nodes but no A nodes before a Ynode
      }
      uncensored <- IsUncensored(inputs$uncensored, inputs$all.nodes$C, cur.node=final.Ynode)
      intervention.match <- InterventionMatch(inputs$intervention.match, inputs$all.nodes$A, cur.node=final.Ynode)
      for (i in 1:num.regimes) {
        if (is.duplicate[i]) {
          msm.weights[i, j] <- 0
        } else {
          msm.weights[i, j] <- sum(uncensored & intervention.match[, i]) / nrow(inputs$data)
        } 
      }
    }
  } else if (is.null(inputs$msm.weights)) {
    msm.weights <- array(1, dim=c(n, num.regimes, num.final.Ynodes))
  } else {
    msm.weights <- inputs$msm.weights
  }
  if (is.equal(dim(msm.weights), c(num.regimes, num.final.Ynodes))) {
    msm.weights <- array(rep(msm.weights, each=n), dim=c(n, num.regimes, num.final.Ynodes))
  } 
  if (anyNA(msm.weights) || any(msm.weights < 0)) stop("all msm.weights must be >= 0 and not NA")
  return(msm.weights)
}

# Converts a general formula to a main terms formula and combine summary measures with baseline covariates
# Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S1 + S3 + S4" where 
# S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
ConvertToMainTerms <- function(data, msm, summary.measures, nodes) {
  baseline.column.names <- names(data)[nodes$baseline]
  summary.column.names <- colnames(summary.measures)
  rhs.vars <- RhsVars(msm)
  if (length(intersect(baseline.column.names, summary.column.names)) > 0) stop("Baseline covariate columns of data and columns of summary.measures may not have the same name")
  if (!all(rhs.vars %in% c(baseline.column.names, summary.column.names))) stop("All right hand side variables in working.msm must be either column names of summary.measures or column names of baseline covariates")
  baseline.column.names <- intersect(baseline.column.names, rhs.vars)
  baseline.data <- data[, baseline.column.names, drop=FALSE]
  num.regimes <- dim(summary.measures)[1]
  num.summary.measures <- dim(summary.measures)[2]
  num.final.Ynodes <- dim(summary.measures)[3]
  n <- nrow(data)
  for (j in 1:num.final.Ynodes) {
    for (i in 1:num.regimes) {
      combined.summary.measures <- model.matrix(as.formula(msm), data.frame(Y=1, baseline.data, matrix(summary.measures[i, , j], nrow=n, ncol=num.summary.measures, byrow=TRUE, dimnames=list(NULL, colnames(summary.measures)))))
      if (i == 1 && j == 1) {
        #initialize here now that we know how many columns there are
        main.terms.summary.measures <- array(dim=c(n, ncol(combined.summary.measures), num.regimes, num.final.Ynodes))
        beta.names <- colnames(combined.summary.measures) #this is the same for all i and j
      }
      main.terms.summary.measures[, , i, j] <- combined.summary.measures
    }
  }
  colnames(main.terms.summary.measures) <- paste("S", 1:ncol(main.terms.summary.measures), sep="") #temp names
  main.terms.msm <- paste("Y ~ -1 +", paste(colnames(main.terms.summary.measures), collapse=" + ")) #formula using temp names 
  return(list(msm=main.terms.msm, summary.measures=main.terms.summary.measures, beta.names=beta.names, baseline.column.names=baseline.column.names))
}

# Convert censoring nodes stored as binaries into factors (factors are recommended but binaries are currently accepted)
ConvertCensoringNodes <- function(data, Cnodes, has.deterministic.functions=FALSE) {
  error.msg <- "in data, all Cnodes should be factors with two levels, 'censored' and 'uncensored'\n See ?BinaryToCensoring \n (binary is also accepted, where 0=censored, 1=uncensored, but this is not recommended)"
  for (i in Cnodes) {
    col <- data[, i]
    if (is.numeric(col)) {
      if (! all(col %in% c(0, 1, NA))) stop(error.msg)
      data[, i] <- BinaryToCensoring(is.uncensored=col)
      if (has.deterministic.functions) warning("Censoring nodes have been converted from binaries to factors - see ?BinaryToCensoring.\n Note that if you are writing your own deterministic.g.function or deterministic.Q.function that censoring nodes are converted to factors\n before these functions are called.")
    } else if (is.factor(col)) {
      if (! all(levels(col) %in% c("censored", "uncensored"))) {
        stop("all levels of data[, Cnodes] should be in censored, uncensored (NA should not be a level)") 
      }
      #no action required
    } else {
      stop(error.msg) 
    } 
  }
  return(data)
}

#Before passing data to SuperLearner, convert factor to binary
ConvertCensoringNodeToBinary <- function(x) {
  stopifnot(is.factor(x) && all(levels(x) %in% c("censored", "uncensored")))
  b <- rep(NA_integer_, length(x))
  b[x == "censored"] <- 0L
  b[x == "uncensored"] <- 1L
  return(b)
}

# We don't want to show all of the warnings 
GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred", "prediction from a rank-deficient fit may be misleading") 
  if (update.step) {
    # It's ok if the updating step doesn't converge, we'll fix it in FixScoreEquation 
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}

SetSeedIfRegressionTesting <- function() {
  #if comparing outputs between different versions of ltmle, we need to sync random numbers 
  #before calling SuperLearner or FixScoreEquation since these use random numbers
  seed <- Sys.getenv("LTMLE.REGRESSION.TESTING.SEED") 
  stopifnot(length(seed) == 1)
  if (seed != "") {
    seed <- as.numeric(seed)
    stopifnot(is.finite(seed))
    set.seed(seed)
  }
  invisible(NULL)
}

ltmle.glm <- function(formula, family, data, weights) {
  try.result <- try(m <- speedglm::speedglm(formula=formula, family=family, data=data, weights=weights, maxit=100), silent = TRUE)
  if (inherits(try.result, "try-error")) { 
    ShowGlmMessage()
    if (is.null(weights)) {
      m <- glm(formula=formula, family=family, data=data, control=glm.control(maxit=100)) 
    } else {
      m <- glm(formula=formula, family=family, data=data.frame(data, weights), weights=weights, control=glm.control(maxit=100))
    }
  } 
  return(m)
} 

ltmle.glm.fit <- function(x, y, weights, family, offset, intercept) {
  try.result <- try({
    m <- speedglm::speedglm.wfit(y=y, X=x, family=family, weights=weights, offset=offset, intercept=intercept, maxit=100)
    class(m) <- c("speedglm", "speedlm")
  }, silent = TRUE)
  if (inherits(try.result, "try-error")) {
    ShowGlmMessage()
    m <- glm.fit(x=x, y=y, family=family, weights=weights, offset=offset, intercept=intercept, control=glm.control(maxit=100)) 
    class(m) <- c("glm", "lm")
  }
  return(m)
}

ShowGlmMessage <- function() {
  message("speedglm failed, using glm instead. If you see a lot of this message and you have large absolute values in data[, Lnodes], you may get a speed performance improvement by rescaling these values.")
}

Default.SL.Library <- list("SL.glm",
                           "SL.stepAIC",
                           "SL.bayesglm", 
                           c("SL.glm", "screen.corP"), 
                           c("SL.step", "screen.corP"), 
                           c("SL.step.forward", "screen.corP"), 
                           c("SL.stepAIC", "screen.corP"), 
                           c("SL.step.interaction", "screen.corP"), 
                           c("SL.bayesglm", "screen.corP")
)  
