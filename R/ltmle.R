 # Longitudinal TMLE to estimate an intervention-specific mean outcome or marginal structural model


#v0.9: 5/3/13 Joshua Schwab  released package

#v0.91 - started fixing det.ac.map bug - lots to clean up:
#rename det.ac.map [done]

#fix det.Q.map [done]
#clean up CheckInputs [done]
#testdetgfun gives diff answers - i think they should be the same but maybe not? fix this and then check that appropriate error messages are thrown if detgfun is inconsistent with data [fixed/done]
# make the "all subs false" error clearer [done]
#SetA should also set all Cnodes to 1 [done]
#return all glm fits [done] error if Anodes and Cnodes both NULL [done] 
#replace all "d" with "data" [done]

#check with maya re ex 4 (says "check this")

#export BinaryToCensoring [in switch functions, just need to export]
#convert binary C to factor in ltmleMSM.private [done]
#document (change existing docs [done, except creating rd file for deterministic templates and example functions], clean g and Q examples [done] - also note that is.deterministic in det.Q.fun will be used by EstimateG unless user checks for this in function [check with maya])
#include example det g and Q functions [done - in switchfunctions]
#write functions to create shell for det.g and det.Q functions [in progress-see LinhBug/ dswitch functions - mostly done, check with maya re det q function]
#update release notes, update version, update authors
#write sam re changing default on survivalFunction to true

# General code flow:
#ltmle -> ltmleMSM.private(pooledMSM=T) -> ...
#ltmleMSM(pooledMSM=T) -> ltmleMSM.private(pooledMSM=T) -> ...
#ltmleMSM(pooledMSM=F) -> ltmle -> ltmleMSM.private(pooledMSM=T) -> ...

#longitudinal targeted maximum liklihood estimation for E[Y_a]
#' @export
ltmle <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=TRUE, Qform=NULL, gform=NULL, 
                  abar, rule=NULL, gbounds=c(0.01, 1), Yrange=NULL, deterministic.g.function=NULL, stratify=FALSE, 
                  SL.library=NULL, estimate.time=nrow(data) > 50, gcomp=FALSE, mhte.iptw=FALSE, 
                  iptw.only=FALSE, deterministic.Q.function=NULL) {
  if (!is.null(rule)) {
    if (!(missing(abar) || is.null(abar))) stop("'abar' should not be specified when using a 'rule' function")
    abar <- t(apply(data, 1, rule))
  }
  if (is.vector(abar)) {
    abar <- matrix(rep(abar, each=nrow(data)), nrow=nrow(data))
  } else if (is.null(abar)) {
    abar <- matrix(nrow=nrow(data), ncol=0)
  }
  regimens <- abar
  dim(regimens) <- c(nrow(regimens), ncol(regimens), 1)
  working.msm <- "Y ~ -1 + S1"
  summary.measures <- array(1, dim=c(1, 1, 1))
  colnames(summary.measures) <- "S1"
  summary.baseline.covariates <- NULL

  temp <- ltmleMSM.private(data=data, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=survivalOutcome, Qform=Qform, gform=gform, Yrange=Yrange, gbounds=gbounds, deterministic.g.function=deterministic.g.function, SL.library=SL.library, regimens=regimens, working.msm=working.msm, summary.measures=summary.measures, summary.baseline.covariates=summary.baseline.covariates, final.Ynodes=NULL, pooledMSM=TRUE, stratify=stratify, weight.msm=FALSE, estimate.time=estimate.time, gcomp=gcomp, normalizeIC=FALSE, mhte.iptw=FALSE, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function) #it doesn't matter whether mhte.iptw is T or F when pooledMSM=T
  nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes)
  data <- ConvertCensoringNodes(data, Cnodes)
  iptw.list <- CalcIPTW(data, nodes, abar, drop3(temp$cum.g[, , 1, drop=F]), mhte.iptw) #get cum.g for regimen 1 (there's only 1 regimen)
  
  r <- list()
  if (iptw.only) {
    tmle <- NA
    tmle.IC <- rep(NA, nrow(data))
  } else {
    names(temp$beta) <- NULL
    tmle <- plogis(temp$beta)
    tmle.IC <- as.numeric(temp$IC)
  }
  r$estimates <- c(tmle=tmle, iptw=iptw.list$iptw.estimate, naive=iptw.list$naive.estimate)
  r$IC <- list(tmle=tmle.IC, iptw=iptw.list$iptw.IC)
  
  if (gcomp) {
    names(r$estimates)[1] <- names(r$IC)[1] <- "gcomp"
  }
  r$cum.g <- temp$cum.g[, , 1] #only one regimen
  r$call <- match.call()
  r$gcomp <- gcomp
  r$fit <- temp$fit
  r$fit$g <- r$fit$g[[1]]  #only one regimen

  r$formulas <- temp$formulas
  r$binaryOutcome <- temp$binaryOutcome
  r$transformOutcome <- temp$transformOutcome==TRUE #Want to store transformOutcome flag without attributes

  if (temp$transformOutcome) {
    Yrange <- attr(temp$transformOutcome, "Yrange")
    #back transform estimate and IC
    r$estimates[1] <- r$estimates[1]*diff(Yrange) + min(Yrange)
    estName <- names(r$estimates)[1]
    r$IC[[estName]] <- r$IC[[estName]]*diff(Yrange)
  }

  class(r) <- "ltmle"
  return(r)
}

#longitudinal targeted maximum likelihood estimation for a marginal structural model
#' @export 
ltmleMSM <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=TRUE, Qform=NULL, gform=NULL, gbounds=c(0.01, 1), Yrange=NULL, deterministic.g.function=NULL, SL.library=NULL, regimens, working.msm, summary.measures, summary.baseline.covariates=NULL, final.Ynodes=NULL, pooledMSM=TRUE, stratify=FALSE, weight.msm=TRUE, estimate.time=nrow(data) > 50, gcomp=FALSE, mhte.iptw=FALSE, iptw.only=FALSE, deterministic.Q.function=NULL, memoize=TRUE) {
  if (memoize && require(memoise)) {
    glm.ltmle.memoized <- memoize(glm.ltmle)
  }
  
  if (is.list(regimens)) {
    if (!all(do.call(c, lapply(regimens, is.function)))) stop("If 'regimens' is a list, then all elements should be functions.")
    regimens <- aperm(simplify2array(lapply(regimens, function(rule) apply(data, 1, rule)), higher=TRUE), c(2, 1, 3)) 
  }
  result <- ltmleMSM.private(data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, estimate.time, gcomp, normalizeIC=TRUE, mhte.iptw, iptw.only, deterministic.Q.function)
  result$call <- match.call()
  return(result) 
}

# This just shields the normalizeIC parameter, which should always be TRUE except when being called by ltmle
ltmleMSM.private <- function(data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, estimate.time, gcomp, normalizeIC, mhte.iptw, iptw.only, deterministic.Q.function) {

  nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes)
  Qform <- CreateLYNodes(data, nodes, check.Qform=TRUE, Qform=Qform)$Qform
  data <- ConvertCensoringNodes(data, Cnodes, has.deterministic.functions=!is.null(deterministic.g.function) && is.null(deterministic.Q.function))
  if (is.null(final.Ynodes)) {
    final.Ynodes <- max(nodes$Y)
  } else {
    final.Ynodes <- NodeToIndex(data, final.Ynodes)
  }
   
  if (identical(SL.library, 'default')) SL.library <- get("Default.SL.Library")
  
  if (is.null(Qform)) Qform <- GetDefaultForm(data, nodes, is.Qform=TRUE, stratify, survivalOutcome)
  if (is.null(gform)) gform <- GetDefaultForm(data, nodes, is.Qform=FALSE, stratify, survivalOutcome)
  
  if (length(dim(summary.measures)) == 2) {
    num.final.Ynodes <- length(final.Ynodes)
    summary.measures <- array(repmat(summary.measures, m=1, n=num.final.Ynodes), dim=c(nrow(summary.measures), ncol(summary.measures), num.final.Ynodes), dimnames=list(rownames(summary.measures), colnames(summary.measures), NULL))
  }
  
  #error checking (also transform Y in data if using continuous Y)
  check.results <- CheckInputs(data, nodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, deterministic.Q.function)
  data <- check.results$data
  
  if (estimate.time) EstimateTime(data, nodes, survivalOutcome, Qform, gform, gbounds, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only)
  
  result <- MainCalcs(data, nodes, survivalOutcome, Qform, gform, gbounds, deterministic.g.function, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, normalizeIC, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.function)
  result$gcomp <- gcomp
  result$formulas <- list(Qform=Qform, gform=gform)
  result$binaryOutcome <- check.results$binaryOutcome
  result$transformOutcome <- check.results$transformOutcome
  class(result) <- "ltmleMSM"
  return(result)
}

# Loop over final Ynodes, run main calculations
MainCalcs <- function(data, nodes, survivalOutcome, Qform, gform, gbounds, deterministic.g.function, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, normalizeIC, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.function) {
  if (! pooledMSM) {
    if (! is.null(summary.baseline.covariates)) stop("summary.baseline.covariates not supported")
    return(NonpooledMSM(data, nodes, survivalOutcome, Qform, gform, gbounds, deterministic.g.function, stratify, SL.library, regimens, working.msm, final.Ynodes, summary.measures, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.function))
  }
  if (iptw.only) {
    final.Ynodes <- final.Ynodes[length(final.Ynodes)]
  }
  # Several functions in the pooled version are only written to accept main terms MSM
  # Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ S1 + S1 + S3 + S4" where 
  # S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
  main.terms <- ConvertToMainTerms(working.msm, summary.measures)
  working.msm <- main.terms$msm
  summary.measures <- main.terms$summary.measures    
  num.final.Ynodes <- length(final.Ynodes)
  
  #summary.measures: num.regimens x num.measures x num.final.ynodes
  n <- nrow(data)
  num.regimens <- dim(regimens)[3]
  Qstar <- array(dim=c(n, num.regimens, num.final.Ynodes))
  weights <- matrix(nrow=num.regimens, ncol=num.final.Ynodes)
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- final.Ynodes[j]
    if (is.matrix(gform)) {
      gform1 <- gform[, nodes$AC < final.Ynode, drop=FALSE]
    } else {
      gform1 <- gform[nodes$AC < final.Ynode]
    }
    
    #It would be better to reuse g instead of calculating the same thing every time final.Ynode varies (note: g does need to be recalculated for each abar/regimen) - memoizing gets around this to some degree but it could be written better
    fixed.tmle <- FixedTimeTMLE(data[, 1:final.Ynode, drop=FALSE], nodes$A[nodes$A <= final.Ynode], nodes$C[nodes$C <= final.Ynode], nodes$L[nodes$L <= final.Ynode], nodes$Y[nodes$Y <= final.Ynode], survivalOutcome, Qform[nodes$LY <= final.Ynode], gform1, gbounds,  deterministic.g.function, SL.library, regimens[, nodes$A <= final.Ynode, , drop=FALSE], working.msm, summary.measures[, , j, drop=FALSE], summary.baseline.covariates, stratify, weight.msm, gcomp, iptw.only, deterministic.Q.function)
    if (iptw.only) return(list(cum.g=fixed.tmle$cum.g))
    if (j == 1) {
      IC <- fixed.tmle$IC
    } else {
      IC <- IC + fixed.tmle$IC
    }
    Qstar[, , j] <- fixed.tmle$Qstar #[n x num.regimens]
    weights[, j] <- fixed.tmle$weights #[num.regimens]  
  }
  fitted.msm <- FitPooledMSM(working.msm, Qstar, summary.measures, weights, summary.baseline.covariates) 
  IC <- FinalizeIC(IC, summary.measures, summary.baseline.covariates, Qstar, fitted.msm$m.beta, weights, normalizeIC)
  beta <- coef(fitted.msm$m)
  names(beta) <- main.terms$beta.names
  return(list(IC=IC, msm=fitted.msm$m, beta=beta, cum.g=fixed.tmle$cum.g, fit=fixed.tmle$fit)) #note: only returns cum.g and fit for the last final.Ynode
}

# Fit the MSM
FitPooledMSM <- function(working.msm, Qstar, summary.measures, weights, summary.baseline.covariates) {
  #Qstar: n x num.regimens x num.final.ynodes
  #summary.measures: num.regimens x num.summary.measures x num.final.ynodes
  #weights: num.regimens x num.final.ynodes
  if (! is.null(summary.baseline.covariates)) stop("need to update for summary.baseline.covariates")
  
  n <- dim(Qstar)[1]
  num.regimens <- dim(Qstar)[2]
  num.final.ynodes <- dim(Qstar)[3]
  num.summary.measures <- dim(summary.measures)[2]
  
  X <- matrix(nrow=n * num.regimens * num.final.ynodes, ncol=num.summary.measures)
  colnames(X) <- colnames(summary.measures)
  cnt <- 1
  for (j in 1:num.final.ynodes) {
    for (i in 1:num.regimens) {
      X[cnt:(cnt+n-1), ] <- matrix(summary.measures[i, , j], nrow=n, ncol=num.summary.measures, byrow=T)
      cnt <- cnt + n
    }
  }
  Y <- as.vector(Qstar)
  weight.vec <- rep(weights, each=n)
  
  m <- glm(as.formula(working.msm), data=data.frame(Y, X), family="quasibinomial", weights=weight.vec, na.action=na.exclude) 
  m.beta <- predict(m, type="response")
  dim(m.beta) <- dim(Qstar)
  return(list(m=m, m.beta=m.beta))
}

# ltmleMSM for a single final.Ynode
FixedTimeTMLE <- function(data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome, Qform, gform, gbounds, deterministic.g.function, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, stratify, weight.msm, gcomp, iptw.only, deterministic.Q.function) {
  #summary.measures: num.regimens x num.summary.measures
  #summary.baseline.covariates: names/indicies: num.summary.baseline.covariates x 1 
  #stacked.summary.measures: (n*num.regimens) x (num.summary.measures + num.summary.baseline.covariates)
  stacked.summary.measures <- GetStackedSummaryMeasures(summary.measures, data[, summary.baseline.covariates, drop=FALSE])
  
  if (identical(SL.library, 'default')) SL.library <- get("Default.SL.library") #Using get to avoid the "no visible binding for global variable" note in R CMD check
  nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes)
  
  SL.library.Q <- GetLibrary(SL.library, "Q")
  SL.library.g <- GetLibrary(SL.library, "g")
  
  num.regimens <- dim(regimens)[3]
  n <- nrow(data)
  num.betas <- ncol(model.matrix(as.formula(working.msm), data=data.frame(Y=1, stacked.summary.measures)))
  tmle <- weights <- rep(NA, num.regimens)
  IC <- matrix(0, nrow=n, ncol=num.betas)
  cum.g <- array(0, dim=c(n, length(nodes$AC), num.regimens))
  if (num.regimens > 1) {
    is.duplicate <- duplicated(regimens, MARGIN=3)
  } else {
    is.duplicate <- FALSE  #without this if statement there were problems when abar=NULL (ncol(regimens)=0)
  }
  
  fit.g <- vector("list", num.regimens)
  for (i in 1:num.regimens) {
    if (is.duplicate[i]) {
      weights[i] <- 0
      g.list <- list(fit=list("no g fits because this is a duplicate regimen"))
    } else {
      abar <- GetABar(regimens, i)
      # estimate each g factor, and cumulative probabilities
      g.list <- EstimateG(data, survivalOutcome, gform, nodes, abar=abar, deterministic.g.function, stratify, gbounds, SL.library.g, deterministic.Q.function)
      cum.g[, , i] <- g.list$cum.g
      weights[i] <- ComputeGA(data, nodes$A, nodes$C, abar, final.Ynode=max(nodes$Y), weight.msm)
    } 
    fit.g[[i]] <- g.list$fit
  }
  if (iptw.only) return(list(cum.g=cum.g))
  Qstar.kplus1 <- matrix(data[, nodes$Y[length(nodes$Y)]], nrow=n, ncol=num.regimens)
  logitQ <- matrix(nrow=n, ncol=num.regimens)
  
  regimens.with.positive.weight <- which(weights > 0)
  fit.Q <- fit.Qstar <- vector("list", length(nodes$LY))
  names(fit.Q) <- names(fit.Qstar) <- names(data)[nodes$LY]
  for (j in length(nodes$LY):1){
    cur.node <- nodes$LY[j]
    deterministic.list.origdata <- IsDeterministic(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g=FALSE, survivalOutcome)
    uncensored <- IsUncensored(data, nodes$C, cur.node)
    intervention.match <- subs <- matrix(nrow=n, ncol=num.regimens)
    for (i in regimens.with.positive.weight) {
      abar <- GetABar(regimens, i)
      intervention.match[, i] <- InterventionMatch(data, abar=abar, nodes$A, cur.node)  
      newdata <- SetA(data, abar=abar, nodes, cur.node)
      deterministic.list.newdata <- IsDeterministic(newdata, cur.node, deterministic.Q.function, nodes, called.from.estimate.g=FALSE, survivalOutcome)
      if (stratify) {
        subs[, i] <- uncensored & intervention.match & !deterministic.list.origdata$is.deterministic
      } else {
        subs[, i] <- uncensored & !deterministic.list.origdata$is.deterministic
      }
      #check that sum(subs[,i]) > 0
      if (any(subs[, i])) {
        Q.est <- Estimate(Qform[j], data=data.frame(Q.kplus1=Qstar.kplus1[, i], data), family="quasibinomial", newdata=newdata, subs=subs[, i], SL.library=SL.library.Q, type="link")
        logitQ[, i] <- Q.est$predicted.values
        #pass back Q.est$fit
      } else {
        if (! all(deterministic.list.newdata$is.deterministic)) {
          msg <- paste0("ltmle failed trying to estimate ", Qform[j], " because there are no observations that are\nuncensored", ifelse(stratify, ", follow abar,", ""), " and are not set deterministically due to death or deterministic.Q.function\n")
          stop(msg)
        }
        Q.est <- list(fit="no estimation of Q at this node because all rows are set deterministically")
      }
    }
    if (all(deterministic.list.newdata$is.deterministic)) {
      #no updating needed if all rows are set deterministically
      Qstar <- matrix(deterministic.list.newdata$Q, nrow=n, ncol=num.regimens)
      if (max(abs(Qstar.kplus1 - Qstar)) > 1e-8) {
        #if Qstar.kplus1 != Qstar when all deterministic score equation will not be solved
        stop("inconsistency in deterministic data - all rows are set deterministically but the deterministically set values are not equal to Qstar.kplus1") 
      }
      Qstar.est <- list(fit="no updating at this node because all rows are set deterministically")
    } else {
      ACnode.index  <- which.max(nodes$AC[nodes$AC < cur.node])
      update.list <- UpdateQ(Qstar.kplus1, logitQ, stacked.summary.measures, subs, cum.g[, ACnode.index, ], working.msm, uncensored, intervention.match, weights, gcomp)
      Qstar <- update.list$Qstar
      
      Qstar[deterministic.list.newdata$is.deterministic, ] <- deterministic.list.newdata$Q
      curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimens.with.positive.weight)
      if (any(abs(colSums(curIC)) > 0.001) && !gcomp) {
        fix.score.list <- FixScoreEquation(Qstar.kplus1, update.list$h.g.ratio, uncensored, intervention.match, deterministic.list.newdata, update.list$off, update.list$X, regimens.with.positive.weight)
        Qstar <- fix.score.list$Qstar
        curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimens.with.positive.weight)
        update.list$fit <- fix.score.list$fit
      } 
      IC <- IC + curIC
    }
    Qstar.kplus1 <- Qstar
    fit.Q[[j]] <- Q.est$fit
    fit.Qstar[[j]] <- update.list$fit
  }
  #tmle <- colMeans(Qstar)
  return(list(IC=IC, Qstar=Qstar, weights=weights, cum.g=cum.g, fit=list(g=fit.g, Q=fit.Q, Qstar=fit.Qstar))) 
}

#final step in calculating TMLE influence curve
FinalizeIC <- function(IC, summary.measures, summary.baseline.covariates, Qstar, m.beta, weights, normalizeIC) {
  #mBeta, Qstar: n x num.regimens x num.final.ynodes
  #summary.measures: num.regimens x num.summary.measures x num.final.ynodes
  #weights: num.regimens x num.final.ynodes
  
  num.betas <- ncol(IC)
  num.regimens <- nrow(summary.measures)
  n <- nrow(Qstar)
  num.final.ynodes <- dim(Qstar)[3]
  
  if (num.betas != ncol(summary.measures)) stop("this will fail if num.betas != num.summary.measures") 
  if (! is.null(summary.baseline.covariates)) stop("need to update for summary.baseline.covariates")
  
  finalIC <- matrix(0, nrow=n, ncol=num.betas)
  for (j in 1:num.final.ynodes) {
    for (i in 1:num.regimens) {
      if (weights[i, j] > 0) {
        m1 <- matrix(Qstar[, i, j] - m.beta[, i, j], ncol=1)   #n x 1
        m2 <- matrix(summary.measures[i, , j], nrow=1)  #1 x num.betas  
        finalIC <- finalIC + weights[i, j] * (m1 %*% m2)
      }
    }  
  }
  
  if (any(abs(colSums(finalIC)) > 0.001 )) {
    msg <- capture.output(cat("final IC problem", colSums(finalIC)))
    warning(paste(msg, collapse="\n"))
  }
  IC <- IC + finalIC
  if (normalizeIC) {
    IC  <- NormalizeIC(IC, summary.measures, m.beta, ignore.bad.ic=FALSE, weights) 
  }
  return(IC) 
}

# Normalize the influence curve matrix
NormalizeIC <- function(IC, summary.measures, m.beta, ignore.bad.ic=FALSE, weights) {    
  #need to update for summary.baseline.covariates
  num.betas <- ncol(IC)
  num.regimens <- nrow(summary.measures)
  num.final.ynodes <- dim(summary.measures)[3]
  
  C <- matrix(0, nrow=num.betas, ncol=num.betas)
  for (j in 1:num.final.ynodes) {
    for (i in 1:num.regimens) {
      if (weights[i, j] > 0) {
        m.beta.temp <- m.beta[, i, j]  #This is a n x 1 vector that should be all the same, but could have NA values
        index <- min(which(!is.na(m.beta.temp)))
        stopifnot(is.finite(index)) #if all NA, this will be Inf
        stopifnot(all(m.beta.temp[!is.na(m.beta.temp)] == m.beta.temp[index]))
        m.beta.temp <- m.beta.temp[index]
        h <- matrix(summary.measures[i, , j], ncol=1) * weights[i, j] 
        C <- C + h %*% t(h) * m.beta.temp * (1 - m.beta.temp) / weights[i, j]
      }
    }
  }
  if (abs(det(C)) < 1e-18 ) {
    normalized.IC <- matrix(nrow=nrow(IC), ncol=num.betas)
    warning("det(C) = 0, normalized IC <- NA")
  } else {
    normalized.IC <- t(solve(C, t(IC))) #IC %*% solve(C) 
    if (!ignore.bad.ic && any(abs(colSums(normalized.IC)) > 0.001 )) {
      msg <- capture.output({
        cat("normalized IC problem", colSums(normalized.IC), "\n")
        cat("inv(C) = \n")
        print(solve(C))
        })
      warning(paste(msg, collapse="\n"))
    }
  }
  return(normalized.IC)
}

# Get a single regimen from the regimens array
GetABar <- function(regimens, i) {
  abar <- regimens[, , i]
  if (dim(regimens)[2] == 1) abar <- matrix(abar, ncol=1) #if there's only 1 Anode, make sure abar comes back as a matrix
  return(abar)
}

# Combine summary.measures and summary.baseline.covariates into one stacked matrix
GetStackedSummaryMeasures <- function(summary.measures, summary.baseline.covariates) {
  #summary.measures: num.regimens x num.summary.measures
  #summary.baseline.covariates: n x num.summary.baseline.covariates
  #stacked.summary.measures: (n*num.regimens) x (num.summary.measures + num.summary.baseline.covariates)
  summary.baseline.covariates <- as.matrix(summary.baseline.covariates)
  num.regimens <- nrow(summary.measures)
  n <- nrow(summary.baseline.covariates)
  stacked.matrix <- cbind(matrix(rep(summary.measures, each=n), nrow=n*num.regimens),
                          matrix(rep(summary.baseline.covariates, times=num.regimens), nrow=n*num.regimens))
  colnames(stacked.matrix) <- c(colnames(summary.measures), colnames(summary.baseline.covariates))
  rownames(stacked.matrix) <- paste0(rep(paste0("n", 1:n), times=num.regimens), rep(paste0("r", 1:num.regimens), each=n))
  return(stacked.matrix)
}

# Targeting step - update the initial fit of Q using clever covariates
UpdateQ <- function(Qstar.kplus1, logitQ, stacked.summary.measures, subs, cum.g, working.msm, uncensored, intervention.match, weights, gcomp) { 
  #Q, Qstar.kplus1: n x num.regimens
  #cum.g: n x num.regimens (already indexed for this node)
  #subs: n x num.regimens
  #uncensored: n x 1
  #intervention.match: n x num.regimens
  #summary.measures: num.regimens x num.summary.measures
  #summary.baseline.covariates: names/indicies: num.summary.baseline.covariates x 1
  #weights: num.regimens x 1
  #stacked.summary.measures: (n*num.regimens) x (num.summary.measures + num.summary.baseline.covariates)
  #h.g.ratio: n x num.regimens x num.betas
  
  n <- nrow(logitQ)
  num.regimens <- ncol(logitQ)
  off <- as.vector(logitQ)
  Y <- as.vector(Qstar.kplus1)
  weight.vec <- rep(weights, each=n)
  X <- stacked.summary.measures / matrix(cum.g, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures))
  indicator <- matrix(uncensored, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) *matrix(intervention.match, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) #I(A=rule and uncensored)
  f <- as.formula(paste(working.msm, "+ offset(off) - 1"))
  data.temp <- data.frame(Y, X * indicator, off)
  if (gcomp) {
    Qstar <- plogis(logitQ)
    m <- "no Qstar fit because gcomp=TRUE (so no updating step)"
  } else {
    ctrl <- glm.control(trace=FALSE, maxit=1000)
    SuppressGivenWarnings(m <- glm(f, data=data.temp, subset=as.vector(subs), family="quasibinomial", weights=weight.vec, control=ctrl), GetWarningsToSuppress(TRUE)) #this should include the indicators
    newdata <- data.frame(Y, X, off)
    SuppressGivenWarnings(Qstar <- matrix(predict(m, newdata=newdata, type="response"), nrow=nrow(logitQ)), GetWarningsToSuppress(TRUE))  #this should NOT include the indicators  #note: could also use plogis(off + X %*% coef(m))
  }
  h.g.ratio <- model.matrix(f, model.frame(f, data=data.temp, na.action=na.pass)) #this should include the indicators
  dim(h.g.ratio) <- c(n, num.regimens, ncol(h.g.ratio))
  for (i in 1:num.regimens) {
    if (weights[i] > 0) {
      h.g.ratio[, i, ] <- h.g.ratio[, i, ] * weights[i] 
    } else {
      h.g.ratio[, i, ] <- 0   #cum.g is 0 so X is NA so h.g.ratio is NA when weight is 0
    }
  }
  if (ncol(stacked.summary.measures) != dim(h.g.ratio)[3]) stop("only works if working.msm has only main terms")
  return(list(Qstar=Qstar, h.g.ratio=h.g.ratio, X=X, off=off, fit=m))
}

# Sometimes GLM doesn't converge and the updating step of TMLE doesn't solve the score equation (sum of TMLE influence curve not equal to zero). This function attempts to solve the score equation directly using various optimizers.
FixScoreEquation <- function(Qstar.kplus1, h.g.ratio, uncensored, intervention.match, deterministic.list, off, X, regimens.with.positive.weight) {
  CalcScore <- function(e) {
    Qstar <- QstarFromE(e)
    ICtemp <- CalcIC(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimens.with.positive.weight)
    return(sum(colSums(ICtemp) ^ 2)) #each column has to add to zero
  }
  
  QstarFromE <- function(e) {
    Qstar <- plogis(off + X %*% e) #X: n x (num.summary.measures + num.summary.baseline.covariates) (which should be num.beta);  e: num.beta x 1 
    Qstar <- matrix(Qstar, nrow=length(uncensored)) #this should NOT include the indicators
    Qstar[deterministic.list$is.deterministic, ] <- deterministic.list$Q
    return(Qstar)
  }
  
  FindMin <- function(minimizer) {
    if (minimizer == "DEoptim") {
      num.tries <- 1
    } else {
      num.tries <- 30 #increase this to solve more problems at the cost of longer run time
    }
    init.e <- numeric(num.betas) #first try an initial estimate of epsilon=0
    for (i in 1:num.tries) {
      if (minimizer == "nlminb") {
        m <- nlminb(start=init.e, objective=CalcScore, control=list(abs.tol=max.objective, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
        e <- m$par
        obj.val <- m$objective
      } else if (minimizer == "optim") {
        m <- optim(par=init.e, fn=CalcScore, control=list(abstol=max.objective, reltol=1e-14, maxit=2000))
        e <- m$par
        obj.val <- m$value
      } else if (minimizer == "nlm") {
        m <- nlm(f=CalcScore, p=init.e)
        e <- m$estimate
        obj.val <- m$minimum
      } else if (minimizer == "DEoptim") {
        m1 <- nlminb(start=init.e, objective=CalcScore, control=list(abs.tol=max.objective, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
        m <- DEoptim(CalcScore, lower=rep(-100,num.betas), upper=rep(100,num.betas), control=DEoptim.control(VTR=max.objective, itermax=20000, strategy=5, initialpop=rbind(init.e, m1$par, matrix(rnorm(10*num.betas^2, sd=3), ncol=num.betas)), NP= 2 + 10*num.betas, trace=500))
        e <- m$optim$bestmem
        obj.val <- m$optim$bestval
      } else {
        stop("bad minimizer")
      }
      if (obj.val < max.objective) {
        m$ltmle.msg <- paste("updating step using glm failed to solve score equation; solved using", minimizer)
        return(list(e=e, solved=TRUE, m=m))
      }
      init.e <- rnorm(num.betas) #if the first try didn't work, try a random initial estimate of epsilon
    }
    return(list(e=numeric(num.betas), solved=FALSE, m="score equation not solved!")) #return Q (not updated)
  }
  max.objective <- 0.0001 ^ 2
  num.betas <- ncol(X)
  
  for (offset.lbound in c(1e-8,0.0001, 0.001, 0.01)) {
    off <- Bound(off, qlogis(c(offset.lbound, 1-offset.lbound)))
    l <- FindMin("nlminb")
    if (! l$solved) l <- FindMin("optim")
    if (! l$solved) l <- FindMin("nlm")
    if (l$solved) break
  }
  if (! l$solved && require(DEoptim)) l <- FindMin("DEoptim")   #comment out this line to turn off DEoptim (it's slow)
  if (! l$solved) {warning("------ all minimizers failed, returning non-updated Q ------")}
  Qstar <- QstarFromE(l$e)
  return(list(Qstar=Qstar, fit=l$m))
}

# Estimate how long it will take to run ltmleMSM
EstimateTime <- function(data, nodes, survivalOutcome, Qform, gform, gbounds, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only) {
  sample.index <- sample(nrow(data), size=50)
  start.time <- Sys.time()
  if (is.matrix(gform)) gform <- gform[sample.index, , drop=F]
  try.result <- try(  MainCalcs(data[sample.index, ], nodes, survivalOutcome, Qform, gform, gbounds, deterministic.g.function=NULL, SL.library, regimens[sample.index, , , drop=F], working.msm, summary.measures, summary.baseline.covariates[sample.index, , drop=F], final.Ynodes, normalizeIC=FALSE, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.function=NULL), silent=TRUE)
  if (inherits(try.result, "try-error")) {
    message("Timing estimate unavailable")
  } else {
    elapsed.time <- Sys.time() - start.time 
    est.time <- round(sqrt(as.double(elapsed.time, units="mins") * nrow(data) / 50), digits=0)
    if (est.time == 0) est.time <- "< 1"
    message("Estimate of time to completion: ", est.time, " minutes ")
  }
  invisible(NULL)
}

# Get summary measures for one or two ltmle objects (standard errors, p-values, confidence intervals)
# If two objects, include effect measures (additive effect, relative risk, odds ratio)
#' @S3method summary ltmle
summary.ltmle <- function(object, control.object=NULL, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  #object is treatment, control.object is control
  if (! is.null(control.object) && class(control.object) != "ltmle") stop("the control.object argument to summary.ltmle must be of class ltmle")
  if (! estimator %in% c("tmle", "iptw", "gcomp", "naive")) stop("estimator should be one of: tmle, iptw, gcomp, naive")
  if (estimator == "tmle" && object$gcomp) stop("estimator 'tmle' is not available because ltmleMSM was called with gcomp=TRUE")
  if (estimator == "gcomp" && !object$gcomp) stop("estimator 'gcomp' is not available because ltmleMSM was called with gcomp=FALSE")
  treatment.summary <- GetSummary(object$estimates[estimator], object$IC[[estimator]], loggedIC=FALSE)
  if (! is.null(control.object)) {
    control.summary <- GetSummary(control.object$estimates[estimator], control.object$IC[[estimator]], loggedIC=FALSE)
    effect.measures <- GetEffectMeasures(est0=control.object$estimates[estimator], IC0=control.object$IC[[estimator]], est1=object$estimates[estimator], IC1=object$IC[[estimator]], binaryOutcome=object$binaryOutcome && control.object$binaryOutcome)
    effect.measures.summary <- lapply(effect.measures, function (x) GetSummary(x$est, x$IC, x$loggedIC))
  } else {
    control.summary <- effect.measures.summary <- NULL
  }
  ans <- list(treatment=treatment.summary, control=control.summary, effect.measures=effect.measures.summary, treatment.call=object$call, control.call=control.object$call, estimator=estimator)
  class(ans) <- "summary.ltmle"
  return(ans)
}

# Get summary measures for MSM parameters (standard errors, p-values, confidence intervals)
#' @S3method summary ltmleMSM
summary.ltmleMSM <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
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
  
  n <- nrow(IC)
  v <- apply(IC, 2, var)
  std.dev <- sqrt(v/n)
  pval <- 2 * pnorm(-abs(estimate / std.dev))
  CI <- GetCI(estimate, std.dev)
  cmat <- cbind(estimate, std.dev, CI, pval)
  dimnames(cmat) <- list(names(estimate), c("Estimate", "Std. Error", "CI 2.5%", "CI 97.5%", "p-value"))
  ans <- list(cmat=cmat, estimator=estimator, transformOutcome=object$transformOutcome)
  class(ans) <- "summary.ltmleMSM"
  return(ans)
}

# Print method for summary.ltmleMSM
#' @S3method print summary.ltmleMSM
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
  invisible(x)
}

# Print method for summary.ltmle
#' @S3method print summary.ltmle
print.summary.ltmle <- function(x, ...) {
  cat("Estimator: ", x$estimator, "\n")
  if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
  if (is.null(x$control)) {
    PrintCall(x$treatment.call)
    PrintSummary(x$treatment)
  } else {
    cat("Treatment ")
    PrintCall(x$treatment.call)
    cat("Control ")
    PrintCall(x$control.call)
    cat("Treatment Estimate:\n")
    PrintSummary(x$treatment)
    cat("\nControl Estimate:\n")
    PrintSummary(x$control)
    cat("\nAdditive Effect:\n")
    PrintSummary(x$effect.measure$ATE)
    if (!is.null(x$effect.measure$RR)) {
      cat("\nRelative Risk:\n")
      PrintSummary(x$effect.measure$RR)
    }
    if (!is.null(x$effect.measure$OR)) {
      cat("\nOdds Ratio:\n")
      PrintSummary(x$effect.measure$OR)
    }
  }
  invisible(x)
}

# Print method for ltmleMSM
#' @S3method print ltmleMSM
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
#' @S3method print ltmle
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
  cat("   Parameter Estimate: ", signif(x$estimate, 5))
  cat("\n    Estimated Std Err: ", signif(x$std.dev, 5))
  cat("\n              p-value: ", ifelse(x$pvalue <= 2*10^-16, "<2e-16",signif(x$pvalue, 5)))
  cat("\n    95% Conf Interval:",paste("(", signif(x$CI[1], 5), ", ", signif(x$CI[2], 5), ")", sep=""),"\n")
  invisible(x)
}

# Calculate estimate, standard deviation, p-value, confidence interval
GetSummary <- function(estimate, IC, loggedIC) {
  if (is.null(IC)) {
    std.dev <- NA
  } else {
    n <- length(IC)
    std.dev <- sqrt(var(IC) / n)
  }
  if (loggedIC) {
    pvalue <- 2 * pnorm(-abs(log(estimate) / std.dev))
    CI <- exp(GetCI(log(estimate), std.dev))
  } else {
    pvalue <- 2 * pnorm(-abs(estimate / std.dev))
    CI <- GetCI(estimate, std.dev)
  }
  
  return(list(estimate=estimate, std.dev=std.dev, pvalue=pvalue, CI=CI))
}

# Calculate 95% confidence interval
GetCI <- function(estimate, std.dev) {
  x <- qnorm(0.975) * std.dev
  CI <- cbind("2.5%"=estimate - x, "97.5%"=estimate + x)
  return(CI)
}

# Calculate Average Treatment Effect, Relative Risk, Odds Ratio
GetEffectMeasures <- function(est0, IC0, est1, IC1, binaryOutcome) {  
  names(est0) <- names(est1) <- NULL
  est.ATE <- est1 - est0

  if (binaryOutcome) {
    est.RR <- est1 / est0
    est.OR <- (est1/(1-est1)) / (est0 / (1-est0))
  }

  if (is.null(IC0)) {
    ATE.IC <- RR.IC <- OR.IC <- NULL
  } else {
    ATE.IC <- -IC0 + IC1   
    if (binaryOutcome) {
      logRR.IC <- -1/est0 * IC0 + 1/est1 * IC1
      logOR.IC <- 1/(est0^2 - est0) * IC0 + 1/(est1 - est1^2) * IC1 
    }
  }

  result <- list(ATE=list(est=est.ATE, IC=ATE.IC, loggedIC=FALSE))
  if (binaryOutcome) {
    result$RR <- list(est=est.RR, IC=logRR.IC, loggedIC=TRUE)
    result$OR <- list(est=est.OR, IC=logOR.IC, loggedIC=TRUE)
  }

  return(result)
}

# Calculate IPTW and naive estimates
CalcIPTW <- function(data, nodes, abar, cum.g, mhte.iptw) {
  n <- nrow(data)
  final.Ynode <- nodes$Y[length(nodes$Y)]
  
  uncensored <- IsUncensored(data, nodes$C, final.Ynode)
  intervention.match <- InterventionMatch(data, abar, nodes$A, final.Ynode)  #A==abar for all A
  index <- uncensored & intervention.match
  
  Y <- data[index, final.Ynode]
  g <- cum.g[index, ncol(cum.g)]
  
  iptw.IC <- rep(0, n)
  if (mhte.iptw) {
    iptw.estimate <- sum( Y / g ) / sum(1 / g) 
    iptw.IC[index] <- ((Y - iptw.estimate) / g) / (1/n * sum (1 / g))
  } else {
    iptw.estimate <- sum(Y / g) / n
    iptw.IC[index] <- Y / g  
    iptw.IC <- iptw.IC - iptw.estimate 
  }  
  
  naive.estimate <- mean( (data[, final.Ynode])[index] )
  return(list(iptw.estimate=iptw.estimate, naive.estimate=naive.estimate, iptw.IC=iptw.IC))
}

# Parametric estimation of each g-factor
EstimateG <- function(data, survivalOutcome, gform, nodes, abar, deterministic.g.function, stratify, gbounds, SL.library, deterministic.Q.function) {
  gmat <- matrix(NaN, nrow=nrow(data), ncol=length(nodes$AC))
  uncensored <- rep(TRUE, nrow(data))
  fit <- vector("list", length(nodes$AC))
  names(fit) <- names(data)[nodes$AC]
  for (i in 1:length(nodes$AC)) {
    cur.node <- nodes$AC[i]
    if (cur.node %in% nodes$A) {
      cur.abar <- abar[, nodes$A == cur.node]
    } else {
      cur.abar <- rep(1, nrow(data))  #if this is a cnode, abar is always 1
    }
    newdata <- SetA(data, abar, nodes, cur.node)
    deterministic.origdata <- IsDeterministic(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g=TRUE, survivalOutcome)$is.deterministic #deterministic due to death or Q.function
    deterministic.newdata <- IsDeterministic(newdata, cur.node, deterministic.Q.function, nodes, called.from.estimate.g=TRUE, survivalOutcome)$is.deterministic #deterministic due to death or Q.function - using data modified so A = abar
    if (is.numeric(gform)) {
      probAis1 <- gform[, i]  #if gform is numeric, it's a matrix of probAis1
      g.est <- list(fit="gform passed as numeric, so no estimation took place")
    } else {
      deterministic.g.list.origdata <- IsDeterministicG(data, cur.node, deterministic.g.function, nodes) #deterministic due to acnode map - using original data
      deterministic.g.list.newdata <- IsDeterministicG(newdata, cur.node, deterministic.g.function, nodes) #deterministic due to acnode map - using data modified so A = abar
      deterministic.g.origdata <- deterministic.g.list.origdata$is.deterministic
      uncensored <- IsUncensored(data, nodes$C, cur.node)
      
      if (stratify) {
        intervention.match <- InterventionMatch(data, abar, nodes$A, nodes$AC[i]) 
        subs <- uncensored & intervention.match & !deterministic.origdata & !deterministic.g.origdata
      } else {
        subs <- uncensored & !deterministic.origdata & !deterministic.g.origdata
      }
      
      if (all(deterministic.g.list.newdata$is.deterministic | deterministic.newdata)) {
        # all rows are set deterministically, no need to estimate
        g.est <- list(fit="all rows are set deterministically, no estimation at this node")
        probAis1 <- rep(NaN, nrow(data)) #this will be filled in below
      } else {
        # not all rows are set deterministically
        if (any(subs)) {
          g.est <- Estimate(gform[i], data=data, subs=subs, family="binomial", newdata=newdata, SL.library=SL.library, type="response")
          probAis1 <- g.est$predicted.values
        } else {
          msg <- paste0("ltmle failed trying to estimate ", gform[i], " because there are no observations that are\nuncensored", ifelse(stratify, ", follow abar,", ""), " and are not set deterministically due to death or deterministic.g.function or deterministic.Q.function\n")
          stop(msg)
        }
      }
      probAis1[deterministic.g.list.newdata$is.deterministic] <- deterministic.g.list.newdata$prob1
    } 
    #probAis1 is prob(a=1), gmat is prob(a=abar)
    #cur.abar can be NA after censoring/death if treatment is dynamic
    gmat[!is.na(cur.abar) & cur.abar == 1, i] <- probAis1[!is.na(cur.abar) & cur.abar == 1]
    gmat[!is.na(cur.abar) & cur.abar == 0, i] <- 1 - probAis1[!is.na(cur.abar) & cur.abar == 0]
    
    gmat[deterministic.newdata, i] <- 1  #a=abar deterministically after death or other deterministic Q
    fit[[i]] <- g.est$fit
  }
  cum.g <- CalcCumG(gmat, gbounds)
  return(list(cum.g=cum.g, fit=fit))
}

# Truncate values within supplied bounds
Bound <- function(x, bounds) {
  x[x < min(bounds)] <- min(bounds)
  x[x > max(bounds)] <- max(bounds)
  return(x)
}

# Convert named nodes to indicies of nodes
NodeToIndex <- function(data, node) {
  if (! is.data.frame(data)) stop("data must be a data frame")
  if (is.numeric(node) || is.null(node)) return(node)
  if (! is.character(node)) stop("nodes must be numeric, character, or NULL")
  index <- match(node, names(data))
  if (any(is.na(index))) {
    stop(paste("\nnamed node(s) not found:", node[is.na(index)]))
  }
  return(index)
}

# Run GLM or SuperLearner
Estimate <- function(form, data, subs, family, newdata, SL.library, type) {
  stopifnot(type %in% c("link", "response"))
  f <- as.formula(form)
  
  if (any(is.na(data[subs, LhsVars(f)]))) stop("NA in Estimate")
  if (is.null(SL.library) || length(RhsVars(f)) == 0) { #in a formula like "Y ~ 1", call glm
    #estimate using GLM
    if (sum(subs) > 1) {
      SuppressGivenWarnings({
        m <- get.stack("glm.ltmle.memoized", mode="function", ifnotfound=glm.ltmle)(form, data=data[subs, all.vars(f), drop=F], family=family, control=glm.control(trace=FALSE, maxit=1000)) #there's probably a better way to do this
        predicted.values <- predict(m, newdata=newdata, type=type)
      }, GetWarningsToSuppress())
    } else {
      #glm breaks when sum(subs) == 1
      predicted.values <- rep(data[subs, LhsVars(f)], nrow(newdata))
      m <- "fit not returned because there was only 1 observation to fit"
    }
  } else {
    #estimate using SuperLearner
    if (family == "quasibinomial") family <- "binomial"
    rhs <- setdiff(RhsVars(f), rownames(alias(f, data=data[subs,])$Complete))  #remove aliased columns from X - these can cause problems if they contain NAs and the user is expecting the column to be dropped
    new.subs <- apply(newdata[, rhs, drop=FALSE], 1, function (x) !any(is.na(x)))  #remove NA values from newdata - these will output to NA anyway and cause errors in SuperLearner
    
    m <- SuperLearner(Y=data[subs, LhsVars(f)], X=data[subs, rhs, drop=FALSE], SL.library=SL.library, verbose=FALSE, family=family, newX=newdata[new.subs, rhs, drop=FALSE])
    predicted.values <- rep(NA, nrow(newdata))
    predicted.values[new.subs] <- m$SL.predict
    if (max(predicted.values, na.rm=T) > 1 || min(predicted.values, na.rm=T) < 0) {
      stop("predicted.values > 1 or < 0")
    }
    if (type == "link") {
      stopifnot(family == "binomial")
      predicted.values <- qlogis(Bound(predicted.values, bounds=c(0.0001, 0.9999)))
    }
  }
  return(list(predicted.values=predicted.values, fit=m))
}

# This is here for memoizing
glm.ltmle <- function(f, data, family, control) {
  return(glm(f, data=data, family=family, control=control))
}

# Calculate bounded cumulative G
CalcCumG <- function(g, gbounds) {
  cum.g <- Bound(t(apply(g, 1, cumprod)), gbounds)
  cum.g <- matrix(cum.g, nrow=nrow(g)) #to fix problems where apply returns a vector
  return(cum.g)
}

# Determine which patients are following specified treatment regimen (abar)
#return vector of [numObservations x 1] I(A==abar) from Anodes[1] to the Anode just before cur.node
InterventionMatch <- function(data, abar, Anodes, cur.node) {
  intervention.match <- XMatch(data, abar, Anodes, cur.node, all, default=TRUE)
  return(intervention.match)
}

# Determine which patients are uncensored
#return vector of [numDataRows x 1] I(C=uncensored) from Cnodes[1] to the Cnode just before cur.node
IsUncensored <- function(data, Cnodes, cur.node) {
  uncensored <- XMatch(data, Xbar="uncensored", Cnodes, cur.node, all, default=TRUE)
  return(uncensored)
}

# Determine which patients have died or have Q set deterministically by user function
# return list:
#    is.deterministic: vector of [numObservations x 1] - true if patient is already dead before cur.node or set by deterministic.Q.function
#    Q.value: vector of [which(is.deterministic) x 1] - value of Q
IsDeterministic <- function(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g, survivalOutcome) {
  #set Q.value to 1 if previous y node is 1
  if (survivalOutcome) {
     is.deterministic <- XMatch(data, Xbar=1, nodes$Y, cur.node, any, default=FALSE) #deterministic if any previous y node is 1
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
  if (any(inconsistent.rows)) stop(paste("At node:",names(data)[cur.node], "deterministic.Q.function is inconsistent with data - Q.value is either 0 or 1 but this does not match the final Y node value\nCheck data rows:", rownames(data)[det.list$is.deterministic][inconsistent.rows]))
  
  #return combined values
  Q.value <- rep(NA, nrow(data))
  Q.value[is.deterministic] <- 1
  Q.value[det.list$is.deterministic] <- det.list$Q.value
  is.deterministic <- is.deterministic | det.list$is.deterministic
  Q.value <- Q.value[is.deterministic]
  if (any(is.na(c(is.deterministic, Q.value)))) stop("NA in is.deterministic or Q.value")
  return(list(is.deterministic=is.deterministic, Q.value=Q.value))
}

# Determine which patients have an Anode value which is deterministic 
# For example, deterministic.g.function may be used to specify that once a patient starts treatment, they stay on treatment and this should be taken into consideration during estimation of G
IsDeterministicG <- function(data, cur.node, deterministic.g.function, nodes) {
  default <- list(is.deterministic=rep(FALSE, nrow(data)), prob1=NULL)
  if (is.null(deterministic.g.function)) return(default)
  #put this in a try-catch?
  det.list <- deterministic.g.function(data=data, current.node=cur.node, nodes=nodes)
  if (is.null(det.list)) return(default)
  if (!is.list(det.list) || !setequal(names(det.list), c("is.deterministic", "prob1")) || length(det.list) != 2) stop("deterministic.g.function should return a list with names: is.deterministic, prob1")
  if (! length(det.list$prob1) %in% c(1, length(which(det.list$is.deterministic)))) stop("the length of the 'prob1' element of deterministic.g.function's return argument should be either 1 or length(which(det.list$is.deterministic))")
  
  inconsistent.rows <- (det.list$prob1 %in% c(0,1)) & (det.list$prob1 != data[det.list$is.deterministic, cur.node]) & !is.na(data[det.list$is.deterministic, cur.node])
  if (any(inconsistent.rows)) stop(paste("At node:",names(data)[cur.node], "deterministic.g.function is inconsistent with data - prob1 is either 0 or 1 but this does not match the node value.\nCheck data rows:", rownames(data)[det.list$is.deterministic][inconsistent.rows]))
  return(det.list)
}

#Utility function called by IsUncensored, IsDeterministic - compares history in d within Xnodes prior to cur.node to Xbar
#
#any.all should be either the function 'any' or the function 'all'
#default: value to return if value is NA or there are no nodes before cur.node (TRUE for InterventionMatch and IsUncensored because NA indicates person matches intervention/is uncensored until they died; FALSE for IsDeterministic because NA indicates person was alive until they were censored)
XMatch <- function(data, Xbar, Xnodes, cur.node, any.all, default) {
  if (!any(Xnodes < cur.node)) return(rep(default, nrow(data)))
  last.Xnode.index <- which.max(Xnodes[Xnodes < cur.node])
  
  Xnodes.subset <- Xnodes[1:last.Xnode.index]
  if (identical(Xbar, 1) || identical(Xbar, "uncensored")) {
    Xbar.subset <- Xbar 
  } else {
    Xbar.subset <- Xbar[, 1:last.Xnode.index]
  } 
  d.subset <- data[, Xnodes.subset, drop=FALSE]
  matches <- apply(d.subset == Xbar.subset, 1, any.all)
  matches[is.na(matches)] <- default
  return(matches)
}

# Calculate the TMLE influence curve for one node
CalcIC <- function(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimens.with.positive.weight) {
  n <- nrow(Qstar)
  num.regimens <- ncol(Qstar)
  num.betas <- dim(h.g.ratio)[3] #h.g.ratio: n x num.regimens x num.betas
  
  IC <- matrix(0, nrow=n, ncol=num.betas)
  for (i in regimens.with.positive.weight) {
    index <- uncensored & intervention.match[, i]
    if (any(h.g.ratio[index, i, ] != 0)) {
      regimenIC <- matrix(0, nrow=n, ncol=num.betas)
      regimenIC[index, ] <- (Qstar.kplus1[index, i] - Qstar[index, i]) * h.g.ratio[index, i, ]
      IC <- IC + regimenIC
    }
  }
  return(IC)
}

#Set the Anodes of d to abar and Cnodes to uncensored (up to and including cur.node - cur.node itself is included for consistency checking in DeterministicG)
SetA <- function(data, abar, nodes, cur.node) {
  Anode.index <- nodes$A <= cur.node
  data[, nodes$A[Anode.index]] <- abar[, Anode.index]
  
  Cnode.index <- nodes$C <= cur.node
  data[, nodes$C[Cnode.index]] <- 1 #uncensored
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
CheckInputs <- function(data, nodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, deterministic.Q.function) {
  if (!all(is.null(GetLibrary(SL.library, "Q")), is.null(GetLibrary(SL.library, "g")))) library("SuperLearner")
  #each set of nodes should be sorted - otherwise causes confusion with gform, Qform, abar
  if (is.unsorted(nodes$A, strictly=TRUE)) stop("Anodes must be in increasing order")
  if (is.unsorted(nodes$C, strictly=TRUE)) stop("Cnodes must be in increasing order")
  if (is.unsorted(nodes$L, strictly=TRUE)) stop("Lnodes must be in increasing order")
  if (is.unsorted(nodes$Y, strictly=TRUE)) stop("Ynodes must be in increasing order")
  if (is.unsorted(final.Ynodes, strictly=TRUE)) stop("final.Ynodes must be in increasing order")

  if (length(nodes$L) > 0) {
    if (min(nodes$L) < min(nodes$AC)) stop("Lnodes are not allowed before A/C nodes. If you want to include baseline nodes, include them in data but not in Lnodes")
    if (max(nodes$L) > max(nodes$Y)) stop("Lnodes are not allowed after the final Y node")
  }
  
  all.nodes <- c(nodes$A, nodes$C, nodes$L, nodes$Y)
  if (length(all.nodes) > length(unique(all.nodes))) stop("A node cannot be listed in more than one of Anodes, Cnodes, Lnodes, Ynodes")
  if (is.null(nodes$Y)) stop("Ynodes cannot be null")
  if (is.null(nodes$AC)) stop("Anodes and Cnodes cannot both be null")

  if (min(all.nodes) < ncol(data)) {
    if (!all((min(all.nodes):ncol(data)) %in% all.nodes)) {
      stop("All nodes after the first of A-, C-, L-, or Ynodes must be in A-, C-, L-, or Ynodes")
    }
  }

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
    }
  } else {
    if (! is.numeric(gform)) stop("gform should be a character vector or numeric")
    g <- as.matrix(gform)
    if (nrow(g) != nrow(data)) stop("if gform is numeric, it should have the same number of rows as data")
    if (ncol(g) != length(nodes$AC)) stop("if gform is numeric, it should have the same number of columns as length(c(Anodes, Cnodes))")
  }
  
  if (! is.character(Qform)) stop("Qform should be a character vector")
  if (length(Qform) != length(nodes$LY)) {
    stop("length of Qform is not equal to number of L/Y nodes")
  }
  for (i in 1:length(Qform)) {
    if (LhsVars(Qform[i]) != "Q.kplus1") stop("LHS of each Qform should be Q.kplus1")
    if (length(names(Qform[i])) == 0) stop("Each element of Qform must be named. The name must match the name of the corresponding L/Y node in data.")
    if (names(Qform[i]) != names(data)[nodes$LY[i]]) stop("The name of each element of Q must match the name of the corresponding L/Y node in data.")
    parents <- names(data)[1:(nodes$LY[i]-1)]
    if (!all(RhsVars(Qform[i]) %in% parents)) {
      stop("Some nodes in Qform[", i, "] are not parents of ", names(Qform[i]))
    }    
  }
  
  if (length(gbounds) != 2) stop("gbounds should have length 2")
  if (! (is.null(deterministic.g.function) || is.function(deterministic.g.function))) {
    stop("deterministic.g.function should be a function or NULL")
  }
  
  if (! all(unlist(data[, nodes$A]) %in% c(0, 1, NA))) stop("in data, all Anodes should be binary")
  #note: Cnodes are checked in ConvertCensoringNodes

  all.Y <- unlist(data[, nodes$Y])

  binaryOutcome <- all(all.Y %in% c(0, 1, NA))

  transformOutcome <- FALSE

  if (!binaryOutcome) {
    if (survivalOutcome) {
      stop("When survivalOutcome is TRUE, all Ynodes should be 0, 1, or NA")
    } else { #not a survival Outcome
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
        attr(transformOutcome, 'Yrange') <- Yrange 
        data[,nodes$Y] <- (data[, nodes$Y]-min(Yrange))/diff(Yrange) 
      } else {
        #if Yrange was not specified, get it
        Yrange <- range(all.Y, na.rm=TRUE)
        if (min(Yrange) < 0 || max(Yrange) > 1) {
          #And see if we need to transform
          transformOutcome <- TRUE
          attr(transformOutcome, 'Yrange') <- Yrange
          message("Some Ynodes are not in [0, 1], and Yrange was NULL, so all Y nodes are being\ntransformed to (Y-min.of.all.Ys)/range.of.all.Ys") 
          data[,nodes$Y] <- (data[, nodes$Y]-min(Yrange))/diff(Yrange)        
        }
      }
    }
  } else { #Is a binary outcome
    if (!is.null(Yrange) && !is.equal(Yrange, c(0L,1L))) {
      stop("All Ynodes are 0, 1, or NA, but Yrange is something other than NULL or c(0,1)")
    }
  }
  

  
  for (i in nodes$Y) {
    deterministic <- IsDeterministic(data, cur.node=i, deterministic.Q.function=NULL, nodes, called.from.estimate.g=FALSE, survivalOutcome)$is.deterministic
    if (survivalOutcome && any(is.na(data[deterministic, i])) || ! all(data[deterministic, i] == 1)) stop("For survival outcomes, once a Ynode jumps to 1 (e.g. death), all subsequent Ynode values should be 1.")    
  }
  
  if (! is.equal(dim(regimens)[1:2], c(nrow(data), length(nodes$A)))) stop("Problem with abar or regimens:\n   In ltmleMSM, regimens should have dimensions n x num.Anodes x num.regimens\n   In ltmle, abar should be a matrix with dimensions n x num.Anodes or a vector with length num.Anodes")
  num.regimens <- dim(regimens)[3]
  stopifnot(num.regimens == nrow(summary.measures))
  if (!all(regimens %in% c(0, 1, NA))) stop("all regimens should be binary")
  for (i in seq_along(nodes$A)) {
    cur.node <- nodes$A[i]
    uncensored <- IsUncensored(data, nodes$C, cur.node)
    deterministic <- IsDeterministic(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g=TRUE, survivalOutcome)$is.deterministic
    if (any(is.na(regimens[uncensored & !deterministic, i, ]))) {
      stop("NA in regimens/abar not allowed (except after censoring/death)")
    }
  }
 
  if (! is.null(summary.baseline.covariates)) stop("summary.baseline.covariates is currently in development and is not yet supported - set summary.baseline.covariates to NULL")
  if (LhsVars(working.msm) != "Y") stop("the left hand side variable of working.msm should always be 'Y' [this may change in future releases]")
  if (! all(RhsVars(working.msm) %in% colnames(summary.measures))) stop("all right hand side variables in working.msm should be found in the column names of summary.measures")
  dynamic.regimens <- !all(duplicated(regimens)[2:nrow(data)])
  #if (dynamic.regimens && weight.msm) stop("dynamic regimens are not currently supported with weight.msm=TRUE [under development]")
    return(list(data=data, binaryOutcome=binaryOutcome, transformOutcome=transformOutcome))
}


# Get the default Q or g formula - each formula consists of all parent nodes except censoring and event nodes [also except A nodes if stratifying]
GetDefaultForm <- function(data, nodes, is.Qform, stratify, survivalOutcome) {
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

  #Prints formulas with automatic wrapping thanks to print.formula
  message(ifelse(is.Qform, "Qform", "gform"),
          " not specified, using defaults:")
  lapply(seq_along(form), function(i, names) {
          message("formula for ", names[i], ":")
          #Using print on a formula because it nicely wraps
          message(capture.output(print(as.formula(form[i]), showEnv=FALSE)))
        }, names=names(form))
  message("")

  return(form)
}

# Organize nodes
CreateNodes <- function(data, Anodes, Cnodes, Lnodes, Ynodes) {  
  Anodes <- NodeToIndex(data, Anodes)
  Cnodes <- NodeToIndex(data, Cnodes)
  Lnodes <- NodeToIndex(data, Lnodes)
  Ynodes <- NodeToIndex(data, Ynodes)
  
  nodes <- list(A=Anodes, C=Cnodes, L=Lnodes, Y=Ynodes, AC=sort(c(Anodes, Cnodes)))
  nodes$LY <- CreateLYNodes(data, nodes, check.Qform=FALSE)
  return(nodes)
}

# Get the LY nodes but don't include "blocks" of L/Y nodes uninterrupted by A/C nodes
CreateLYNodes <- function(data, nodes, check.Qform, Qform) {
  LYnodes <- sort(c(nodes$L, nodes$Y))
  #if there are no A/C nodes between two or more LY nodes, only the first LY node in the block is considered an LY node
  nodes.to.remove <- NULL
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
  if (is.null(names(SL.library))) return(SL.library)
  if (! identical(sort(names(SL.library)), sort(c("Q", "g")))) stop("If SL.library has names, it must have two names: Q and g")
  if (! estimate.type %in% c("Q", "g")) stop("bad estimate.type")
  return(SL.library[[estimate.type]])
}

# The non-pooled version of the ltmleMSM 
NonpooledMSM <- function(data, nodes, survivalOutcome, Qform, gform, gbounds, deterministic.g.function, stratify, SL.library, regimens, working.msm, final.Ynodes, summary.measures, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.function) {  
  tmle.index <- ifelse(gcomp, "gcomp", "tmle")
  num.regimens <- dim(regimens)[3]
  num.final.Ynodes <- length(final.Ynodes)
  
  num.ACnodes <- sum(nodes$AC < max(final.Ynodes))
  tmle <- iptw <- weights <- matrix(nrow=num.regimens, ncol=num.final.Ynodes)
  IC <- IC.iptw <- array(dim=c(num.regimens, num.final.Ynodes, nrow(data)))
  cum.g <- array(dim=c(nrow(data), num.ACnodes, num.regimens))
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- final.Ynodes[j]
    if (is.matrix(gform)) {
      gform1 <- gform[, nodes$AC < final.Ynode, drop=FALSE]
    } else {
      gform1 <- gform[nodes$AC < final.Ynode]
    }
    is.duplicate <- duplicated(regimens[, nodes$A < final.Ynode, , drop=FALSE], MARGIN=3)
    
    #It would be better to reuse g instead of calculating the same thing every time final.Ynode varies (note: g does need to be recalculated for each abar/regimen) - memoizing gets around this to some degree but it could be written better
    for (i in which(! is.duplicate)) {
      abar <- drop3(regimens[, , i, drop=F])
      abar <- abar[, nodes$A < final.Ynode, drop=FALSE]
      
      weights[i, j] <- ComputeGA(data[, 1:final.Ynode, drop=FALSE], nodes$A[nodes$A <= final.Ynode], nodes$C[nodes$C <= final.Ynode], abar, final.Ynode, weight.msm)
      
      if (weights[i, j] > 0) {
        result <- ltmle(data=data[, 1:final.Ynode, drop=FALSE], Anodes=nodes$A[nodes$A <= final.Ynode], Cnodes=nodes$C[nodes$C <= final.Ynode], Lnodes=nodes$L[nodes$L <= final.Ynode], Ynodes=nodes$Y[nodes$Y <= final.Ynode], survivalOutcome=survivalOutcome, Qform=Qform[nodes$LY <= final.Ynode], gform=gform1, abar=abar, gbounds=gbounds, deterministic.g.function=deterministic.g.function, stratify=stratify, SL.library=SL.library, estimate.time=FALSE, gcomp=gcomp, mhte.iptw=mhte.iptw, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function)
        tmle[i, j] <- result$estimates[tmle.index]
        iptw[i, j] <- min(1, result$estimates["iptw"])
        IC[i, j, ] <- result$IC[[tmle.index]]
        IC.iptw[i, j, ] <- result$IC[["iptw"]]
      }
      if (j == num.final.Ynodes) {
        if (weights[i, j] == 0) {
          #we didn't calculate cum.g because weight was 0 but we need to return it
          result <- ltmle(data=data[, 1:final.Ynode, drop=FALSE], Anodes=nodes$A[nodes$A <= final.Ynode], Cnodes=nodes$C[nodes$C <= final.Ynode], Lnodes=nodes$L[nodes$L <= final.Ynode], Ynodes=nodes$Y[nodes$Y <= final.Ynode], survivalOutcome=survivalOutcome, Qform=Qform[nodes$LY <= final.Ynode], gform=gform1, abar=abar, gbounds=gbounds, deterministic.g.function=deterministic.g.function, stratify=stratify, SL.library=SL.library, estimate.time=FALSE, gcomp=gcomp, mhte.iptw=mhte.iptw, iptw.only=TRUE, deterministic.Q.function=deterministic.Q.function)
        }
        cum.g[, , i] <- result$cum.g      
      }
    }
  }
  if (iptw.only) {
    m <- list()
  } else {
    m <- FitMSM(tmle, summary.measures, working.msm, IC, weights)
  }
  m.iptw <- FitMSM(iptw, summary.measures, working.msm, IC.iptw, weights)
  
  return(list(IC=m$beta.IC, msm=m$msm, beta=m$beta, cum.g=cum.g, beta.iptw=m.iptw$beta, IC.iptw=m.iptw$beta.IC))
}

# Called by NonpooledMSM to fit the MSM
FitMSM <- function(tmle, summary.measures, working.msm, IC, weights) {
  #summary.measures: num.regimens x num.measures x num.final.Ynodes
  #IC is num.regimens x num.final.Ynodes x n
  #tmle, weights:  num.regimens x num.final.Ynodes
  
  num.regimens <- nrow(tmle)
  num.final.Ynodes <- ncol(tmle)
  num.summary.measures <- ncol(summary.measures)
  
  n <- dim(IC)[3]
  Y <- as.vector(tmle)  
  weight.vec <- as.vector(weights)
  X <- apply(summary.measures, 2, rbind)
  summary.data <- data.frame(Y, X)
  m <- glm(formula=as.formula(working.msm), family="quasibinomial", data=summary.data, x=TRUE, na.action=na.exclude, weights=weight.vec)
  model.mat <- model.matrix.NA(as.formula(working.msm), summary.data) #model matrix (includes intercept and interactions)
  if (! is.equal(names(coef(m)), colnames(model.mat))) stop("re-ordering error - expecting same order of betas and model.mat columns")
  
  if (is.null(IC)) return(coef(m))
  num.coef <- length(coef(m))
  C <- matrix(0, nrow=num.coef, ncol=num.coef)
  for (j in 1:num.final.Ynodes) {
    for (i in 1:num.regimens) {
      if (! is.na(tmle[i, j])) {
        index <- sub2ind(row=i, col=j, num.rows=num.regimens)
        if (weights[i, j] > 0) { #if weight is 0, amount that would be added is zero (but divides by zero and causes NaN)
          h <- matrix(model.mat[index, ], ncol=1) * weights[i, j]  #index needs to pick up regimen i, time j
          m.beta <- predict(m, type="response")[index] 
          C <- C + h %*% t(h) * m.beta * (1 - m.beta) / weights[i, j]
          if (any(is.na(C))) {stop("NA in C")}
        }
      }
    }
  }
  beta.IC <- matrix(0, n, num.coef)
  for (k in 1:num.coef) {
    for (j in 1:num.final.Ynodes) {
      for (i in 1:num.regimens) {
        if (! is.na(tmle[i, j])) {
          index <- sub2ind(row=i, col=j, num.rows=num.regimens)
          h <- matrix(model.mat[index, ], ncol=1) * weights[i, j]   #index needs to pick up regimen i, time j
          if (abs(det(C)) < 1e-17) {
            W <- rep(NA, num.coef)
          } else {
            W <- solve(C, h)  #finds inv(C) * h
          }
          beta.IC[, k] <- beta.IC[, k] + W[k] * IC[i, j, ]
        }
      }
    }
  }
  return(list(beta=coef(m), msm=m, beta.IC=beta.IC))
}


#compute G - prob of following abar given alive, uncensored; conditioning on past treatment/no censoring, but not L, W
ComputeGA <- function(data, Anodes, Cnodes, abar, final.Ynode, weight.msm) {
  if (!weight.msm) return(1)
  uncensored <- IsUncensored(data, Cnodes, cur.node=final.Ynode)
  intervention.match <- InterventionMatch(data, abar, Anodes, cur.node=final.Ynode)
  gA <- sum(uncensored & intervention.match) / nrow(data)
  return(gA)
}

# Converts a general formula to a main terms formula
# Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ S1 + S1 + S3 + S4" where 
# S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
ConvertToMainTerms <- function(msm, summary.measures) {
  #need to update for summary.baseline.covariates
  num.regimens <- dim(summary.measures)[1]
  num.final.Ynodes <- dim(summary.measures)[3]
  if (num.final.Ynodes > 1) {
    stacked.summary.measures <- apply(summary.measures, 2, rbind)
  } else {
    stacked.summary.measures <- summary.measures #without this "if" there were problems with 1 regimen, 1 final [there's probably a better way to deal with this]
    dim(stacked.summary.measures) <- dim(summary.measures)[1:2]
    dimnames(stacked.summary.measures) <- dimnames(summary.measures)[1:2]
  }
  temp.model.matrix <- model.matrix(as.formula(msm), data=data.frame(Y=1, stacked.summary.measures))
  num.betas <- ncol(temp.model.matrix)
  main.terms.summary.measures <- array(dim=c(num.regimens, num.betas, num.final.Ynodes))
  cnt <- 1
  for (i in 1:num.final.Ynodes) {
    main.terms.summary.measures[, , i] <- temp.model.matrix[cnt:(cnt + num.regimens - 1), ]
    cnt <- cnt + num.regimens
  }
  colnames(main.terms.summary.measures) <- paste("S", 1:ncol(main.terms.summary.measures), sep="") #temp names
  main.terms.msm <- paste("Y ~ -1 +", paste(colnames(main.terms.summary.measures), collapse=" + ")) #formula using temp names 
  return(list(msm=main.terms.msm, summary.measures=main.terms.summary.measures, beta.names=colnames(temp.model.matrix)))
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

# We don't want to show all of the warnings 
GetWarningsToSuppress <- function(update.step=FALSE) {
  warnings.to.suppress <- c("glm.fit: fitted probabilities numerically 0 or 1 occurred", "prediction from a rank-deficient fit may be misleading")
  if (update.step) {
    # It's ok if the updating step doesn't converge, we'll fix it in FixScoreEquation 
    warnings.to.suppress <- c(warnings.to.suppress, "glm.fit: algorithm did not converge")
  }
  return(warnings.to.suppress)
}
