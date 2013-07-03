# Longitudinal TMLE to estimate an intervention-specific mean outcome or marginal structural model

#v0.9: 5/3/13 Joshua Schwab  released package

# General code flow:
#ltmle -> ltmleMSM.private(pooledMSM=T) -> ...
#ltmleMSM(pooledMSM=T) -> ltmleMSM.private(pooledMSM=T) -> ...
#ltmleMSM(pooledMSM=F) -> ltmle -> ltmleMSM.private(pooledMSM=T) -> ...

#longitudinal targeted maximum liklihood estimation for E[Y_a]
#' @export
ltmle <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, Qform=NULL, gform=NULL, 
                  abar, rule=NULL, gbounds=c(0.01, 1), deterministic.acnode.map=NULL, stratify=FALSE, 
                  SL.library=NULL, estimate.time=nrow(data) > 50, gcomp=FALSE, mhte.iptw=FALSE, 
                  iptw.only=FALSE, deterministic.Q.map=NULL) {
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

  temp <- ltmleMSM.private(data=data, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, Qform=Qform, gform=gform, gbounds=gbounds, deterministic.acnode.map=deterministic.acnode.map, SL.library=SL.library, regimens=regimens, working.msm=working.msm, summary.measures=summary.measures, summary.baseline.covariates=summary.baseline.covariates, final.Ynodes=NULL, pooledMSM=TRUE, stratify=stratify, weight.msm=FALSE, estimate.time=estimate.time, gcomp=gcomp, normalizeIC=FALSE, mhte.iptw=FALSE, iptw.only=iptw.only, deterministic.Q.map=deterministic.Q.map) #it doesn't matter whether mhte.iptw is T or F when pooledMSM=T
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
  r$cum.g <- temp$cum.g
  r$call <- match.call()
  r$gcomp <- gcomp
  class(r) <- "ltmle"
  return(r)
}

#longitudinal targeted maximum likelihood estimation for a marginal structural model
#' @export 
ltmleMSM <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, Qform=NULL, gform=NULL, gbounds=c(0.01, 1), deterministic.acnode.map=NULL, SL.library=NULL, regimens, working.msm, summary.measures, summary.baseline.covariates=NULL, final.Ynodes=NULL, pooledMSM=TRUE, stratify=FALSE, weight.msm=TRUE, estimate.time=nrow(data) > 50, gcomp=FALSE, mhte.iptw=FALSE, iptw.only=FALSE, deterministic.Q.map=NULL, memoize=TRUE) {
  if (memoize && require(memoise)) {
    glm.ltmle.memoized <- memoize(glm.ltmle)
  }
  
  if (is.list(regimens)) {
    if (!all(do.call(c, lapply(regimens, is.function)))) stop("If 'regimens' is a list, then all elements should be functions.")
    regimens <- aperm(simplify2array(lapply(regimens, function(rule) apply(data, 1, rule)), higher=TRUE), c(2, 1, 3)) 
  }
  result <- ltmleMSM.private(data, Anodes, Cnodes, Lnodes, Ynodes, Qform, gform, gbounds, deterministic.acnode.map, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, estimate.time, gcomp, normalizeIC=TRUE, mhte.iptw, iptw.only, deterministic.Q.map)
  result$call <- match.call()
  return(result) 
}

# This just shields the normalizeIC parameter, which should always be TRUE except when being called by ltmle
ltmleMSM.private <- function(data, Anodes, Cnodes, Lnodes, Ynodes, Qform, gform, gbounds, deterministic.acnode.map, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, estimate.time, gcomp, normalizeIC, mhte.iptw, iptw.only, deterministic.Q.map) {
  
  nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes)
  Qform <- CreateLYNodes(data, nodes, check.Qform=TRUE, Qform=Qform)$Qform
  data <- ConvertCensoringNodes(data, Cnodes)
  if (is.null(final.Ynodes)) {
    final.Ynodes <- max(nodes$Y)
  } else {
    final.Ynodes <- NodeToIndex(data, final.Ynodes)
  }
   
  if (identical(SL.library, 'default')) SL.library <- list("SL.glm", "SL.glmnet", "SL.stepAIC", "SL.bayesglm", c("SL.glm", "screen.corP"), c("SL.glmnet", "screen.corP"), c("SL.step", "screen.corP"), c("SL.step.forward", "screen.corP"), c("SL.stepAIC", "screen.corP"), c("SL.step.interaction", "screen.corP"), c("SL.bayesglm", "screen.corP"))
  
  PrintForms <- function(t) {
    #Prints formulas with automatic wrapping thanks to print.formula
    invisible(lapply(seq_along(t), function(i, names) {
    cat("\t", names[i], ":\n")
    print(as.formula(t[i]), showEnv=FALSE)
    }, names=names(t)))
    cat("\n")
  }
  if (is.null(Qform)) {
    Qform <- GetDefaultForm(data, nodes, is.Qform=TRUE, stratify)
    cat("Qform not specified, using defaults:\n")
    PrintForms(Qform)
  }
  if (is.null(gform)) {
    gform <- GetDefaultForm(data, nodes, is.Qform=FALSE, stratify)
    cat("gform not specified, using defaults:\n")
    PrintForms(gform)
  }


  
  if (length(dim(summary.measures)) == 2) {
    num.final.Ynodes <- length(final.Ynodes)
    summary.measures <- array(repmat(summary.measures, m=1, n=num.final.Ynodes), dim=c(nrow(summary.measures), ncol(summary.measures), num.final.Ynodes), dimnames=list(rownames(summary.measures), colnames(summary.measures), NULL))
  }
  
  # error checking (also convert deterministic.acnode.map and deterministic.Q.map to index)
  deterministic.maps <- CheckInputs(data, nodes, Qform, gform, gbounds, deterministic.acnode.map, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, deterministic.Q.map)
  deterministic.acnode.map <- deterministic.maps$deterministic.acnode.map
  deterministic.Q.map <- deterministic.maps$deterministic.Q.map
  
  if (estimate.time) EstimateTime(data, nodes, Qform, gform, gbounds, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only)
  
  result <- MainCalcs(data, nodes, Qform, gform, gbounds, deterministic.acnode.map, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, normalizeIC, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.map)
  result$gcomp <- gcomp
  class(result) <- "ltmleMSM"
  return(result)
}

# Loop over final Ynodes, run main calculations
MainCalcs <- function(data, nodes, Qform, gform, gbounds, deterministic.acnode.map, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, normalizeIC, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.map) {
  if (! pooledMSM) {
    if (! is.null(summary.baseline.covariates)) stop("summary.baseline.covariates not supported")
    return(NonpooledMSM(data, nodes$A, nodes$C, nodes$L, nodes$Y, Qform, gform, gbounds, deterministic.acnode.map, stratify, SL.library, regimens, working.msm, final.Ynodes, summary.measures, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.map))
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
    det.ac.map <- TruncateDeterministicNodeMap(data, deterministic.acnode.map, final.Ynode)
    det.Q.map <- TruncateDeterministicNodeMap(data, deterministic.Q.map, final.Ynode)
    fixed.tmle <- FixedTimeTMLE(data[, 1:final.Ynode, drop=FALSE], nodes$A[nodes$A <= final.Ynode], nodes$C[nodes$C <= final.Ynode], nodes$L[nodes$L <= final.Ynode], nodes$Y[nodes$Y <= final.Ynode], Qform[nodes$LY <= final.Ynode], gform1, gbounds, det.ac.map, SL.library, regimens[, nodes$A <= final.Ynode, , drop=FALSE], working.msm, summary.measures[, , j, drop=FALSE], summary.baseline.covariates, stratify, weight.msm, gcomp, iptw.only, det.Q.map)
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
  return(list(IC=IC, msm=fitted.msm$m, beta=beta, cum.g=fixed.tmle$cum.g)) #note: only returns cum.g for the last final.Ynode
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
FixedTimeTMLE <- function(data, Anodes, Cnodes, Lnodes, Ynodes, Qform, gform, gbounds, deterministic.acnode.map, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, stratify, weight.msm, gcomp, iptw.only, deterministic.Q.map) {
  #summary.measures: num.regimens x num.summary.measures
  #summary.baseline.covariates: names/indicies: num.summary.baseline.covariates x 1 
  #stacked.summary.measures: (n*num.regimens) x (num.summary.measures + num.summary.baseline.covariates)
  stacked.summary.measures <- GetStackedSummaryMeasures(summary.measures, data[, summary.baseline.covariates, drop=FALSE])
  
  if (identical(SL.library, 'default')) SL.library <- list("SL.glm", "SL.glmnet", "SL.stepAIC", "SL.bayesglm", c("SL.glm", "screen.corP"), c("SL.glmnet", "screen.corP"), c("SL.step", "screen.corP"), c("SL.step.forward", "screen.corP"), c("SL.stepAIC", "screen.corP"), c("SL.step.interaction", "screen.corP"), c("SL.bayesglm", "screen.corP"))
  nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes)
  
  if (is.null(Qform)) Qform <- GetDefaultForm(data, nodes, is.Qform=TRUE, stratify)
  if (is.null(gform)) gform <- GetDefaultForm(data, nodes, is.Qform=FALSE, stratify)
   
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
  
  for (i in 1:num.regimens) {
    if (is.duplicate[i]) {
      weights[i] <- 0
    } else {
      abar <- GetABar(regimens, i)
      # estimate each g factor, and cumulative probabilities
      gList <- EstimateG(data, gform, nodes, abar=abar, deterministic.acnode.map, stratify, gbounds, SL.library.g, deterministic.Q.map)
      
      cum.g[, , i] <- gList$cum.g
      weights[i] <- ComputeGA(data, nodes$A, nodes$C, abar, final.Ynode=max(nodes$Y), weight.msm)
    } 
  }
  if (iptw.only) return(list(cum.g=cum.g))
  Qstar.kplus1 <- matrix(data[, nodes$Y[length(nodes$Y)]], nrow=n, ncol=num.regimens)
  Q <- matrix(nrow=n, ncol=num.regimens)
  
  regimens.with.positive.weight <- which(weights > 0)
  for (j in length(nodes$LY):1){
    cur.node <- nodes$LY[j]
    deterministic.list <- IsDeterministic(data, nodes$Y, cur.node, deterministic.Q.map, called.from.estimate.g=FALSE)
    deterministic <- deterministic.list$is.deterministic
    uncensored <- IsUncensored(data, nodes$C, cur.node)
    intervention.match <- subs <- matrix(nrow=n, ncol=num.regimens)
    for (i in regimens.with.positive.weight) {
      abar <- GetABar(regimens, i)
      intervention.match[, i] <- InterventionMatch(data, abar=abar, nodes$A, cur.node)  
      if (stratify) {
        newdata <- data
        subs[, i] <- uncensored & intervention.match & !deterministic
      } else {
        newdata <- SetA(data, abar=abar, nodes$A, cur.node)
        subs[, i] <- uncensored & !deterministic
      }
      Q[, i] <- Estimate(Qform[j], d=data.frame(Q.kplus1=Qstar.kplus1[, i], data), family="quasibinomial", newdata=newdata, subs=subs[, i], SL.library=SL.library.Q)
    }
    ACnode.index  <- which.max(nodes$AC[nodes$AC < cur.node])
    update.list <- UpdateQ(Qstar.kplus1, Q, stacked.summary.measures, subs, cum.g[, ACnode.index, ], working.msm, uncensored, intervention.match, weights, gcomp)
    Qstar <- update.list$Qstar
    Qstar[deterministic, ] <- deterministic.list$Q
    
    curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimens.with.positive.weight)
    if (any(abs(colSums(curIC)) > 0.001) && !gcomp) {
      Qstar <- FixScoreEquation(Qstar.kplus1, update.list$h.g.ratio, uncensored, intervention.match, deterministic.list, update.list$off, update.list$X, regimens.with.positive.weight)
      curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimens.with.positive.weight)
    }
    
    IC <- IC + curIC
    Qstar.kplus1 <- Qstar
  }
  #tmle <- colMeans(Qstar)
  return(list(IC=IC, Qstar=Qstar, weights=weights, cum.g=cum.g)) 
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
  
  if (any(abs(colSums(finalIC)) > 0.001 )) {cat("final IC problem", colSums(finalIC), "\n")}
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
    cat("det(C) = 0, normalized IC <- NA\n")
  } else {
    normalized.IC <- t(solve(C, t(IC))) #IC %*% solve(C) 
    if (!ignore.bad.ic && any(abs(colSums(normalized.IC)) > 0.001 )) {
      cat("normalized IC problem", colSums(normalized.IC), "\n")
      cat("inv(C) = \n")
      print(solve(C))
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
UpdateQ <- function(Qstar.kplus1, Q, stacked.summary.measures, subs, cum.g, working.msm, uncensored, intervention.match, weights, gcomp) { 
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
  
  n <- nrow(Q)
  num.regimens <- ncol(Q)
  off <- qlogis(Bound(as.vector(Q), c(.0001, .9999)))
  Y <- as.vector(Qstar.kplus1)
  weight.vec <- rep(weights, each=n)
  X <- stacked.summary.measures / matrix(cum.g, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures))
  indicator <- matrix(uncensored, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) *matrix(intervention.match, nrow=nrow(stacked.summary.measures), ncol=ncol(stacked.summary.measures)) #I(A=rule and uncensored)
  f <- as.formula(paste(working.msm, "+ offset(off) - 1"))
  d <- data.frame(Y, X * indicator, off)
  if (gcomp) {
    Qstar <- Q
  } else {
    ctrl <- glm.control(trace=FALSE, maxit=1000)
    SuppressGivenWarnings(m <- glm(f, data=d, subset=as.vector(subs), family="quasibinomial", weights=weight.vec, control=ctrl), GetWarningsToSuppress(TRUE)) #this should include the indicators
    
    newdata <- data.frame(Y, X, off)
    SuppressGivenWarnings(Qstar <- matrix(predict(m, newdata=newdata, type="response"), nrow=nrow(Q)), GetWarningsToSuppress(TRUE))  #this should NOT include the indicators  #note: could also use plogis(off + X %*% coef(m))
  }
  h.g.ratio <- model.matrix(f, model.frame(f, data=d, na.action=na.pass)) #this should include the indicators
  dim(h.g.ratio) <- c(n, num.regimens, ncol(h.g.ratio))
  for (i in 1:num.regimens) {
    if (weights[i] > 0) {
      h.g.ratio[, i, ] <- h.g.ratio[, i, ] * weights[i] 
    } else {
      h.g.ratio[, i, ] <- 0   #cum.g is 0 so X is NA so h.g.ratio is NA when weight is 0
    }
  }
  if (ncol(stacked.summary.measures) != dim(h.g.ratio)[3]) stop("only works if working.msm has only main terms")
  return(list(Qstar=Qstar, h.g.ratio=h.g.ratio, X=X, off=off))
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
    init.e <- numeric(num.betas)
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
        m <- DEoptim(CalcScore, lower=rep(-100,num.betas), upper=rep(100,num.betas), control=DEoptim.control(VTR=max.objective, itermax=20000, strategy=5, initialpop=rbind(init.e, m1$par, matrix(rnorm(10*num.betas^2, sd=3), ncol=num.betas)), NP= 2 + 10*num.betas))
        e <- m$optim$bestmem
        obj.val <- m$optim$bestval
      } else {
        stop("bad minimizer")
      }
      if (obj.val < max.objective) {
        return(list(e=e, solved=TRUE))
      }
      init.e <- rnorm(num.betas)
    }
    return(list(e=numeric(num.betas), solved=FALSE)) #return Q (not updated)
  }
  
  max.objective <- 0.0001 ^ 2
  num.betas <- ncol(X)
  l <- FindMin("nlminb")
  if (! l$solved) l <- FindMin("optim")
  if (! l$solved) l <- FindMin("nlm")
  if (! l$solved && require(DEoptim)) l <- FindMin("DEoptim")   #comment out this line to turn off DEoptim (it's slow)
  if (! l$solved) {cat("\n\n -------- all minimizers failed, returning non-updated Q ------ \n\n\n")}
  Qstar <- QstarFromE(l$e)
  return(Qstar)
}

# Estimate how long it will take to run ltmleMSM
EstimateTime <- function(data, nodes, Qform, gform, gbounds, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only) {
  sample.index <- sample(nrow(data), size=50)
  start.time <- Sys.time()
  if (is.matrix(gform)) gform <- gform[sample.index, , drop=F]
  try.result <- try(  MainCalcs(data[sample.index, ], nodes, Qform, gform, gbounds, deterministic.acnode.map=NULL, SL.library, regimens[sample.index, , , drop=F], working.msm, summary.measures, summary.baseline.covariates[sample.index, , drop=F], final.Ynodes, normalizeIC=FALSE, pooledMSM, stratify, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.map=NULL), silent=TRUE)
  if (inherits(try.result, "try-error")) {
    cat("Timing estimate unavailable\n")
  } else {
    elapsed.time <- Sys.time() - start.time 
    est.time <- round(sqrt(as.double(elapsed.time, units="mins") * nrow(data) / 50), digits=0)
    if (est.time == 0) est.time <- "< 1"
    cat("Estimate of time to completion:", est.time, "minutes \n")
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
    effect.measures <- GetEffectMeasures(est0=control.object$estimates[estimator], IC0=control.object$IC[[estimator]], est1=object$estimates[estimator], IC1=object$IC[[estimator]])
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
  ans <- list(cmat=cmat, estimator=estimator)
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
    cat("\nControl ")
    PrintCall(x$control.call)
    cat("Treatment Estimate:\n")
    PrintSummary(x$treatment)
    cat("\nControl Estimate:\n")
    PrintSummary(x$control)
    cat("\nAdditive Effect:\n")
    PrintSummary(x$effect.measure$ATE)
    cat("\nRelative Risk:\n")
    PrintSummary(x$effect.measure$RR)
    cat("\nOdds Ratio:\n")
    PrintSummary(x$effect.measure$OR)
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
  cat("Call:  ", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
}

# Print estimate, standard error, p-value, confidence interval
PrintSummary <- function(x) {
  cat("   Parameter Estimate: ", signif(x$estimate, 5))
  cat("\n    Estimated Std Err: ", signif(x$std.dev^2, 5))
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
GetEffectMeasures <- function(est0, IC0, est1, IC1) {  
  names(est0) <- names(est1) <- NULL
  ATE <- est1 - est0
  RR <- est1 / est0
  OR <- (est1/(1-est1)) / (est0 / (1-est0))
  
  if (is.null(IC0)) {
    ATE.IC <- RR.IC <- OR.IC <- NULL
  } else {
    ATE.IC <- -IC0 + IC1   
    logRR.IC <- -1/est0 * IC0 + 1/est1 * IC1
    logOR.IC <- 1/(est0^2 - est0) * IC0 + 1/(est1 - est1^2) * IC1 
  }
  return(list(ATE=list(est=ATE, IC=ATE.IC, loggedIC=FALSE), RR=list(est=RR, IC=logRR.IC, loggedIC=TRUE), OR=list(est=OR, IC=logOR.IC, loggedIC=TRUE)))
}

# Calculate IPTW and naive estimates
CalcIPTW <- function(d, nodes, abar, cum.g, mhte.iptw) {
  n <- nrow(d)
  final.Ynode <- nodes$Y[length(nodes$Y)]
  
  uncensored <- IsUncensored(d, nodes$C, final.Ynode)
  intervention.match <- InterventionMatch(d, abar, nodes$A, final.Ynode)  #A==abar for all A
  index <- uncensored & intervention.match
  
  Y <- d[index, final.Ynode]
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
  
  naive.estimate <- mean( (d[, final.Ynode])[index] )
  return(list(iptw.estimate=iptw.estimate, naive.estimate=naive.estimate, iptw.IC=iptw.IC))
}

# Parametric estimation of each g-factor
EstimateG <- function(d, gform, nodes, abar, deterministic.acnode.map, stratify, gbounds, SL.library, deterministic.Q.map) {
  gmat <- matrix(NaN, nrow=nrow(d), ncol=length(nodes$AC))
  uncensored <- rep(TRUE, nrow(d))
  for (i in 1:length(nodes$AC)) {
    cur.node <- nodes$AC[i]
    if (cur.node %in% nodes$A) {
      cur.abar <- abar[, nodes$A == cur.node]
    } else {
      cur.abar <- rep(1, nrow(d))  #if this is a cnode, abar is always 1
    }
    deterministic <- IsDeterministic(d, nodes$Y, cur.node, deterministic.Q.map, called.from.estimate.g=TRUE)$is.deterministic #deterministic due to death or Q.map
    if (! is.numeric(gform)) {
      deterministic.g.list <- IsDeterministicG(d, cur.node, deterministic.acnode.map) #deterministic due to acnode map
      deterministic.g <- deterministic.g.list$is.deterministic
      
      uncensored <- IsUncensored(d, nodes$C, cur.node)
      
      if (stratify) {
        intervention.match <- InterventionMatch(d, abar, nodes$A, nodes$AC[i]) 
        subs <- uncensored & intervention.match & !deterministic & !deterministic.g
        newdata <- d
      } else {
        subs <- uncensored & !deterministic & !deterministic.g
        newdata <- SetA(d, abar, nodes$A, cur.node)
      }
      
      probAis1 <- rep(NaN, nrow(d))
      if (any(subs)) {
        probAis1 <- Estimate(gform[i], d=d, subs=subs, family="binomial", newdata=newdata, SL.library=SL.library)
      } else {
        stop("in EstimateG, all subs are FALSE")
        #warning("in EstimateG, all subs are FALSE") #not clear what to do here
        #probAis1 <- rep(0.5, nrow(d))
      }
      
      probAis1[deterministic.g] <- deterministic.g.list$prob1
    } else {
      probAis1 <- gform[, i]  #if gform is numeric, it's a matrix of probAis1
    }
    #probAis1 is prob(a=1), gmat is prob(a=abar)
    #cur.abar can be NA after censoring/death if treatment is dynamic
    gmat[!is.na(cur.abar) & cur.abar == 1, i] <- probAis1[!is.na(cur.abar) & cur.abar == 1]
    gmat[!is.na(cur.abar) & cur.abar == 0, i] <- 1 - probAis1[!is.na(cur.abar) & cur.abar == 0]
    
    gmat[deterministic, i] <- 1  #a=abar deterministically after death or other deterministic Q
  }
  cum.g <- CalcCumG(gmat, gbounds)
  return(list(cum.g=cum.g))
}

# Truncate values within supplied bounds
Bound <- function(x, bounds) {
  x[x < min(bounds)] <- min(bounds)
  x[x > max(bounds)] <- max(bounds)
  return(x)
}

# Convert named nodes to indicies of nodes
NodeToIndex <- function(d, node) {
  if (! is.data.frame(d)) stop("data must be a data frame")
  if (is.numeric(node) || is.null(node)) return(node)
  if (! is.character(node)) stop("nodes must be numeric, character, or NULL")
  index <- match(node, names(d))
  if (any(is.na(index))) {
    stop(paste("\nnamed node(s) not found:", node[is.na(index)]))
  }
  return(index)
}

# Run GLM or SuperLearner
Estimate <- function(form, d, subs, family, newdata, SL.library) {
  f <- as.formula(form)
  
  if (any(is.na(d[subs, LhsVars(f)]))) stop("NA in Estimate")
  if (is.null(SL.library) || length(RhsVars(f)) == 0) { #in a formula like "Y ~ 1", call glm
    #estimate using GLM
    if (sum(subs) > 1) {
      SuppressGivenWarnings({
        m <- get.stack("glm.ltmle.memoized", mode="function", ifnotfound=glm.ltmle)(form, data=d[subs, all.vars(f), drop=F], family=family, control=glm.control(trace=FALSE, maxit=1000)) #there's probably a better way to do this
        predicted.values <- predict(m, newdata=newdata, type="response")
      }, GetWarningsToSuppress())
    } else {
      #glm breaks when sum(subs) == 1
      predicted.values <- rep(d[subs, LhsVars(f)], nrow(newdata))
    }
  } else {
    #estimate using SuperLearner
    if (family == "quasibinomial") family <- "binomial"
    rhs <- setdiff(RhsVars(f), rownames(alias(f, data=d[subs,])$Complete))  #remove aliased columns from X - these can cause problems if they contain NAs and the user is expecting the column to be dropped
    new.subs <- apply(newdata[, rhs, drop=FALSE], 1, function (x) !any(is.na(x)))  #remove NA values from newdata - these will output to NA anyway and cause errors in SuperLearner
    
    m <- SuperLearner(Y=d[subs, LhsVars(f)], X=d[subs, rhs, drop=FALSE], SL.library=SL.library, verbose=FALSE, family=family, newX=newdata[new.subs, rhs, drop=FALSE])
    predicted.values <- rep(NA, nrow(newdata))
    predicted.values[new.subs] <- m$SL.predict
    if (max(predicted.values, na.rm=T) > 1 || min(predicted.values, na.rm=T) < 0) {
      stop("predicted.values > 1 or < 0")
    }
  }
  return(predicted.values)
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
InterventionMatch <- function(d, abar, Anodes, cur.node) {
  intervention.match <- XMatch(d, abar, Anodes, cur.node, all, default=TRUE)
  return(intervention.match)
}

# Determine which patients are uncensored
#return vector of [numDataRows x 1] I(C=1) from Cnodes[1] to the Cnode just before cur.node
IsUncensored <- function(d, Cnodes, cur.node) {
  uncensored <- XMatch(d, Xbar=1, Cnodes, cur.node, all, default=TRUE)
  return(uncensored)
}

# Determine which patients have died
# return list:
#    is.deterministic: vector of [numObservations x 1] - true if patient is already dead before cur.node or set by deterministic.Q.map
#    Q.value: vector of [which(is.deterministic) x 1] - value of Q
IsDeterministic <- function(d, Ynodes, cur.node, deterministic.Q.map, called.from.estimate.g) {
  is.deterministic <- XMatch(d, Xbar=1, Ynodes, cur.node, any, default=FALSE) #deterministic if any previous y node is 1
  Q.value <- rep(NA, nrow(d))
  Q.value[is.deterministic] <- 1
  
  for (i in seq_along(deterministic.Q.map)) {
    if (deterministic.Q.map[[i]]$node < cur.node) {
      if (!called.from.estimate.g || deterministic.Q.map[[i]]$implies.deterministic.g) {
        is.det <- deterministic.Q.map[[i]]$is.deterministic
        val <- deterministic.Q.map[[i]]$Q.value
        
        prev.Q <- Q.value[is.deterministic & is.det]
        temp <- numeric(nrow(d))
        temp[is.det] <- val
        new.Q <- temp[is.deterministic & is.det]
        
        if (! is.equal(prev.Q, new.Q)) {
          stop(paste("inconsistent deterministic Q at node:", names(d)[deterministic.Q.map[[i]]$node]))
        }
        Q.value[is.det] <- val
        is.deterministic <- is.deterministic | is.det
      }
    }
  } 
  Q.value <- Q.value[is.deterministic]
  if (any(is.na(c(is.deterministic, Q.value)))) stop("NA in is.deterministic or Q.value")
  return(list(is.deterministic=is.deterministic, Q.value=Q.value))
}

#Utility function called by IsUncensored, IsDeterministic - compares history in d within Xnodes prior to cur.node to Xbar
#
#any.all should be either the function 'any' or the function 'all'
#default: value to return if value is NA or there are no nodes before cur.node (TRUE for InterventionMatch and IsUncensored because NA indicates person matches intervention/is uncensored until they died; FALSE for IsDeterministic because NA indicates person was alive until they were censored)
XMatch <- function(d, Xbar, Xnodes, cur.node, any.all, default) {
  if (!any(Xnodes < cur.node)) return(rep(default, nrow(d)))
  last.Xnode.index <- which.max(Xnodes[Xnodes < cur.node])
  
  Xnodes.subset <- Xnodes[1:last.Xnode.index]
  if (identical(Xbar, 1)) {
    Xbar.subset <- 1 
  } else {
    Xbar.subset <- Xbar[, 1:last.Xnode.index]
  } 
  d.subset <- d[, Xnodes.subset, drop=FALSE]
  matches <- apply(d.subset == Xbar.subset, 1, any.all)
  matches[is.na(matches)] <- default
  return(matches)
}

# Determine which patients have an Anode value which is deterministic 
# For example, deterministic.acnode.map may be used to specify that once a patient starts treatment, they stay on treatment and this should be taken into consideration during estimation of G
IsDeterministicG <- function(d, cur.node, deterministic.acnode.map) {
  default <- list(is.deterministic=rep(FALSE, nrow(d)), prob1=NULL)
  if (is.null(deterministic.acnode.map)) return(default)
  index <- which(sapply(deterministic.acnode.map, function (x) x$node) == cur.node)
  if (length(index) == 0) return(default)
  if (length(index) > 1) stop("bad acnode.map")
  is.deterministic <- deterministic.acnode.map[[index]]$is.deterministic
  return(list(is.deterministic=is.deterministic, prob1=deterministic.acnode.map[[index]]$prob1))
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

#Set the Anodes of d to abar
SetA <- function(d, abar, Anodes, cur.node) {
  if (!any(Anodes < cur.node)) return(d)
  last.Anode.index  <- which.max(Anodes[Anodes < cur.node]) 
  d[, Anodes[1:last.Anode.index]] <- abar[, 1:last.Anode.index]
  return(d)
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
CheckInputs <- function(data, nodes, Qform, gform, gbounds, deterministic.acnode.map, SL.library, regimens, working.msm, summary.measures, summary.baseline.covariates, final.Ynodes, pooledMSM, stratify, weight.msm, deterministic.Q.map) {
  if (!all(is.null(GetLibrary(SL.library, "Q")), is.null(GetLibrary(SL.library, "g")))) library("SuperLearner")
  #each set of nodes should be sorted - otherwise causes confusion with gform, Qform, abar
  if (is.unsorted(nodes$A, strictly=TRUE)) stop("Anodes must be in increasing order")
  if (is.unsorted(nodes$C, strictly=TRUE)) stop("Cnodes must be in increasing order")
  if (is.unsorted(nodes$L, strictly=TRUE)) stop("Lnodes must be in increasing order")
  if (is.unsorted(nodes$Y, strictly=TRUE)) stop("Ynodes must be in increasing order")
  if (is.unsorted(final.Ynodes, strictly=TRUE)) stop("final.Ynodes must be in increasing order")

  if (! is.null(nodes$L)) {
    if (min(nodes$L) < min(nodes$AC)) stop("Lnodes are not allowed before A/C nodes. If you want to include baseline nodes, include them in data but not in Lnodes")
    if (max(nodes$L) > max(nodes$Y)) stop("Lnodes are not allowed after the final Y node")
  }
  
  all.nodes <- c(nodes$A, nodes$C, nodes$L, nodes$Y)
  if (length(all.nodes) > length(unique(all.nodes))) stop("A node cannot be listed in more than one of Anodes, Cnodes, Lnodes, Ynodes")
  if (is.null(nodes$Y)) stop("Ynodes cannot be null")
  if (is.null(nodes$AC)) stop("Anodes and Cnodes cannot both be null")
  
  if (is.character(gform)) {
    if (length(gform) != length(nodes$AC)) stop("length(gform) != length(c(Anodes, Cnodes))")
    for (i in 1:length(gform)) {
      if (LhsVars(gform[i]) != names(data)[nodes$AC[i]]) {
        cat("gform[", i,"] = ")
        print(gform[i])
        cat("names(data)[ACnodes[", i, "]] = ", names(data)[nodes$AC[i]], "\n")
        cat("note: ACnodes = sort(c(Anodes, Cnodes))")
        stop("The LHS variable of gform[i] should match names(data)[ACnodes[i]]")
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
    cat("L/Y nodes: ", names(data)[nodes$LY], "\n")
    cat("length(Qform) = ", length(Qform), "\n")
    stop("length of Qform is not equal to number of L/Y nodes")
  }
  for (i in 1:length(Qform)) {
    if (LhsVars(Qform[i]) != "Q.kplus1") stop("LHS of each Qform should be Q.kplus1")
    if (length(names(Qform[i])) == 0) stop("Each element of Qform must be named. The name must match the name of the corresponding L/Y node in data.")
    if (names(Qform[i]) != names(data)[nodes$LY[i]]) stop("The name of each element of Q must match the name of the corresponding L/Y node in data.")
  }
  
  if (length(gbounds) != 2) stop("gbounds should have length 2")
  if (! is.null(deterministic.acnode.map)) {
    for (i in 1:length(deterministic.acnode.map)) {
      map <- deterministic.acnode.map[[i]]
      if (! setequal(names(map), c("node", "is.deterministic", "prob1"))) stop("each element of deterministic.acnode.map should have names: node, is.deterministic, prob1")
      deterministic.acnode.map[[i]]$node <- NodeToIndex(data, deterministic.acnode.map[[i]]$node)
      if (! deterministic.acnode.map[[i]]$node %in% nodes$AC) stop("each node element of deterministic.acnode.map should be in Anodes or Cnodes")
      if (! length(map$prob1) %in% c(1, length(which(map$is.deterministic)))) stop(paste("length of deterministic.acnode.map[[", i, "]]$prob1  should be either 1 or length(which(deterministic.acnode.map[[", i, "]]$is.deterministic))", sep=""))
    }
  }
  
  if (! is.null(deterministic.Q.map)) {
    for (i in 1:length(deterministic.Q.map)) {
      map <- deterministic.Q.map[[i]]
      if (! setequal(names(map), c("node", "is.deterministic", "Q.value", "implies.deterministic.g"))) stop("each element of deterministic.Q.map should have names: node, is.deterministic, Q.value, implies.deterministic.g")
      if (! map$implies.deterministic.g %in% c(0,1)) stop(paste0("length of deterministic.Q.map[[", i, "]]$implies.deterministic.g should be binary"))
      deterministic.Q.map[[i]]$node <- NodeToIndex(data, map$node)
      if (! length(map$Q.value) %in% c(1, length(which(map$is.deterministic)))) stop(paste0("length of deterministic.Q.map[[", i, "]]$Q.value  should be either 1 or length(which(deterministic.Q.map[[", i, "]]$is.deterministic))"))
    }
  }

  CheckDeterministicACNodeMap(data, deterministic.acnode.map)
  
  if (! is.null(deterministic.Q.map)) {
    finalY <- data[, max(final.Ynodes)]
    for (i in nodes$LY[nodes$LY <= max(final.Ynodes)]) {
      deterministic.list <- IsDeterministic(data, nodes$Y, cur.node=i, deterministic.Q.map, called.from.estimate.g=FALSE)
      if (any((deterministic.list$Q.value %in% c(0,1)) & (deterministic.list$Q.value != finalY[deterministic.list$is.deterministic]))) stop("deterministic.Q.map is inconsistent with data")
    }
  }
  
  if (! all(unlist(data[, nodes$A]) %in% c(0, 1, NA))) stop("in data, all Anodes should be binary")
  #note: Cnodes are checked in ConvertCensoringNodes
  
  for (i in nodes$Y) {
    if (! all(data[, i] %in% c(0, 1, NA))) {
      index <- ! is.na(data[, i])
      if (any(index) && any(data[index, i] >= 1 | data[index, i] < 0)) stop("in data, all Ynodes should either be binary or in [0, 1)")
    }
    deterministic <- IsDeterministic(data, nodes$Y, cur.node=i, deterministic.Q.map=NULL, called.from.estimate.g=FALSE)$is.deterministic
    if (any(is.na(data[deterministic, i])) || ! all(data[deterministic, i] == 1)) stop("This function assumes that once a Ynode jumps to 1 (e.g. death), all subsequent Ynode values will also be 1. Your data does not follow this assumption.")    
  }
  
  if (! is.equal(dim(regimens)[1:2], c(nrow(data), length(nodes$A)))) stop("Problem with abar or regimens:\n   In ltmleMSM, regimens should have dimensions n x num.Anodes x num.regimens\n   In ltmle, abar should be a matrix with dimensions n x num.Anodes or a vector with length num.Anodes")
  num.regimens <- dim(regimens)[3]
  stopifnot(num.regimens == nrow(summary.measures))
  if (!all(regimens %in% c(0, 1, NA))) stop("all regimens should be binary")
  for (i in seq_along(nodes$A)) {
    cur.node <- nodes$A[i]
    uncensored <- IsUncensored(data, nodes$C, cur.node)
    deterministic <- IsDeterministic(data, nodes$Y, cur.node, deterministic.Q.map, called.from.estimate.g=TRUE)$is.deterministic
    if (any(is.na(regimens[uncensored & !deterministic, i, ]))) {
      stop("NA in regimens/abar not allowed (except after censoring/death)")
    }
  }
 
  if (! is.null(summary.baseline.covariates)) stop("summary.baseline.covariates is currently in development and is not yet supported - set summary.baseline.covariates to NULL")
  if (LhsVars(working.msm) != "Y") stop("the left hand side variable of working.msm should always be 'Y' [this may change in future releases]")
  if (! all(RhsVars(working.msm) %in% colnames(summary.measures))) stop("all right hand side variables in working.msm should be found in the column names of summary.measures")
  dynamic.regimens <- !all(duplicated(regimens)[2:nrow(data)])
  if (dynamic.regimens && weight.msm) stop("dynamic regimens are not currently supported with weight.msm=TRUE [under development]")
  
  return(list(deterministic.acnode.map=deterministic.acnode.map, deterministic.Q.map=deterministic.Q.map))
}

# Check deterministic.acnode.map for errors
CheckDeterministicACNodeMap <- function(data, deterministic.acnode.map) {
  ok <- TRUE
  if (is.null(deterministic.acnode.map)) return(NULL)
  for (i in 1:length(deterministic.acnode.map)) {
    isdet <- deterministic.acnode.map[[i]]$is.deterministic
    node <- deterministic.acnode.map[[i]]$node
    prob <- deterministic.acnode.map[[i]]$prob
    d <- data[isdet, node]
    index <- !is.na(d) & d == 1 & prob == 0
    if (is.character(node)) {
      name <- node
    } else {
      name <- names(data)[node]
    }
    if (any(index)) {
      cat("deterministic.acnode.map indicates Prob(", name, ") = 1 is 0, but these nodes are 1:\n", sep="")
      print(head(data.frame(data[isdet,node,drop=F][index,,drop=F], prob=prob[index])))
      ok <- FALSE
    }
    index <- !is.na(d) & d == 0 & prob == 1
    if (any(index)) {
      cat("deterministic.acnode.map indicates Prob(", name, ") = 1 is 1, but these nodes are 0:\n", sep="")
      print(head(data.frame(data[isdet,node,drop=F][index,,drop=F], prob=prob[index])))
      ok <- FALSE
    }
  }
  if (! ok) stop("Inconsistent data in deterministic.acnode.map - see diagnostic output above.")
  return(NULL)
}

# Get the default Q or g formula - each formula consists of all parent nodes except censoring and event nodes [also except A nodes if stratifying]
GetDefaultForm <- function(data, nodes, is.Qform, stratify) {
  if (is.Qform) {
    lhs <- rep("Q.kplus1", length(nodes$LY))
    node.set <- nodes$LY
  } else {
    lhs <- names(data)[nodes$AC]
    node.set <- nodes$AC
  }
  if (stratify) {
    stratify.nodes <- c(nodes$C, nodes$Y, nodes$A)
  } else {
    stratify.nodes <- c(nodes$C, nodes$Y)
  }
  form <- NULL
  for (i in 1:length(node.set)) {
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
      cat("L/Y nodes (after removing blocks)  : ", names(data)[new.LYnodes], "\n")
      cat("Qform names                        : ", names(Qform), "\n")
      cat(paste("The following nodes are not being considered as L/Y nodes because they are part of a block of L/Y nodes. They are being dropped from Qform:\n"), paste(names(Qform)[removed.Qform.index], "\n", collapse=" "))
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
NonpooledMSM <- function(data, Anodes, Cnodes, Lnodes, Ynodes, Qform, gform, gbounds, deterministic.acnode.map, stratify, SL.library, regimens, working.msm, final.Ynodes, summary.measures, weight.msm, gcomp, mhte.iptw, iptw.only, deterministic.Q.map) {  
  tmle.index <- ifelse(gcomp, "gcomp", "tmle")
  num.regimens <- dim(regimens)[3]
  num.final.Ynodes <- length(final.Ynodes)
  
  LYnodes <- sort(c(Lnodes, Ynodes))
  ACnodes <- sort(c(Anodes, Cnodes))
  num.ACnodes <- sum(ACnodes < max(final.Ynodes))
  tmle <- iptw <- weights <- matrix(nrow=num.regimens, ncol=num.final.Ynodes)
  IC <- IC.iptw <- array(dim=c(num.regimens, num.final.Ynodes, nrow(data)))
  cum.g <- array(dim=c(nrow(data), num.ACnodes, num.regimens))
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- final.Ynodes[j]
    if (is.matrix(gform)) {
      gform1 <- gform[, ACnodes < final.Ynode, drop=FALSE]
    } else {
      gform1 <- gform[ACnodes < final.Ynode]
    }
    is.duplicate <- duplicated(regimens[, Anodes < final.Ynode, , drop=FALSE], MARGIN=3)
    
    #It would be better to reuse g instead of calculating the same thing every time final.Ynode varies (note: g does need to be recalculated for each abar/regimen) - memoizing gets around this to some degree but it could be written better
    for (i in which(! is.duplicate)) {
      abar <- drop3(regimens[, , i, drop=F])
      abar <- abar[, Anodes < final.Ynode, drop=FALSE]
      det.ac.map <- TruncateDeterministicNodeMap(data, deterministic.acnode.map, final.Ynode)
      det.Q.map <- TruncateDeterministicNodeMap(data, deterministic.Q.map, final.Ynode)
      
      weights[i, j] <- ComputeGA(data[, 1:final.Ynode, drop=FALSE], Anodes[Anodes <= final.Ynode], Cnodes[Cnodes <= final.Ynode], abar, final.Ynode, weight.msm)
      
      if (weights[i, j] > 0) {
        result <- ltmle(data[, 1:final.Ynode, drop=FALSE], Anodes[Anodes <= final.Ynode], Cnodes[Cnodes <= final.Ynode], Lnodes[Lnodes <= final.Ynode], Ynodes[Ynodes <= final.Ynode], Qform[LYnodes <= final.Ynode], gform1, abar, gbounds, det.ac.map, stratify, SL.library, estimate.time=FALSE, gcomp=gcomp, mhte.iptw=mhte.iptw, iptw.only=iptw.only, deterministic.Q.map=det.Q.map)
        tmle[i, j] <- result$estimates[tmle.index]
        iptw[i, j] <- min(1, result$estimates["iptw"])
        IC[i, j, ] <- result$IC[[tmle.index]]
        IC.iptw[i, j, ] <- result$IC[["iptw"]]
      }
      if (j == num.final.Ynodes) {
        if (weights[i, j] == 0) {
          #we didn't calculate cum.g because weight was 0 but we need to return it
          result <- ltmle(data[, 1:final.Ynode, drop=FALSE], Anodes[Anodes <= final.Ynode], Cnodes[Cnodes <= final.Ynode], Lnodes[Lnodes <= final.Ynode], Ynodes[Ynodes <= final.Ynode], Qform[LYnodes <= final.Ynode], gform1, abar, gbounds, det.ac.map, stratify, SL.library, estimate.time=FALSE, gcomp=gcomp, mhte.iptw=mhte.iptw, iptw.only=TRUE, deterministic.Q.map=deterministic.Q.map)
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
          if (any(is.na(C))) {print("NA in C"); browser()}
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
          if (abs(det(C)) < 1e-18) {
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

# Truncate the parts of deterministic.acnode.map after final.Ynode
# Most of the code should work ok without this but it helps to pass CheckInputs
TruncateDeterministicNodeMap <- function(data, deterministic.node.map, final.Ynode) {
  det.map <- deterministic.node.map
  if (! is.null(det.map)) {
    index <- NULL
    for (k in 1:length(det.map)) {
      node.index <- NodeToIndex(data, det.map[[k]]$node)
      if (node.index <= final.Ynode) {
        index <- c(index, k)
      }
    }
    det.map <- det.map[index]
    if (length(det.map) == 0) det.map <- NULL
  }
  return(det.map)
}

# Convert censoring nodes stored as factors into binaries (the ltmle code was originally written for binaries but factors confusion in inputs)
ConvertCensoringNodes <- function(data, Cnodes) {
  CensoringToBinary <- function(x) {
    if (! all(levels(x) %in% c("censored", "uncensored"))) {
      stop("all levels of data[, Cnodes] should be in censored, uncensored (NA should not be a level)")
    }
    b <- rep(NA_integer_, length(x))
    b[x == "censored"] <- 0L
    b[x == "uncensored"] <- 1L
    return(b)
  }
  
  error.msg <- "in data, all Cnodes should be factors with two levels, 'censored' and 'uncensored' \n (binary is also accepted, where 0=censored, 1=uncensored, but is not recommended)"
  for (i in Cnodes) {
    col <- data[, i]
    if (is.numeric(col)) {
      if (! all(col %in% c(0, 1, NA))) stop(error.msg)
    } else if (is.factor(col)) {
      data[, i] <- CensoringToBinary(col)
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
