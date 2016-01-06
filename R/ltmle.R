# Longitudinal TMLE to estimate an intervention-specific mean outcome or marginal structural model

# General code flow:
#  ltmle -> CreateInputs -> LtmleFromInputs -> LtmleMSMFromInputs -> ...
#  ltmleMSM -> CreateInputs -> LtmleMSMFromInputs -> ...

#longitudinal targeted maximum liklihood estimation for E[Y_a]
#' @export
ltmle <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=NULL, Qform=NULL, gform=NULL, 
                  abar, rule=NULL, gbounds=c(0.01, 1), Yrange=NULL, deterministic.g.function=NULL, stratify=FALSE, 
                  SL.library=NULL, estimate.time=TRUE, gcomp=FALSE, 
                  iptw.only=FALSE, deterministic.Q.function=NULL, IC.variance.only=FALSE, observation.weights=NULL) {
  msm.inputs <- GetMSMInputsForLtmle(data, abar, rule, gform)
  inputs <- CreateInputs(data=data, Anodes=Anodes, Cnodes=Cnodes, Lnodes=Lnodes, Ynodes=Ynodes, survivalOutcome=survivalOutcome, Qform=Qform, gform=msm.inputs$gform, Yrange=Yrange, gbounds=gbounds, deterministic.g.function=deterministic.g.function, SL.library=SL.library, regimes=msm.inputs$regimes, working.msm=msm.inputs$working.msm, summary.measures=msm.inputs$summary.measures, final.Ynodes=msm.inputs$final.Ynodes, stratify=stratify, msm.weights=msm.inputs$msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, IC.variance.only=IC.variance.only, observation.weights=observation.weights, msm.family=msm.inputs$msm.family) 
  result <- LtmleFromInputs(inputs)
  result$call <- match.call()
  return(result)
}

RegimesFromAbar <- function(data, abar, rule) {
  if (!is.null(rule)) {
    if (!(missing(abar) || is.null(abar))) stop("'abar' should not be specified when using a 'rule' function")
    abar <- t(apply(data, 1, rule))
  }
  if (is.vector(abar)) {
    abar <- matrix(rep(abar, each=nrow(data)), nrow=nrow(data))
  } else if (is.null(abar)) {
    abar <- matrix(nrow=nrow(data), ncol=0)
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
  if (is.numeric(gform)) {
    stopifnot(is.matrix(gform))
    dim(gform) <- c(nrow(gform), ncol(gform), 1)
  }
  msm.inputs <- list(regimes=regimes, working.msm=working.msm, summary.measures=summary.measures, gform=gform, final.Ynodes=NULL, msm.weights=msm.weights, msm.family=binomial()) 
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

#longitudinal targeted maximum likelihood estimation for a marginal structural model
#' @export 
ltmleMSM <- function(data, Anodes, Cnodes=NULL, Lnodes=NULL, Ynodes, survivalOutcome=NULL, Qform=NULL, gform=NULL, gbounds=c(0.01, 1), Yrange=NULL, deterministic.g.function=NULL, SL.library=NULL, regimes, working.msm, summary.measures, final.Ynodes=NULL, stratify=FALSE, msm.weights="empirical", estimate.time=TRUE, gcomp=FALSE, iptw.only=FALSE, deterministic.Q.function=NULL, memoize=TRUE, IC.variance.only=FALSE, observation.weights=NULL, msm.family=binomial()) {
  if (memoize && requireNamespace("memoise")) {
    glm.ltmle.memoized <- memoise::memoize(glm.ltmle)
  }
  
  inputs <- CreateInputs(data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, estimate.time, gcomp, iptw.only, deterministic.Q.function, IC.variance.only, observation.weights, msm.family)
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
CreateInputs <- function(data, Anodes, Cnodes, Lnodes, Ynodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, estimate.time, gcomp, iptw.only, deterministic.Q.function, IC.variance.only, observation.weights, msm.family) {
  if (is.list(regimes)) {
    if (!all(do.call(c, lapply(regimes, is.function)))) stop("If 'regimes' is a list, then all elements should be functions.")
    regimes <- aperm(simplify2array(lapply(regimes, function(rule) apply(data, 1, rule)), higher=TRUE), c(2, 1, 3)) 
  }
  if (!(is.null(regimes) || dim(regimes) != 3)) {
    stop("regimes must be an array with 3 dimensions (unless Anodes is NULL, in which case regimes can be NULL)")
  }
  if (is.null(regimes) || dim(regimes)[3]==0) {
    if (length(Anodes) != 0) {
      stop("regimes must not be NULL (or have dim(regimes)[3]==0) unless Anodes is also NULL")
    }
    regimes <- array(dim=c(nrow(data), 0, 1))
  }
  if (is.logical(regimes)) {
    regimes <- regimes * 1
    message("abar or regimes was passed as logical and was converted to numeric")
  }
  nodes <- CreateNodes(data, Anodes, Cnodes, Lnodes, Ynodes)
  Qform <- CreateLYNodes(data, nodes, check.Qform=TRUE, Qform=Qform)$Qform
  data <- ConvertCensoringNodes(data, Cnodes, has.deterministic.functions=!is.null(deterministic.g.function) && is.null(deterministic.Q.function))
  if (is.null(final.Ynodes)) {
    final.Ynodes <- max(nodes$Y)
  } else {
    final.Ynodes <- NodeToIndex(data, final.Ynodes)
  }
  
  #Using get to avoid the "no visible binding for global variable" note in R CMD check
  if (identical(SL.library, 'default')) SL.library <- get("Default.SL.Library")
  SL.library.Q <- GetLibrary(SL.library, "Q")
  SL.library.g <- GetLibrary(SL.library, "g")
    
  if (is.null(summary.measures)) {
    summary.measures <- matrix(nrow=dim(regimes)[3], ncol=0)
  }
  if (length(dim(summary.measures)) == 2) {
    num.final.Ynodes <- length(final.Ynodes)
    summary.measures <- array(repmat(summary.measures, m=1, n=num.final.Ynodes), dim=c(nrow(summary.measures), ncol(summary.measures), num.final.Ynodes), dimnames=list(rownames(summary.measures), colnames(summary.measures), NULL))
  }
  if (is.null(observation.weights)) observation.weights <- rep(1, nrow(data))
  
  #error checking (also get value for survivalOutcome if NULL)
  check.results <- CheckInputs(data, nodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, deterministic.Q.function, observation.weights, msm.family) 
  survivalOutcome <- check.results$survivalOutcome
  
  if (!isTRUE(attr(data, "called.from.estimate.variance", exact=TRUE))) { 
    data <- CleanData(data, nodes, deterministic.Q.function, survivalOutcome)
  }
  untransformed.data <- data
  transform.list <- TransformOutcomes(data, nodes, Yrange)
  data <- transform.list$data
  transformOutcome <- transform.list$transformOutcome
  binaryOutcome <- check.results$binaryOutcome
  
  if (is.null(Qform)) Qform <- GetDefaultForm(data, nodes, is.Qform=TRUE, stratify, survivalOutcome, showMessage=TRUE)
  if (is.null(gform)) gform <- GetDefaultForm(data, nodes, is.Qform=FALSE, stratify, survivalOutcome, showMessage=TRUE)
  
  inputs <- list(data=data, untransformed.data=untransformed.data, nodes=nodes, survivalOutcome=survivalOutcome, Qform=Qform, gform=gform, gbounds=gbounds, Yrange=Yrange, deterministic.g.function=deterministic.g.function, SL.library.Q=SL.library.Q, SL.library.g=SL.library.g, regimes=regimes, working.msm=working.msm, summary.measures=summary.measures, final.Ynodes=final.Ynodes, stratify=stratify, msm.weights=msm.weights, estimate.time=estimate.time, gcomp=gcomp, iptw.only=iptw.only, deterministic.Q.function=deterministic.Q.function, binaryOutcome=binaryOutcome, transformOutcome=transformOutcome, IC.variance.only=IC.variance.only, observation.weights=observation.weights, msm.family=msm.family)
  class(inputs) <- "ltmleInputs"
  return(inputs)
}

# Prevent misspelled argument after $ from returning NULL
`$.ltmleInputs` <- function(x, name) {
  if (! (name %in% names(x))) stop(paste(name, "is not an element of x"))
  value <- x[[name]]
  if (identical(value, "ltmleInputs-NULL")) {
    return(NULL)
  } else {
    return(value)
  }
}

# Prevent misspelled argument after $<- from adding new element
`$<-.ltmleInputs` <- function(x, name, value) {
  if (! (name %in% names(x))) stop(paste(name, "is not an element of x"))
  if (is.null(value)) {
    value <- "ltmleInputs-NULL"
  }
  x[[name]] <- value
  return(x)
}

# Loop over final Ynodes, run main calculations
MainCalcs <- function(inputs) {
  #if (inputs$iptw.only) { #could use something like this to save time, but we want to pass all final ynodes to CalcIPTW below
  #  inputs$final.Ynodes <- inputs$final.Ynodes[length(inputs$final.Ynodes)]
  #}
  # Several functions in the pooled version are only written to accept main terms MSM
  # Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S1 + S3 + S4" where 
  # S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
  main.terms <- ConvertToMainTerms(inputs$data, inputs$working.msm, inputs$summary.measures, inputs$nodes)
  inputs$working.msm <- main.terms$msm
  combined.summary.measures <- main.terms$summary.measures   
  baseline.column.names <- main.terms$baseline.column.names
  num.final.Ynodes <- length(inputs$final.Ynodes)
  num.betas <- length(main.terms$beta.names)
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes       note: num.measures is summary measures and baseline covariates, converted to main terms
  n <- nrow(inputs$data)
  num.regimes <- dim(inputs$regimes)[3]
  Qstar <- g.ratio <- array(dim=c(n, num.regimes, num.final.Ynodes))
  msm.weights <- GetMsmWeights(inputs) #n x num.regimes x num.final.Ynodes
  new.var.y <- array(dim=c(num.betas, num.betas, num.final.Ynodes))
  IC <- matrix(0, n, num.betas)
  #store IC for each final Ynode, compare var(IC) to sum(var(IC.ynode))
  IC.y <- array(dim=c(n, num.betas, num.final.Ynodes))
  for (j in 1:num.final.Ynodes) {
    #It would be better to reuse g instead of calculating the same thing every time final.Ynode varies (note: g does need to be recalculated for each abar/regime) - memoizing gets around this to some degree but it could be written better
    fixed.tmle <- FixedTimeTMLE(SubsetInputs(inputs, final.Ynode=inputs$final.Ynodes[j]), drop3(msm.weights[, , j, drop=FALSE]), dropn(combined.summary.measures[, , , j, drop=FALSE], n=4), baseline.column.names)
    IC <- IC + fixed.tmle$IC
    IC.y[, , j] <- fixed.tmle$IC
    Qstar[, , j] <- fixed.tmle$Qstar # n x num.regimes
    new.var.y[, , j] <- fixed.tmle$est.var 
    g.ratio[, , j] <- fixed.tmle$g.ratio
  }

  iptw <- CalcIPTW(inputs$data, inputs$nodes, inputs$working.msm, inputs$regimes, combined.summary.measures, inputs$final.Ynodes, fixed.tmle$cum.g, msm.weights, inputs$observation.weights, inputs$msm.family)
  names(iptw$beta) <- main.terms$beta.names

  if (inputs$iptw.only) {
    beta <- rep(NA, length(iptw$beta))
    fitted.msm <- NULL
    variance.estimate <- NULL
  } else {
    fitted.msm <- FitPooledMSM(inputs$working.msm, Qstar, combined.summary.measures, msm.weights * inputs$observation.weights, inputs$msm.family)
    if (isTRUE(attr(inputs$data, "called.from.estimate.variance", exact=TRUE))) { 
      #variance estimate is not needed, this avoids some warnings
      IC <- matrix(NA, n, num.betas)
      C.old <- matrix(NA, num.betas, num.betas)
    } else {
      IC <- FinalizeIC(IC, combined.summary.measures, Qstar, fitted.msm$m.beta, msm.weights, g.ratio, inputs$observation.weights) #n x num.betas
      C.old <- NormalizeIC(IC, combined.summary.measures, fitted.msm$m.beta, msm.weights, g.ratio = array(1, dim=c(n, num.regimes, num.final.Ynodes)), inputs$observation.weights, inputs$msm.family) #C without using g.ratio (setting g.ratio to 1)
    }
 
    if (inputs$IC.variance.only) {   
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
      C <- NormalizeIC(IC, combined.summary.measures, fitted.msm$m.beta, msm.weights, g.ratio, inputs$observation.weights, inputs$msm.family)
      variance.estimate <- safe.solve(C) %*% new.var %*% t(safe.solve(C))
    }
    
    IC <- t(safe.solve(C.old, t(IC))) #IC %*% solve(C) 
    beta <- coef(fitted.msm$m)
    names(beta) <- main.terms$beta.names
  }
  return(list(IC=IC, msm=fitted.msm$m, beta=beta, cum.g=fixed.tmle$cum.g, cum.g.unbounded=fixed.tmle$cum.g.unbounded, fit=fixed.tmle$fit, variance.estimate=variance.estimate, beta.iptw=iptw$beta, IC.iptw=iptw$IC, Qstar=Qstar, g.ratio=g.ratio)) #note: only returns cum.g and fit for the last final.Ynode
}


CalcIPTW <- function(data, nodes, working.msm, regimes, combined.summary.measures, final.Ynodes, cum.g, msm.weights, observation.weights, msm.family) {
  if (isTRUE(attr(data, "called.from.estimate.variance", exact=TRUE))) { 
    return(list(beta=NA, IC=matrix(NA, 1, 1)))
  }
  n <- nrow(data)
  num.regimes <- dim(regimes)[3]
  num.final.Ynodes <- length(final.Ynodes)
  Y.vec <- X.mat <- weight.vec <- NULL
  save.xy <- list()
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- final.Ynodes[j]
    for (i in 1:num.regimes) {
      abar <- drop3(regimes[, , i, drop=F])
      abar <- abar[, nodes$A < final.Ynode, drop=FALSE]
      uncensored <- IsUncensored(data, nodes$C, final.Ynode)
      intervention.match <- InterventionMatch(data, abar, nodes$A, final.Ynode)  
      index <- uncensored & intervention.match
      col.index <- which.max(nodes$AC[nodes$AC < final.Ynode]) 
      
      Y <- data[index, final.Ynode]
      g <- cum.g[index, col.index, i] 
      X <- combined.summary.measures[index, , i, j]
      if (is.vector(X)) { #if only one summary.measure or sum(index)==1, X is dropped to vector
        X <- matrix(X, nrow=sum(index), ncol=ncol(combined.summary.measures))
      }
      weight <- msm.weights[index, i, j] * observation.weights[index] / g
      weight[msm.weights[index, i, j] == 0 | observation.weights[index] == 0] <- 0 #avoid problems where weight and g are both 0
      
      save.xy[[length(save.xy) + 1]] <- list(X=X, Y=Y, weight=weight, index=index)
      Y.vec <- c(Y.vec, Y)
      X.mat <- rbind(X.mat, X)
      weight.vec <- c(weight.vec, weight) 
    }
  }
  colnames(X.mat) <- colnames(combined.summary.measures)
  
  if (nrow(X.mat) == 0) {
    #this happens if there are no rows uncensored and intervention.match
    warning("no rows uncensored and matching regimes/abar - IPTW returns NA")
    num.beta <- ncol(combined.summary.measures)
    return(list(beta=rep(NA, num.beta), IC=matrix(nrow=n, ncol=num.beta)))
  }
  SuppressGivenWarnings(m.glm <- glm(formula(working.msm), family=msm.family, data=data.frame(Y=Y.vec, X.mat, weight.vec), weights=scale(weight.vec, center=FALSE)), "non-integer #successes in a binomial glm!") #note: scale weights because there were rare problems where large weights caused convergence problems
  beta <- coef(m.glm)
  IC <- matrix(0, nrow=n, ncol=length(beta))  #n x num.betas
  m.beta <- array(dim=c(n, num.regimes, num.final.Ynodes)) 
  cnt <- 0
  for (j in 1:num.final.Ynodes) {
    final.Ynode <- final.Ynodes[j]
    for (i in 1:num.regimes) {
      newdata <- data.frame(combined.summary.measures[, , i, j])
      colnames(newdata) <- colnames(combined.summary.measures) #needed if only one summary measure
      SuppressGivenWarnings(m.beta[, i, j] <- predict(m.glm, newdata=newdata, type="response"), "prediction from a rank-deficient fit may be misleading")
      
      cnt <- cnt + 1
      XY.list <- save.xy[[cnt]]
      IC[XY.list$index, ] <- IC[XY.list$index, ] + XY.list$weight * XY.list$X * (XY.list$Y - m.beta[XY.list$index, i, j]) #recycles weight, Y, m.beta
    }
  }

  C <- NormalizeIC(IC, combined.summary.measures, m.beta, msm.weights, g.ratio=array(1, dim=c(n, num.regimes, num.final.Ynodes)), observation.weights=observation.weights, msm.family=msm.family) 
  normalized.IC <- t(safe.solve(C, t(IC)))  
  return(list(beta=beta, IC=normalized.IC))
}

# ltmleMSM for a single final.Ynode
FixedTimeTMLE <- function(inputs, msm.weights, combined.summary.measures, baseline.column.names) {
  inputs$summary.measures <- NULL #just to make sure it isn't used - should only use combined.summary.measures 
  #combined.summary.measures: n x num.measures x num.regimes   (num.measures=num.summary.measures + num.baseline.covariates)
  
  data <- inputs$data
  nodes <- inputs$nodes
  
  num.regimes <- dim(inputs$regimes)[3]
  n <- nrow(data)
  num.betas <- ncol(combined.summary.measures)
  tmle <- rep(NA, num.regimes)
  IC <- matrix(0, nrow=n, ncol=num.betas)
  cum.g <- cum.g.unbounded <- prob.A.is.1 <- array(0, dim=c(n, length(nodes$AC), num.regimes))
  cum.g.meanL <- cum.g.meanL.unbounded <- array(0, dim=c(n, length(nodes$AC), num.regimes, length(nodes$LY)-1))
  fit.g <- vector("list", num.regimes)
  for (i in 1:num.regimes) {
    if (all(msm.weights[, i] == 0)) {
      g.list <- list(fit=list("no g fit because regime weight is 0"))
    } else {
      # estimate each g factor, and cumulative probabilities
      g.list <- EstimateG(inputs, regime.index=i)
      cum.g[, , i] <- g.list$cum.g
      cum.g.unbounded[, , i] <- g.list$cum.g.unbounded
      cum.g.meanL[, , i, ] <- g.list$cum.g.meanL
      cum.g.meanL.unbounded[, , i, ] <- g.list$cum.g.meanL.unbounded
      prob.A.is.1[, , i] <- g.list$prob.A.is.1
    } 
    fit.g[[i]] <- g.list$fit
  }
  if (inputs$iptw.only) return(list(cum.g=cum.g, cum.g.unbounded=cum.g.unbounded, IC=NA, Qstar=NA, est.var=NA, g.ratio=NA))
  
  est.var <- matrix(0, num.betas, num.betas)  
  logitQ <- matrix(nrow=n, ncol=num.regimes)
  regimes.with.positive.weight <- which(apply(msm.weights > 0, 2, any))
  if (length(regimes.with.positive.weight) == 0) stop("All regimes have weight 0 (one possible reason is that msm.weights=NULL and no data rows match any of the regimes and are uncensored)")
  fit.Qstar <- vector("list", length(nodes$LY))
  names(fit.Qstar) <- names(data)[nodes$LY]
  fit.Q <- vector("list", length(regimes.with.positive.weight))
  for (i in regimes.with.positive.weight) fit.Q[[i]] <- fit.Qstar
  Qstar.kplus1 <- matrix(data[, max(nodes$Y)], nrow=n, ncol=num.regimes)
  
  for (LYnode.index in length(nodes$LY):1) {
    cur.node <- nodes$LY[LYnode.index]
    deterministic.list.origdata <- IsDeterministic(data, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=FALSE, inputs$survivalOutcome)
    uncensored <- IsUncensored(data, nodes$C, cur.node)
    intervention.match <- subs <- matrix(nrow=n, ncol=num.regimes)
    for (i in regimes.with.positive.weight) {
      abar <- GetABar(inputs$regimes, i)
      intervention.match[, i] <- InterventionMatch(data, abar=abar, nodes$A, cur.node)  
      newdata <- SetA(data, abar=abar, nodes, cur.node)
      deterministic.list.newdata <- IsDeterministic(newdata, cur.node, inputs$deterministic.Q.function, nodes, called.from.estimate.g=FALSE, inputs$survivalOutcome)
      if (inputs$stratify) {
        subs[, i] <- uncensored & intervention.match[, i] & !deterministic.list.origdata$is.deterministic
      } else {
        subs[, i] <- uncensored & !deterministic.list.origdata$is.deterministic
      }
      if (any(subs[, i])) {
        Q.est <- Estimate(inputs$Qform[LYnode.index], data=data.frame(data, Q.kplus1=Qstar.kplus1[, i]), family="quasibinomial", newdata=newdata, subs=subs[, i], SL.library=inputs$SL.library.Q, type="link", nodes=nodes, observation.weights=inputs$observation.weights)
        logitQ[, i] <- Q.est$predicted.values
      } else {
        if (! all(deterministic.list.newdata$is.deterministic)) {
          msg <- paste0("ltmle failed trying to estimate ", inputs$Qform[LYnode.index], " because there are no observations that are\nuncensored", ifelse(inputs$stratify, ", follow abar,", ""), " and are not set deterministically due to death or deterministic.Q.function\n")
          stop(msg)
        }
        Q.est <- list(fit="no estimation of Q at this node because all rows are set deterministically")
      }
      fit.Q[[i]][[LYnode.index]] <- Q.est$fit 
    }
    if (all(deterministic.list.newdata$is.deterministic)) {
      #no updating needed if all rows are set deterministically
      Qstar <- matrix(deterministic.list.newdata$Q, nrow=n, ncol=num.regimes)
      if (max(abs(Qstar.kplus1 - Qstar)) > 1e-8) {
        #if Qstar.kplus1 != Qstar when all deterministic score equation will not be solved
        stop("inconsistency in deterministic data - all rows are set deterministically but the deterministically set values are not equal to Qstar.kplus1") 
      }
      Qstar.est <- list(fit="no updating at this node because all rows are set deterministically")
    } else {
      ACnode.index  <- which.max(nodes$AC[nodes$AC < cur.node])
      update.list <- UpdateQ(Qstar.kplus1, logitQ, combined.summary.measures, subs, cum.g[, ACnode.index, ], inputs$working.msm, uncensored, intervention.match, msm.weights, inputs$gcomp, inputs$observation.weights)
      Qstar <- update.list$Qstar
      Qstar[deterministic.list.newdata$is.deterministic, ] <- deterministic.list.newdata$Q
      curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
      curIC.relative.error <- abs(colSums(curIC)) / apply(abs(combined.summary.measures), 2, mean)
      if (any(curIC.relative.error > 0.001) && !inputs$gcomp) {
        fix.score.list <- FixScoreEquation(Qstar.kplus1, update.list$h.g.ratio, uncensored, intervention.match, deterministic.list.newdata, update.list$off, update.list$X, regimes.with.positive.weight)
        Qstar <- fix.score.list$Qstar
        curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
        update.list$fit <- fix.score.list$fit      
      }
      est.var <- est.var + EstimateVariance(inputs, combined.summary.measures, regimes.with.positive.weight, uncensored, deterministic.list.newdata, Qstar, Qstar.kplus1, cur.node, msm.weights, LYnode.index, ACnode.index, cum.g, prob.A.is.1, baseline.column.names, cum.g.meanL, cum.g.unbounded, cum.g.meanL.unbounded, inputs$observation.weights, is.last.LYnode=(LYnode.index==length(nodes$LY)))
    }
    IC <- IC + curIC 
    Qstar.kplus1 <- Qstar
    fit.Qstar[[LYnode.index]] <- update.list$fit
  }
  g.ratio <- CalcGUnboundedToBoundedRatio(inputs, cum.g, cum.g.meanL, cum.g.unbounded, cum.g.meanL.unbounded)
  #tmle <- colMeans(Qstar)
  return(list(IC=IC, Qstar=Qstar, cum.g=cum.g, cum.g.unbounded=cum.g.unbounded, g.ratio=g.ratio, est.var=est.var, fit=list(g=fit.g, Q=fit.Q, Qstar=fit.Qstar))) 
}

EstimateVariance <- function(inputs, combined.summary.measures, regimes.with.positive.weight, uncensored, deterministic.list.newdata, Qstar, Qstar.kplus1, cur.node, msm.weights, LYnode.index, ACnode.index, cum.g, prob.A.is.1, baseline.column.names, cum.g.meanL, cum.g.unbounded, cum.g.meanL.unbounded, observation.weights, is.last.LYnode) {
  if (inputs$IC.variance.only) return(NA)
  if (!inputs$binaryOutcome) stop("IC.variance.only=FALSE not currently compatible with non binary outcomes")
  if (!is.null(inputs$deterministic.Q.function)) stop("IC.variance.only=FALSE not currently compatible with deterministic.Q.function")
  if (inputs$gcomp) stop("IC.variance.only=FALSE not currently compatible with gcomp")
  if (inputs$stratify) stop("IC.variance.only=FALSE not currently compatible with stratify=TRUE")
  TmleOfVariance <- function(Z, Z.meanL) {
    if (all(is.na(Z))) stop("all Z are NA in EstimateVariance")
    #length(Z.meanL) == 0 --> point treatment case, return mean(Z); all Z can be zero when summary.measure is 0
    if (length(Z.meanL) == 0 || all(Z==0 | is.na(Z))) {
      Qstar <- Scale(Z, 0, 1)
      return(list(EZd1 = mean(Z, na.rm=T), Qstar = Qstar))
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
    var.tmle <- ltmle(sparsity.data, Anodes=temp.nodes$A, Cnodes=temp.nodes$C, Lnodes=temp.nodes$L, Ynodes=temp.nodes$Y, survivalOutcome=FALSE, Qform=Qform, gform=drop3(prob.A.is.1[, 1:ACnode.index, d1, drop=FALSE]), abar=drop3(inputs$regimes[, nodes$A <= cur.node, d1, drop=FALSE]), gbounds=inputs$gbounds, stratify=stratify, estimate.time=FALSE, deterministic.Q.function=det.q.function, IC.variance.only=TRUE, observation.weights=observation.weights) 
  
    EZd1 <- var.tmle$estimates["tmle"] * diff(range(Z, na.rm=T)) + min(Z, na.rm=T)
    return(list(EZd1 = EZd1, Qstar=var.tmle$Qstar))
  }
  
  EqualRegimesIndex <- function(dd1, dd2) {
    #index of each observation where regime d1 matches regime d2 up to cur.node
    return(apply(drop3(inputs$regimes[, inputs$nodes$A <= cur.node, dd1, drop=F]) == drop3(inputs$regimes[, inputs$nodes$A <= cur.node, dd2, drop=F]), 1, all)) 
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
  
  nodes <- inputs$nodes  
  num.regimes <- dim(inputs$regimes)[3]
  num.betas <- ncol(combined.summary.measures)
  n <- nrow(inputs$data)
  
  #used in ltmle call below
  if (inputs$survivalOutcome) {
    det.q.function <- function(data, current.node, nodes, called.from.estimate.g) {
      if (!any(nodes$Y < current.node)) return(NULL)
      prev.Y <- data[, nodes$Y[nodes$Y < current.node], drop=F]
      prev.Y[is.na(prev.Y)] <- 0
      is.deterministic <- apply(prev.Y == 1, 1, any)
      Q.value <- data[is.deterministic, max(nodes$Y)] #this is 0 before scaling but may be nonzero after scaling
      return(list(is.deterministic=is.deterministic, Q.value=Q.value))   
    }
  } else {
    det.q.function <- NULL
  }
  static.treatment <- IsStaticTreatment()
  variance.estimate <- matrix(0, num.betas, num.betas)
  alive <- !deterministic.list.newdata$is.deterministic
  Sigma <- array(dim=c(n, num.regimes, num.regimes))
  #fixme - Sigma is not exactly symmetric due to SetA on d1, but could probably save time with approximate symmetry
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
        Q.data <- inputs$data[alive, 1:cur.node, drop=F]
        resid.sq <- (Qstar.kplus1[alive, d1] - Qstar[alive, d1]) * (Qstar.kplus1[alive, d2] - Qstar[alive, d2]) 
        resid.sq.range <- range(resid.sq, na.rm=T)
        if (diff(resid.sq.range) > 0) {
          Q.data[, cur.node] <- (resid.sq - resid.sq.range[1]) / diff(resid.sq.range)
          names(Q.data)[cur.node]  <- "Q.kplus1" #ugly - to match with Qform
          m <- glm(formula = inputs$Qform[LYnode.index], family = "quasibinomial", data = Q.data, control=glm.control(trace=FALSE, maxit=1000)) 
          Q.newdata <- SetA(data = Q.data, abar = GetABar(regimes = inputs$regimes, d1)[alive, , drop=F], nodes = nodes, cur.node = cur.node)
          SuppressGivenWarnings(Q.resid.sq.pred <- predict(m, newdata = Q.newdata, type = "response"), "prediction from a rank-deficient fit may be misleading")
          Sigma[alive, d1, d2] <- Q.resid.sq.pred * diff(resid.sq.range) + resid.sq.range[1]
        } else {
          resid.sq.value <- min(resid.sq, na.rm = T) #all values are the same, just get one non-NA
          Sigma[alive, d1, d2] <- resid.sq.value
        }
        Sigma[!alive, d1, d2] <- 0
      }
    }
  }
  

  if (static.treatment) {
    for (d1 in regimes.with.positive.weight) {   
      #Z.without.h1h1 <- Sigma[, d1, d1] / cum.g[, ACnode.index, d1] #without h1*h1'
      #Z.without.h1h1.meanL <- 1 / cum.g.meanL[, ACnode.index, d1, ]
      Z.without.sum.meas <- Sigma[, d1, d1] / cum.g[, ACnode.index, d1] * cum.g.unbounded[, ACnode.index, d1] / cum.g[, ACnode.index, d1] * msm.weights[, d1]^2 * observation.weights^2
      Z.without.sum.meas.meanL <- 1 / cum.g.meanL[, ACnode.index, d1, ] * cum.g.meanL.unbounded[, ACnode.index, d1, ] / cum.g.meanL[, ACnode.index, d1, ] * msm.weights[, d1]^2 * observation.weights^2
      var.tmle <- TmleOfVariance(Z.without.sum.meas, Z.without.sum.meas.meanL)
      no.V <- all(combined.summary.measures[1, , d1] == combined.summary.measures[, , d1]) 
      if (no.V) {
        #h1 <- combined.summary.measures[1, , d1] * msm.weights[1, d1]  #num.betas x 1
        #variance.estimate <- variance.estimate + (h1 %*% t(h1)) * var.tmle$EZd1
        variance.estimate <- variance.estimate + (combined.summary.measures[1, , d1] %*% t(combined.summary.measures[1, , d1])) * var.tmle$EZd1
      } else {
        #has V (so combined.summary.measures varies)
        baseline.msm <- "Qstar ~ 1"
        if (length(baseline.column.names) > 0) {
          baseline.msm <- paste(baseline.msm, "+", paste(baseline.column.names, collapse=" + "), "+", paste0("I(", baseline.column.names, "^2)", collapse=" + "))
        }            
        m <- glm(baseline.msm, family = "quasibinomial", data=data.frame(Qstar=var.tmle$Qstar, inputs$data[, baseline.column.names, drop=FALSE]), control=glm.control(trace=FALSE, maxit=1000)) 
        SuppressGivenWarnings(pred.Qstar <- predict(m, type = "response") * diff(range(Z.without.sum.meas, na.rm=T)) + min(Z.without.sum.meas, na.rm=T), "prediction from a rank-deficient fit may be misleading")  #n x 1
        variance.estimate.sum <- matrix(0, num.betas, num.betas)
        for (i in 1:n) {
          #h1 <- combined.summary.measures[i, , d1] * msm.weights[i, d1]
          #variance.estimate.sum <- variance.estimate.sum + (h1 %*% t(h1)) * pred.Qstar[i]
          variance.estimate.sum <- variance.estimate.sum + (combined.summary.measures[i, , d1] %*% t(combined.summary.measures[i, , d1])) * pred.Qstar[i]
        }
        variance.estimate <- variance.estimate + variance.estimate.sum / n
      }
    }
  } else {
    for (beta.index2 in 1:num.betas) {
      for (d1 in regimes.with.positive.weight) {          
        Z.base <- rep(0, n)  #Z without h1(d1, V, beta.index1)
        Z.base.meanL <- matrix(0, n, dim(cum.g.meanL)[4])
        for (d2 in regimes.with.positive.weight) {
          equal.regimes.index <- EqualRegimesIndex(d1, d2) #index of each observation where regime d1 matches regime d2 
          h1 <- combined.summary.measures[, beta.index2, d2] * msm.weights[, d2]
          Z.base[equal.regimes.index] <- Z.base[equal.regimes.index] + h1[equal.regimes.index] * Sigma[equal.regimes.index, d1, d2] / cum.g[equal.regimes.index, ACnode.index, d1] * observation.weights[equal.regimes.index] #this is equivalent to using cum.g.unbounded in the denominator and multiplying by phi=cum.g.unbounded/cum.g.bounded  
          Z.base.meanL[equal.regimes.index, ] <- Z.base.meanL[equal.regimes.index, ] + h1[equal.regimes.index] * 1 / cum.g.meanL[equal.regimes.index, ACnode.index, d1, ] * observation.weights[equal.regimes.index] #recycles
        }        
        for (beta.index1 in 1:num.betas) {  
          if (beta.index1 >= beta.index2) {
            Z <- combined.summary.measures[, beta.index1, d1] * msm.weights[, d1] * cum.g.unbounded[, ACnode.index, d1] / cum.g[, ACnode.index, d1] * observation.weights * Z.base 
            Z.meanL <- combined.summary.measures[, beta.index1, d1] * msm.weights[, d1] * cum.g.meanL.unbounded[, ACnode.index, d1, ] / cum.g.meanL[, ACnode.index, d1, ] * observation.weights * Z.base.meanL
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
    orig.variance.estimate <- variance.estimate
    try.result <- try({
      near.pd <- Matrix::nearPD(variance.estimate) #may cause an error if variance.estimate is negative definite
      variance.estimate <- as.matrix(near.pd$mat)
    })
    if (inherits(try.result, "try-error") || !near.pd$converged || any((abs(orig.variance.estimate - variance.estimate) / orig.variance.estimate) > 0.1)) {
      warning("Covariance matrix from EstimateVariance not positive definite, unable to compute standard errors. You may want to try IC.variance.only=TRUE.")
      variance.estimate <- matrix(nrow=num.betas, ncol=num.betas)
    }
  }
  return(variance.estimate)
}

CalcGUnboundedToBoundedRatio <- function(inputs, cum.g, cum.g.meanL, cum.g.unbounded, cum.g.meanL.unbounded) {
  n <- dim(cum.g)[1]
  num.AC.nodes <- dim(cum.g)[2]
  num.regimes <- dim(cum.g)[3]
  if (inputs$IC.variance.only) return(matrix(1, n, num.regimes))
  if (! any(is.na(cum.g))) return(AsMatrix(cum.g.unbounded[, num.AC.nodes, ] / cum.g[, num.AC.nodes, ]))
  #cum.g is NA after censoring - for censored observations use cum.g.meanL  
  #If censored at node j, set all nodes > j to meanl. 
  #[,,k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L (na.rm=T)
  g.ratio <- matrix(NA, n, num.regimes)
  for (i in 1:num.regimes) {
    g.ratio.temp <- cbind(cum.g.meanL.unbounded[, num.AC.nodes, i, ] / cum.g.meanL[, num.AC.nodes, i, ], cum.g.unbounded[, num.AC.nodes, i] / cum.g[, num.AC.nodes, i])
    index <- max.col(!is.na(g.ratio.temp), "last")
    g.ratio[, i] <- g.ratio.temp[sub2ind(1:n, col = index, num.rows = n)]
  }
  if (any(is.na(g.ratio))) stop("NA in g.ratio")
  return(g.ratio)
}

# remove any information in ltmleInputs after final.Ynode
SubsetInputs <- function(inputs, final.Ynode) {
  if (is.numeric(inputs$gform)) {
    stopifnot(length(dim(inputs$gform)) == 3)
    inputs$gform <- inputs$gform[, inputs$nodes$AC < final.Ynode, , drop=FALSE]
  } else {
    inputs$gform <- inputs$gform[inputs$nodes$AC < final.Ynode]
  }
  inputs$Qform <- inputs$Qform[inputs$nodes$LY <= final.Ynode]
  inputs$data <- inputs$data[, 1:final.Ynode, drop=FALSE]
  inputs$untransformed.data <- inputs$untransformed.data[, 1:final.Ynode, drop=FALSE]
  inputs$regimes <- inputs$regimes[, inputs$nodes$A <= final.Ynode, , drop=FALSE]
  inputs$summary.measures <- drop3(inputs$summary.measures[, , inputs$final.Ynodes == final.Ynode, drop=FALSE])
  inputs$final.Ynodes <- inputs$final.Ynodes[inputs$final.Ynodes <= final.Ynode]
  inputs$nodes <- lapply(inputs$nodes, function (x) x[x <= final.Ynode])
  
  return(inputs)
}

# Fit the MSM
FitPooledMSM <- function(working.msm, Qstar, combined.summary.measures, msm.weights, msm.family) {
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
  
  m <- glm(as.formula(working.msm), data=data.frame(Y, X), family=msm.family, weights=scale(weight.vec, center = FALSE), na.action=na.exclude, control=glm.control(maxit=1000)) 
  SuppressGivenWarnings(m.beta <- predict(m, type="response"), "prediction from a rank-deficient fit may be misleading")
  dim(m.beta) <- dim(Qstar)
  return(list(m=m, m.beta=m.beta))
}

#final step in calculating TMLE influence curve
FinalizeIC <- function(IC, combined.summary.measures, Qstar, m.beta, msm.weights, g.ratio, observation.weights) {
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
  if (any(abs(colSums(finalIC)) > 0.001 )) {
    msg <- capture.output(cat("final IC problem", colSums(finalIC)))
    warning(paste(msg, collapse="\n"))
  }
  IC <- IC + finalIC
  return(IC) 
}

# Normalize the influence curve matrix
NormalizeIC <- function(IC, combined.summary.measures, m.beta, msm.weights, g.ratio, observation.weights, msm.family) {    
  #combined.summary.measures: n x num.measures x num.regimes x num.final.Ynodes   (num.measures=num.summary.measures + num.baseline.covariates)
  #g.ratio = g.unbounded / g.bounded : n x num.regimes x num.final.Ynodes
  n <- nrow(IC)
  num.betas <- ncol(IC)
  num.regimes <- dim(combined.summary.measures)[3]
  num.final.Ynodes <- dim(combined.summary.measures)[4]
  
  if (is.equal(msm.family, binomial())) {
    is.binomial <- TRUE
    if (max(m.beta) > (1 - 1e-6) || min(m.beta) < 1e-6) {
      warning("All predicted probabilities are all very close to 0 or 1. Unable to compute standard errors.")
      return(matrix(NA, nrow=num.betas, ncol=num.betas))
    }
  } else if (is.equal(msm.family, gaussian())) {
    is.binomial <- FALSE
  } else {
    stop("unexpected msm.family")  #if other families are possible, need to adjust scalar below
  }
  C <- array(0, dim=c(num.betas, num.betas, n))
  for (j in 1:num.final.Ynodes) {
    for (i in 1:num.regimes) {
      positive.msm.weights <- which(msm.weights[, i, j] > 0)
      tempC <- array(0, dim=c(num.betas, num.betas, n))
      for (k in positive.msm.weights) {
        h <- matrix(combined.summary.measures[k, , i, j], ncol=1) * g.ratio[k, i, j]
        if (is.binomial) {
          scalar <- m.beta[k, i, j] * (1 - m.beta[k, i, j])
        } else {
          scalar <- 1
        }
        tempC[, , k] <- h %*% t(h) * scalar * observation.weights[k] * msm.weights[k, i, j] 
      }
      if (any(is.na(tempC))) stop("NA in tempC")
      C <- C + tempC
    }
  }
  C <- apply(C, c(1, 2), mean)
  if (rcond(C) < 1e-12) {
    C <- matrix(NA, nrow=num.betas, ncol=num.betas)
    warning("rcond(C) near 0, standard errors not available")
  } else {
    normalized.IC <- t(safe.solve(C, t(IC))) #IC %*% solve(C) 
    if (!any(abs(colSums(IC)) > 0.001) && !any(is.na(normalized.IC)) && any(abs(colSums(normalized.IC)) > 0.001)) {
      msg <- capture.output({
        cat("normalized IC problem", colSums(normalized.IC), "\n")
        cat("inv(C) = \n")
        print(safe.solve(C))
      })
      warning(paste(msg, collapse="\n"))
    }
  }
  return(C)
}

# Get a single regime from the regimes array
GetABar <- function(regimes, i) {
  abar <- AsMatrix(regimes[, , i]) #if there's only 1 Anode, make sure abar comes back as a matrix
  return(abar)
}

# Targeting step - update the initial fit of Q using clever covariates
UpdateQ <- function(Qstar.kplus1, logitQ, combined.summary.measures, subs, cum.g, working.msm, uncensored, intervention.match, msm.weights, gcomp, observation.weights) { 
  #logitQ, Qstar.kplus1: n x num.regimes
  #cum.g: n x num.regimes (already indexed for this node)
  #subs: n x num.regimes
  #uncensored: n x 1
  #intervention.match: n x num.regimes
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
  
  stacked.summary.measures <- apply(combined.summary.measures, 2, rbind)
  
  weight.vec <- uncensored * observation.weights * as.vector(intervention.match) / as.vector(cum.g) * as.vector(msm.weights) #recycles uncensored and observation.weights
  f <- as.formula(paste(working.msm, "+ offset(off)"))
  data.temp <- data.frame(Y, stacked.summary.measures, off)
  newdata <- data.temp  
  if (gcomp) {
    Qstar <- plogis(logitQ)
    m <- "no Qstar fit because gcomp=TRUE (so no updating step)"
  } else {
    ctrl <- glm.control(trace=FALSE, maxit=1000)
    
    if (any(subs & weight.vec>0)) {
      SuppressGivenWarnings(m <- glm(f, data=data.temp, subset=as.vector(subs) & (weight.vec > 0), family="quasibinomial", weights=scale(weight.vec, center=FALSE), control=ctrl), GetWarningsToSuppress(TRUE)) #this should include the indicators; only include weight.vec>0 because others have NAs
      SuppressGivenWarnings(Qstar <- matrix(predict(m, newdata=newdata, type="response"), nrow=nrow(logitQ)), GetWarningsToSuppress(TRUE))  #this should NOT include the indicators  #note: could also use plogis(off + X %*% coef(m)) [but this has problems with NAs in coef(m)?]
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
  if (ncol(stacked.summary.measures) != dim(h.g.ratio)[3]) stop("only works if working.msm has only main terms")
  return(list(Qstar=Qstar, h.g.ratio=h.g.ratio, X=stacked.summary.measures, off=off, fit=m)) 
}

# Sometimes GLM doesn't converge and the updating step of TMLE doesn't solve the score equation (sum of TMLE influence curve not equal to zero). This function attempts to solve the score equation directly using various optimizers. [Note: this was needed more often before we made changes to the updating algorithm.]
FixScoreEquation <- function(Qstar.kplus1, h.g.ratio, uncensored, intervention.match, deterministic.list, off, X, regimes.with.positive.weight) {
  CalcScore <- function(e) {
    Qstar <- QstarFromE(e)
    ICtemp <- CalcIC(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight)
    return(sum(colSums(ICtemp) ^ 2)) #each column has to add to zero
  }
  
  QstarFromE <- function(e) {
    Qstar <- plogis(off + X %*% e) #X: n x (num.summary.measures + num.baseline.covariates) (which should be num.beta);  e: num.beta x 1 
    dim(Qstar) <- dim(Qstar.kplus1)
    Qstar[deterministic.list$is.deterministic, ] <- deterministic.list$Q
    return(Qstar)
  }
  
  FindMin <- function(minimizer) {
    num.tries <- 30 #increase this to solve more problems at the cost of longer run time
    init.e <- numeric(num.betas) #first try an initial estimate of epsilon=0
    for (i in 1:num.tries) {
      if (minimizer == "nlminb") {
        m <- nlminb(start=init.e, objective=CalcScore, control=list(abs.tol=max.objective, eval.max=500, iter.max=500, x.tol=1e-14, rel.tol=1e-14))
        e <- m$par
        obj.val <- m$objective
      } else if (minimizer == "optim") {
        if (num.betas == 1) {
          m <- optimize(f = CalcScore, interval=c(-1, 1) * 10^i)
          e <- m$minimum
          obj.val <- m$objective
        } else {
          m <- optim(par=init.e, fn=CalcScore, control=list(abstol=max.objective, reltol=1e-14, maxit=2000))
          e <- m$par
          obj.val <- m$value
        }
      } else if (minimizer == "nlm") {
        m <- nlm(f=CalcScore, p=init.e)
        e <- m$estimate
        obj.val <- m$minimum
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
  
  for (offset.lbound in c(1e-8, 0.0001, 0.001, 0.01)) {
    off <- Bound(off, qlogis(c(offset.lbound, 1-offset.lbound)))
    l <- FindMin("nlminb")
    if (! l$solved) l <- FindMin("optim")
    if (! l$solved) l <- FindMin("nlm")
    if (l$solved) break
  }
  if (! l$solved) {
    stop("all minimizers failed")
  }
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
  small.inputs$data <- inputs$data[sample.index, ]
  small.inputs$regimes <- inputs$regimes[sample.index, , , drop=F]
  small.inputs$observation.weights <- inputs$observation.weights[sample.index]
  if (is.numeric(inputs$gform)) small.inputs$gform <- inputs$gform[sample.index, , , drop=F]
  if (length(dim(inputs$msm.weights)) == 3) small.inputs$msm.weights <- inputs$msm.weights[sample.index, , , drop=F]
  start.time <- Sys.time()
  #try.result <- MainCalcs(small.inputs)
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

# Get standard error, p-value, and confidence interval for one ltmle object 
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
#' @export 
summary.ltmleEffectMeasures <- function(object, estimator=ifelse(object$gcomp, "gcomp", "tmle"), ...) {
  info <- GetSummaryLtmleMSMInfo(object, estimator)
  beta <- info$estimate
  IC <- info$IC
  y0 <- plogis(beta[1])
  y1 <- plogis(beta[1] + beta[2])
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
    measures.variance.estimate <- NULL #if IC.variance.only=T
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
  
  ans <- list(call=object$call, effect.measures=measures.max, variance.estimate.ratio=info$variance.estimate.ratio, transformOutcome=object$transformOutcome, estimator=estimator)
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
#' @export 
print.ltmleEffectMeasures <- function(x, ...) {
  PrintCall(x$call)
  cat("Use summary(...) to get estimates, standard errors, p-values, and confidence intervals for treatment EYd, control EYd, additive effect, relative risk, and odds ratio.\n")
  invisible(x)
}

# Print method for summary.ltmleEffectMeasures
#' @export 
print.summary.ltmleEffectMeasures <- function(x, ...) {
  cat("Estimator: ", x$estimator, "\n")
  if (x$estimator=="gcomp") {cat("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.\n")}
  PrintCall(x$call)
  lapply(x$effect.measures, PrintSummary)
  CheckVarianceEstimateRatio(x)
  if (x$transformOutcome) {
    Yrange <- attr(x$transformOutcome, "Yrange")
    cat("NOTE: All parameters are based on the transformed outcome ( Y -", min(Yrange),
        ")/(", max(Yrange),"-", min(Yrange),")")
  }
  invisible(x)
}

# Print a warning message if the TMLE based variance estimate is much greater than the IC based variance estimate 
CheckVarianceEstimateRatio <- function(summary.obj) {
  if (any(is.na(summary.obj$variance.estimate.ratio))) {
    warning("Unable to compute standard errors.")
    return(NULL)
  }
  if (any(summary.obj$variance.estimate.ratio > 100)) {
    warning.msg <- paste0("max(TMLE based variance estimate / IC based variance estimate) = ", floor(max(summary.obj$variance.estimate.ratio)), ".\nWhen this ratio is greater than 100, both variance estimates are less likely to be accurate.")
    warning(warning.msg)
  }
}

# Print method for ltmleMSM
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
      stop("unexpected x$long.name") 
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

# Parametric estimation of each g-factor
EstimateG <- function(inputs, regime.index) {
  abar <- GetABar(inputs$regimes, regime.index)
  gmat <- prob.A.is.1 <- matrix(NaN, nrow=nrow(inputs$data), ncol=length(inputs$nodes$AC))
  gmat.meanL <- prob.A.is.1.meanL <- cum.g.meanL <- cum.g.meanL.unbounded <- array(NaN, dim=c(nrow(inputs$data), length(inputs$nodes$AC), length(inputs$nodes$LY) - 1))
  uncensored <- rep(TRUE, nrow(inputs$data))
  fit <- vector("list", length(inputs$nodes$AC))
  names(fit) <- names(inputs$data)[inputs$nodes$AC]
  if (!inputs$IC.variance.only) { 
    abar.meanL <- abar
    for (i in sseq(1, length(inputs$nodes$A))) {
      if (any(is.na(abar.meanL[, i]))) { #abar.meanL needs to be nonNA
        abar.meanL[is.na(abar.meanL[, i]), i] <- Mode(abar.meanL[, i], na.rm = TRUE)
      }
    }
  }
  for (i in 1:length(inputs$nodes$AC)) {
    cur.node <- inputs$nodes$AC[i]
    uncensored <- IsUncensored(inputs$data, inputs$nodes$C, cur.node)
    newdata <- SetA(inputs$data, abar, inputs$nodes, cur.node)
    if (!inputs$IC.variance.only) newdata.meanL <- SetA(inputs$data, abar.meanL, inputs$nodes, cur.node)
    deterministic.origdata <- IsDeterministic(inputs$data, cur.node, inputs$deterministic.Q.function, inputs$nodes, called.from.estimate.g=TRUE, inputs$survivalOutcome)$is.deterministic #deterministic due to death or Q.function
    deterministic.newdata <- IsDeterministic(newdata, cur.node, inputs$deterministic.Q.function, inputs$nodes, called.from.estimate.g=TRUE, inputs$survivalOutcome)$is.deterministic #deterministic due to death or Q.function - using data modified so A = abar
    if (is.numeric(inputs$gform)) {
      if (!inputs$IC.variance.only) stop("IC.variance.only=FALSE not currently compatible with numeric gform")
      prob.A.is.1[, i] <- inputs$gform[, i, regime.index]  #if gform is numeric, it's a matrix of prob.A.is.1
      g.est <- list(fit="gform passed as numeric, so no estimation took place")
    } else {
      deterministic.g.list.origdata <- IsDeterministicG(inputs$data, cur.node, inputs$deterministic.g.function, inputs$nodes, using.newdata=F) #deterministic due to acnode map - using original data
      deterministic.g.list.newdata <- IsDeterministicG(newdata, cur.node, inputs$deterministic.g.function, inputs$nodes, using.newdata=T) #deterministic due to acnode map - using data modified so A = abar
      deterministic.g.origdata <- deterministic.g.list.origdata$is.deterministic
      if (inputs$stratify) {
        intervention.match <- InterventionMatch(inputs$data, abar, inputs$nodes$A, inputs$nodes$AC[i]) 
        subs <- uncensored & intervention.match & !deterministic.origdata & !deterministic.g.origdata
      } else {
        subs <- uncensored & !deterministic.origdata & !deterministic.g.origdata
      }
      
      if (all(deterministic.g.list.newdata$is.deterministic | deterministic.newdata)) {
        # all rows are set deterministically, no need to estimate
        g.est <- list(fit="all rows are set deterministically, no estimation at this node")
        #fill in prob.A.is.1 below
      } else {
        # not all rows are set deterministically
        if (any(subs)) {
          g.est <- Estimate(inputs$gform[i], data=inputs$data, subs=subs, family="quasibinomial", newdata=newdata, SL.library=inputs$SL.library.g, type="response", nodes=inputs$nodes, observation.weights=inputs$observation.weights)
          prob.A.is.1[, i] <- g.est$predicted.values
          if (!inputs$IC.variance.only) {           
            #n x numACnodes x (numLYnodes - 1)
            #[,,k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L (na.rm=T)
            prob.A.is.1.meanL[, i, ] <- PredictProbAMeanL(newdata.meanL, inputs$nodes, subs, g.est$fit, SL.XY=g.est$SL.XY)
          }
        } else {
          msg <- paste0("ltmle failed trying to estimate ", inputs$gform[i], " because there are no observations that are\nuncensored", ifelse(inputs$stratify, ", follow abar,", ""), " and are not set deterministically due to death or deterministic.g.function or deterministic.Q.function\n")
          stop(msg)
        }
      }
      prob.A.is.1[deterministic.g.list.newdata$is.deterministic, i] <- prob.A.is.1.meanL[deterministic.g.list.newdata$is.deterministic, i, ] <- deterministic.g.list.newdata$prob1
    } 
    #prob.A.is.1 is prob(a=1), gmat is prob(a=abar)
    #cur.abar can be NA after censoring/death if treatment is dynamic
    if (cur.node %in% inputs$nodes$A) {
      cur.abar <- abar[, inputs$nodes$A == cur.node]
      if (!inputs$IC.variance.only) cur.abar.meanL <- abar.meanL[, inputs$nodes$A == cur.node]
    } else {
      cur.abar <- cur.abar.meanL <- rep(1, nrow(inputs$data))  #if this is a cnode, abar is always 1 (uncensored)
    }
    gmat[, i] <- CalcGVec(prob.A.is.1[, i], cur.abar, deterministic.newdata)
    if (!inputs$IC.variance.only) gmat.meanL[, i, ] <- apply(AsMatrix(prob.A.is.1.meanL[, i, ]), 2, CalcGVec, cur.abar.meanL, deterministic.newdata)
    if (any(is.na(gmat[uncensored, i]))) stop("Error - NA in g. g should only be NA after censoring. If you passed numeric gform, make sure there are no NA values except after censoring. Otherwise something has gone wrong.")
    fit[[i]] <- g.est$fit
  }
  cum.g.unbounded <- CalcCumG(gmat, c(0, 1))
  cum.g <- CalcCumG(gmat, inputs$gbounds)
  if (!inputs$IC.variance.only) {
    for (i in sseq(1, dim(gmat.meanL)[3])) {
      cum.g.meanL[, , i] <- CalcCumG(drop3(gmat.meanL[, , i, drop=F]), inputs$gbounds)
      cum.g.meanL.unbounded[, , i] <- CalcCumG(drop3(gmat.meanL[, , i, drop=F]), c(0,1)) 
    }
  }
  return(list(cum.g=cum.g, cum.g.unbounded=cum.g.unbounded, cum.g.meanL=cum.g.meanL, fit=fit, prob.A.is.1=prob.A.is.1, cum.g.meanL.unbounded=cum.g.meanL.unbounded))
}

PredictProbAMeanL <- function(data, nodes, subs, fit, SL.XY) {
  #Predict the probability that A=1 if L and Y nodes are set to their mean (or median) values
  
  #probAis1.meanL is n x num.LYnodes - 1
  #probAis1.meanL[, k] is prob.A.is.1 with all L and Y nodes after and including LYnodes[k] set to mean of L
  
  #somewhat inefficient - for W A.1 L.2 A.2 L.3 A.3 Y, does P(A.1=1) setting L.3 to mean and then L.2 and L.3 to mean, but none of these can be used in P(A.1=1) because they're after A.1
  
  #A is already set to abar in data
  probAis1.meanL <- matrix(NaN, nrow(data), length(nodes$LY) - 1)
  if (ncol(probAis1.meanL) == 0) return(probAis1.meanL)
  all.LY.nodes <- sort(union(nodes$L, nodes$Y)) #not the same as nodes$LY, which removes blocks
  newdata <- data
  LYindex <- length(nodes$LY)
  for (i in length(all.LY.nodes):1) { 
    regression.node <- all.LY.nodes[i]
    L <- data[subs, regression.node]
    if (is.numeric(L) && !IsBinary(L)) {
      meanL <- mean(L, na.rm = TRUE)
    } else {
      meanL <- Mode(L, na.rm = TRUE) #for factors and binaries
    }
    newdata[, regression.node] <- meanL
    if (regression.node %in% nodes$LY[1:length(nodes$LY)-1]) {
      LYindex <- LYindex - 1
      SuppressGivenWarnings({
        if ("SuperLearner" %in% class(fit)) {
          probAis1.meanL[, LYindex] <- predict(fit, newdata = cbind(newdata, added.constant=1), X=SL.XY$X, Y=SL.XY$Y, onlySL=TRUE)$pred
        } else {
          probAis1.meanL[, LYindex] <- predict(fit, newdata = newdata, type = "response")
        }
      }, "prediction from a rank-deficient fit may be misleading")
    }
  }
  if (any(is.na(probAis1.meanL[, 1]))) stop("NA in probAis1.meanL[, 1]")
  return(probAis1.meanL)
}



CalcGVec <- function(prob.A.is.1, cur.abar, deterministic.newdata) {
  g <- rep(NA, length(prob.A.is.1))
  g[!is.na(cur.abar) & cur.abar == 1] <- prob.A.is.1[!is.na(cur.abar) & cur.abar == 1]
  g[!is.na(cur.abar) & cur.abar == 0] <- 1 - prob.A.is.1[!is.na(cur.abar) & cur.abar == 0]    
  g[deterministic.newdata] <- 1  #a=abar deterministically after death or other deterministic Q
  return(g)
}

# Truncate values within supplied bounds
Bound <- function(x, bounds) {
  stopifnot(length(bounds) == 2 && !any(is.na(bounds)))
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
Estimate <- function(form, data, subs, family, newdata, SL.library, type, nodes, observation.weights) {
  stopifnot(type %in% c("link", "response"))
  if (form == "IDENTITY") {
    predicted.values <- data[, "Q.kplus1"]
    if (type == "link") {
      stopifnot(family %in% c("binomial", "quasibinomial"))
      predicted.values <- qlogis(Bound(predicted.values, bounds=c(0.0001, 0.9999)))
    }
    m <- "no fit because form == IDENTITY"
    return(list(predicted.values=predicted.values, fit=m))
  }
  data <- ConvertCensoringNodesToBinary(data, nodes$C) #convert factors to binaries for compatability with glm and some SL libraries
  f <- as.formula(form)
  if (any(is.na(data[subs, LhsVars(f)]))) stop("NA in Estimate")
  observation.weights <- observation.weights[subs]
  if (is.null(SL.library) || length(RhsVars(f)) == 0) { #in a formula like "Y ~ 1", call glm
    #estimate using GLM
    if (sum(subs) > 1) {
      SuppressGivenWarnings({
        m <- get.stack("glm.ltmle.memoized", mode="function", ifnotfound=glm.ltmle)(f, data=data[subs, all.vars(f), drop=F], observation.weights=observation.weights, family=family, control=glm.control(trace=FALSE, maxit=1000)) #there's probably a better way to do this
        predicted.values <- predict(m, newdata=newdata, type=type)
      }, GetWarningsToSuppress())
    } else {
      #glm breaks when sum(subs) == 1
      predicted.values <- rep(data[subs, LhsVars(f)], nrow(newdata))
      m <- "fit not returned because there was only 1 observation to fit"
    }
    SL.XY <- NULL
  } else {
    #estimate using SuperLearner
    if (family == "quasibinomial") family <- "binomial"
    
    rhs <- setdiff(RhsVars(f), rownames(alias(f, data=data[subs,])$Complete))  #remove aliased columns from X - these can cause problems if they contain NAs and the user is expecting the column to be dropped
    new.subs <- apply(newdata[, rhs, drop=FALSE], 1, function (x) !any(is.na(x)))  #remove NA values from newdata - these will output to NA anyway and cause errors in SuperLearner
    Y <- data[subs, LhsVars(f)]
    X <- data[subs, rhs, drop=FALSE]
    newX <- newdata[new.subs, rhs, drop=FALSE]
    if (ncol(X) == 1) {
      #SuperLearner crashes if there are screening algorithms and only one column - add a constant
      X <- cbind(X, added.constant=1)
      newX <- cbind(newX, added.constant=1)
    }
    try.result <- try({
      SuppressGivenWarnings(m <- SuperLearner::SuperLearner(Y=Y, X=X, SL.library=SL.library, verbose=FALSE, family=family, newX=newX, obsWeights=observation.weights), c("non-integer #successes in a binomial glm!", "prediction from a rank-deficient fit may be misleading")) 
    })
    SL.XY <- list(X=X, Y=Y)
    GetSLStopMsg <- function(Y) ifelse(all(Y %in% c(0, 1, NA)), "", "\n Note that many SuperLeaner libraries crash when called with continuous dependent variables, as in the case of initial Q regressions with continuous Y or subsequent Q regressions even if Y is binary.")
    if (inherits(try.result, "try-error")) {
      stop(paste("\n\nError occured during call to SuperLearner:\n", form, GetSLStopMsg(Y), "\n The error reported is:\n", try.result))
    }
    if (all(is.na(m$SL.predict))) {
      stop(paste("\n\nSuperLearner returned all NAs during regression:\n", form, GetSLStopMsg(Y)))
    }
    predicted.values <- rep(NA, nrow(newdata))
    predicted.values[new.subs] <- m$SL.predict
    if (max(predicted.values, na.rm=T) > 1 || min(predicted.values, na.rm=T) < 0) {
      msg <- paste("SuperLearner returned predicted.values > 1 or < 0: [min, max] = [", min(predicted.values, na.rm=T), ",", max(predicted.values, na.rm=T), "]. Bounding to [0,1]")
      warning(msg)
      predicted.values <- Bound(predicted.values, bounds=c(0, 1))
    }
    if (type == "link") {
      stopifnot(family == "binomial")
      predicted.values <- qlogis(Bound(predicted.values, bounds=c(0.0001, 0.9999)))
    }
  }
  return(list(predicted.values=predicted.values, fit=m, SL.XY=SL.XY))
}

# This is here for memoizing
glm.ltmle <- function(f, data, family, control, observation.weights) {
  #note: observation.weights is in Estimate (the environment of f)
  return(glm(f, data=data.frame(data[, all.vars(f), drop=F], observation.weights), family=family, control=control, weights=scale(observation.weights, center=FALSE)))
}

# Calculate bounded cumulative G
CalcCumG <- function(g, gbounds) {
  cum.g <- AsMatrix(Bound(t(apply(g, 1, cumprod)), gbounds)) #AsMatrix to fix problems where apply returns a vector
  return(cum.g)
}

# Determine which patients are following specified treatment regime (abar)
#return vector of [numObservations x 1] I(A==abar) from Anodes[1] to the Anode just before cur.node
# note: if calling from outside ltmle:::, cur.node needs to be the node index, not a string!
InterventionMatch <- function(data, abar, Anodes, cur.node) {
  intervention.match <- XMatch(data, abar, Anodes, cur.node, all, default=TRUE)
  return(intervention.match)
}

# Determine which patients are uncensored
#return vector of [numDataRows x 1] I(C=uncensored) from Cnodes[1] to the Cnode just before cur.node
# note: if calling from outside ltmle:::, cur.node needs to be the node index, not a string!
IsUncensored <- function(data, Cnodes, cur.node) {
  if (! all(sapply(data[, Cnodes], is.factor))) stop("something has gone wrong in ltmle:::IsUncensored - all Cnodes should have been converted to factors")
  uncensored <- XMatch(data, Xbar="uncensored", Cnodes, cur.node, all, default=TRUE)
  return(uncensored)
}

# Determine which patients have died or have Q set deterministically by user function before cur.node
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
  if (any(inconsistent.rows)) stop(paste("At node:",names(data)[cur.node], "deterministic.Q.function is inconsistent with data - Q.value is either 0 or 1 but this does not match the final Y node value\nCheck data rows:", paste(head(rownames(data)[det.list$is.deterministic][inconsistent.rows]), collapse=" ")))
  
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
IsDeterministicG <- function(data, cur.node, deterministic.g.function, nodes, using.newdata) {
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
CalcIC <- function(Qstar.kplus1, Qstar, h.g.ratio, uncensored, intervention.match, regimes.with.positive.weight) {
  n <- nrow(Qstar)
  num.regimes <- ncol(Qstar)
  num.betas <- dim(h.g.ratio)[3] #h.g.ratio: n x num.regimes x num.betas
  
  IC <- matrix(0, nrow=n, ncol=num.betas)
  for (i in regimes.with.positive.weight) {
    index <- uncensored & intervention.match[, i]
    if (any(h.g.ratio[index, i, ] != 0)) {
      regimeIC <- matrix(0, nrow=n, ncol=num.betas)
      regimeIC[index, ] <- (Qstar.kplus1[index, i] - Qstar[index, i]) * h.g.ratio[index, i, ]
      IC <- IC + regimeIC
    }
  }
  return(IC)
}

#Set the Anodes of d to abar and Cnodes to uncensored (up to and including cur.node - cur.node itself is included for consistency checking in DeterministicG)
SetA <- function(data, abar, nodes, cur.node) {
  Anode.index <- nodes$A <= cur.node
  data[, nodes$A[Anode.index]] <- abar[, Anode.index]
  
  Cnode.index <- nodes$C <= cur.node
  data[, nodes$C[Cnode.index]] <- factor(rep("uncensored", nrow(data))) #recycled
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
CheckInputs <- function(data, nodes, survivalOutcome, Qform, gform, gbounds, Yrange, deterministic.g.function, SL.library, regimes, working.msm, summary.measures, final.Ynodes, stratify, msm.weights, deterministic.Q.function, observation.weights, msm.family) {
  stopifnot(length(dim(regimes)) == 3)
  num.regimes <- dim(regimes)[3]
  if (!all(is.null(GetLibrary(SL.library, "Q")), is.null(GetLibrary(SL.library, "g")))) {
    if (!requireNamespace("SuperLearner")) stop("SuperLearner package is required if SL.library is not NULL")
  } 
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
  if (min(nodes$Y) < min(nodes$AC)) stop("Ynodes are not currently allowed before A/C nodes.")
  
  all.nodes <- c(nodes$A, nodes$C, nodes$L, nodes$Y)
  if (length(all.nodes) > length(unique(all.nodes))) stop("A node cannot be listed in more than one of Anodes, Cnodes, Lnodes, Ynodes")
  if (is.null(nodes$Y)) stop("Ynodes cannot be null")
  if (is.null(nodes$AC)) stop("Anodes and Cnodes cannot both be null")
  
  if (min(all.nodes) < ncol(data)) {
    if (!all((min(all.nodes):ncol(data)) %in% all.nodes)) {
      stop("All nodes after the first of A-, C-, L-, or Ynodes must be in A-, C-, L-, or Ynodes")
    }
  }
  for (reserved.name in c("observation.weights", "Q.kplus1", "Qstar")) {
    #these might cause conflicts
    if (reserved.name %in% names(data)) stop(paste(reserved.name, "is reserved and may not be used as a column name of data"))
  }
  
  #If gform is NULL, it will be set by GetDefaultForm; no need to check here
  if (!is.null(gform)) {
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
      if (nrow(gform) != nrow(data)) stop("if gform is numeric, it should have the same number of rows as data")
      if (ncol(gform) != length(nodes$AC)) stop("if gform is numeric, it should have the same number of columns as length(c(Anodes, Cnodes))")
      if (length(dim(gform)) != 3 || dim(gform)[3] != num.regimes) stop("if gform is numeric, dim[3] should be num.regimes")
      if (!is.null(deterministic.g.function)) stop("if gform is numeric, deterministic.g.function must be NULL")
      if (max(gform, na.rm=T) > 1 || min(gform, na.rm=T) < 0) stop("if gform is numeric, all values should be probabilities")
    }
  }
  
  #If Qform is NULL, it will be set by GetDefaultForm; no need to check here
  if (!is.null(Qform)) {
    if (! is.character(Qform)) stop("Qform should be a character vector")
    if (length(Qform) != length(nodes$LY)) {
      stop("length of Qform is not equal to number of L/Y nodes")
    }
    for (i in 1:length(Qform)) {
      if (length(names(Qform[i])) == 0) stop("Each element of Qform must be named. The name must match the name of the corresponding L/Y node in data.")
      if (names(Qform[i]) != names(data)[nodes$LY[i]]) stop("The name of each element of Q must match the name of the corresponding L/Y node in data.")
      if (Qform[i] != "IDENTITY") {
        #This is only meant to be used in a call by ltmle:::EstimateVariance
        if (LhsVars(Qform[i]) != "Q.kplus1") stop("LHS of each Qform should be Q.kplus1")
        parents <- names(data)[1:(nodes$LY[i]-1)]
        if (!all(RhsVars(Qform[i]) %in% parents)) {
          stop("Some nodes in Qform[", i, "] are not parents of ", names(Qform[i]))
        }    
      }
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
    if (is.null(survivalOutcome)) {
      survivalOutcome <- FALSE
    }    
    if (survivalOutcome) {
      stop("When survivalOutcome is TRUE, all Ynodes should be 0, 1, or NA")
    } 
  }
  
  for (i in nodes$Y) {
    uncensored <- IsUncensored(data, nodes$C, cur.node=i)
    deterministic <- IsDeterministic(data, cur.node=i, deterministic.Q.function=NULL, nodes, called.from.estimate.g=FALSE, survivalOutcome)$is.deterministic #pass deterministic.Q.function=NULL so we're only picking up deaths (if surivalOutcome=FALSE, deterministic will all be FALSE)
    if (any(is.na(data[deterministic, i])) || ! all(data[deterministic, i] == 1)) stop("For survival outcomes, once a Ynode jumps to 1 (e.g. death), all subsequent Ynode values should be 1.")    
    if (any(is.na(data[uncensored, i]))) stop("Ynodes may not be NA except after censoring")
  }
  
  if (! is.equal(dim(regimes)[1:2], c(nrow(data), length(nodes$A)))) stop("Problem with abar or regimes:\n   In ltmleMSM, regimes should have dimensions n x num.Anodes x num.regimes\n   In ltmle, abar should be a matrix with dimensions n x num.Anodes or a vector with length num.Anodes")
  stopifnot(num.regimes == nrow(summary.measures))
  if (!all(regimes %in% c(0, 1, NA))) stop("all regimes should be binary")
  for (i in seq_along(nodes$A)) {
    cur.node <- nodes$A[i]
    uncensored <- IsUncensored(data, nodes$C, cur.node)
    deterministic <- IsDeterministic(data, cur.node, deterministic.Q.function, nodes, called.from.estimate.g=TRUE, survivalOutcome)$is.deterministic
    if (any(is.na(regimes[uncensored & !deterministic, i, ]))) {
      stop("NA in regimes/abar not allowed (except after censoring/death)")
    }
  }
  
  if ((length(dim(summary.measures)) != 3) || ! is.equal(dim(summary.measures)[c(1, 3)], c(num.regimes, length(final.Ynodes)))) stop("summary.measures should be an array with dimensions num.regimes x num.summary.measures x num.final.Ynodes")
  if (class(working.msm) != "character") stop("class(working.msm) must be 'character'")
  if (LhsVars(working.msm) != "Y") stop("the left hand side variable of working.msm should always be 'Y' [this may change in future releases]")
  if (!is.vector(observation.weights) || length(observation.weights) != nrow(data) || any(is.na(observation.weights)) || any(observation.weights < 0) || max(observation.weights) == 0) stop("observation.weights must be NULL or a vector of length nrow(data) with no NAs, no negative values, and at least one positive value")
  
  if (!(is.equal(msm.family, binomial()) || is.equal(msm.family, gaussian()))) stop("msm.family must be either binomial() or gaussian() - an object of class 'family', not a string or function")
  return(list(survivalOutcome=survivalOutcome, binaryOutcome=binaryOutcome))
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
  deterministic.Q.function.depends.on.called.from.estimate.g <- length(grep("called.from.estimate.g", as.character(body(deterministic.Q.function)))) > 0
  for (i in 1:(ncol(data)-1)) {
    if (any(is.na(data[ua, 1:i]))) stop("NA values are not permitted in data except after censoring or a survival event")
    is.deterministic <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=TRUE, survivalOutcome=survivalOutcome)$is.deterministic #check determinisitic including node i 
    
    if (deterministic.Q.function.depends.on.called.from.estimate.g) {
      is.deterministic.Q <- ua & IsDeterministic(data, cur.node=i + 1, deterministic.Q.function=deterministic.Q.function, nodes=nodes, called.from.estimate.g=FALSE, survivalOutcome=survivalOutcome)$is.deterministic 
      if (any(is.deterministic[ua] & !is.deterministic.Q[ua])) stop("Any row set deterministic by deterministic.Q.function(..., called.from.estimate.g=TRUE) must imply that the row is also set deterministic by deterministic.Q.function(..., called.from.estimate.g=FALSE)") #det.Q.fun(T) should imply det.Q.fun(F)
    }
    
    ua[ua] <- !is.deterministic[ua]
    if (any(is.na(ua))) stop("internal ltmle error - ua should not be NA in CleanData")
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
      if (any(is.na(ua))) stop("internal ltmle error - ua should not be NA in CleanData")
    } 
  }
  if (changed && showMessage) {
    message("Note: for internal purposes, all nodes after a censoring event are set to NA and \n all nodes (except Ynodes) are set to NA after Y=1 if survivalFunction is TRUE (or if specified by deterministic.Q.function).\n Your data did not conform and has been adjusted. This may be relevant if you are \n writing your own deterministic function(s) or debugging ltmle.")
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
      inputs.subset <- SubsetInputs(inputs, final.Ynode)
      uncensored <- IsUncensored(inputs.subset$data, inputs.subset$nodes$C, cur.node=final.Ynode)
      if (dim(inputs.subset$regimes)[2] > 0) {
        is.duplicate <- duplicated(inputs.subset$regimes, MARGIN=3)
      } else {
        is.duplicate <- c(FALSE, rep(TRUE, num.regimes - 1))  #in case there are C nodes but no A nodes before a Ynode
      }
      for (i in 1:num.regimes) {
        if (is.duplicate[i]) {
          msm.weights[i, j] <- 0
        } else {
          intervention.match <- InterventionMatch(inputs.subset$data, abar=GetABar(inputs.subset$regimes, i), inputs.subset$nodes$A, cur.node=final.Ynode)
          msm.weights[i, j] <- sum(uncensored & intervention.match) / nrow(inputs.subset$data)
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
  } else if (! is.equal(dim(msm.weights), c(n, num.regimes, num.final.Ynodes))) {
    stop("dim(msm.weights) should be c(n, num.regimes, num.final.Ynodes) or c(num.regimes, num.final.Ynodes)")
  }
  if (any(is.na(msm.weights)) || any(msm.weights < 0)) stop("all msm.weights must be >= 0 and not NA")
  return(msm.weights)
}

# Converts a general formula to a main terms formula and combine summary measures with baseline covariates
# Ex: If working.msm is "Y ~ X1*X2", convert to "Y ~ -1 + S1 + S1 + S3 + S4" where 
# S1 is 1 (intercept), S2 is X1, S3 is X2, S4 is X1:X2
ConvertToMainTerms <- function(data, msm, summary.measures, nodes) {
  baseline.column.names <- names(data)[seq(1, min(c(nodes$A, nodes$L, nodes$C, nodes$Y)) - 1)]
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

#Before passing data to SuperLearner, convert factors to binary
ConvertCensoringNodesToBinary <- function(data, Cnodes) {
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
