## ----setup, include = FALSE---------------------------------------------------
library(ltmle)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
set.seed(123)
rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
n <- 1000
ua <- rep(TRUE, n)   #ua = uncensored and alive
L1 <- A1 <- Y1 <- C2.binary <- L2 <- A2 <- Y2 <- rep(NA_real_, n)
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

## -----------------------------------------------------------------------------
result <- ltmle(data, Anodes=c("A1","A2"), Cnodes=c("C1", "C2"), 
                Lnodes=c("L1", "L2"), Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
                survivalOutcome=TRUE)
summary(result) 

## -----------------------------------------------------------------------------
deterministic.g.function <- function(data, current.node, nodes) {
  if (names(data)[current.node] == "A1") {
    det <- (data$L1 < -2 | data$L1 > 2) & !is.na(data$L1)
    prob1 <- ((data$L1 < -2) * 1 + (data$L1 > 2) * 0.1)[det]
  } else if (names(data)[current.node] == "A2") {
    det <- data$A1 == 1 & !is.na(data$A1)
    prob1 <- 1
  } else if (names(data[current.node]) %in% c("C1", "C2")){
    return(NULL)  #this returns the default of no deterministic links 
    #note that it is not necessary to specify that prior censoring indicates future censoring
  } else {
    stop("unexpected current.node")
  }
  return(list(is.deterministic=det, prob1=prob1))  
}
result <- ltmle(data, Anodes=c("A1","A2"), Cnodes=c("C1", "C2"), 
                Lnodes=c("L1", "L2"), Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
                survivalOutcome=TRUE, deterministic.g.function = deterministic.g.function)
summary(result) 

## -----------------------------------------------------------------------------
n <- 1000
L2 <- A2 <- Y2 <- as.numeric(rep(NA, n))
W <- rnorm(n)
A1 <- rexpit(W)
Y1 <- rexpit(W - A1)
alive <- Y1 == 0
L2[alive] <- (0.5 * W - 0.9 * A1 + rnorm(n))[alive]
completed.study <- alive & L2 > 0

## -----------------------------------------------------------------------------
A2[alive & !completed.study] <- rexpit((0.5 * W + 0.8 * L2)[alive & !completed.study])

Y2[alive & !completed.study] <- rexpit((L2 - 0.8 * A1 - A2)[alive & !completed.study])
Y2[alive & completed.study] <- 0
Y2[!alive] <- 1  # if a patient dies at time 1, record death at time 2 as well
data <- data.frame(W, A1, Y1, L2, A2, Y2)

## -----------------------------------------------------------------------------
det.Q.fun.1 <- function(data, current.node, nodes, called.from.estimate.g) {
  L2.index <- which(names(data) == "L2")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data$L2 > 0 & !is.na(data$L2)
  return(list(is.deterministic=is.deterministic, Q.value=0))
}

result.scenario1 <- ltmle(data, Anodes=c("A1","A2"), Lnodes="L2", Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
  SL.library=NULL, estimate.time=FALSE, deterministic.Q.function=det.Q.fun.1, survivalOutcome=TRUE)
summary(result.scenario1)

## -----------------------------------------------------------------------------
A2[alive] <- rexpit((0.5 * W + 0.8 * L2)[alive])  #patients can change treatment after leaving study
Y2[alive & !completed.study] <- rexpit((L2 - 0.8 * A1 - A2)[alive & !completed.study])
Y2[alive & completed.study] <- 0
Y2[!alive] <- 1  # if a patient dies at time 1, record death at time 2 as well
data <- data.frame(W, A1, Y1, L2, A2, Y2)

det.Q.fun.2 <- function(data, current.node, nodes, called.from.estimate.g) {
  #there is no deterministic information when calculating g - treatment may still change
  if (called.from.estimate.g) return(NULL)  
  
  L2.index <- which(names(data) == "L2")
  stopifnot(length(L2.index) == 1)
  L2.in.history <- L2.index < current.node
  if (! L2.in.history) return(NULL)
  
  is.deterministic <- data$L2 > 0 & !is.na(data$L2)
  return(list(is.deterministic=is.deterministic, Q.value=0))
}

result.scenario2 <- ltmle(data, Anodes=c("A1","A2"), Lnodes="L2", Ynodes=c("Y1", "Y2"), abar=c(1, 1), 
 SL.library=NULL, estimate.time=FALSE, deterministic.Q.function=det.Q.fun.2, survivalOutcome=TRUE)
summary(result.scenario2)

