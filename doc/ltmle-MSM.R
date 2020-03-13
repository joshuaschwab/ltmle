## ----setup, include = FALSE---------------------------------------------------
library(ltmle)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
n <- 5000
time.points <- 5
prev.L <- rnorm(n)
prev.A <- rep(0, n)
sum.A <- rep(0, n)
data <- data.frame(L_0 = prev.L)
for (t in 1:time.points) {
  L <- 0.1 * prev.L + 0.3 * prev.A + rnorm(n)
  A <- rexpit(L)
  
  data1 <- data.frame(L, A)
  names(data1) <- paste0(c("L_", "A_"), t)
  data <- cbind(data, data1)
  
  prev.A <- A
  prev.L <- L
  
  sum.A <- sum.A + A
}
data$Y <- rexpit(sum.A / time.points + L)
head(data)

## -----------------------------------------------------------------------------
regime.matrix <- as.matrix(expand.grid(rep(list(0:1), time.points)))
dim(regime.matrix)
head(regime.matrix, 20)
num.regimes <- 2^time.points
regimes <- array(dim = c(n, time.points, num.regimes)) #n x numAnodes x numRegimes = n x time.points x 2^time.points
summary.measures <- array(dim = c(num.regimes, 1, 1)) #numRegimes x num.summary.measures x num.final.Ynodes = 2^time.points x 1 x 1
for (i in 1:num.regimes) {
  regimes[, , i] <- matrix(regime.matrix[i, ], byrow = TRUE, nrow = n, ncol = time.points)
  summary.measures[i, 1, 1] <- sum(regime.matrix[i, ])
}
colnames(summary.measures) <- "time.on.treatment"
regimes[1:3, , 1:3]
summary.measures

## -----------------------------------------------------------------------------
result1 <- ltmleMSM(data, Anodes = paste0("A_", 1:time.points), 
                    Lnodes = paste0("L_", 0:time.points), Ynodes = "Y", 
                    regimes = regimes, summary.measures = summary.measures, 
                    working.msm = "Y ~ time.on.treatment + I(time.on.treatment^2)",
                    variance.method = "ic")
result2 <- ltmleMSM(data, Anodes = paste0("A_", 1:time.points), 
                    Lnodes = paste0("L_", 0:time.points), Ynodes = "Y", 
                    regimes = regimes, summary.measures = summary.measures, 
                    working.msm = "Y ~ time.on.treatment + I(time.on.treatment^2)")
summary(result1)
summary(result2)

## -----------------------------------------------------------------------------
at.least.3 <- summary.measures[, "time.on.treatment", 1] >= 3
regimes.3 <- regimes[, , at.least.3]
summary.measures.3 <- summary.measures[at.least.3, , , drop = F]
dim(regimes.3)
summary.measures.3

result <- ltmleMSM(data, Anodes = paste0("A_", 1:time.points), 
                   Lnodes = paste0("L_", 0:time.points), 
                   Ynodes = "Y", regimes = regimes.3, 
                   summary.measures = summary.measures.3, 
                   working.msm = "Y ~ time.on.treatment + I(time.on.treatment^2)")
summary(result)

## -----------------------------------------------------------------------------
data(sampleDataForLtmleMSM)
head(sampleDataForLtmleMSM$data, 20)
dim(sampleDataForLtmleMSM$regimes)
sampleDataForLtmleMSM$regimes[1:5, , ]
sampleDataForLtmleMSM$summary.measures

## -----------------------------------------------------------------------------
Anodes <- c("A0", "A1", "A2")
Lnodes <- c("CD4_1", "CD4_2")
Ynodes <- c("Y1", "Y2", "Y3")

## -----------------------------------------------------------------------------
msm.weights <- matrix(1:12, nrow=4, ncol=3)  

## -----------------------------------------------------------------------------
result.regimes <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, 
                   Lnodes=Lnodes, Ynodes=Ynodes, 
                   survivalOutcome=TRUE,
                   regimes=sampleDataForLtmleMSM$regimes, 
                   summary.measures=sampleDataForLtmleMSM$summary.measures,
                   final.Ynodes=Ynodes, 
                   working.msm="Y ~ male + time + pmax(time - switch.time, 0)", 
                   msm.weights=msm.weights, estimate.time=FALSE)
print(summary(result.regimes))

## -----------------------------------------------------------------------------
regimesList <- list(function(row) c(1,1,1),
                     function(row) c(0,1,1),
                     function(row) c(0,0,1),
                     function(row) c(0,0,0))
result.regList <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, 
                   Lnodes=Lnodes, Ynodes=Ynodes, 
                   survivalOutcome=TRUE, regimes=regimesList, 
                   summary.measures=sampleDataForLtmleMSM$summary.measures,
                   final.Ynodes=Ynodes, 
                   working.msm="Y ~ male + time + pmax(time - switch.time, 0)", 
                   msm.weights=msm.weights, estimate.time=FALSE)

## -----------------------------------------------------------------------------
print(summary(result.regList))         

## -----------------------------------------------------------------------------
result <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, 
                   Lnodes=Lnodes, Ynodes=Ynodes, 
                   survivalOutcome=TRUE,
                   regimes=sampleDataForLtmleMSM$regimes, 
                   summary.measures=sampleDataForLtmleMSM$summary.measures[, , c(1, 3)],
                   final.Ynodes=c("Y1", "Y3"), 
                   working.msm="Y ~ male + time + pmax(time - switch.time, 0)", 
                   estimate.time=FALSE)
summary(result)

## -----------------------------------------------------------------------------
result <- ltmleMSM(sampleDataForLtmleMSM$data, Anodes=Anodes, 
                    Lnodes=Lnodes, Ynodes=Ynodes, 
                    survivalOutcome=TRUE,
                    regimes=sampleDataForLtmleMSM$regimes, 
                    summary.measures=sampleDataForLtmleMSM$summary.measures[, , 3],
                    final.Ynodes="Y3", 
                    working.msm="Y ~ male + switch.time", 
                    estimate.time=FALSE)
summary(result)

## -----------------------------------------------------------------------------
GenerateData <- function(n, abar = NULL) {
  W <- rnorm(n)
  if (is.null(abar)) {
    A1 <- rexpit(W)
  } else {
    A1 <- abar$a1
  }
  L <- plogis(rnorm(n) + 0.3 * W + 0.5 * A1)
  if (is.null(abar)) {
    A2 <- rexpit(-0.5 * W + A1 - 0.6 * L)
  } else {
    A2 <- as.integer(L > abar$theta)
  }
  Y <- rexpit(-2 + W + A1 + L + 2 * A2)
  if (is.null(abar)) {
    return(data.frame(W, A1, L, A2, Y))
  } else {
    return(mean(Y))
  }
}

## -----------------------------------------------------------------------------
set.seed(11)
n <- 10000
data <- GenerateData(n)
regimes <- array(dim = c(n, 2, 10)) #n x num.Anodes x num.regimes
theta.set <- seq(0, 1, length.out = 5)
summary.measures <- array(theta.set, dim = c(10, 2, 1))
colnames(summary.measures) <- c("a1", "theta")
cnt <- 0
for (a1 in 0:1) {
  for (theta.index in 1:5) {
    cnt <- cnt + 1
    regimes[, 1, cnt] <- a1
    regimes[, 2, cnt] <- data$L > theta.set[theta.index]
    summary.measures[cnt, , 1] <- c(a1, theta.set[theta.index])
  }
}
summary.measures
head(data, 3)
regimes[1:3, , ]

## -----------------------------------------------------------------------------
working.msm <- "Y ~ a1*theta"
summary(ltmleMSM(data, Anodes = c("A1", "A2"), Lnodes = "L", Ynodes = "Y", 
              regimes = regimes, summary.measures = summary.measures, 
              working.msm = working.msm))

## -----------------------------------------------------------------------------
truth <- rep(NA_real_, 10)
cnt <- 0
for (a1 in 0:1) {
  for (theta.index in 1:5) {
    cnt <- cnt + 1
    truth[cnt] <- GenerateData(n = 1e6, 
                    abar = list(a1 = a1, theta = theta.set[theta.index]))
  }
}

## -----------------------------------------------------------------------------
m.true <- glm(working.msm, 
              data = data.frame(Y = truth, summary.measures[, , 1]), 
              family = "quasibinomial")
m.true

