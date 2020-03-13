## ----setup, include = FALSE---------------------------------------------------
library(ltmle)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
n <- 10000
W <- rnorm(n)
A <- rexpit(W)
Y <- rexpit(W + A)
data <- data.frame(W, A, Y)
head(data)

## -----------------------------------------------------------------------------
mean(Y)
mean(Y[A == 1])

## -----------------------------------------------------------------------------
result <- ltmle(data, Anodes = "A", Ynodes = "Y", abar = 1)
result

## -----------------------------------------------------------------------------
n.large <- 1e6
W <- rnorm(n.large)
A <- 1
Y <- rexpit(W + A)
mean(Y)

## -----------------------------------------------------------------------------
n <- 1000
W1 <- rnorm(n)
W2 <- rbinom(n, size=1, prob=0.3)   
W3 <- rnorm(n)
A <- rexpit(-1 + 2 * W1 + W3)
Y <- rexpit(-0.5 + 2 * W1^2 + 0.5 * W2 - 0.5 * A + 0.2 * W3 * A - 1.1 * W3)
data <- data.frame(W1, W2, W3, A, Y)

## -----------------------------------------------------------------------------
result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", abar=1, SL.library="default")

## -----------------------------------------------------------------------------
summary(result)

## -----------------------------------------------------------------------------
summary(result, estimator="iptw")

## -----------------------------------------------------------------------------
result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
                Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), 
                gform="A ~ W1 + W3", abar=1, SL.library="default")
summary(result)

## -----------------------------------------------------------------------------
result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
 Qform=c(Y="Q.kplus1 ~ I(W1^2) + W2 + W3*A"), gform="A ~ W1 + W3", 
 abar=1, SL.library=NULL)
summary(result)

## -----------------------------------------------------------------------------
result <- ltmle(data, Anodes="A", Lnodes=NULL, Ynodes="Y", 
                      abar=list(1, 0), SL.library="default")
summary(result)

## -----------------------------------------------------------------------------
n <- 100000
W <- rnorm(n)
C <- BinaryToCensoring(is.censored = rexpit(W))
summary(C)
Y <- rep(NA, n)
Y[C == "uncensored"] <- rexpit(W[C == "uncensored"])
data <- data.frame(W, C, Y)
head(data, 20)
result <- ltmle(data, Anodes = NULL, Cnodes = "C", Ynodes = "Y", abar = NULL)
summary(result)

## -----------------------------------------------------------------------------
mean(data$Y)
mean(data$Y, na.rm = T)

## -----------------------------------------------------------------------------
n <- 1000
W <- rnorm(n)
A1 <- rexpit(W)
L <- 0.3 * W + 0.2 * A1 + rnorm(n)
A2 <- rexpit(W + A1 + L)
Y <- rexpit(W - 0.6 * A1 + L - 0.8 * A2)
data <- data.frame(W, A1, L, A2, Y)
head(data)

## -----------------------------------------------------------------------------
ltmle(data, Anodes=c("A1", "A2"), Lnodes="L", Ynodes="Y", abar=c(0, 0))

## -----------------------------------------------------------------------------
n <- 1000
W <- rnorm(n)
A1 <- rexpit(W)
C <- BinaryToCensoring(is.censored = rexpit(0.6 * W - 0.5 * A1))
uncensored <- C == "uncensored"
L <- A2 <- Y <- rep(NA, n)
L[uncensored] <- (0.3 * W[uncensored] + 0.2 * A1[uncensored] + rnorm(sum(uncensored)))
A2[uncensored] <- rexpit(W[uncensored] + A1[uncensored] + L[uncensored])
Y[uncensored] <- rexpit(W[uncensored] - 0.6 * A1[uncensored] + L[uncensored] - 0.8 * A2[uncensored])
data <- data.frame(W, A1, C, L, A2, Y)
head(data)

## -----------------------------------------------------------------------------
ltmle(data, Anodes=c("A1", "A2"), Cnodes = "C", Lnodes="L", Ynodes="Y", abar=c(1, 0))

## -----------------------------------------------------------------------------
abar <- matrix(nrow=n, ncol=2)
abar[, 1] <- 1
abar[, 2] <- L > 0

result.abar <- ltmle(data, Anodes=c("A1", "A2"), Cnodes = "C", Lnodes="L", Ynodes="Y", abar=abar)
result.abar

## -----------------------------------------------------------------------------
rule <- function(row) c(1, row["L"] > 0)

result.rule <- ltmle(data, Anodes=c("A1", "A2"), Cnodes = "C", Lnodes="L", Ynodes="Y", rule=rule)
result.rule

## -----------------------------------------------------------------------------
summary(result.abar)
summary(result.rule)

## -----------------------------------------------------------------------------
n <- 1000
W <- rnorm(n)
A <- rexpit(4 * W)
Y <- rexpit(W + A)
df <- data.frame(W, A, Y)

## -----------------------------------------------------------------------------
r1 <- ltmle(df, Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE)
print(summary(r1))

## -----------------------------------------------------------------------------
r2 <- ltmle(df, Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, 
 variance.method="ic")
print(summary(r2))

## -----------------------------------------------------------------------------
r3 <- ltmle(df, Anodes="A", Ynodes="Y", abar = 1, estimate.time=FALSE, 
 variance.method="iptw")
print(summary(r3))

## -----------------------------------------------------------------------------
print(summary(r1, estimator="iptw"))
print(summary(r2, estimator="iptw")) #the same - variance.method only affects TMLE
print(summary(r3, estimator="iptw")) #the same - variance.method only affects TMLE

## -----------------------------------------------------------------------------
summary(r1$cum.g)
summary(r1$cum.g.unbounded)
head(data.frame(df, g = r1$cum.g, g.unbounded = r1$cum.g.unbounded), 20)

## -----------------------------------------------------------------------------
num.households <- 500
people.in.household <- round(runif(num.households, min = 1, max = 10))
length(people.in.household)
n <- sum(people.in.household) 
n
W.household <- rnorm(num.households)
length(W.household)
W.household.expanded <- rep(W.household, times = people.in.household)
W.indiv <- rnorm(n)
length(W.indiv)
A <- rexpit(1.5 * W.household.expanded + 0.4 * W.indiv)
Y <- rexpit(-1 + 2.3 * W.household.expanded - 0.6 * W.indiv + 1.2 * A)

## -----------------------------------------------------------------------------
id <- 1:num.households 
id.expanded <- rep(id, times = people.in.household)
data <- data.frame(W.household.expanded, W.indiv, A, Y)
head(cbind(id.expanded, data), 20)

## -----------------------------------------------------------------------------
result.without.id <- ltmle(data, Anodes = "A", Ynodes = "Y", abar = 0)
result.with.id <- ltmle(data, Anodes = "A", Ynodes = "Y", abar = 0, id = id.expanded)

## -----------------------------------------------------------------------------
summary(result.without.id)
summary(result.with.id)

## -----------------------------------------------------------------------------
length(result.without.id$IC$tmle)
length(result.with.id$IC$tmle)

## -----------------------------------------------------------------------------
n <- 1000
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

## -----------------------------------------------------------------------------
result <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
 Ynodes=grep("^Y", names(data)), abar=c(1, 0)) 
summary(result)

## -----------------------------------------------------------------------------
result <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
 Ynodes=grep("^Y", names(data)), abar=c(1, 0), 
 Qform=c(L1a="Q.kplus1 ~ 1", L2a="Q.kplus1 ~ 1"))
summary(result)

## -----------------------------------------------------------------------------
result <- ltmle(data, Anodes=c(3, 7), Lnodes=c("L1a", "L1b", "L2a", "L2b"), 
 Ynodes=grep("^Y", names(data)), abar=c(1, 0), 
 Qform=c(L1a="Q.kplus1 ~ 1", L1b="Q.klus1~A1", 
 Y1="Q.kplus1~L1a", L2a="Q.kplus1 ~ 1", L2b="Q.klus1~A1", Y2="Q.kplus1~A2 + gender"))
summary(result)

