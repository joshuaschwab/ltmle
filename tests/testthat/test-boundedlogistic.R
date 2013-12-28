context("Test boundedlogistic")

test_that("boundedlogistic isn't too far off glm", {
  set.seed(88)

  n <- 100
  p <- 5

  tol <- 1e-3 #this isn't that good. need to improve the optimization method

  niter <- 10
  for (i in 1:niter) {
   X <- matrix(runif(n*p), n)
   beta <- runif(p+1) * 3 - 1.5
   pY <- plogis(cbind(1,X) %*% matrix(beta))
   Y <- rbinom(n,1,pY)
   dat <- data.frame(Y=Y, X=X)

   wt <-runif(n)*10
   rounded.wt <- round(wt)

   expect_equal(coef(boundedlogistic(Y~., dat)),
                coef(glm(Y~.,quasibinomial, dat)),
                tolerance=tol, check.attributes=FALSE)


   expect_equal(coef(boundedlogistic(Y~., dat, weight=wt)), 
                coef(glm(Y~.,quasibinomial, dat, weight=wt)),
                tolerance=tol, check.attributes=FALSE)

   expect_equal(coef(boundedlogistic(Y~., dat, weight=rounded.wt)), 
                coef(glm(Y~.,quasibinomial, dat, weight=rounded.wt)),
                tolerance=tol, check.attributes=FALSE)

 }

 })