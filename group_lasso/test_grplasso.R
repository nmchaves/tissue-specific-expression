rm(list=ls())
library(grplasso)
## ISSUE NO CV EASILY USABLE
## Perform the Logistic Group Lasso on a random dataset

cv.grplasso <- function(x,y,index,nfolds = 5) {
  # select the parameters
  lambda <- lambdamax(x, y = y, index = index, penscale = sqrt, model = LogReg()) * 0.1^(0:10)
  N <- nrow(x)
  y <- drop(y)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=5 recommended")

  # fit the model and store the outputs
  foldid <- sample(rep(seq(nfolds), length = N))
  outlist <- as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    which <- foldid == i
    y_sub <- y[!which]
    x_sub <- x[!which, , drop = FALSE]
    fit <- grplasso(x = x_sub, y = y_sub, index = index, lambda = lambda, model = LogReg(), penscale = sqrt,
                    control = grpl.control(update.hess = "lambda", trace = 0))
    coeffs <- fit$coefficients
    # compute the auc for each of the lambda parameters
  }
  print(outlist)
}

set.seed(1)
n <- 50  ## observations
p <- 4   ## variables
## First variable (intercept) not penalized, two groups having 2 degrees
## of freedom each
index <- c(NA, 1, 1, 3, 3)

## Create a random design matrix, including the intercept (first column)
x <- cbind(1, matrix(rnorm(p * n), nrow = n))
colnames(x) <- c("Intercept", paste("X", 1:4, sep = ""))

par <- c(0, 2.1, -1.8, 0, 0)
prob <- 1 / (1 + exp(-x %*% par))
mean(pmin(prob, 1 - prob)) ## Bayes risk
y <- rbinom(n, size = 1, prob = prob) ## binary response vector

## Use a multiplicative grid for the penalty parameter lambda, starting
## at the maximal lambda value
lambda <- lambdamax(x, y = y, index = index, penscale = sqrt, model = LogReg()) * 0.1^(0:10)


## Fit the solution path on the lambda grid
fit <- grplasso(x, y = y, index = index, lambda = lambda, model = LogReg(), penscale = sqrt,
                control = grpl.control(update.hess = "lambda", trace = 3))
## Plot coefficient paths
plot(fit)

pred <- predict(fit)

plot(pred,pred.resp)

pred <- predict(fit, x,type='response')
