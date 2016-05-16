rm(list=ls())

library(grplasso)
library(AUC)
library(ggplot2)
## ISSUE NO CV EASILY USABLE
## Perform the Logistic Group Lasso on a random dataset

cv.grplasso <- function(x, y, index,nfolds = 5, nlamdas = 20, plot.error=FALSE) {
  # select the parameters
  lambda <- lambdamax(x, y = y, index = index, penscale = sqrt, model = LogReg()) * 0.5^seq(0,4,len=nlamdas)
  N <- nrow(x)
  y <- drop(y)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=5 recommended")

  # fit the model and store the outputs
  foldid <- sample(rep(seq(nfolds), length = N))
  aucs <- matrix(,nrow=nfolds,ncol=length(lambda))
  for (i in seq(nfolds)) {
    which <- foldid == i
    y_train <- y[!which]
    x_train <- x[!which, , drop = FALSE]
    y_test <- y[which]
    x_test <- x[which, , drop = FALSE]
    fit <- grplasso(x = x_train, y = y_train, index = index, lambda = lambda, model = LogReg(), penscale = sqrt,
                    control = grpl.control(update.hess = "lambda", trace = 0))
    coeffs <- fit$coefficients
    pred.resp <- predict(fit, x_test, type = "response")
    # compute the auc for each of the lambda parameters
    # auc_score = accuracy(y,coeffs)
    for (j in seq(length(lambda))) {
      predictions <- pred.resp[,j]
      aucs[i,j] <- auc(accuracy(predictions,as.factor(y_test)))
    }
  }
  error.mean = apply(aucs,2,mean)
  error.sd = apply(aucs,2,sd)
  # plot the c.v. error 
  if (plot.error) {
    error.high = error.mean + error.sd
    error.low = error.mean - error.sd
    plot(lambda,error.mean,
       xlab="Lambda (log)",ylab="CV Error", 
       log="x",xlim=rev(range(lambda)),ylim = c(0,1))
    arrows(lambda,error.high,lambda,error.low,col=2,angle=90,length=0.05,code=3)
  }
  # select the variable that minmizes the average CV
  max.idx <- which.max(error.mean);
  lambda.opt <- lambda[max.idx]
  cv.error <- error.mean[max.idx]
  return(list(lambda = lambda.opt, error=cv.error))
}

split_data <- function(gene_features, labels, gene_ids_ordered, train_set_size=0.7) {
  # train_set_size: Fraction of genes used for training set
  num_samples = length(labels)
  num_features = dim(gene_features)[1]
  num_train_samples = ceiling(train_set_size*num_samples)
  sample_order = sample(1:num_samples, num_samples) # permute the samples
  train_idx = sample_order[1:num_train_samples]
  test_idx =  sample_order[(num_train_samples+1):num_samples]

  train = list(x=gene_features[train_idx,],y=labels[train_idx],ids=gene_ids_ordered[train_idx])
  test = list(x=gene_features[test_idx,],y=labels[test_idx],ids=gene_ids_ordered[test_idx])
  return(list(train=train, test=test))
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

data <- split_data(x,y,y)
cat('----------------------------------------\n')
cat('dimensionality of data: ',dim(data$train$x)[2],'\n')
cat('# of training samples:  ',length(data$train$y),'\n')
cat('# of test samples:      ',length(data$test$y),'\n')

# apply cv for grplasso
cat('----------------------------------------\n')
cat('Training group lasso classifier...')
x <- data$train$x
y <- data$train$y
cv.result <- cv.grplasso(x,y,index,plot.error=TRUE)
## Re-fit the model with the best tuning paramter from cross-validation
fit <- grplasso(x, y = y, index = index, lambda = cv.result$lambda, model = LogReg(), penscale = sqrt,
                control = grpl.control(update.hess = "lambda", trace = 0))
cat('done\n')
cat('cross-validaiton error:',cv.result$error,'\n')
cat('parameter (lambda) tuned as:',cv.result$lambda,'\n')

## compute the test error
cat('----------------------------------------\n')
prediction <- predict(fit, data$test$x, type = "response")
auc.val <- auc(accuracy(prediction,as.factor(data$test$y)))
cat('Test Error:',auc.val,'\n')



