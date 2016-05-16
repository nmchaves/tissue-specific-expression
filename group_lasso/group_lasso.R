rm(list=ls())

# load libraries
library(grplasso)
library(AUC)
library(ggplot2)

cv.grplasso <- function(x, y, index, nfolds = 5, nlamdas = 20, plot.error=FALSE) {
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

split_data <- function(gene_features, labels, train_set_size=0.67) {
  # train_set_size: Fraction of genes used for training set
  num_samples = length(labels)
  print(num_samples)
  num_features = dim(gene_features)[1]
  num_train_samples = ceiling(train_set_size*num_samples)
  sample_order = sample(1:num_samples, num_samples) # permute the samples
  train_idx = sample_order[1:num_train_samples]
  test_idx =  sample_order[(num_train_samples+1):num_samples]

  train = list(x=gene_features[train_idx,],y=labels[train_idx])
  test = list(x=gene_features[test_idx,],y=labels[test_idx])
  return(list(train=train, test=test))
}

load.pos.neg.sets <- function(dir.name,pos.name,neg.name,grp.name) {
  pos.data <- read.table(paste(dir.name,pos.name,sep=''),sep='\t',row.names = 1,skip=2)
  neg.data <- read.table(paste(dir.name,neg.name,sep=''),sep='\t',row.names = 1,skip=2)
  response = c(rep(1,dim(pos.data)[1]),rep(0,dim(neg.data)[1]))
  # TODO: CHANGE GROUP
  group = sample(54,dim(pos.data)[2],replace = TRUE) # random group selection
  # print(response)
  # print(dim(pos.data))
  # print(rownames(pos.data))
  # print(rownames(neg.data))
  return(list(x=rbind(pos.data,neg.data),y=response,group=group))
}


set.seed(1)
n <- 50  ## observations
p <- 4   ## variables
## First variable (intercept) not penalized, two groups having 2 degrees
## of freedom each

## Create a random design matrix, including the intercept (first column)
x <- cbind(1, matrix(rnorm(p * n), nrow = n))
colnames(x) <- c("Intercept", paste("X", 1:4, sep = ""))
par <- c(0, 2.1, -1.8, 0, 0)
prob <- 1 / (1 + exp(-x %*% par))
mean(pmin(prob, 1 - prob)) ## Bayes risk
y <- rbinom(n, size = 1, prob = prob) ## binary response vector

## load all data
cat('----------------------------------------\n')
cat('Loading data...\n')
dir.name <- '/Users/jjzhu/Documents/GTEx/CS341_Code/aws_mock/experiment_inputs_subset/'
pos.name <- 'GO:0000578_pos.txt'  # filename for positive set
neg.name <- 'GO:0000578_neg_0.txt' # filename for negative set
grp.name <- '
full.data <- load.pos.neg.sets(dir.name,pos.name,neg.name)
index <- c(NA, full.data$group)
full.x <- cbind(1, full.data$x)  # add intercept
full.y <- full.data$y

## split data
cat('----------------------------------------\n')
data <- split_data(full.x,full.y)
cat('dimensionality of data: ',dim(data$train$x)[2],'\n')
cat('# of training samples:  ',length(data$train$y),'\n')
cat('# of test samples:      ',length(data$test$y),'\n')

## apply cv for grplasso
cat('----------------------------------------\n')
cat('Training group lasso classifier...')
x <- data$train$x
y <- data$train$y
cv.result <- cv.grplasso(x,y,index,nfolds=3,plot.error=TRUE)
## Re-fit the model with the best tuning paramter from cross-validation
fit <- grplasso(x, y = y, index = index, lambda = cv.result$lambda, model = LogReg(), 
                penscale = sqrt,control = grpl.control(update.hess = "lambda", trace = 0))
cat('done\n')
cat('cross-validaiton error:',cv.result$error,'\n')
cat('parameter (lambda) tuned as:',cv.result$lambda,'\n')

## compute the test error
cat('----------------------------------------\n')
prediction <- predict(fit, data$test$x, type = "response")
auc.val <- auc(accuracy(prediction,as.factor(data$test$y)))
cat('Test Error:',auc.val,'\n')



