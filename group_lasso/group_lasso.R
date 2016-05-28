rm(list=ls())

# load libraries
library(grplasso)
library(AUC)
library(ggplot2)
library(caret)

cv.grplasso <- function(x, y, index, nfolds = 5, nlamdas = 20, plot.error=FALSE) {
  # select the parameters
  lambda <- lambdamax(x, y=y, index=index, penscale = sqrt, model = LogReg()) * 0.5^seq(0,5,len=nlamdas)
  N <- nrow(x)
  y <- drop(y)
  if (nfolds < 3) 
    stop("nfolds must be bigger than 3; nfolds=5 recommended")
  aucs <- matrix(,nrow=nfolds,ncol=length(lambda))
  flds <- createFolds(y,nfolds,list=TRUE,returnTrain=FALSE)
  for (i in seq(nfolds)) {
    test.index <- flds[[i]]
    y_test <- y[test.index]
    x_test <- x[test.index, , drop = FALSE]
    y_train <- y[-test.index]
    x_train <- x[-test.index, , drop = FALSE]
    cat('- fold number:',i,'\n')
    fit <- grplasso(x = x_train, y = y_train, index = index, lambda = lambda, model = LogReg(), penscale = sqrt,
                    control = grpl.control(update.hess = "lambda", trace = 0))
    coeffs <- fit$coefficients
    pred.resp <- predict(fit, x_test, type = "response")
    for (j in seq(length(lambda))) {
      predictions <- pred.resp[,j]
      aucs[i,j] <- auc(accuracy(predictions,as.factor(y_test)))
    }
  }
  aucs = aucs[complete.cases(aucs),]
  error.mean = apply(aucs,2,mean)
  error.sd = apply(aucs,2,sd)
  # plot the c.v. error 
  if (plot.error) {
    error.high = error.mean + error.sd
    error.low = error.mean - error.sd
    plot(lambda,error.mean,
         xlab="Lambda (log)",ylab="CV AUC Score", 
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
  set.seed(1)
  # train_set_size: Fraction of genes used for training set
  train.index <- createDataPartition(labels, p = train_set_size, list = FALSE)
  train = list(x=gene_features[train.index,],y=labels[train.index])
  test = list(x=gene_features[-train.index,],y=labels[-train.index])
  return(list(train=train, test=test))
}

group.to.int <- function(group) {
  type.names = unique(group)
  type.counts = rep(0,length(type.names))
  group.idx = rep(0,length(group))
  for (i in 1:length(group.idx)) {
    type.idx = which(type.names == group[i])
    group.idx[i] = type.idx
    type.counts[type.idx] = type.counts[type.idx] + 1
  }
  names(type.counts) = type.names
  return(list(types=type.counts,idx=group.idx))
}

load.pos.neg.sets <- function(pos.name,neg.name,group.info) {
  pos.data <- read.table(pos.name,sep='\t',row.names=1,skip=2)
  neg.data <- read.table(neg.name,sep='\t',row.names=1,skip=2)
  
  # group: tissue type, tissue specific type
  data <- rbind(pos.data,neg.data)
  response <- c(rep(1,dim(pos.data)[1]),rep(0,dim(neg.data)[1]))
  return(list(x=data,
              y=response,
              group=group.info$idx,
              types=group.info$types))
}

plot.coefficient.path <- function(x,y,index) {
  lambda <- lambdamax(x, y = y, index = index, penscale = sqrt, model = LogReg()) * 0.5^seq(0,5,len=10)
  fit <- grplasso(x = x, y = y, index = index, lambda = lambda, model = LogReg(), penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0))
  plot(fit,log='x')
}

############
### MAIN ###
############
main <- '/Users/jasonzhu/Documents/CS341_Code/'
dir.name <- paste(main,'data/pca_experiment_inputs/',sep='')

## load all go names
all.go.name <- paste(main,'data/GO_terms_final_gene_counts.txt',sep='')
go.names <- read.table(all.go.name)$V1

## load group information
origin.filename <- paste(main,'data/local_large/log_norm_pca_transcript_rpkm_in_go_nonzero_exp.txt',sep='')
title.fields <- strsplit(readLines(origin.filename, n=1) ,'\t')[[1]]
group.info <- group.to.int(title.fields[5:length(title.fields)])

for (go.idx in 1:length(go.names)){
  set.seed(1)
  go.term <- as.character(go.names[go.idx])
  for (neg_idx in 2:9) {
    
    # neg_idx <- 1
    neg_pfx <- paste(paste('_neg_',neg_idx,sep=''),'.txt.txt',sep='')
    out_pfx <- paste(paste(go.term,'_',sep=''),neg_idx,sep='')
    out.name <- paste(paste(paste(main,'data/grplasso_results/grplasso_',sep=''),out_pfx,sep=''),'.txt',sep='')
    cat('========================================\n')
    cat('Target output: ', out.name,'\n')
    
    ## load all data
    cat('----------------------------------------\n')
    cat('Loading raw data from',go.term,'...\n')
    pos.name <- paste(dir.name,paste(go.term,'_pos.txt.txt',sep=''),sep='')   # filename for positive set
    neg.name <- paste(dir.name,paste(go.term,neg_pfx,sep=''),sep='') # filename for negative set
    full.data <- load.pos.neg.sets(pos.name,neg.name,group.info)
    
    ## feature extraction for each tissue type
    # cat('----------------------------------------\n')
    cat('NOT Reducing dimension of group features...\n')
    # dim.red <- reduce.features(full.data,ndim=ndim)
    dim.red <- full.data
    index <- c(NA, dim.red$group)
    full.x <- cbind(1, dim.red$x)  # add intercept
    full.y <- dim.red$y
    
    ## fit the data with coefficient path
    # plot.coefficient.path(full.x,full.y,index)
    
    ## split data
    # cat('----------------------------------------\n')
    data <- split_data(full.x,full.y)
    cat('dimensionality of data: ',dim(data$train$x)[2],'\n')
    cat('# of training samples:  ',length(data$train$y),'\n')
    cat('# of test samples:      ',length(data$test$y),'\n')
    
    ## apply cv for grplasso
    # cat('----------------------------------------\n')
    cat('Training group lasso classifier...\n')
    x <- as.matrix(data$train$x)
    y <- data$train$y
    cv.result <- cv.grplasso(x,y,index,nfolds=5,plot.error=FALSE)
    ## Re-fit the model with the best tuning paramter from cross-validation
    cat('Fitting model with lambda: ', cv.result$lambda)
    fit <- grplasso(x, y = y, index = index, lambda = cv.result$lambda, model = LogReg(), 
                    penscale = sqrt,control = grpl.control(update.hess = "lambda", trace = 0))
    cat('done\n')
    cat('cross-validaiton error:',cv.result$error,'\n')
    cat('parameter (lambda) tuned as:',cv.result$lambda,'\n')
    
    ## compute the test error
    # cat('----------------------------------------\n')
    prediction <- predict(fit, data$test$x, type = "response")
    auc.val <- auc(accuracy(prediction,as.factor(data$test$y)))
    cat('Test Error:',auc.val,'\n')
    
    ## store the coefficients of the fit
    sink(out.name)
    model <- 'group_lasso'
    cat('# Prediction results for:\t',go.term,'\n')
    cat('# Model used:\t',model,'\n')
    cat('# ROC AUC score:\t',auc.val,'\n')
    # cat('# Dimension per tissue:\t',ndim,'\n')
    if (length(full.data$types) > length(dim.red$types)) {
      cat('# Tissues used: ',length(dim.red$types),'\n')
    } else {
      cat('# All tissues were included\n')
    }
    cat('# Best penalty parameter (CV):', cv.result$lambda,'\n')
    grplasso.ceoff <- fit$coefficients
    for (i in 1:length(dim.red$types)) {
      cat('# tissue\t',i,'\t', dim.red$types[i],'\n')
    }
    cat('# Coefficients:\n')
    for (i in 1:length(dim.red$group)) {
      cat(dim.red$group[i],'\t',grplasso.ceoff[i+1],'\n')
    }
    cat('# Gene ID\tLabel\tPrediction\n')
    for (i in 1:length(prediction)) {
      cat(rownames(data$test$x)[i],'\t',data$test$y[i],'\t',prediction[i],'\n')
    }
    sink()
    # file.show(out.name)
    cat('saved result to output\n')
    # cat('----------------------------------------\n')
    
  }
}
