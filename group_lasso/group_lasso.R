rm(list=ls())

# load libraries
library(grplasso)
library(AUC)
library(ggplot2)

cv.grplasso <- function(x, y, index, nfolds = 5, nlamdas = 20, plot.error=FALSE) {
  # select the parameters
  lambda <- lambdamax(x, y = y, index = index, penscale = sqrt, model = LogReg()) * 0.5^seq(0,5,len=nlamdas)
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
  # train_set_size: Fraction of genes used for training set
  num_samples = length(labels)
  num_features = dim(gene_features)[1]
  num_train_samples = ceiling(train_set_size*num_samples)
  sample_order = sample(1:num_samples, num_samples) # permute the samples
  train_idx = sample_order[1:num_train_samples]
  test_idx =  sample_order[(num_train_samples+1):num_samples]
  
  train = list(x=gene_features[train_idx,],y=labels[train_idx])
  test = list(x=gene_features[test_idx,],y=labels[test_idx])
  return(list(train=train, test=test))
}

group.to.int <- function(group2col, specific=TRUE) {
  if (specific) {
    group = group2col[,2]
  } else {
    group = group2col[,1]
  }
  type.names = unique(group)
  type.counts = rep(0,length(type.names))
  group.idx = rep(0,dim(group2col)[1])
  for (i in 1:length(group.idx)) {
    type.idx = which(type.names == group[i])
    group.idx[i] = type.idx
    type.counts[type.idx] = type.counts[type.idx] + 1
  }
  names(type.counts) = type.names
  return(list(types=type.counts,idx=group.idx))
}

load.pos.neg.sets <- function(pos.name,neg.name,grp.name,specific=TRUE,transform=FALSE) {
  pos.data <- read.table(pos.name,sep='\t',row.names=1,skip=2)
  neg.data <- read.table(neg.name,sep='\t',row.names=1,skip=2)
  group    <- read.table(grp.name,sep='\t',row.names=1,skip=1) 
  # group: tissue type, tissue specific type
  data <- rbind(pos.data,neg.data)
  response <- c(rep(1,dim(pos.data)[1]),rep(0,dim(neg.data)[1]))
  group.info <- group.to.int(group,specific)
  
  if (transform) {
    # log10(x+1) and standardize each feature
    # iterate through each column
    for (i in 1:dim(data)[2]) {
      feature <- data[,i] 
      log.feat <- apply(t(feature),1,function(x) log10(x+1))
      log.feat <- log.feat - mean(log.feat)  # centering
      if (sd(log.feat) > 1e-10)  {           # scaling
        log.feat <- log.feat / sd(log.feat)
      }
      data[,i] <- log.feat  # update the feature
    }
  }
  
  return(list(x=data,
              y=response,
              group=group.info$idx,
              types=group.info$types))
}

reduce.features <- function(grouped.data, ndim=3) {
  # reduce the dimension of each feature to ndim, and remove features that have less than ndim
  type.counts <- grouped.data$types
  x <- grouped.data$x
  if (dim(x)[1] < ndim) {
    stop('Error: # of samples is less than # of subtype dimentions: ', dim(x)[1],'<',ndim)
  }
  nusable <- 0
  for (i in 1:length(type.counts)) {
    if (type.counts[i] < ndim) {
      cat('Warning: ',names(type.counts)[i],' has too few dimensions and not included in the model\n')
    } else {
      nusable = nusable + 1
    }
  }
  new.x <- matrix(nrow=dim(x)[1], ncol=ndim*nusable)
  new.groups <- numeric(ndim*nusable)
  new.types <- numeric(0)
  for (i in 1:length(type.counts)) {
    if (type.counts[i] < ndim) {
      next
    } else {
      new.types <- c(names(type.counts)[i],new.types)
    }
    i.type <- length(new.types)
    sel <- which(grouped.data$group == i)
    # check if any columns are zero
    remove <- numeric(0)
    for (j in 1:length(sel)) {
      if ( sd(x[,sel[j]]) < 0.001) {
        # print(x[,sel[j]])
        remove <- c(remove,j) # remove the entry
      }
    }
    if (length(remove)>0) { 
      sel <- sel[-remove] 
      # cat('Warning: removed columns: ',remove,' with zero variance\n')
    } 
    pca <- prcomp(x[ ,sel],center = TRUE, scale. = TRUE) 
    new.x[ ,((i.type-1)*ndim+1):(i.type*ndim)] <- pca$x[ ,1:ndim]
    new.groups[((i.type-1)*ndim+1):(i.type*ndim)] <- rep(i.type,ndim)
  }
  # print(new.types)
  # print(new.groups)
  # print(dim(new.x))
  # stop('lol')
  rownames(new.x) <- rownames(grouped.data$x)
  return(list(x=new.x,
              y=grouped.data$y,
              group=new.groups,
              types=new.types))
}

plot.coefficient.path <- function(x,y,index) {
  lambda <- lambdamax(x, y = y, index = index, penscale = sqrt, model = LogReg()) * 0.5^seq(0,5,len=10)
  fit <- grplasso(x = x, y = y, index = index, lambda = lambda, model = LogReg(), penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0))
  plot(fit,log='x')
}

## set seed for reproducability
set.seed(2)

main <- '/Users/jasonzhu/Documents/CS341_Code/'
dir.name <- paste(main,'data/experiment_inputs/',sep='')
grp.name <- paste(main,'data/samples_to_tissues_map.txt',sep='')
all.go.name <- paste(main,'data/GO_terms_final_gene_counts.txt',sep='')
go.names <- read.table(all.go.name)$V1
ndim <- 5

for (go.idx in 529:length(go.names)){
  go.term <- as.character(go.names[go.idx])
  # go.term <- 'GO:0000578'
  neg_idx <- 0
  neg_pfx <- paste(paste('_neg_',neg_idx,sep=''),'.txt.txt',sep='')
  out_pfx <- paste(paste(go.term,'_',sep=''),neg_idx,sep='')
  out.name <- paste(paste(paste(main,'data/grplasso_results/grplasso_',sep=''),out_pfx,sep=''),'.txt',sep='')

  ## load all data
  cat('----------------------------------------\n')
  cat('Loading raw data from',go.term,'...\n')
  pos.name <- paste(dir.name,paste(go.term,'_pos.txt.txt',sep=''),sep='')   # filename for positive set
  neg.name <- paste(dir.name,paste(go.term,neg_pfx,sep=''),sep='') # filename for negative set
  full.data <- load.pos.neg.sets(pos.name,neg.name,grp.name,specific=TRUE,transform=TRUE)
  
  ## feature extraction for each tissue type
  # cat('----------------------------------------\n')
  cat('Reducing dimension of group features...\n')
  dim.red <- reduce.features(full.data,ndim=ndim)
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
  cat('Training group lasso classifier...')
  x <- data$train$x
  y <- data$train$y
  cv.result <- cv.grplasso(x,y,index,nfolds=3,plot.error=FALSE)
  ## Re-fit the model with the best tuning paramter from cross-validation
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
  cat('# Dimension per tissue:\t',ndim,'\n')
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
