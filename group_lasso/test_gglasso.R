# install.packages("gglasso")

# load gglasso library
library(gglasso)
# load data set
data(colon)

# define group index
group <- rep(1:20,each=5)
# plot coefficient path from gglasso fit
m2 <- gglasso(x=colon$x, y=colon$y, group=group, loss="logit")
par(mfrow=c(1,3)) 
plot(m2) # plots the coefficients against the log-lambda sequence 
plot(m2,group=TRUE) # plots group norm against the log-lambda sequence 
plot(m2,log.l=FALSE) # plots against the lambda sequence


# 5-fold cross validation using group lasso
# penalized logisitic regression
cv <- cv.gglasso(x=colon$x, y=colon$y, group=group, loss="logit",
                 pred.loss="misclass", lambda.factor=0.05, nfolds=5)
# the coefficients at lambda = lambda.min, newx = x[1,]
pre = predict(cv$gglasso.fit, newx = colon$x[1:10,],
              s = cv$lambda.min, type = "class")

plot(cv)


# make plots
par(mfrow=c(1,3))
plot(m1) # plots the coefficients against the log-lambda sequence
plot(m1,group=TRUE) # plots group norm against the log-lambda sequence
plot(m1,log.l=FALSE) # plots against the lambda sequence
