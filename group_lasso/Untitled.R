rm(list=ls())
library(grplasso)
## Use the Logistic Group Lasso on the splice data set
data(splice)
## Define a list with the contrasts of the factors
contr <- rep(list("contr.sum"), ncol(splice) - 1)
names(contr) <- names(splice)[-1]