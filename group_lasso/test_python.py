import numpy as np
np.random.seed(42)
# A three-dimensional matrix is defined as:
shape = (4, 4, 4)
# The number of samples is defined as:
num_samples = 50
# The number of features per sample is defined as:
num_ft = shape[0] * shape[1] * shape[2]
# Define X randomly as simulated data
X = np.random.rand(num_samples, num_ft)
# Define y as zeros or ones
y = np.random.randint(0, 2, (num_samples, 1))

import parsimony.estimators as estimators
import parsimony.algorithms as algorithms
import parsimony.functions.nesterov.gl as gl
k = 0.0  # l2 ridge regression coefficient
l = 0.1  # l1 lasso coefficient

groups = [range(0, 2 * num_ft / 3), range(num_ft/ 3, num_ft)]
print groups
A = gl.linear_operator_from_groups(num_ft, groups)

lambdas = [1e-8, 1e-4, 1, 1e3, 1e10];
for g in lambdas:
    print g
    # g = 0.1  # group lasso coefficient
    estimator = estimators.LogisticRegressionL1L2GL(k, l, g, A=A, algorithm=algorithms.proximal.FISTA(), algorithm_params=dict(max_iter=1000))
    print estimator
    res = estimator.fit(X, y)
    # print "Estimated prediction rate =", estimator.score(X, y)
    print "Prediction error = ", estimator.score(X, y)

# print estimator.beta
