# rome(RObust M-Estimator)

We provide efficient procedures for fitting the robust M-estimators with Lasso penalty. The idea is to bring the
exact coordinate descent algorithm ([2010](#ref-glmnet)) into robust loss functions, such as Huber loss function.

Currently, the beta version (0.1.0) is released, and the version
provides the basic functions, including fit, cross-validation, and other
functions for illustration. Please email Younghoon Kim
<yk748@cornell.edu> if any bugs/errors have been discovered.

## References

<div id="refs" class="references">

<div id="ref-glmnet">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software, Articles* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.
