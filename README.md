# rome(RObust M-Estimator)

We provide efficient procedures for fitting the robust M-estimators with Lasso penalty. The idea is to bring the
exact coordinate descent algorithm ([2010](#ref-glmnet)) into robust loss functions, such as Huber loss function.

Currently, the beta version (0.8.0) is released, and the version
provides the basic functions, including fit, cross-validation, and other
functions for illustration. Please email Younghoon Kim
<yk748@cornell.edu> if any bugs/errors have been discovered.

## References

<div id="ref-rome">

Kim, Younghoon, Basu, Sumanta, and Loh, Po-Ling. 2025.
“Exact Coordinate Descent for High-Dimensional Regularized Robust M-Estimators.” *In progress*.

<div id="ref-glmnet">

Friedman, Jerome, Hastie, Trevor, and Tibshirani, Robert. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software, Articles* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.
