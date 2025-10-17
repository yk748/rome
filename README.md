# rome(RObust M-Estimator)

We provide efficient procedures for fitting robust M-estimators with a Lasso penalty. The main idea is to extend the exact coordinate descent algorithm ([2010](#ref-glmnet)) to robust loss functions, such as the Huber loss.

The current version (1.0.0) provides basic functionalities, including model fitting, cross-validation, and illustrative examples. This version includes only the $L_1$-penalized Huber regression but will be continuously updated to incorporate additional loss and penalty functions. The algorithm is decribed in Kim, Loh, and Basu ([2025](#ref-rome)).

Please email Younghoon Kim
<yk748@cornell.edu> if you encounter any bugs or errors. All errors are my own.

## References

<div id="ref-rome">

Kim, Younghoon, Loh, Po-Ling, and Basu, Sumanta. 2025.
“Exact Coordinate Descent for High-Dimensional Regularized Regularized Huber Regression.” <https://arxiv.org/abs/2510.13715>.

<div id="ref-glmnet">

Friedman, Jerome, Hastie, Trevor, and Tibshirani, Robert. 2010.
“Regularization Paths for Generalized Linear Models via Coordinate
Descent.” *Journal of Statistical Software, Articles* 33 (1): 1–22.
<https://doi.org/10.18637/jss.v033.i01>.
