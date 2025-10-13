#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

/* Huber loss gradient */
double huber_grad(double v, double thresh) {
    if (fabs(v) <= thresh) return v;
    return thresh * (v > 0 ? 1.0 : -1.0);
}

/* Soft-thresholding operator */
double soft_threshold(double v, double thres) {
    if (v > thres) return v - thres;
    else if (v < -thres) return v + thres;
    else return 0.0;
}

/* Matrix-vector product: y = X * beta (column-major) */
void matvec_mult(const double* X, const double* beta, double* y, int n, int p) {
    for (int i = 0; i < n; i++) y[i] = 0.0;
    for (int j = 0; j < p; j++)
        for (int i = 0; i < n; i++)
            y[i] += X[i + j * n] * beta[j];
}

/* X^T * v (crossprod, column-major) */
void crossprod(const double* X, const double* v, double* out, int n, int p) {
    for (int j = 0; j < p; j++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++)
            sum += X[i + j * n] * v[i];
        out[j] = sum;
    }
}

/* Euclidean norm of vector difference */
double diff_norm(const double* a, const double* b, int p) {
    double sum = 0.0;
    for (int j = 0; j < p; j++) {
        double diff = a[j] - b[j];
        sum += diff * diff;
    }
    return sqrt(sum);
}

/* Power iteration to approximate largest eigenvalue of X^T X / n */
double max_eigen_XtX(const double* X, int n, int p) {
    double* b = calloc(p, sizeof(double));
    double* Xb = calloc(n, sizeof(double));
    for (int j = 0; j < p; j++) b[j] = 1.0; // initial guess

    double lambda = 0.0;
    for (int iter = 0; iter < 100; iter++) {
        /* Xb = X * b */
        for (int i = 0; i < n; i++) Xb[i] = 0.0;
        for (int j = 0; j < p; j++)
            for (int i = 0; i < n; i++)
                Xb[i] += X[i + j * n] * b[j];

        /* b_new = X^T * Xb / n */
        double norm = 0.0;
        for (int j = 0; j < p; j++) {
            double sum = 0.0;
            for (int i = 0; i < n; i++)
                sum += X[i + j * n] * Xb[i];
            b[j] = sum / n;
            norm += b[j] * b[j];
        }

        norm = sqrt(norm);
        for (int j = 0; j < p; j++) b[j] /= norm;  // normalize

        lambda = norm;  // approximate largest eigenvalue
    }

    free(b); free(Xb);
    return lambda;
}

/* Core composite gradient descent function */
void comp_grad_huber_l1_(
    const double* X, const double* y, int n, int p,
    double delta, double lambda, double step_size,
    int max_iter, double tol, int verbose,
    double* beta_out
) {
    double* beta = calloc(p, sizeof(double));
    double* beta_new = calloc(p, sizeof(double));
    double* r = malloc(n * sizeof(double));
    double* grad_r = malloc(n * sizeof(double));
    double* grad_beta = malloc(p * sizeof(double));

    /* Step size computation */
    if (step_size <= 0) {
        double L = max_eigen_XtX(X, n, p);
        step_size = 1.0 / L;
        if (verbose) Rprintf("Computed step size: %g\n", step_size);
    }

    for (int iter = 0; iter < max_iter; iter++) {
        /* r = y - X * beta */
        matvec_mult(X, beta, r, n, p);
        for (int i = 0; i < n; i++) r[i] = y[i] - r[i];

        /* grad_r = huber_grad(r, delta) */
        for (int i = 0; i < n; i++) grad_r[i] = huber_grad(r[i], delta);

        /* grad_beta = - X^T * grad_r / n */
        crossprod(X, grad_r, grad_beta, n, p);
        for (int j = 0; j < p; j++) grad_beta[j] = -grad_beta[j] / n;

        /* beta_new = soft_threshold(beta - step_size * grad_beta, lambda * step_size) */
        for (int j = 0; j < p; j++) {
            double z = beta[j] - step_size * grad_beta[j];
            beta_new[j] = soft_threshold(z, lambda * step_size);
        }

        /* Check convergence */
        if (diff_norm(beta_new, beta, p) < tol) {
            if (verbose) Rprintf("Converged at iteration %d\n", iter + 1);
            break;
        }

        /* Update beta */
        for (int j = 0; j < p; j++) beta[j] = beta_new[j];
    }

    /* Copy result */
    for (int j = 0; j < p; j++) beta_out[j] = beta[j];

    free(beta); free(beta_new); free(r); free(grad_r); free(grad_beta);
}

/* .Call wrapper for R */
SEXP comp_grad_huber_l1(SEXP X_, SEXP y_, SEXP delta_, SEXP lambda_,
    SEXP step_size_, SEXP max_iter_, SEXP tol_, SEXP verbose_) {

    /* Ensure numeric vectors/matrices */
    PROTECT(X_ = coerceVector(X_, REALSXP));
    PROTECT(y_ = coerceVector(y_, REALSXP));

    /* Get matrix dimensions */
    SEXP dim = getAttrib(X_, R_DimSymbol);
    int n = INTEGER(dim)[0];
    int p = INTEGER(dim)[1];

    /* Extract parameters from SEXP */
    double delta = asReal(delta_);
    double lambda = asReal(lambda_);
    double step_size = asReal(step_size_);
    int max_iter = asInteger(max_iter_);
    double tol = asReal(tol_);
    int verbose = asInteger(verbose_);

    /* Allocate output vector */
    SEXP beta_out = PROTECT(allocVector(REALSXP, p));

    /* Call the computational function */
    comp_grad_huber_l1_(REAL(X_), REAL(y_), n, p,
        delta, lambda, step_size,
        max_iter, tol, verbose,
        REAL(beta_out));

    UNPROTECT(3);  /* X_, y_, beta_out */
    return beta_out;
}