#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>

// Utility functions ----------------------------------------------------
static double huber_loss_(double v, double delta) {
    if (fabs(v) <= delta) return 0.5 * v * v;
    return delta * (fabs(v) - 0.5 * delta);
}

static double huber_grad_(double v, double delta) {
    if (fabs(v) <= delta) return v;
    return delta * (v > 0 ? 1.0 : -1.0);
}

static double soft_thresh_(double z, double thresh) {
    if (z > thresh) return z - thresh;
    if (z < -thresh) return z + thresh;
    return 0.0;
}

// Objective function ----------------------------------------------------
static double objective_(const double* beta, const double* X,
    const double* y, int n, int p,
    double delta, double lambda) {
    double loss = 0.0;
    for (int i = 0; i < n; i++) {
        double r = y[i];
        for (int j = 0; j < p; j++) r -= X[i + j * n] * beta[j];
        loss += huber_loss_(r, delta);
    }
    loss /= n;
    double pen = 0.0;
    for (int j = 0; j < p; j++) pen += fabs(beta[j]);
    return loss + lambda * pen;
}

// Power iteration: lambda_max(X'X / n)  ----------------------------------------------------
static double max_eigen_(const double* X, int n, int p) {
    double* b = calloc(p, sizeof(double));
    double* Xb = calloc(n, sizeof(double));
    for (int j = 0; j < p; j++) b[j] = 1.0;
    double lam = 0.0;
    for (int it = 0; it < 200; it++) {
        for (int i = 0; i < n; i++) {
            Xb[i] = 0.0;
            for (int j = 0; j < p; j++) Xb[i] += X[i + j * n] * b[j];
        }
        double norm = 0.0;
        for (int j = 0; j < p; j++) {
            double s = 0.0;
            for (int i = 0; i < n; i++) s += X[i + j * n] * Xb[i];
            b[j] = s / n;
            norm += b[j] * b[j];
        }
        norm = sqrt(norm);
        for (int j = 0; j < p; j++) b[j] /= norm;
        lam = norm;
    }
    free(b); free(Xb);
    return lam;
}

// Exact 1-D coordinate solver (fun_cd)  ----------------------------------------------------
typedef struct { double kink; double B; } KinkPair;

static int cmp_kink(const void* a, const void* b) {
    double d = ((KinkPair*)a)->kink - ((KinkPair*)b)->kink;
    return (d > 0) - (d < 0);
}

static double fun_cd_(const double* xj, const double* rj,
    const double* thres_j, int n,
    double delta, double lambda) {
    KinkPair* kp = malloc(2 * n * sizeof(KinkPair));
    for (int i = 0; i < n; i++) {
        kp[i].kink = rj[i] - thres_j[i];
        kp[i].B = xj[i] * xj[i] / n;
        kp[i + n].kink = rj[i] + thres_j[i];
        kp[i + n].B = -xj[i] * xj[i] / n;
    }
    qsort(kp, 2 * n, sizeof(KinkPair), cmp_kink);

    double* slope = malloc(2 * n * sizeof(double));
    double cusum = 0.0;
    for (int i = 0; i < 2 * n; i++) {
        cusum += kp[i].B;
        slope[i] = cusum;
    }

    double Sj = 0.0;
    for (int i = 0; i < n; i++) Sj -= delta * fabs(xj[i]) / n;
    double deriv = Sj
        + (kp[0].kink >= 0 ? lambda : 0.0)
        + (kp[0].kink <= 0 ? -lambda : 0.0);

    if (deriv >= 0.0) { free(kp); free(slope); return 0.0; }

    double out = R_NaN;
    for (int i = 0; i < 2 * n - 1; i++) {
        Sj += slope[i] * (kp[i + 1].kink - kp[i].kink);
        double deriv_next = Sj
            + (kp[i + 1].kink >= 0 ? lambda : 0.0)
            + (kp[i + 1].kink <= 0 ? -lambda : 0.0);
        if (deriv < 0.0 && deriv_next >= 0.0) {
            out = kp[i].kink - deriv * (kp[i + 1].kink - kp[i].kink)
                / (deriv_next - deriv);
            break;
        }
        deriv = deriv_next;
    }
    free(kp); free(slope);
    return out;
}


// Gradient Descent (GD) ----------------------------------------------------
static void run_gd_(const double* X, const double* y, int n, int p,
    double delta, double lambda, int max_iter,
    double* obj_path) {

    double* beta = calloc(p, sizeof(double));
    double* r = malloc(n * sizeof(double));
    double* grad_r = malloc(n * sizeof(double));
    double* grad_b = malloc(p * sizeof(double));
    double* beta_new = malloc(p * sizeof(double));

    double L = max_eigen_(X, n, p);
    double eta = 1.0 / L;

    for (int iter = 0; iter < max_iter; iter++) {

        for (int i = 0; i < n; i++) {
            r[i] = y[i];
            for (int j = 0; j < p; j++) r[i] -= X[i + j * n] * beta[j];
        }

        for (int i = 0; i < n; i++) grad_r[i] = huber_grad_(r[i], delta);
        for (int j = 0; j < p; j++) {
            double s = 0.0;
            for (int i = 0; i < n; i++) s += X[i + j * n] * grad_r[i];
            grad_b[j] = -s / n;
        }

        for (int j = 0; j < p; j++)
            beta[j] = soft_thresh_(beta[j] - eta * grad_b[j], lambda * eta);

        obj_path[iter] = objective_(beta, X, y, n, p, delta, lambda);
    }
    free(beta); free(r); free(grad_r); free(grad_b); free(beta_new);
}

// Coordinate Gradient Descent (CGD) ----------------------------------------------------
static void run_cgd_(const double* X, const double* y, int n, int p,
    double delta, double lambda, int max_iter,
    double* obj_path) {

    double* beta = calloc(p, sizeof(double));
    double* r = malloc(n * sizeof(double));
    double* a = malloc(p * sizeof(double));  /* a_j = (1/n) sum X_ij^2 */

    for (int j = 0; j < p; j++) {
        double s = 0.0;
        for (int i = 0; i < n; i++) s += X[i + j * n] * X[i + j * n];
        a[j] = s / n;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        /* Refresh residual at start of sweep */
        for (int i = 0; i < n; i++) {
            r[i] = y[i];
            for (int j = 0; j < p; j++) {
                r[i] -= X[i + j * n] * beta[j];
            }
        }
        for (int j = 0; j < p; j++) {
            double grad_j = 0.0;
            for (int i = 0; i < n; i++) {
                grad_j -= huber_grad_(r[i], delta) * X[i + j * n];
            }
            grad_j /= n;

            double mu = 1.0 / a[j];
            double beta_j_new = soft_thresh_(beta[j] - mu * grad_j, lambda * mu);
            double change = beta_j_new - beta[j];

            /* Rank-1 residual update */
            if (change != 0.0)
                for (int i = 0; i < n; i++) r[i] -= X[i + j * n] * change;

            beta[j] = beta_j_new;
        }
        obj_path[iter] = objective_(beta, X, y, n, p, delta, lambda);
    }
    free(beta); free(r); free(a);
}

// Exact Coordinate Descent (CD, proposed) ----------------------------------------------------
static void run_cd_(const double* X, const double* y, int n, int p,
    double delta, double lambda, int max_iter,
    double* obj_path) {

    double* beta = calloc(p, sizeof(double));
    double* r = malloc(n * sizeof(double));
    double* thres = malloc(n * p * sizeof(double));  /* delta/|x_ij| */
    double* xj_buf = malloc(n * sizeof(double));
    double* rj_buf = malloc(n * sizeof(double));
    double* thres_j = malloc(n * sizeof(double));

    for (int j = 0; j < p; j++)
        for (int i = 0; i < n; i++)
            thres[i + j * n] = delta / fabs(X[i + j * n]);

    for (int iter = 0; iter < max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            r[i] = y[i];
            for (int j = 0; j < p; j++) {
                r[i] -= X[i + j * n] * beta[j];
            }
        }

        for (int j = 0; j < p; j++) {
            for (int i = 0; i < n; i++) {
                xj_buf[i] = X[i + j * n];
                thres_j[i] = thres[i + j * n];
                rj_buf[i] = (r[i] + xj_buf[i] * beta[j]) / xj_buf[i];
            }

            double sum_j = 0.0; double sign_multiplier = -1.0;
            for (int i = 0; i < n; i++) {
                double r_ij = rj_buf[i] - beta[j];
                if (fabs(r_ij) <= thres_j[i]) {
                    sum_j += xj_buf[i] * xj_buf[i] * r_ij;
                }
                else {
                    if (r_ij > 0) {
                        sign_multiplier = 1.0;
                    }
                    sum_j += fabs(xj_buf[i]) * delta * sign_multiplier;
                }
                    
            }
            sum_j /= n;

            double beta_j_new;
            if (fabs(sum_j) <= lambda) {
                beta_j_new = beta[j];
            }
            else {
                double fit = fun_cd_(xj_buf, rj_buf, thres_j, n, delta, lambda);
                if (R_IsNaN(fit)) {
                    beta_j_new = beta[j];
                }
                else {
                    beta_j_new = fit;
                }
            }

            double change = beta_j_new - beta[j];
            if (change != 0.0)
                for (int i = 0; i < n; i++) {
                    r[i] -= xj_buf[i] * change;
                }
            beta[j] = beta_j_new;
        }
        obj_path[iter] = objective_(beta, X, y, n, p, delta, lambda);
    }
    free(beta); free(r); free(thres); free(xj_buf); free(rj_buf); free(thres_j);
}

// .Call wrappers ----------------------------------------------------
SEXP run_gd(SEXP X_, SEXP y_, SEXP delta_, SEXP lambda_, SEXP max_iter_) {
    PROTECT(X_ = coerceVector(X_, REALSXP));
    PROTECT(y_ = coerceVector(y_, REALSXP));
    SEXP dim = getAttrib(X_, R_DimSymbol);
    int n = INTEGER(dim)[0], p = INTEGER(dim)[1];
    int max_iter = asInteger(max_iter_);
    SEXP out = PROTECT(allocVector(REALSXP, max_iter));
    run_gd_(REAL(X_), REAL(y_), n, p,
        asReal(delta_), asReal(lambda_), max_iter, REAL(out));
    UNPROTECT(3); return out;
}

SEXP run_cgd(SEXP X_, SEXP y_, SEXP delta_, SEXP lambda_, SEXP max_iter_) {
    PROTECT(X_ = coerceVector(X_, REALSXP));
    PROTECT(y_ = coerceVector(y_, REALSXP));
    SEXP dim = getAttrib(X_, R_DimSymbol);
    int n = INTEGER(dim)[0], p = INTEGER(dim)[1];
    int max_iter = asInteger(max_iter_);
    SEXP out = PROTECT(allocVector(REALSXP, max_iter));
    run_cgd_(REAL(X_), REAL(y_), n, p,
        asReal(delta_), asReal(lambda_), max_iter, REAL(out));
    UNPROTECT(3); return out;
}

SEXP run_cd(SEXP X_, SEXP y_, SEXP delta_, SEXP lambda_, SEXP max_iter_) {
    PROTECT(X_ = coerceVector(X_, REALSXP));
    PROTECT(y_ = coerceVector(y_, REALSXP));
    SEXP dim = getAttrib(X_, R_DimSymbol);
    int n = INTEGER(dim)[0], p = INTEGER(dim)[1];
    int max_iter = asInteger(max_iter_);
    SEXP out = PROTECT(allocVector(REALSXP, max_iter));
    run_cd_(REAL(X_), REAL(y_), n, p,
        asReal(delta_), asReal(lambda_), max_iter, REAL(out));
    UNPROTECT(3); return out;
}