#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <R_ext/Applic.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"


// Compute ell2 norm -------------------------------------------
double l2_norm(double* a, double* b, int p) {
    double diff, sum = 0.0;
    for (int j = 0; j < p; j++) {
        diff = a[j] - b[j];  
        sum += diff * diff;         
    }
    return sqrt(sum);
}

// Data frame for sorting -------------------------------------------
typedef struct {
    double kinks;
    double slopes;
    double cumulative_sum;
    double value;
} Pair;

// Quick sort -------------------------------------------
int comparePairs(const void* a, const void* b) {
    const Pair* pairA = (const Pair*)a;
    const Pair* pairB = (const Pair*)b;
    if (pairA->kinks < pairB->kinks) return -1;
    if (pairA->kinks > pairB->kinks) return 1;
    return 0;
}

// Penalized Huber loss minimization for a fixed variable & fixed lambda -------------------------------------------
void huber_fit(double* out, double* rj, double* wxj, double* w2x2j, double* threshold_j, double lambda_fixed, double delta, int n) {

    // Set up the table -------------------------------------------
    Pair* pairs = malloc((2 * n) * sizeof(Pair));
    for (int i = 0; i < n; i++) {
        pairs[i].kinks = rj[i] - threshold_j[i];
        pairs[i].slopes = w2x2j[i] / (double) n;
    }

    for (int i = n; i < (2 * n); i++) {
        pairs[i].kinks = rj[i - n] + threshold_j[i - n];
        pairs[i].slopes = -w2x2j[i - n] / (double) n;
    }

    qsort(pairs, 2 * n, sizeof(Pair), comparePairs);
    double cusum = 0.0;
    for (int i = 0; i < (2 * n); i++) {
        cusum += pairs[i].slopes;
        pairs[i].cumulative_sum = cusum;
    }

    // Initialization -------------------------------------------
    double S = 0.0;
    for (int i = 0; i < n;i++) {
        S += -delta * fabs(wxj[i]) / (double)n;
    }

    if (pairs[0].kinks >= 0) {
        pairs[0].value = S + lambda_fixed;
    }
    else {
        pairs[0].value = S - lambda_fixed;
    }

    // Main loop -------------------------------------------
    double sol = 0.0;
    for (int i = 0; i < (2 * n - 1); i++) {
        S = S + pairs[i].cumulative_sum * (pairs[i + 1].kinks - pairs[i].kinks);
        if (pairs[i + 1].kinks >= 0) {
            pairs[i + 1].value = S + lambda_fixed;
        }
        else {
            pairs[i + 1].value = S - lambda_fixed;
        }

        if (i > 1) {
            if (pairs[i - 1].value < 0 && pairs[i].value >= 0) {
                sol = pairs[i - 1].kinks - pairs[i - 1].value * ((pairs[i].kinks - pairs[i - 1].kinks) / (pairs[i].value - pairs[i - 1].value));
                break;
            }
        }
    }
    *out = sol;
    free(pairs);
}

// Compute partial residuals (rj) -------------------------------------------
void partial_residuals(double* rj, double* y, double* x, double* beta_new, int n, int p, int j) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int k = 0; k < p; k++) {
            if (k != j) {
                sum += x[k* p + i] * beta_new[k];
            }
        }
        rj[i] = (y[i] - sum) / x[j * n + i];
    }
}

// Simple processing with assignment of nonconst -------------------------------------------
void simple_process(double* wx, double* w2x2, double* threshold, double* x, double* w, int n, int p, double delta){
    int i, j, jn; 
    for (j = 0; j < p; j++) {
        jn = j * n;
        for (i = 0; i < n; i++) {
            wx[jn + i] = w[i] * x[jn + i];
            threshold[jn + i] = delta / fabs(wx[jn + i]);
            w2x2[jn + i] = wx[jn + i]* wx[jn + i];
        }
    }
}

// Get cross-product -------------------------------------------
double crossprod(double* v, double* wx, int n, int j) {
    int jn = j * n, i; 
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        sum += wx[jn + i] * v[i];
    }
    return(sum);
}

// Get maximum of cross-product -------------------------------------------
double maxprod(double* v, double* wx, int n, int p) {
    int j; 
    double z, tmp, max = 0.0;
    for (j = 0; j < p; j++) {
        tmp = crossprod(v, wx, n, j);
        z = fabs(tmp);
        if (z > max) {
            max = z;
        }
    }
    return(max);
}

// get the function of Huber  -------------------------------------------
void function_huber(double* u, double* v, double delta, int n) {
    int i; 
    double abs = 0.0;
    for (i = 0; i < n; i++) {
        abs = fabs(v[i]);
        if (abs <= delta) {
            u[i] = 0.5 * abs * abs;
        }
        else {
            u[i] = delta * abs - 0.5 * delta * delta;
        }
    }
}

// get the derivative of Huber  -------------------------------------------
void derivative_huber(double* u, double* v, double delta, int n) {
    int i; 
    double tmp = 0.0;
    for (i = 0; i < n; i++) {
        tmp = fabs(v[i]);
        if (tmp <= delta) {
            u[i] = v[i];
        }
        else {
            u[i] = delta * ((v[i] > 0) - (v[i] < 0));
        }
    }
}

// Exact Coordinate Descent (ECD) -------------------------------------------
static void ecd_huber(double* beta, 
                int* iter, 
                double* lambda, 
                double* x, 
                double* y,
                double* w,
                double* delta_,
                double* tol_,
                double* lambda_min_, 
                int* nlam_,
                int* n_, 
                int* p_, 
                int* ppflag_, 
                int* max_iter_,
                int* user_)
{
    // Declarations -------------------------------------------
    double tol = tol_[0]; 
    double gap_tol = 1e-6;
    double delta = delta_[0];
    double lambda_min = lambda_min_[0];
    int nlam = nlam_[0]; 
    int n = n_[0]; 
    int p = p_[0]; 
    int ppflag = ppflag_[0]; 
    // int scrflag = scrflag_[0]; 
    int max_iter = max_iter_[0]; 
    int user = user_[0];


    int i, j, l, jn, ll;
    double lstep, bj, kkt_j, r_ij, lambda_fixed;
    
    double* gap = malloc(max_iter * sizeof(double));
    double* wx = malloc(n * p * sizeof(double)); // wx
    double* w2x2 = malloc(n * p * sizeof(double)); // w^2x^2
    double* threshold = malloc(n * p * sizeof(double)); // delta/|wx|
    double* threshold_j = malloc(n * sizeof(double)); // delta/|wx_j|
    double* beta_old = malloc(p * sizeof(double));
    double* beta_lag = malloc(p * sizeof(double));
    double* beta_new = malloc(p * sizeof(double));
    double* r = malloc(n * sizeof(double));
    double* rj = malloc(n * sizeof(double));
    double* xj = malloc(n * sizeof(double));
    double* wxj = malloc(n * sizeof(double));
    double* w2x2j = malloc(n * sizeof(double));
    // int* include = R_Calloc(p, int);

    // Preprocessing -------------------------------------------
    // This needs to be added.
    if (ppflag == 0) {
        simple_process(wx, w2x2, threshold, x, w, n, p, delta);
    }

    // Screening rule -------------------------------------------
    // This needs to be added (not used now).
    //if (scrflag == 0) {
    //    for (j = 0; j < p; j++) {
    //        include[j] = 1;
    //    }
    //}

    // Set up lambda -------------------------------------------
   derivative_huber(r, y, delta, n);
    if (user == 0) {
        lambda[0] = maxprod(r, wx, n, p) / (double) n;
        if (lambda_min == 0.0) {
            lambda_min = 0.001;
        }
        lstep = log(lambda_min) / ( (double) nlam - 1.0);
        for (l = 1; l < nlam; l++) {
            lambda[l] = lambda[l - 1] * exp(lstep);
        }
    }


    // Solution path -------------------------------------------
    for (l = 0; l < nlam; l++) {

        lambda_fixed = lambda[l];
        // Call up previous beta -------------------------------------------
        if (l == 0) {
            for (j = 0; j < p; j++) {
                beta_old[j] = 0;
                beta_lag[j] = beta_old[j];
                beta_new[j] = beta_old[j];
            }
        }
        else {
            for (j = 0; j < p; j++) {
                beta_old[j] = beta[(l - 1) * p + j];
                beta_lag[j] = beta_old[j];
                beta_new[j] = beta_old[j];
            }
        }

        // For each fixed lambda, -------------------------------------------
        ll = 0;
        while (ll < max_iter) {

            for (j = 0; j < p; j++) {

                // Get bj, threshold_j & xj -------------------------------------------
                bj = beta_old[j];
                jn = j * n;
                for (i = 0; i < n; i++) {
                    threshold_j[i] = threshold[jn + i];
                    xj[i] = x[jn + i];
                    wxj[i] = wx[jn + i];
                    w2x2j[i] = w2x2[jn + i];
                }
                // compute rj -------------------------------------------
                partial_residuals(rj, y, x, beta_new, n, p, j);

                // Check KKT conditions -------------------------------------------
                kkt_j = 0.0;
                for (i = 0;i < n;i++) {
                    r_ij = rj[i] - bj;
                    if (fabs(r_ij) <= threshold_j[i]) {
                        kkt_j += w2x2j[i] * r_ij / (double)n;
                    }
                    else {
                        kkt_j += fabs(wxj[i]) * delta * ((r_ij > 0) - (r_ij < 0)) / (double)n;
                    }
                }
                if (fabs(kkt_j) <= lambda_fixed) {
                    beta_new[j] = beta_old[j];
                }
                else {
                    double out;
                    huber_fit(&out, rj, wxj, w2x2j, threshold_j, lambda_fixed, delta, n);
                    beta_new[j] = out;
                }
            } // coordinate cycle ends
            ll += 1;

            // Check the solution updates -------------------------------------------
            gap[ll] = l2_norm(beta_new, beta_old, p);
            if (gap[ll] < tol || ll == max_iter) {
                for (j = 0;j < p;j++) {
                    beta[l * p + j] = beta_new[j];
                }
                iter[l] = ll;
                break;
            }
            else {
                if (ll >= 2) {
                    if (fabs(gap[ll] - gap[(ll-2)]) < gap_tol) {
                        for (j = 0;j < p;j++) {
                            beta[l * p + j] = beta_new[j];
                        }
                        iter[l] = ll;
                        break;
                    }
                    else {
                        for (j = 0;j < p;j++) {
                            beta_old[j] = beta_new[j];
                        }
                    }
                }
                else {
                    for (j = 0;j < p;j++) {
                        beta_old[j] = beta_new[j];
                    }
                }
            }

        }// fixed lambda cycle ends
    }

    // Free allocated memories -------------------------------------------
    free(gap);
    free(wx);
    free(w2x2);
    free(threshold);
    free(threshold_j);
    free(beta_old);
    free(beta_lag);
    free(beta_new);
    free(r);
    free(rj);
    free(xj);
    free(wxj);
    free(w2x2j);
    //R_Free(include);
}


static const R_CMethodDef cMethods[] = {
    // name        pointer         Num args
    {"ecd_huber", (DL_FUNC)&ecd_huber, 15},
    {NULL, NULL, 0}   // Placeholder to indicate last one.
};


void R_init_rome(DllInfo* info)
{
	R_registerRoutines(
        info,      // DllInfo
        cMethods,      // .C
        NULL,  // .Call
        NULL,      // Fortran
        NULL       // External
    );
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}