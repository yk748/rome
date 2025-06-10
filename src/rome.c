#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <R_ext/Applic.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"

// From util.c -------------------------------------------
double huber_loss(double v, double thresh);
double huber_grad(double v, double thresh);
void simple_process(double* v, double* x2, double* threshold, double* y, double* x, int n, int p, double delta);
void simple_process_adaptive(double* v, double* wx, double* w2x2, double* threshold, double* y, double* x, double* w, int n, int p, double delta);

// From rome_fit.c -------------------------------------------
int comparePairs(const void* a, const void* b);
void huber_fit(double* out, double* rj, double* first_j, double* second_j, double* threshold_j, double lambda_fixed, double delta, int n, double N);

// Exact Coordinate Descent (ECD) for adaptive penalized Huber regression -------------------------------------------
void ecd_huber_adaptive_active_(double* beta, int* iter, double* lambda, int* saturated, int* numv, double* x, double* y, double* w, double* delta_, double* eps_, double* lambda_min_,
    int* nlam_, int* n_, int* p_, int* ppflag_, int* dfmax_, int* max_iter_, int* user_, int* scrflag_, int* trace_)
{
    // Declarations from R -------------------------------------------
    double eps = eps_[0];
    double delta = delta_[0];
    double lambda_min = lambda_min_[0];
    int nlam = nlam_[0];
    int n = n_[0];
    int p = p_[0];
    int ppflag = ppflag_[0];
    int scrflag = scrflag_[0];
    int max_iter = max_iter_[0];
    int dfmax = dfmax_[0];
    int user = user_[0];
    int trace = trace_[0];

    // Declarations from C -------------------------------------------
    int i, j, l, ll, jn, lp;
    int nnzero = 0; int violations = 0; int nv = 0;
    double N = (double)n;
    double grad_sum, screen_factor, null_dev, tol, lstep, lambda_diff, cut_off, lambda_fixed, 
        max_update, change, upd, sum_sq, first_sum, second_sum, abs_diff, max_diff = 0.0;

    int* include = R_Calloc(p, int);
    double* c_old = R_Calloc(p, double);
    double* c_new = R_Calloc(p, double);
    double* v = R_Calloc(n, double); // v = y -x beta
    double* uj = R_Calloc(n, double); // u = v + x_j beta_j
    double* wx = R_Calloc(n * p, double); // wx
    double* wxj = R_Calloc(n, double); // wx_j
    double* w2x2 = R_Calloc(n * p, double); // w^2x^2
    double* w2x2j = R_Calloc(n, double);// w^2x^2_j
    double* threshold = R_Calloc(n * p, double); // delta/|wx|
    double* threshold_j = R_Calloc(n, double); // delta/|wx_j|

    double* rj = R_Calloc(n, double);
    
    double* beta_old = R_Calloc(p, double);

    // Preprocessing -------------------------------------------
    // This needs to be added.
    if (ppflag == 0) {
        simple_process_adaptive(v, wx, w2x2, threshold, y, x, w, n, p, delta);
    }

    // Set up initial solutions -------------------------------------------
    for (j = 0; j < p; j++) {
        include[j] = 0;
        grad_sum = 0.0;
        jn = j * n;
        for (i = 0; i < n; i++) {
            grad_sum += huber_grad(y[i], delta) * wx[i + jn];
        }
        c_old[j] = fabs(grad_sum) / N;
        c_new[j] = c_old[j];
    }
    screen_factor = 1.0;

    // Compute null deviance and tolerance
    null_dev = 0.0;
    for (i = 0; i < n; i++) {
        null_dev += w[i] * huber_loss(y[i], delta);
    }
    tol = eps * null_dev;
    if (trace) {
        Rprintf("Tolerance = %f\n", tol);
    }

    // Set up lambda -------------------------------------------
    if (user == 0) {
        double max_lambda = 0.0;
        for (j = 0; j < p; j++) {
            grad_sum = 0.0;
            jn = j * n;
            for (i = 0; i < n; i++) {
                grad_sum += x[i + jn] * huber_grad(y[i], delta);
            }
            double abs_val = fabs(grad_sum);
            if (abs_val > max_lambda) {
                max_lambda = abs_val;
            }
        }
        lambda[0] = max_lambda / N;

        if (lambda_min == 0.0) {
            lambda_min = 0.001;
        }
        lstep = log(lambda_min) / ((double)nlam - 1.0);
        for (l = 1; l < nlam; l++) {
            lambda[l] = lambda[l - 1] * exp(lstep);
        }
    }


    // Solution path -------------------------------------------
    for (l = 0; l < nlam; l++) {

        lp = l * p;
        lambda_fixed = lambda[l];
        if (trace) {
            Rprintf("Lambda %d\n", l + 1);
        }

        // Construct the set of predictors that need to be updated -------------------------------------------
        if (screen_factor > 3.0) {
            screen_factor = 3.0;
        }

        if (l == 0) {
            for (j = 0; j < p; j++) {
                beta_old[j] = 0.0;
            }
            lambda_diff = 1.0;
            cut_off = lambda_fixed;
        }
        else {
            for (j = 0; j < p; j++) {
                beta_old[j] = beta[(l - 1) * p + j];
            }
            lambda_diff = lambda[(l - 1)] - lambda_fixed;
            cut_off = (1.0 + screen_factor) * lambda_fixed - screen_factor * lambda_diff;
        }

        // Construct the set of predictors that need to be updated -------------------------------------------
        for (int j = 0; j < p; j++) {
            if (include[j] == 0 && c_old[j] > cut_off) {
                include[j] = 1;
            }
        }
        if (scrflag == 1) {
            screen_factor = 0.0;
        }


        // For each fixed lambda, solve the problem -------------------------------------------
        ll = 0;
        while (ll < max_iter) {
            // Check dfmax
            if (nnzero > dfmax) {
                for (ll = l; ll < nlam; ll++) {
                    iter[ll] = NA_INTEGER;
                }
                saturated[0] = 1;
                break;
            }

            // Update the predictors in the set -------------------------------------------
            while (ll < max_iter) {
                ll += 1;
                max_update = 0.0;
                for (j = 0; j < p; j++) {
                    if (include[j]) {
                        // Get partial residuals -------------------------------------------
                        jn = j * n;
                        for (i = 0;i < n;i++) {
                            uj[i] = v[i] + x[i + jn] * beta_old[j];
                            rj[i] = w[i] * uj[i] / x[i + jn];
                            threshold_j[i] = threshold[i + jn];
                            wxj[i] = wx[i + jn];
                            w2x2j[i] = w2x2[i + jn];
                        }

                        // Get partial residuals -------------------------------------------
                        double out;
                        huber_fit(&out, rj, wxj, w2x2j, threshold_j, lambda_fixed, delta, n, N);
                        beta[lp + j] = out;

                        // Update residuals -------------------------------------------
                        change = beta[lp + j] - beta_old[j];
                        for (i = 0;i < n;i++) {
                            v[i] += -x[i + jn] * change;
                        }
                        beta_old[j] = beta[lp + j];

                        // Update max_update -------------------------------------------
                        sum_sq = 0.0;
                        for (i = 0; i < n; i++) {
                            sum_sq += w2x2j[i];
                        }
                        upd = (sum_sq / N) * (change * change);
                        if (upd > max_update) {
                            max_update = upd;
                        }
                    }
                }
                // Convergence check -------------------------------------------
                if (max_update < tol) {
                    break;
                }
            }// Update ends

            // Check KKT conditions of the predictors not in the set, update E, & count nonzero variables -------------------------------------------
            violations = 0; nnzero = 0;
            for (j = 0;j < p;j++) {
                if (!include[j]) {
                    // Get partial residuals -------------------------------------------
                    jn = j * n;
                    for (i = 0;i < n;i++) {
                        uj[i] = v[i] + x[i + jn] * beta[lp + j];
                        wxj[i] = wx[i + jn];
                    }

                    first_sum = 0.0; second_sum = 0.0;
                    for (i = 0; i < n; i++) {
                        if (fabs(w[i]*uj[i]) > delta) {
                            second_sum += fabs(wxj[i]) * ((wxj[i] > 0) - (wxj[i] < 0)) * ((w[i] * uj[i] > 0) - (w[i] * uj[i] < 0));
                        }
                        else {
                            first_sum += w[i] * uj[i] * wxj[i];
                        }
                    }
                    c_new[j] = -first_sum / N + -delta * second_sum / N;
                    if (fabs(c_new[j]) > lambda_fixed) {
                        include[j] = 1;
                        violations += 1;
                    }
                    c_old[j] = c_new[j];
                }
                if (beta_old[j] != 0) {
                    nnzero++;
                }
            }// KKT check ends
            if (scrflag == 1) {
                for (j = 0; j < p; j++) {
                    abs_diff = fabs(c_old[j] - c_new[j]);
                    if (abs_diff > max_diff) {
                        max_diff = abs_diff;
                    }
                }
                screen_factor = max_diff / lambda_diff;
            }
            if (violations == 0) {
                break;
            }
            nv += violations;
        }// iteration ends
        iter[l] = ll;
        if (trace) {
            Rprintf("# iterations = %d\n", ll);
        }
    }// fixed lambda cycle ends
    numv[0] = nv;

    // Free allocated memories -------------------------------------------
    R_Free(include);
    R_Free(c_old);
    R_Free(c_new);
    R_Free(v);
    R_Free(uj);
    R_Free(wx);
    R_Free(wxj);
    R_Free(w2x2);
    R_Free(w2x2j);
    R_Free(threshold);
    R_Free(threshold_j);
    R_Free(rj);
    R_Free(beta_old);
}

// Exact Coordinate Descent (ECD) for penalized Huber regression -------------------------------------------
void ecd_huber_active_(double* beta, int* iter, double* lambda, int* saturated, int* numv, double* x, double* y, double* delta_, double* eps_, double* lambda_min_,
    int* nlam_, int* n_, int* p_, int* ppflag_, int* dfmax_, int* max_iter_, int* user_, int* scrflag_, int* trace_)
{
    // Declarations from R -------------------------------------------
    double eps = eps_[0];
    double delta = delta_[0];
    double lambda_min = lambda_min_[0];
    int nlam = nlam_[0];
    int n = n_[0];
    int p = p_[0];
    int ppflag = ppflag_[0];
    int scrflag = scrflag_[0];
    int max_iter = max_iter_[0];
    int dfmax = dfmax_[0];
    int user = user_[0];
    int trace = trace_[0];

    // Declarations from C -------------------------------------------
    int i, j, l, ll, jn, lp;
    int nnzero = 0; int violations = 0; int nv = 0;
    double N = (double)n;
    double grad_sum, screen_factor, null_dev, tol, lstep, lambda_diff, cut_off, lambda_fixed,
        max_update, change, upd, sum_sq, first_sum, second_sum, abs_diff, max_diff = 0.0;

    int* include = R_Calloc(p, int);
    double* c_old = R_Calloc(p, double);
    double* c_new = R_Calloc(p, double);
    double* v = R_Calloc(n, double); // v = y -x beta
    double* uj = R_Calloc(n, double); // u = v + x_j beta_j
    double* xj = R_Calloc(n, double); // x_j
    double* x2 = R_Calloc(n * p, double); // x^2
    double* x2j = R_Calloc(n, double);// x^2_j
    double* threshold = R_Calloc(n * p, double); // delta/|x|
    double* threshold_j = R_Calloc(n, double); // delta/|x_j|

    double* rj = R_Calloc(n, double);

    double* beta_old = R_Calloc(p, double);

    // Preprocessing -------------------------------------------
    // This needs to be added.
    if (ppflag == 0) {
        simple_process(v, x2, threshold, y, x, n, p, delta);
    }

    // Set up initial solutions -------------------------------------------
    for (j = 0; j < p; j++) {
        include[j] = 0;
        grad_sum = 0.0;
        jn = j * n;
        for (i = 0; i < n; i++) {
            grad_sum += huber_grad(y[i], delta) * x[i + jn];
        }
        c_old[j] = fabs(grad_sum) / N;
        c_new[j] = c_old[j];
    }
    screen_factor = 1.0;

    // Compute null deviance and tolerance
    null_dev = 0.0;
    for (i = 0; i < n; i++) {
        null_dev += huber_loss(y[i], delta);
    }
    tol = eps * null_dev;
    if (trace) {
        Rprintf("Tolerance = %f\n", tol);
    }


    // Set up lambda -------------------------------------------
    if (user == 0) {
        double max_lambda = 0.0;
        for (j = 0; j < p; j++) {
            grad_sum = 0.0;
            jn = j * n;
            for (i = 0; i < n; i++) {
                grad_sum += x[i + jn] * huber_grad(y[i], delta);
            }
            double abs_val = fabs(grad_sum);
            if (abs_val > max_lambda) {
                max_lambda = abs_val;
            }
        }
        lambda[0] = max_lambda / N;

        if (lambda_min == 0.0) {
            lambda_min = 0.001;
        }
        lstep = log(lambda_min) / ((double)nlam - 1.0);
        for (l = 1; l < nlam; l++) {
            lambda[l] = lambda[l - 1] * exp(lstep);
        }
    }


    // Solution path -------------------------------------------
    for (l = 0; l < nlam; l++) {

        lp = l * p;
        lambda_fixed = lambda[l];
        if (trace) {
            Rprintf("Lambda %d\n", l + 1);
        }

        // Construct the set of predictors that need to be updated -------------------------------------------
        if (screen_factor > 3.0) {
            screen_factor = 3.0;
        }

        if (l == 0) {
            for (j = 0; j < p; j++) {
                beta_old[j] = 0.0;
            }
            lambda_diff = 1.0;
            cut_off = lambda_fixed;
        }
        else {
            for (j = 0; j < p; j++) {
                beta_old[j] = beta[(l - 1) * p + j];
            }
            lambda_diff = lambda[(l - 1)] - lambda_fixed;
            cut_off = (1.0 + screen_factor) * lambda_fixed - screen_factor * lambda_diff;
        }

        // Construct the set of predictors that need to be updated -------------------------------------------
        for (int j = 0; j < p; j++) {
            if (include[j] == 0 && c_old[j] > cut_off) {
                include[j] = 1;
            }
        }

        if (scrflag == 1) {
            screen_factor = 0.0;
        }


        // For each fixed lambda, solve the problem -------------------------------------------
        ll = 0;
        while (ll < max_iter) {
            // Check dfmax
            if (nnzero > dfmax) {
                for (ll = l; ll < nlam; ll++) {
                    iter[ll] = NA_INTEGER;
                }
                saturated[0] = 1;
                break;
            }

            // Update the predictors in the set -------------------------------------------
            while (ll < max_iter) {
                ll += 1;
                max_update = 0.0;
                for (j = 0; j < p; j++) {
                    if (include[j]) {
                        upd = 0.0;

                        // Get partial residuals -------------------------------------------
                        jn = j * n;
                        for (i = 0;i < n;i++) {
                            uj[i] = v[i] + x[i + jn] * beta_old[j];
                            rj[i] = uj[i] / x[i + jn];
                            threshold_j[i] = threshold[i + jn];
                            xj[i] = x[i + jn];
                            x2j[i] = x2[i + jn];
                        }

                        // Get partial residuals -------------------------------------------
                        double out;
                        huber_fit(&out, rj, xj, x2j, threshold_j, lambda_fixed, delta, n, N);
                        beta[lp + j] = out;

                        // Update residuals -------------------------------------------
                        change = beta[lp + j] - beta_old[j];
                        for (i = 0;i < n;i++) {
                            v[i] += -x[i + jn] * change;
                        }

                        // Update max_update -------------------------------------------
                        sum_sq = 0.0;
                        for (i = 0; i < n; i++) {
                            sum_sq += x2j[i];
                        }
                        upd = (sum_sq / N) * (change * change);
                        if (upd > max_update) {
                            max_update = upd;
                        }
                        beta_old[j] = beta[lp + j];
                    }
                }
                // Convergence check -------------------------------------------
                if (max_update < tol) {
                    break;
                }
            }// Update ends

            // Check KKT conditions of the predictors not in the set & update E -------------------------------------------
            violations = 0; nnzero = 0;
            for (j = 0;j < p;j++) {
                if (!include[j]) {

                    // Get partial residuals -------------------------------------------
                    jn = j * n;
                    for (i = 0;i < n;i++) {
                        uj[i] = v[i] + x[i + jn] * beta[lp + j];
                        xj[i] = x[i + jn];
                    }

                    first_sum = 0.0; second_sum = 0.0;
                    for (i = 0; i < n; i++) {
                        if (fabs(uj[i]) > delta) {
                            second_sum += fabs(xj[i]) * ((xj[i] > 0) - (xj[i] < 0)) * ((uj[i] > 0) - (uj[i] < 0));
                        }
                        else {
                            first_sum += uj[i] * xj[i];
                        }
                    }
                    c_new[j] = -first_sum / N + -delta * second_sum / N;
                    if (fabs(c_new[j]) > lambda_fixed) {
                        include[j] = 1;
                        violations += 1;
                    }
                    c_old[j] = c_new[j];
                }
            }// KKT check ends
            if (scrflag == 1) {
                for (j = 0; j < p; j++) {
                    abs_diff = fabs(c_old[j] - c_new[j]);
                    if (abs_diff > max_diff) {
                        max_diff = abs_diff;
                    }
                }
                screen_factor = max_diff / lambda_diff;
            }
            if (violations == 0) {
                break;
            }
            nv += violations;
        }// iteration ends
        iter[l] = ll;
        if (trace) {
            Rprintf("# iterations = %d\n", ll);
        }
    }// fixed lambda cycle ends
    numv[0] = nv;

    // Free allocated memories -------------------------------------------
    R_Free(include);
    R_Free(c_old);
    R_Free(c_new);
    R_Free(v);
    R_Free(uj);
    R_Free(xj);
    R_Free(x2);
    R_Free(x2j);
    R_Free(threshold);
    R_Free(threshold_j);
    R_Free(rj);
    R_Free(beta_old);
}

// Exact Coordinate Descent (ECD) for adaptive penalized Huber regression without screening (for test purpose) -------------------------------------------
void ecd_huber_noscreen_(double* beta, int* iter, double* lambda, double* x, double* y, double* w, double* delta_, double* eps_, double* lambda_min_, 
    int* nlam_, int* n_, int* p_, int* ppflag_, int* max_iter_, int* user_, int* trace_)
{
    // Declarations from R -------------------------------------------
    double eps = eps_[0];
    double delta = delta_[0];
    double lambda_min = lambda_min_[0];
    int nlam = nlam_[0];
    int n = n_[0];
    int p = p_[0];
    int ppflag = ppflag_[0];
    int max_iter = max_iter_[0];
    int user = user_[0];
    int trace = trace_[0];

    // Declarations from C -------------------------------------------
    int i, j, l, ll, jn, lp;
    double N = (double)n;
    double tol, null_dev, grad_sum, lstep, lambda_fixed, KKT_j, change_sum, loss_sum, first_sum, second_sum;

    double* change = R_Calloc(p, double); // beta - beta_old
    double* loss = R_Calloc(max_iter, double);
    double* v = R_Calloc(n, double); // v = y -x beta
    double* uj = R_Calloc(n, double); // u = v + x_j beta_j
    double* wx = R_Calloc(n * p, double); // wx
    double* wxj = R_Calloc(n, double); // wx_j
    double* w2x2 = R_Calloc(n * p, double); // w^2x^2
    double* w2x2j = R_Calloc(n, double);// w^2x^2_j
    double* threshold = R_Calloc(n * p, double); // delta/|wx|
    double* threshold_j = R_Calloc(n, double); // delta/|wx_j|

    double* rj = R_Calloc(n, double);

    double* beta_old = R_Calloc(p, double);

    // Preprocessing -------------------------------------------
    // This needs to be added.
    if (ppflag == 0) {
        simple_process_adaptive(v, wx, w2x2, threshold, y, x, w, n, p, delta);
    }

    // Compute null deviance
    null_dev = 0.0;
    for (i = 0; i < n; i++) {
        null_dev += huber_loss(y[i], delta);
    }
    tol = eps * null_dev;

    // Set up lambda -------------------------------------------
    if (user == 0) {
        double max_lambda = 0.0;
        for (j = 0; j < p; j++) {
            grad_sum = 0.0;
            jn = j * n;
            for (i = 0; i < n; i++) {
                grad_sum += x[i + jn] * huber_grad(y[i], delta);
            }
            double abs_val = fabs(grad_sum);
            if (abs_val > max_lambda) {
                max_lambda = abs_val;
            }
        }
        lambda[0] = max_lambda / N;

        if (lambda_min == 0.0) {
            lambda_min = 0.001;
        }
        lstep = log(lambda_min) / ((double)nlam - 1.0);
        for (l = 1; l < nlam; l++) {
            lambda[l] = lambda[l - 1] * exp(lstep);
        }
    }


    // Solution path -------------------------------------------
    for (l = 0; l < nlam; l++) {

        lp = l * p;
        lambda_fixed = lambda[l];
        if (trace) {
            Rprintf("Lambda %d\n", l + 1);
        }

        // For each fixed lambda, solve the problem -------------------------------------------
        for (ll = 0; ll < max_iter; ll++) {
            for (j = 0; j < p; j++) {
                jn = j * n;
                for (i = 0;i < n;i++) {
                    uj[i] = v[i] + x[i + jn] * beta_old[j];
                    wxj[i] = wx[i + jn];
                }

                // Check KKT conditions -------------------------------------------
                first_sum = 0.0; second_sum = 0.0;
                for (i = 0; i < n; i++) {
                    if (fabs(w[i] * uj[i]) > delta) {
                        second_sum += fabs(wxj[i]) * ((wxj[i] > 0) - (wxj[i] < 0)) * ((w[i] * uj[i] > 0) - (w[i] * uj[i] < 0));
                    }
                    else {
                        first_sum += w[i] * uj[i] * wxj[i];
                    }
                }
                KKT_j = fabs(-first_sum / N + -delta * second_sum / N);

                if (KKT_j <= lambda_fixed) {
                    beta[lp + j] = 0;
                    change[j] = 0;
                }
                else {
                    // Get partial residuals & update -------------------------------------------
                    for (i = 0;i < n;i++) {
                        rj[i] = w[i] * uj[i] / x[i + jn];
                        threshold_j[i] = threshold[i + jn];
                        w2x2j[i] = w2x2[i + jn];
                    }
                    double out;
                    huber_fit(&out, rj, wxj, w2x2j, threshold_j, lambda_fixed, delta, n, N);
                    beta[lp + j] = out;

                    // Update residuals -------------------------------------------
                    change[j] = beta[lp + j] - beta_old[j];
                    for (i = 0;i < n;i++) {
                        v[i] += -x[i + jn] * change[j];
                    }
                    beta_old[j] = beta[lp + j];
                }
            }
            // Convergence check -------------------------------------------
            change_sum = 0.0;
            for (j = 0;j < p;j++) {
                change_sum += change[j] * change[j];
            }
            if (change_sum < eps || ll == max_iter) {
                iter[l] = ll;
                break;
            }
            else {
                loss_sum = 0.0;
                for (i = 0; i < n; i++) {
                    loss_sum += w[i] * huber_loss(v[i], delta);
                }
                loss[ll] = loss_sum / N;
                if (ll >= 4) {
                    if (fabs(loss[ll] - loss[(ll - 2)]) < tol || fabs(loss[ll] - loss[(ll - 1)]) < tol) {
                        iter[l] = ll;
                        break;
                    }
                }
            }
        }// iteration ends
        if (trace) {
            Rprintf("# iterations = %d\n", ll);
        }
    }// fixed lambda cycle ends

    // Free allocated memories -------------------------------------------
    R_Free(loss);
    R_Free(v);
    R_Free(uj);
    R_Free(wx);
    R_Free(wxj);
    R_Free(w2x2);
    R_Free(w2x2j);
    R_Free(threshold);
    R_Free(threshold_j);
    R_Free(rj);
    R_Free(beta_old);
}

static const R_CMethodDef cMethods[] = {
    // name        pointer         Num args
    {"ecd_huber_adaptive_active_", (DL_FUNC)&ecd_huber_adaptive_active_, 20},
    {"ecd_huber_active_", (DL_FUNC)&ecd_huber_active_, 19},
    {"ecd_huber_noscreen_", (DL_FUNC)&ecd_huber_noscreen_, 16},
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
    R_forceSymbols(info, FALSE);
}