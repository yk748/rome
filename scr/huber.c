#include <R.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <R_ext/Applic.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"

// Miscellaneous
double yhat(double* x, double* v, int p, int i) {
    int start_index = i * p;
    double sum = 0.0;
    for (int j = 0; j < p; j++) {
        sum += x[start_index + j] * v[j];
    }
    return sum;
}

// Huber loss gradient :
double huber_loss_grad(double v, double delta) {
    double val = 0.0;
    if (fabs(v) <= delta)
        val = v;
    else
        val = delta * ((v > 0) - (v < 0));
    return val;
}

// Huber loss:
double huber_loss(double v, double delta) {
    double val = 0.0;
    if (fabs(v) <= delta)
        val = v * v / 2;
    else
        val = delta * fabs(v) - delta * delta / 2;
    return val;
}

// Comparison function for qsort
typedef struct {
    double kinks;
    double slopes;
    double cumulative_sum;
    double value;
} Pair;

int comparePairs(const void* a, const void* b) {
    const Pair* pairA = (const Pair*)a;
    const Pair* pairB = (const Pair*)b;
    if (pairA->kinks < pairB->kinks) return -1;
    if (pairA->kinks > pairB->kinks) return 1;
    return 0;
}

// partial gradient of penalized Huber regression:
static void pg_Huber_reg(double *y, double *x, double *w, double *beta, int n, int p, double delta) {

    double val = 0.0; double tmp = 0.0; double res = 0.0;
    for (int i = 0; i < n; i++) {
        res = w[i] * (y[i] - yhat(x, beta, p, i));
        tmp = tmp + huber_loss_grad(res, delta);
    }
    val = tmp / n;
}

// Exact coordinate descent for Huber, univariate:
static void cd_Huber(double *x, double *w, double *r, int n, double delta, double lambda) {

    // set up the table
    double wght_x[n];
    double incr[n];
    double thres[n];
    double tmp = 0.0;
    for (int i = 0; i < n;i++) {
        wght_x[i] = fabs(w[i] * x[i]);
        incr[i] = wght_x[i] * wght_x[i] / n;
        thres[i] = delta / wght_x[i];
        tmp = tmp + wght_x[i];
    }
    double S = -delta * tmp / n;

    // Allocate memory for pairs
    Pair* pairs = malloc((2*n) * sizeof(Pair));
    for (int i = 0; i < n; i++) {
        pairs[i].kinks = r[i] - thres[i];
        pairs[i].slopes = incr[i];
    }

    for (int i = n; i < (2*n); i++) {
        pairs[i].kinks = r[i - n + 1] - thres[i - n + 1];
        pairs[i].slopes = -incr[i - n + 1];
    }
    qsort(pairs, 2*n, sizeof(Pair), comparePairs);

    // Calculate cumulative sums
    double cusum = 0;
    for (int i = 0; i < n; i++) {
        cusum += pairs[i].slopes;
        pairs[i].cumulative_sum = cusum;
    }

    // Initialization
    if (pairs[0].kinks < 0) {
        pairs[0].value = S - lambda;
    }
    else {
        pairs[0].value = S + lambda;
    }

    double output = 0.0;
    for (int i = 0; i < (2 * n-1);i++) {
        S = S + pairs[i].cumulative_sum * (pairs[i + 1].kinks - pairs[i].kinks);
        if (pairs[0].kinks < 0) {
            pairs[i + 1].value = S - lambda;
        }
        else {
            pairs[i + 1].value = S + lambda;
        }

        if (i > 1) {
            if (pairs[i - 1].value < 0 && pairs[i].value >= 0) {
                output = pairs[i - 1].kinks - pairs[i - 1].value * (pairs[i].kinks - pairs[i - 1].kinks) / (pairs[i].value - pairs[i - 1].value);
                break;
            }
        }
    }


}



// Define the name and address and number of arguments
static const R_CMethodDef cMethods[] = {
    // name          pointer         Num args
  {"pg_Huber_reg", (DL_FUNC) &pg_Huber_reg, 7},
  {"cd_Huber", (DL_FUNC) &cd_Huber, 6},
  {NULL,     NULL, 0}   // Placeholder to indicate last one.
};

//  Initialize the DLL
void R_init_rome(DllInfo *dll) {
    R_registerRoutines(
        dll,     // DllInfo
        cMethods, // .C
        NULL,     // .Call
        NULL,     // Fortran
        NULL      // External
    );
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}