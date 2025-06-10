#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <R_ext/Applic.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"


// Huber loss -------------------------------------------
double huber_loss(double v, double thresh) {
    if (fabs(v) <= thresh) {
        return 0.5 * v * v;
    }
    else {
        return thresh * fabs(v) - 0.5 * thresh * thresh;
    }
}

// Huber gradient -------------------------------------------
double huber_grad(double v, double thresh) {
    if (fabs(v) <= thresh) {
        return v;
    }
    else {
        return thresh * ((v > 0) - (v < 0)); // sign(v)
    }
}

// Simple processing with assignment of nonconst -------------------------------------------
void simple_process(double* v, double* x2, double* threshold, double* y, double* x, int n, int p, double delta) {
    int i, j, jn;
    for (i = 0; i < n; i++) {
        v[i] = y[i];
        for (j = 0; j < p; j++) {
            jn = j * n;
            threshold[jn + i] = delta / fabs(x[jn + i]);
            x2[jn + i] = x[jn + i] * x[jn + i];
        }
    }
}

// Adaptive simple processing with assignment of nonconst -------------------------------------------
void simple_process_adaptive(double* v, double* wx, double* w2x2, double* threshold, double* y, double* x, double* w, int n, int p, double delta) {
    int i, j, jn;
    for (i = 0; i < n; i++) {
        v[i] = y[i];
        for (j = 0; j < p; j++) {
            jn = j * n;
            wx[jn + i] = w[i] * x[jn + i];
            threshold[jn + i] = delta / fabs(wx[jn + i]);
            w2x2[jn + i] = wx[jn + i] * wx[jn + i];
        }
    }
}


