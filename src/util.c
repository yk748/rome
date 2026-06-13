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

// Postprocessing -------------------------------------------
void postprocess(double* beta, double* shift, double* scale, int* nonconst, int nlam, int p) {
    int l, j, lp;
    for (l = 0; l < nlam; l++) {
        lp = l * p;
        for (j = 0; j < p; j++) {
            if (nonconst[j]) {
                beta[lp + j] /= scale[j];
            }
        }
    }
}

// Standardize: center by mean, scale by sd -------------------------------------------
void standardize(double* v, double* x2, double* threshold,
    double* shift, double* scale, int* nonconst,
    double* y, double* x, int n, int p, double delta)
{
    int i, j, jn;
    double xm, xvar, xsd;

    for (i = 0; i < n; i++) {
        v[i] = y[i];
    }

    for (j = 0; j < p; j++) {
        jn = j * n;
        nonconst[j] = 0;
        shift[j] = 0.0;
        scale[j] = 1.0;

        // Mean
        xm = 0.0;
        for (i = 0; i < n; i++) {
            xm += x[jn + i];
        }
        xm /= n;

        // Center and variance
        xvar = 0.0;
        for (i = 0; i < n; i++) {
            x[jn + i] -= xm;
            xvar += x[jn + i] * x[jn + i];
        }
        xvar /= n;
        xsd = sqrt(xvar);

        if (xsd > 1e-6) {
            nonconst[j] = 1;
            shift[j] = xm;
            scale[j] = xsd;
            for (i = 0; i < n; i++) {
                x[jn + i] /= xsd;
                x2[jn + i] = x[jn + i] * x[jn + i];
                threshold[jn + i] = delta / fabs(x[jn + i]);
            }
        }
        else {
            // Constant column: zero out
            for (i = 0; i < n; i++) {
                x[jn + i] = 0.0;
                x2[jn + i] = 0.0;
                threshold[jn + i] = R_PosInf;
            }
        }
    }
}

// Standardize (adaptive) -------------------------------------------
void standardize_adaptive(double* v, double* wx, double* w2x2, double* threshold,
    double* shift, double* scale, int* nonconst,
    double* y, double* x, double* w, int n, int p, double delta)
{
    int i, j, jn;
    double xm, xvar, xsd, wxi;

    for (i = 0; i < n; i++) {
        v[i] = y[i];
    }

    for (j = 0; j < p; j++) {
        jn = j * n;
        nonconst[j] = 0;
        shift[j] = 0.0;
        scale[j] = 1.0;

        xm = 0.0;
        for (i = 0; i < n; i++) {
            xm += x[jn + i];
        }
        xm /= n;

        xvar = 0.0;
        for (i = 0; i < n; i++) {
            x[jn + i] -= xm;
            xvar += x[jn + i] * x[jn + i];
        }
        xvar /= n;
        xsd = sqrt(xvar);

        if (xsd > 1e-6) {
            nonconst[j] = 1;
            shift[j] = xm;
            scale[j] = xsd;
            for (i = 0; i < n; i++) {
                x[jn + i] /= xsd;
                wxi = w[i] * x[jn + i];
                wx[jn + i] = wxi;
                w2x2[jn + i] = wxi * wxi;
                threshold[jn + i] = delta / fabs(wxi);
            }
        }
        else {
            for (i = 0; i < n; i++) {
                x[jn + i] = 0.0;
                wx[jn + i] = 0.0;
                w2x2[jn + i] = 0.0;
                threshold[jn + i] = R_PosInf;
            }
        }
    }
}

// Rescale: shift by min, scale by range -------------------------------------------
void rescale(double* v, double* x2, double* threshold,
    double* shift, double* scale, int* nonconst,
    double* y, double* x, int n, int p, double delta)
{
    int i, j, jn;
    double cmin, cmax, crange;

    for (i = 0; i < n; i++) {
        v[i] = y[i];
    }

    for (j = 0; j < p; j++) {
        jn = j * n;
        nonconst[j] = 0;
        shift[j] = 0.0;
        scale[j] = 1.0;

        cmin = x[jn]; cmax = x[jn];
        for (i = 1; i < n; i++) {
            if (x[jn + i] < cmin) {
                cmin = x[jn + i];
            }
            else if (x[jn + i] > cmax) {
                cmax = x[jn + i];
            }
        }
        crange = cmax - cmin;

        if (crange > 1e-6) {
            nonconst[j] = 1;
            shift[j] = cmin;
            scale[j] = crange;
            for (i = 0; i < n; i++) {
                x[jn + i] = (x[jn + i] - cmin) / crange;
                x2[jn + i] = x[jn + i] * x[jn + i];
                threshold[jn + i] = delta / fabs(x[jn + i]);
            }
        }
        else {
            for (i = 0; i < n; i++) {
                x[jn + i] = 0.0;
                x2[jn + i] = 0.0;
                threshold[jn + i] = R_PosInf;
            }
        }
    }
}

// Rescale (adaptive) -------------------------------------------
void rescale_adaptive(double* v, double* wx, double* w2x2, double* threshold,
    double* shift, double* scale, int* nonconst,
    double* y, double* x, double* w, int n, int p, double delta)
{
    int i, j, jn;
    double cmin, cmax, crange, wxi;

    for (i = 0; i < n; i++) {
        v[i] = y[i];
    }

    for (j = 0; j < p; j++) {
        jn = j * n;
        nonconst[j] = 0;
        shift[j] = 0.0;
        scale[j] = 1.0;

        cmin = x[jn]; cmax = x[jn];
        for (i = 1; i < n; i++) {
            if (x[jn + i] < cmin) {
                cmin = x[jn + i];
            }
            else if (x[jn + i] > cmax) {
                cmax = x[jn + i];
            }
        }
        crange = cmax - cmin;

        if (crange > 1e-6) {
            nonconst[j] = 1;
            shift[j] = cmin;
            scale[j] = crange;
            for (i = 0; i < n; i++) {
                x[jn + i] = (x[jn + i] - cmin) / crange;
                wxi = w[i] * x[jn + i];
                wx[jn + i] = wxi;
                w2x2[jn + i] = wxi * wxi;
                threshold[jn + i] = delta / fabs(wxi);
            }
        }
        else {
            for (i = 0; i < n; i++) {
                x[jn + i] = 0.0;
                wx[jn + i] = 0.0;
                w2x2[jn + i] = 0.0;
                threshold[jn + i] = R_PosInf;
            }
        }
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


