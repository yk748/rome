#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <R_ext/Applic.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"

// Data frame for sorting -------------------------------------------
typedef struct {
    double kinks;
    double slopes;
    double cumulative_sum;
    double value;
} Pair;

// Quick sort -------------------------------------------
int comparePairs(const void* a, const void* b) {
    double diff = ((Pair*)a)->kinks - ((Pair*)b)->kinks;
    return (diff > 0) - (diff < 0);  // Returns 1, -1, or 0
}

// Penalized Huber loss minimization for a fixed variable & fixed lambda -------------------------------------------
void huber_fit(double* out, double* rj, double* first_j, double* second_j, double* threshold_j, double lambda_fixed, double delta, int n, double N) {

    int i;
    // Set up the table -------------------------------------------
    Pair* pairs = malloc((2 * n) * sizeof(Pair));
    for (i = 0; i < n; i++) {
        pairs[i].kinks = rj[i] - threshold_j[i];
        pairs[i].slopes = second_j[i] / N;
    }

    for (i = n; i < (2 * n); i++) {
        pairs[i].kinks = rj[i - n] + threshold_j[i - n];
        pairs[i].slopes = -second_j[i - n] / N;
    }

    qsort(pairs, 2 * n, sizeof(Pair), comparePairs);
    double cusum = 0.0;
    for (i = 0; i < (2 * n); i++) {
        cusum += pairs[i].slopes;
        pairs[i].cumulative_sum = cusum;
    }

    // Initialization -------------------------------------------
    double S = 0.0;
    for (i = 0; i < n;i++) {
        S += -delta * fabs(first_j[i]) / N;
    }

    if (pairs[0].kinks >= 0) {
        pairs[0].value = S + lambda_fixed;
    }
    else {
        pairs[0].value = S - lambda_fixed;
    }

    // Main loop -------------------------------------------
    double sol = 0.0;
    for (i = 0; i < (2 * n - 1); i++) {
        S = S + pairs[i].cumulative_sum * (pairs[i + 1].kinks - pairs[i].kinks);
        if (pairs[i + 1].kinks > 0) {
            pairs[i + 1].value = S + lambda_fixed;
        }
        else if (pairs[i + 1].kinks < 0) {
            pairs[i + 1].value = S - lambda_fixed;
        }
        else {
            pairs[i + 1].value = S;
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
