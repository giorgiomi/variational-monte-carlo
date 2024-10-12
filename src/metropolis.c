#include "functions.h"
#include <math.h>
#include <stdio.h>

// acceptance
double acceptance(double *r, double *r_proposed, double *var_param, int N) {
    double alpha = var_param[0];
    double beta[2] = {var_param[1], var_param[2]};

    double psi_curr = psi(r, var_param, N);
    double psi_prop = psi(r_proposed, var_param, N);

    // printf("psi_curr: %f, psi_prop: %f\n", psi_curr, psi_prop);

    double P_proposed = psi_prop * psi_prop;
    double P = psi_curr * psi_curr;

    double factor = P_proposed / P;
    if (factor > 1.) {
        return 1.;
    } else {
        return factor;
    }
}