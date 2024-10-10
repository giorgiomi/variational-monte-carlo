#include "functions.h"
#include <math.h>

// acceptance
double acceptance(double *r, double *r_proposed, double *var_param, int N) {
    double alpha = var_param[0];
    double beta[2] = {var_param[1], var_param[2]};

    double P_proposed = psi(r_proposed, var_param, N) * psi(r_proposed, var_param, N);
    double P = psi(r, var_param, N) * psi(r, var_param, N);

    double factor = P_proposed / P;
    if (factor > 1.) {
        return 1.;
    } else {
        return factor;
    }
}