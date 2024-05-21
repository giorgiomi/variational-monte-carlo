#include <math.h>

// Scalar product
double scalar_product(double *a, double *b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// U function
double u(double r, double *beta) {
    double beta1 = beta[0];
    double beta2 = beta[1];
    return pow(beta1 / r, beta2);
}