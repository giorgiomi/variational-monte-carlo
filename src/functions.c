#include <math.h>

// scalar product
double scalar_product(double *a, double *b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// u function
double u(double r, double *beta) {
    double beta1 = beta[0];
    double beta2 = beta[1];
    return pow(beta1 / r, beta2);
}

// first derivative of u function
double u_prime(double r, double *beta) {
    double beta1 = beta[0];
    double beta2 = beta[1];
    return -beta2 * pow(beta1 / r, beta2) / r;
}

// second derivative of u function
double u_doubleprime(double r, double *beta) {
    double beta1 = beta[0];
    double beta2 = beta[1];
    return beta2 * (1. + beta2) * pow(beta1 / r, beta2) / (r * r);
}

// LJ
double lennard_jones(double r, double *param) {
    double eps = param[0];
    double sigma = param[1];
    return 4 * eps * (pow(sigma / r, 12) - pow(sigma / r, 6));
}

// kinetic energy
double kinetic_energy(double *r, double *param, int N, double h2_2m) {
    double alpha = param[0];
    double beta[2] = {param[1], param[2]};
    double res = 0.;

    // cycle through k to calculate each kinetic contribution
    for (int k = 0; k < N; k++) {
        double rk[3] = {r[3 * k], r[3 * k + 1], r[3 * k + 2]};
        double rk_mod2 = scalar_product(rk, rk);
        res += 3. / alpha;
        res += -rk_mod2 / (alpha * alpha);

        // cycle through j != k
        for (int j = 0; j < N; j++) {
            if (j != k) {
                double rkj[3] = {rk[0] - r[3 * j], rk[1] - r[3 * j + 1], rk[2] - r[3 * j + 2]};
                double rkj_mod = sqrt(scalar_product(rkj, rkj));
                res += 0.5 * u_doubleprime(rkj_mod, beta);
                res += u_prime(rkj_mod, beta) / rkj_mod;
                res += -u_prime(rkj_mod, beta) * scalar_product(rk, rkj) / (rkj_mod * alpha);

                // cycle through l != k
                for (int l = 0; l < N; l++) {
                    if (l != k) {
                        double rkl[3] = {rk[0] - r[3 * l], rk[1] - r[3 * l + 1], rk[2] - r[3 * l + 2]};
                        double rkl_mod = sqrt(scalar_product(rkl, rkl));
                        res += u_prime(rkj_mod, beta) * u_prime(rkl_mod, beta) * scalar_product(rkj, rkl) / (rkj_mod * rkl_mod);
                    }
                }
            }
        }
    }

    res *= h2_2m;
    return res;
}