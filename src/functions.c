#include <math.h>
#include <stdio.h>

#define H2_2M 6.0596 // h bar squared over 2 time mass, [Å^2 K]
#define EPS 10.22   //[K]
#define SIGMA 2.556 //[Å]
#define A0 5.       //[Å]

// copy array
void copy_array(double *a, double *b, int N){
    for (int i = 0; i < N; i++){
        b[i] = a[i];
    }
    return;
}

// print array
void print_array(double *a, int N){
    for (int i = 0; i < N; i++){
        printf("%f ", a[i]);
    }
    printf("\n");
    return;
}

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
double lennard_jones(double r) {
    return 4 * EPS * (pow(SIGMA / r, 12) - pow(SIGMA / r, 6));
}

// psi
// double psi(double *r, double *param, int N) {
//     double alpha = param[0];
//     double beta[2] = {param[1], param[2]};
// }

// kinetic energy
double kinetic_energy(double *r, double *param, int N) {
    double alpha = param[0];
    double beta[2] = {param[1], param[2]};
    double res = 0.;

    // cycle through i to calculate each kinetic contribution
    for (int i = 0; i < N; i++) {
        double ri[3] = {r[3 * i], r[3 * i + 1], r[3 * i + 2]};
        double ri_mod2 = scalar_product(ri, ri);
        res += 3. / alpha;
        res -= ri_mod2 / (alpha * alpha);

        // cycle through j != i
        for (int j = 0; j < N; j++) {
            if (j != i) {
                double rij[3] = {ri[0] - r[3 * j], ri[1] - r[3 * j + 1], ri[2] - r[3 * j + 2]};
                double rij_mod = sqrt(scalar_product(rij, rij));
                res += 0.5 * u_doubleprime(rij_mod, beta);
                res += u_prime(rij_mod, beta) / rij_mod;
                res -= u_prime(rij_mod, beta) * scalar_product(ri, rij) / (rij_mod * alpha);

                // cycle through l != i
                for (int l = 0; l < N; l++) {
                    if (l != i) {
                        double ril[3] = {ri[0] - r[3 * l], ri[1] - r[3 * l + 1], ri[2] - r[3 * l + 2]};
                        double ril_mod = sqrt(scalar_product(ril, ril));
                        res -= 0.25 * u_prime(rij_mod, beta) * u_prime(ril_mod, beta) * scalar_product(rij, ril) / (rij_mod * ril_mod);
                    }
                }
            }
        }
    }

    res *= H2_2M;
    return res;
}

// harmonic potential after 1 move
double harmonic_potential(double *r_old, double *r_new, double VH_old, int part_index, int N) {
    double m_omega2 = H2_2M * 2. / pow(A0, 4);
    double r2_old = scalar_product(r_old + 3 * part_index, r_old + 3 * part_index);
    double r2_new = scalar_product(r_new + 3 * part_index, r_new + 3 * part_index);
    
    return VH_old - 0.5 * m_omega2 * (r2_old - r2_new);
}

// LJ potential single contribution
double LJ_potential_single(double *r, int part_index, int N) {
    double VLJ_new = 0.;
    for (int i = 0; i < N; i++) {
        if (i != part_index) {
            double rij[3] = {r[3 * i] - r[3 * part_index], r[3 * i + 1] - r[3 * part_index + 1], r[3 * i + 2] - r[3 * part_index + 2]};
            double rij_mod = sqrt(scalar_product(rij, rij));
            VLJ_new += lennard_jones(rij_mod);
        }
    }
    return VLJ_new;
}

// LJ potential after 1 move
double LJ_potential(double *r_old, double *r_new, double VLJ_old, int part_index, int N) {
    return VLJ_old - LJ_potential_single(r_old, part_index, N) + LJ_potential_single(r_new, part_index, N);
}

// total potential (brute force)
double potential_bruteforce(double *r, double *param, int N) {
    double alpha = param[0];
    double beta[2] = {param[1], param[2]};
    double res = 0.;
    double m_omega2 = H2_2M * 2. / pow(A0, 4);

    // LJ
    double VLJ = 0.;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double rij[3] = {r[3 * i] - r[3 * j], r[3 * i + 1] - r[3 * j + 1], r[3 * i + 2] - r[3 * j + 2]};
            double rij_mod = sqrt(scalar_product(rij, rij));
            VLJ += lennard_jones(rij_mod);
        }
    }

    // harmonic
    double VH = 0.;
    for (int i = 0; i < N; i++) {
        double ri[3] = {r[3 * i], r[3 * i + 1], r[3 * i + 2]};
        double ri_mod2 = scalar_product(ri, ri);
        VH += 0.5 * m_omega2 * ri_mod2;
    }

    return VLJ + VH;
}