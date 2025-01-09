#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const double EPS = 10.22;   //[K]
const double SIGMA = 2.556; //[Å]
const double A0 = 5.;       //[Å]
const double BETA2 = 5.;    // adimensional, 5 is the correct value to oppose the LJ divergence
const double DELTA = 5.;    // nice value
const double HBAR2_2M = 6.0596;

// copies array a into b
void copy_array(double *a, double *b, int N){
    for (int i = 0; i < N; i++){
        b[i] = a[i];
    }
    return;
}

// scalar product
double scalar_product(double *a, double *b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// u function
double u(double r, double beta) {
    return pow(beta / r, BETA2);
}

// first derivative of u function
double u_prime(double r, double beta) {
    return -BETA2 * u(r, beta) / r;
}

// second derivative of u function
double u_doubleprime(double r, double beta) {
    return BETA2 * (1. + BETA2) * u(r, beta) / (r * r);
}

double psi(double* r, double alpha, double beta) {
    double r1[3] = {r[0], r[1], r[2]};
    double r1_sq = scalar_product(r1, r1);
    double r2[3] = {r[3], r[4], r[5]};
    double r2_sq = scalar_product(r2, r2);
    double r12[3] = {r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2]};
    double r12_norm = sqrt(scalar_product(r12, r12));
    return exp(-0.5 * (1 / alpha) * (r1_sq + r2_sq) - 0.5 * u(r12_norm, beta));
}

double local_energy(double* r, double alpha, double beta) {
    double r1[3] = {r[0], r[1], r[2]};
    double r1_sq = scalar_product(r1, r1);
    double r2[3] = {r[3], r[4], r[5]};
    double r2_sq = scalar_product(r2, r2);
    double r12[3] = {r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2]};
    double r12_norm = sqrt(scalar_product(r12, r12));

    double kinetic = HBAR2_2M*(6/alpha + u_doubleprime(r12_norm, beta) + 2*u_prime(r12_norm, beta)/r12_norm - (r1_sq+r2_sq)/pow(alpha, 2) - u_prime(r12_norm, beta)*r12_norm/alpha - 0.5*pow(u_prime(r12_norm, beta), 2));
    double momega2 = HBAR2_2M*2/pow(A0, 4);
    double harmonic = 0.5 * momega2 * (r1_sq + r2_sq);
    double lj = 4 * EPS * (pow(SIGMA/r12_norm, 12) - pow(SIGMA/r12_norm, 6));

    return kinetic + harmonic + lj;
}

int main(int argc, char **argv) {
    srand(2);

    // parameters
    int N = 2;
    int n_steps = 100000;
    double alpha = 20.0;
    double alpha_step = 1.0;
    double beta = 1.9;
    double beta_step = 0.05;
    int max_var = 15;

    // arrays
    double alpha_values[max_var];
    double beta_values[max_var];
    double energies[max_var][max_var];

    // file
    FILE *f_energies = fopen("data/test.csv", "w");
    fprintf(f_energies, "alpha,beta,energies\n");

    // print on terminal simulation parameters
    printf("==================================================================================\n");    
    printf("Running simulation with N = %d, n_steps = %d\n\n", N, n_steps);

    // start MC
    double* r_old = malloc(3 * N * sizeof(double)); 
    double* r_new = malloc(3 * N * sizeof(double));

    for (int ia = 0; ia < max_var; ia++) {
        alpha += alpha_step;
        alpha_values[ia] = alpha;
        beta = 1.9;
        for (int jb = 0; jb < max_var; jb++) {
            // printf("ia = %d, jb = %d\n", ia, jb);
            beta += beta_step;
            beta_values[jb] = beta;
            printf("\rVariational parameters: alpha = %.1f, beta1 = %.2f", alpha, beta);
            double energy = 0.0;
            double delta_energy = 0.0;

            // initial positions and WF
            for (int i = 0; i < 3*N; i++) {
                r_old[i] = A0 * 2 * (rand() / (1.0 + RAND_MAX) - 0.5);
            }
            double wfold = psi(r_old, alpha, beta);

            // MC cycle
            for (int m = 0; m < n_steps; m++) {
                // trial position
                for (int i = 0; i < 3*N; i++) {
                    r_new[i] = r_old[i] + A0 * (rand() / (1.0 + RAND_MAX) - 0.5);
                }
                double wfnew = psi(r_new, alpha, beta);

                // Metropolis to accept
                double a = rand() / (1.0 + RAND_MAX);
                if (a < (wfnew * wfnew) / (wfold * wfold)) {
                    copy_array(r_new, r_old, 3*N);
                    wfold = wfnew;
                    delta_energy = local_energy(r_old, alpha, beta);
                }
                energy += delta_energy;

            }
            energy /= n_steps;
            // printf("E = %.3f\n", energy);
            energies[ia][jb] = energy;
            fprintf(f_energies, "%.4f,%.4f,%.10f\n", alpha, beta, energy);
        }
    }

    printf("\n\nSimulation completed 🎉🎊\n");
    printf("==================================================================================\n");

    free(r_old);
    free(r_new);
    return 0;
}