#include "functions.h"
#include "metropolis.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define EPS 10.22   //[K]
#define SIGMA 2.556 //[Å]
#define A0 5.       //[Å]
#define BETA2 5.    // adimensional, 5 is the correct value to oppose the LJ divergence

double potential_energy(double *r, double *var_param, int N) {
    return potential_bruteforce(r, var_param, N);
}

int main(int argc, char **argv) {
    srand(1);
    // parameters
    if (argc < 3) {
        printf("Usage: %s N n_steps [alpha_saved]\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    int n_steps = (int)strtod(argv[2], NULL);
    double alpha_saved = (argc == 4) ? atof(argv[3]) : -1.;
    double delta = 5.; // nice value for acceptance

    // variational parameters, to vary
    double alpha_start = A0 * A0;
    double alpha_end = A0 * A0;
    double alpha_step = 1.;
    double beta1 = 2.5; // 0 to remove interaction

    // Create directory name
    char dir_name[50];
    snprintf(dir_name, sizeof(dir_name), "data/INT_%d_%d", N, n_steps);

    // Create the directory
    struct stat st = {0};
    if (stat(dir_name, &st) == -1) {
        int status = mkdir(dir_name, 0777);
        if (status != 0) {
            perror("Error creating directory");
            return 1;
        }
    } else {
        // Directory exists, remove its contents
        char command[100];
        snprintf(command, sizeof(command), "rm -rf %s/*", dir_name);
        int status = system(command);
        if (status != 0) {
            perror("Error clearing directory");
            return 1;
        }
    }

    // files
    char variational_path[100];
    snprintf(variational_path, sizeof(variational_path), "%s/variational.csv", dir_name);
    FILE *f_variational = fopen(variational_path, "w");
    fprintf(f_variational, "alpha,T,T_std,T_lap,T_lap_std,T_grad,T_grad_std,V,V_std,E,E_std\n");

    FILE *f_energy = NULL;
    FILE *f_kinetic_avg = NULL;
    FILE *f_acceptance = NULL;
    FILE *f_psi = NULL;
    FILE *f_pos = NULL;

    if (argc >= 3 && alpha_saved >= 0.) {
        char energy_path[100];
        snprintf(energy_path, sizeof(energy_path), "%s/energy_%.1f.csv", dir_name, alpha_saved);
        f_energy = fopen(energy_path, "w");

        char kinetic_avg_path[100];
        snprintf(kinetic_avg_path, sizeof(kinetic_avg_path), "%s/kinetic_avg_%.1f.csv", dir_name, alpha_saved);
        f_kinetic_avg = fopen(kinetic_avg_path, "w");

        char acceptance_path[100];
        snprintf(acceptance_path, sizeof(acceptance_path), "%s/acceptance_%.1f.csv", dir_name, alpha_saved);
        f_acceptance = fopen(acceptance_path, "w");

        char psi_path[100];
        snprintf(psi_path, sizeof(psi_path), "%s/psi_%.1f.csv", dir_name, alpha_saved);
        f_psi = fopen(psi_path, "w");

        char pos_path[100];
        snprintf(pos_path, sizeof(pos_path), "%s/pos_%.1f.csv", dir_name, alpha_saved);
        f_pos = fopen(pos_path, "w");

        fprintf(f_energy, "i,T,V,E\n");
        fprintf(f_kinetic_avg, "i,T_lap,T_grad\n");
        fprintf(f_acceptance, "i,a\n");
        fprintf(f_psi, "i,psi\n");
        
        fprintf(f_pos, "i");
        for (int i = 0; i < N; i++) {
            char str[100];
            snprintf(str, sizeof(str), ",x%d,y%d,z%d", i, i, i);
            fprintf(f_pos, "%s", str);
        }
        fprintf(f_pos, "\n");
    }

    // print on terminal simulation parameters
    printf("==================================================================================\n");    
    printf("Running simulation with N = %d, n_steps = %d, delta = %.2f, alpha_saved = %.1f\n\n", N, n_steps, delta, alpha_saved);

    // starting simul with variational alpha
    for (double alpha = alpha_start; alpha <= alpha_end; alpha += alpha_step) {
        double var_param[3] = {alpha, beta1, BETA2};
        printf("\rVariational parameters: alpha = %.1f, beta1 = %.2f, beta2 = %.2f", alpha, beta1, BETA2);
        fflush(stdout);

        // positions
        double *r = malloc(3 * N * sizeof(double));

        // initial configuration (positions)
        for (int i = 0; i < 3 * N; i++) {
            double csi = 2. * (rand() / (1.0 + RAND_MAX)) - 1.;
            r[i] = A0 * csi;
        }

        // initial observables
        double T = kinetic_energy(r, var_param, N);
        double V = potential_energy(r, var_param, N);
        double E = T + V;

        // observables avg and std
        double T_avg = 0.;
        double T2_avg = 0.;
        double V_avg = 0.;
        double V2_avg = 0.;
        double E_avg = 0.;
        double E2_avg = 0.;

        // initial observables, kinetic estimators
        double T_lap = kinetic_estimator_laplacian(r, var_param, N);
        double T_grad = kinetic_estimator_gradient(r, var_param, N);

        // kinetic estimatos avg and std
        double T_lap_avg = 0.;
        double T_lap2_avg = 0.;
        double T_grad_avg = 0.;
        double T_grad2_avg = 0.;

        // initial acceptance
        double rej_rate = 0.;

        // initial wavefunction
        double psi_curr = psi(r, var_param, N);

        if (argc >= 3 && alpha == alpha_saved) {
            fprintf(f_energy, "0,%.10e,%.10e,%.10e\n", T, V, E);
            fprintf(f_kinetic_avg, "0,%.10e,%.10e\n", T_lap, T_grad);
            fprintf(f_acceptance, "0,%.10e\n", 1. - rej_rate);
            fprintf(f_psi, "0,%.10e\n", psi_curr);

            fprintf(f_pos, "0");
            for (int i = 0; i < N; i++) {
                fprintf(f_pos, ",%.10e,%.10e,%.10e", r[3 * i], r[3 * i + 1], r[3 * i + 2]);
            }
            fprintf(f_pos, "\n");
        }
        
        // MC simul with fixed alpha
        double *r_old = malloc(3 * N * sizeof(double));
        for (int i = 1; i <= n_steps; i++) {
            // update configuration with M(RT)^2
            int part_index = i % N;
            copy_array(r, r_old, 3 * N);
            
            // update positions with T function (uniform)
            for (int j = 0; j < 3; j++) {
                double csi = 2. * (rand() / (1.0 + RAND_MAX)) - 1.;
                double x_test = csi * delta;
                r[3 * part_index + j] += x_test;
            }
            
            // accepting the proposed step
            double a = acceptance(r_old, r, var_param, N);
            double a_rand = rand() / (1.0 + RAND_MAX);
            if (a < a_rand) {
                copy_array(r_old, r, 3 * N);
                rej_rate += 1.;
            } 
            
            // calculate observables
            T = kinetic_energy(r, var_param, N);
            V = potential_energy(r, var_param, N);
            E = T + V;

            // observables avg and std
            T_avg += T;
            V_avg += V;
            E_avg += E;
            T2_avg += T * T;
            V2_avg += V * V;
            E2_avg += E * E;

            // calculate kinetic estimators
            T_lap = kinetic_estimator_laplacian(r, var_param, N);
            T_grad = kinetic_estimator_gradient(r, var_param, N);

            // kinetic estimators avg and std
            T_lap_avg += T_lap;
            T_lap2_avg += T_lap * T_lap;
            T_grad_avg += T_grad;
            T_grad2_avg += T_grad * T_grad;

            // print on file
            if (argc >= 3 && alpha == alpha_saved) {
                fprintf(f_energy, "%d,%.10e,%.10e,%.10e\n", i, T, V, E);
                fprintf(f_kinetic_avg, "%d,%.10e,%.10e\n", i, T_lap, T_grad);
                fprintf(f_acceptance, "%d,%.10e\n", i, 1. - rej_rate / i);

                double psi_curr = psi(r, var_param, N);
                fprintf(f_psi, "%d,%.10e\n", i, psi_curr);

                fprintf(f_pos, "%d", i);
                for (int j = 0; j < N; j++) {
                    fprintf(f_pos, ",%.10e,%.10e,%.10e", r[3 * j], r[3 * j + 1], r[3 * j + 2]);
                }
                fprintf(f_pos, "\n");
            }
        }
        
        // print on file alpha and observables
        T_avg /= n_steps;
        V_avg /= n_steps;
        E_avg /= n_steps;
        T_lap_avg /= n_steps;
        T_grad_avg /= n_steps;
        double T_std = sqrt(T2_avg / n_steps - T_avg * T_avg);
        double V_std = sqrt(V2_avg / n_steps - V_avg * V_avg);
        double E_std = sqrt(E2_avg / n_steps - E_avg * E_avg);
        double T_lap_std = sqrt(T_lap2_avg / n_steps - T_lap_avg * T_lap_avg);
        double T_grad_std = sqrt(T_grad2_avg / n_steps - T_grad_avg * T_grad_avg);
        fprintf(f_variational, "%.2f,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n", alpha, T_avg, T_std/sqrt(n_steps), T_lap_avg, T_lap_std/sqrt(n_steps), T_grad_avg, T_grad_std/sqrt(n_steps), V_avg, V_std/sqrt(n_steps), E_avg, E_std/sqrt(n_steps));

        free(r);
        free(r_old);

    }

    printf("\n\nSimulation completed 🎉🎊\n");
    printf("==================================================================================\n");
    
    return 0;
}
