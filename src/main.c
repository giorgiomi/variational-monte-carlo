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
    srand(2);
    // parameters
    if (argc < 9) {
        printf("Usage: %s N n_steps alpha_start alpha_end alpha_step beta1_start beta1_end beta1_step [alpha_saved]\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    int n_steps = (int)strtod(argv[2], NULL);
    double alpha_start = atof(argv[3]);
    double alpha_end = atof(argv[4]);
    double alpha_step = atof(argv[5]);
    double beta1_start = atof(argv[6]);
    double beta1_end = atof(argv[7]);
    double beta1_step = atof(argv[8]);
    double alpha_saved = (argc == 10) ? atof(argv[9]) : -1.;
    double delta = 5.; // nice value for acceptance

    if (alpha_saved != -1. && (alpha_saved < alpha_start || alpha_saved > alpha_end)) {
        fprintf(stderr, "Error: alpha_saved (%.2f) must be between alpha_start (%.2f) and alpha_end (%.2f)\n", alpha_saved, alpha_start, alpha_end);
        return 1;
    }

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
    fprintf(f_variational, "alpha,beta1,T,T_std,T_lap,T_lap_std,T_grad,T_grad_std,VHO,VHO_std,VLJ,VLJ_std,E,E_std\n");

    FILE *f_energy = NULL;
    FILE *f_kinetic_avg = NULL;
    FILE *f_acceptance = NULL;
    FILE *f_psi = NULL;
    FILE *f_pos = NULL;
    FILE *f_check = NULL;

    if (argc >= 10 && alpha_saved >= 0.) {
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

        char check_path[100];
        snprintf(check_path, sizeof(check_path), "%s/check_%.1f.csv", dir_name, alpha_saved);
        f_check = fopen(check_path, "w");

        fprintf(f_energy, "i,T,VHO,VLJ,E\n");
        fprintf(f_kinetic_avg, "i,T_lap,T_grad\n");
        fprintf(f_acceptance, "i,a,a_rand,rate\n");
        fprintf(f_psi, "i,psi\n");
        fprintf(f_check, "i,a,T_old,VHO_old,VLJ_old,T_new,VHO_new,VLJ_new\n");
        
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

    // starting simul with variational alpha and beta1
    for (double alpha = alpha_start; alpha <= alpha_end; alpha += alpha_step) {
        for (double beta1 = beta1_start; beta1 <= beta1_end; beta1 += beta1_step) {
            double var_param[3] = {alpha, beta1, BETA2};
            printf("\rVariational parameters: alpha = %.1f, beta1 = %.2f, beta2 = %.2f", alpha, beta1, BETA2);
            fflush(stdout);

            // positions
            double *r_old = malloc(3 * N * sizeof(double));

            // initial configuration (positions)
            for (int i = 0; i < 3 * N; i++) {
                double csi = 2. * (rand() / (1.0 + RAND_MAX)) - 1.;
                r_old[i] = A0 * csi;
            }

            // observables avg and std
            double T_avg = 0.;
            double T2_avg = 0.;
            double VLJ_avg = 0.;
            double VLJ2_avg = 0.;
            double VHO_avg = 0.;
            double VHO2_avg = 0.;
            double E_avg = 0.;
            double E2_avg = 0.;
            double T_lap_avg = 0.;
            double T_lap2_avg = 0.;
            double T_grad_avg = 0.;
            double T_grad2_avg = 0.;
            double rej_rate = 0.;
            
            // observables
            double *r_new = malloc(3 * N * sizeof(double));
            double T_new, VLJ_new, VHO_new, E_new, T_lap_new, T_grad_new;

            // calculate initial observables
            double T_old = kinetic_energy(r_old, var_param, N);
            double VHO_old = HO_potential(r_old, N);
            double VLJ_old = LJ_potential(r_old, N);
            double E_old = T_old + VHO_old + VLJ_old;
            double T_lap_old = kinetic_estimator_laplacian(r_old, var_param, N);
            double T_grad_old = kinetic_estimator_gradient(r_old, var_param, N);

            // MC simul with fixed alpha and beta1
            for (int i = 0; i < n_steps; i++) {
                // update configuration with M(RT)^2
                int part_index = i % N;
                
                // update positions with T function (uniform)
                for (int j = 0; j < 3; j++) {
                    double csi = 2. * (rand() / (1.0 + RAND_MAX)) - 1.;
                    double x_test = csi * delta;
                    r_new[3 * part_index + j] += x_test;
                }
                
                // accepting the proposed step
                double a = acceptance(r_old, r_new, var_param, N);
                double a_rand = rand() / (1.0 + RAND_MAX);
                if (a < a_rand) {
                    // if rejected, restore the old configuration
                    copy_array(r_old, r_new, 3 * N);
                    rej_rate += 1.;
                    
                    // keep the old observables
                    T_new = T_old;
                    VHO_new = VHO_old;
                    VLJ_new = VLJ_old;
                    E_new = E_old;
                    T_lap_new = T_lap_old;
                    T_grad_new = T_grad_old;
                } else {
                    // if accepted, update the observables
                    T_new = kinetic_energy(r_new, var_param, N);
                    VHO_new = HO_potential(r_new, N);
                    VLJ_new = LJ_potential(r_new, N);
                    E_new = T_new + VHO_new + VLJ_new;
                    T_lap_new = kinetic_estimator_laplacian(r_new, var_param, N);
                    T_grad_new = kinetic_estimator_gradient(r_new, var_param, N);  
                }        

                // observables avg and std
                T_avg += T_new;
                VHO_avg += VHO_new;
                VLJ_avg += VLJ_new;
                E_avg += E_new;
                T2_avg += T_new * T_new;
                VHO2_avg += VHO_new * VHO_new;
                VLJ2_avg += VLJ_new * VLJ_new;
                E2_avg += E_new * E_new;

                // kinetic estimators avg and std
                T_lap_avg += T_lap_new;
                T_lap2_avg += T_lap_new * T_lap_new;
                T_grad_avg += T_grad_new;
                T_grad2_avg += T_grad_new * T_grad_new;

                // print on file
                if (argc >= 10 && alpha == alpha_saved) {
                    fprintf(f_energy, "%d,%.10e,%.10e,%.10e,%.10e\n", i, T_new, VHO_new, VLJ_new, E_new);
                    fprintf(f_kinetic_avg, "%d,%.10e,%.10e\n", i, T_lap_new, T_grad_new);
                    fprintf(f_acceptance, "%d,%.10e,%.10e,%.10e\n", i, a, a_rand, 1. - rej_rate / i);
                    fprintf(f_check, "%d,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n", i, a, T_old, VHO_old, VLJ_old, T_new, VHO_new, VLJ_new);

                    double psi_curr = psi(r_new, var_param, N);
                    fprintf(f_psi, "%d,%.10e\n", i, psi_curr);

                    fprintf(f_pos, "%d", i);
                    for (int j = 0; j < N; j++) {
                        fprintf(f_pos, ",%.10e,%.10e,%.10e", r_new[3 * j], r_new[3 * j + 1], r_new[3 * j + 2]);
                    }
                    fprintf(f_pos, "\n");
                }

                // update old observables
                T_old = T_new;
                VHO_old = VHO_new;
                VLJ_old = VLJ_new;
                E_old = E_new;
                T_lap_old = T_lap_new;
                T_grad_old = T_grad_new;

                // update old configuration
                copy_array(r_new, r_old, 3 * N);
            }
            
            // print on file alpha, beta1 and observables
            T_avg /= n_steps;
            VHO_avg /= n_steps;
            VLJ_avg /= n_steps;
            E_avg /= n_steps;
            T_lap_avg /= n_steps;
            T_grad_avg /= n_steps;
            double T_std = sqrt(T2_avg / n_steps - T_avg * T_avg);
            double VHO_std = sqrt(VHO2_avg / n_steps - VHO_avg * VHO_avg);
            double VLJ_std = sqrt(VLJ2_avg / n_steps - VLJ_avg * VLJ_avg);
            double E_std = sqrt(E2_avg / n_steps - E_avg * E_avg);
            double T_lap_std = sqrt(T_lap2_avg / n_steps - T_lap_avg * T_lap_avg);
            double T_grad_std = sqrt(T_grad2_avg / n_steps - T_grad_avg * T_grad_avg);
            fprintf(f_variational, "%.2f,%.2f,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n", alpha, beta1, T_avg, T_std/sqrt(n_steps), T_lap_avg, T_lap_std/sqrt(n_steps), T_grad_avg, T_grad_std/sqrt(n_steps), VHO_avg, VHO_std/sqrt(n_steps), VLJ_avg, VLJ_std/sqrt(n_steps), E_avg, E_std/sqrt(n_steps));

            free(r_new);
            free(r_old);
        }
    }

    printf("\n\nSimulation completed 🎉🎊\n");
    printf("==================================================================================\n");
    
    return 0;
}
