#include "functions.h"
#include "metropolis.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
    switch (argc) {
        case 1:
            printf("Usage: %s N n_steps delta\n", argv[0]);
            return 1;
        case 2:
            printf("Please provide the number of steps and the delta\n");
            return 1;
        case 3:
            printf("Please provide the delta\n");
            return 1;
        case 4:
            break;
        default:
            return 1;
    }
    int N = atoi(argv[1]);
    int n_steps = atoi(argv[2]);
    double delta = atof(argv[3]);

    // variational parameters, to vary
    // double alpha = A0*A0; // 1 first try
    double beta1 = 0.; // 0 to remove interaction

    // files
    FILE *f_param = fopen("data/param.csv", "w");
    FILE *f_variational = fopen("data/variational.csv", "w");
    fprintf(f_variational, "alpha,T,T_std,T_lap,T_lap_std,T_grad,T_grad_std,V,V_std,E,E_std\n");
    fprintf(f_param, "N,n_steps,delta\n");
    fprintf(f_param, "%d,%d,%.2f\n", N, n_steps, delta);

    // starting simul with variational alpha
    double alpha_start = A0 * A0 - 15.;
    double alpha_end = A0 * A0 + 25.;
    for (double alpha = alpha_start; alpha <= alpha_end; alpha += 1.) {
        double var_param[3] = {alpha, beta1, BETA2};

        // print on terminal simulation parameters
        printf("------------------------------------------------------------------------\n");
        printf("Running simulation with N = %d, n_steps = %d, delta = %.2f\n\n", N, n_steps, delta);
        printf("Variational parameters: alpha = %.2f, beta1 = %.2f, beta2 = %.2f\n\n", alpha, beta1, BETA2);

        // positions
        double *r = malloc(3 * N * sizeof(double));

        // files
        // FILE *f_energy = fopen("data/energy.csv", "w");
        // FILE *f_kinetic_avg = fopen("data/kinetic_avg.csv", "w");
        // FILE *f_acceptance = fopen("data/acceptance.csv", "w");
        // fprintf(f_energy, "i,T,V,E\n");
        // fprintf(f_kinetic_avg, "i,T_lap,T_grad\n");
        // fprintf(f_acceptance, "i,a\n");

        // initial configuration (positions)
        for (int i = 0; i < 3 * N; i++) {
            double csi = 2. * (rand() / (1.0 + RAND_MAX)) - 1.;
            r[i] = A0 * csi;
        }

        // initial observables
        double T = kinetic_energy(r, var_param, N);
        double V = potential_energy(r, var_param, N);
        double E = T + V;
        // fprintf(f_energy, "0,%.10e,%.10e,%.10e\n", T, V, E);

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
        // fprintf(f_kinetic_avg, "0,%.10e,%.10e\n", T_lap, T_grad);

        // kinetic estimatos avg and std
        double T_lap_avg = 0.;
        double T_lap2_avg = 0.;
        double T_grad_avg = 0.;
        double T_grad2_avg = 0.;

        // initial acceptance
        double rej_rate = 0.;
        // fprintf(f_acceptance, "0,%.10e\n", 1. - rej_rate);
        
        // MC simulation
        double *r_old = malloc(3 * N * sizeof(double));
        for (int i = 1; i <= n_steps; i++) {
            // print iteration on terminal
            printf("\rIteration: %d/%d", i, n_steps);
            fflush(stdout);

            // update configuration with M(RT)^2
            int part_index = i % N;
            copy_array(r, r_old, 3 * N);
            //print_array(r_old, 3 * N);
            
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
                // printf(" Rej\n");
                rej_rate += 1.;
            } else {
                // printf(" Acc\n");
            }
            // print_array(r, 3 * N);

            // calculate observables
            T = kinetic_energy(r, var_param, N);
            V = potential_energy(r, var_param, N);
            E = T + V;
            // fprintf(f_energy, "%d,%.10e,%.10e,%.10e\n", i, T, V, E);

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
            // fprintf(f_kinetic_avg, "%d,%.10e,%.10e\n", i, T_lap, T_grad);

            // kinetic estimators avg and std
            T_lap_avg += T_lap;
            T_lap2_avg += T_lap * T_lap;
            T_grad_avg += T_grad;
            T_grad2_avg += T_grad * T_grad;

            // print acceptance rate
            // fprintf(f_acceptance, "%d,%.10e\n", i, 1. - rej_rate / i);
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

        printf("\n\nSimulation completed\n");
        printf("------------------------------------------------------------------------\n");
        free(r);
        free(r_old);

    }
    

    
    return 0;
}