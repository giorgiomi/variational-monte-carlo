#include "functions.h"
#include "metropolis.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EPS 10.22   //[K]
#define SIGMA 2.556 //[Å]
#define A0 5.       //[Å]
#define BETA2 4.    // adimensional, 4 first try

double potential_energy(double *r, double *var_param, int N) {
    return potential_bruteforce(r, var_param, N);
}

int main(int argc, char **argv) {
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
    double alpha = 1.; // 1 first try
    double beta1 = 0.; // 0 to remove interaction
    double var_param[3] = {alpha, beta1, BETA2};

    // print on terminal simulation parameters
    printf("------------------------------------------------------------------------\n");
    printf("Running simulation with N = %d, n_steps = %d, delta = %.2f\n\n", N, n_steps, delta);
    printf("Variational parameters: alpha = %.2f, beta1 = %.2f, beta2 = %.2f\n\n", alpha, beta1, BETA2);

    // positions
    double *r = malloc(3 * N * sizeof(double));

    // files
    FILE *f_energy = fopen("data/energy.csv", "w");
    FILE *f_kinetic_avg = fopen("data/kinetic_avg.csv", "w");
    FILE *f_acceptance = fopen("data/acceptance.csv", "w");
    fprintf(f_energy, "i,T,V,E\n");
    fprintf(f_kinetic_avg, "i,T_lap,T_grad\n");
    fprintf(f_acceptance, "i,a\n");


    // initial configuration (positions)
    for (int i = 0; i < 3 * N; i++) {
        double csi = 2. * (rand() / (1.0 + RAND_MAX)) - 1.;
        r[i] = delta * csi;
    }

    // initial observables
    double T = kinetic_energy(r, var_param, N);
    double V = potential_energy(r, var_param, N);
    double E = T + V;
    fprintf(f_energy, "0,%.10e,%.10e,%.10e\n", T, V, E);

    // initial observables, kinetic averages
    double T_lap = kinetic_average_laplacian(r, var_param, N);
    double T_grad = kinetic_average_gradient(r, var_param, N);
    fprintf(f_kinetic_avg, "0,%.10e,%.10e\n", T_lap, T_grad);

    // initial acceptance
    double acc_rate = 0.;
    fprintf(f_acceptance, "0,%.10e\n", acc_rate);
     
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
        //print_array(r, 3 * N);

        //double E_proposed = energy(r, var_param, N);
        
        // accepting the proposed step
        double a = acceptance(r_old, r, part_index, var_param, N);
        double a_rand = rand() / (1.0 + RAND_MAX);
        if (a < a_rand) {
            copy_array(r_old, r, 3 * N);
            acc_rate += 1.;
        }
        //printf("Step %d\n", i);

        // calculate observables
        T = kinetic_energy(r, var_param, N);
        V = potential_energy(r, var_param, N);
        E = T + V;
        fprintf(f_energy, "%d,%.10e,%.10e,%.10e\n", i, T, V, E);

        // calculate kinetic averages
        T_lap = kinetic_average_laplacian(r, var_param, N);
        T_grad = kinetic_average_gradient(r, var_param, N);
        fprintf(f_kinetic_avg, "%d,%.10e,%.10e\n", i, T_lap, T_grad);

        // print acceptance rate
        fprintf(f_acceptance, "%d,%.10e\n", i, acc_rate / i);
    }

    printf("\n\nSimulation completed\n");
    printf("------------------------------------------------------------------------\n");
    free(r);
    free(r_old);

    return 0;
}