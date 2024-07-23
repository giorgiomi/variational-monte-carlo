#include "functions.h"
#include "metropolis.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EPS 10.22   //[K]
#define SIGMA 2.556 //[Å]
#define A0 5.       //[Å]

double energy(double *r, double *var_param, int N) {
    return kinetic_energy(r, var_param, N);
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
    int N = atoi(argv[1]); //2
    int n_steps = atoi(argv[2]); //200
    double delta = atof(argv[3]); //0.1

    // variational parameters, to vary
    double alpha = 1.;
    double beta1 = 1.;
    double var_param[2] = {alpha, beta1};

    // positions
    double *r = malloc(3 * N * sizeof(double));

    // file
    FILE *f_energy = fopen("data/energy.csv", "w");
    fprintf(f_energy, "i,E\n");

    // initial configuration

    double *old_r = malloc(3 * N * sizeof(double));
    for (int i = 0; i < n_steps; i++) {
        // calculate observables
        double E = energy(r, var_param, N);
        fprintf(f_energy, "%d,%.10e\n", i, E);
        
        // update configuration with M(RT)^2
        int part_index = i % N;
        copy_array(r, old_r, 3 * N);
        //print_array(old_r, 3 * N);
        
        // update positions with T function (uniform)
        for (int j = 0; j < 3; j++) {
            double csi = 2. * (rand() / (1.0 + RAND_MAX)) - 1.;
            double x_test = csi * delta;
            r[3 * part_index + j] += x_test;
        }
        //print_array(r, 3 * N);

        double E_proposed = energy(r, var_param, N);
        
        // accepting the proposed step
        double a = acceptance(old_r, r, part_index, var_param, N);
        double a_rand = rand() / (1.0 + RAND_MAX);
        if (a < a_rand) {
            copy_array(old_r, r, 3 * N);
        }
        //printf("Step %d\n", i);
    }

    return 0;
}