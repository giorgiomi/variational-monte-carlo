#include "functions.h"
#include <math.h>

// acceptance
double acceptance(double *r, double *r_proposed, int part_index, double *var_param, int N) {
    double alpha = var_param[0];
    double beta[2] = {var_param[1], var_param[2]};

    double part_r[3] = {r[3 * part_index], r[3 * part_index + 1], r[3 * part_index + 2]};
    double part_r_proposed[3] = {r_proposed[3 * part_index], r_proposed[3 * part_index + 1], r_proposed[3 * part_index + 2]};
    double kinetic = -(scalar_product(part_r_proposed, part_r_proposed) - scalar_product(part_r, part_r))/alpha;
    double potential = 0.;
    for (int i = 0; i < N; i++) {
        if (i != part_index) {
            double r_i[3] = {r[3 * i], r[3 * i + 1], r[3 * i + 2]};
            double r_i_proposed[3] = {r_proposed[3 * i], r_proposed[3 * i + 1], r_proposed[3 * i + 2]};
            double r_diff[3] = {r_i[0] - part_r[0], r_i[1] - part_r[1], r_i[2] - part_r[2]};
            double r_diff_proposed[3] = {r_i_proposed[0] - part_r_proposed[0], r_i_proposed[1] - part_r_proposed[1], r_i_proposed[2] - part_r_proposed[2]};

            potential += u(scalar_product(r_diff_proposed, r_diff_proposed), beta) - u(scalar_product(r_diff, r_diff), beta);
        }
    }

    double factor = exp(kinetic + potential);
    if (factor > 1.) {
        return 1.;
    } else {
        return factor;
    }
}