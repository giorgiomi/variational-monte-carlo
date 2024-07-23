#ifndef METROPOLIS_H
#define METROPOLIS_H

/*
Returns the acceptance value for the M(RT)^2 step
r_old: previous configuration
r: current configuration
part_index: index of the particle to update
var_param: variational parameters [alpha, beta1, beta2]
N: number of particles
*/
double acceptance(double *r_old, double *r, int part_index, double *var_param, int N);

#endif