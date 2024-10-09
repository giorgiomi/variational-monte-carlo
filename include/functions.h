#ifndef FUNCTIONS_H
#define FUNCTIONS_H

/*
Copy array a to array b
*/
void copy_array(double *a, double *b, int N);

/*
Print array a as such: a[0] a[1] ... a[N-1]
*/
void print_array(double *a, int N);

/*
Returns the scalar product of two 3D vectors a and b (pointers)
*/
double scalar_product(double *a, double *b);

/*
Returns the function u of a radius r and beta (pointer to beta1 and beta2)
*/
double u(double r, double *beta);

/*
Returns the first derivative of the function u of a radius r and beta (pointer to beta1 and beta2)
*/
double u_prime(double r, double *beta);

/*
Returns the second derivative of the function u of a radius r and beta (pointer to beta1 and beta2)
*/
double u_doubleprime(double r, double *beta);

/*
Returns the wavefunction for N particles
*/
double psi(double *r, double *param, int N);

/*
Returns the Lennard-Jones potential V(r) for two particles at distance r
[EPS, SIGMA] constants
*/
double lennard_jones(double r);

/*
Calculates the kinetic energy for N particles
*/
double kinetic_energy(double *r, double *param, int N);

/*
Calculates the average of the laplactian terms of the kinetic energy for N particles
*/
double kinetic_average_laplacian(double *r, double *param, int N);

/*
Calculates the average of the square gradient terms of the kinetic energy for N particles
*/
double kinetic_average_gradient(double *r, double *param, int N);

/*
Calculates the harmonic potential for N particles after a single particle moved
r_old: old positions
r_new: new positions
VH_old: old harmonic potential
part_index: index of the particle that moved
N: number of particles
*/
double harmonic_potential(double *r_old, double *r_new, double VH_old, int part_index, int N);

/*
Calculates the Lennard-Jones potential contribution for a single particle
r: positions
part_index: index of the particle
N: number of particles
*/
double LJ_potential_single(double *r, int part_index, int N);

/*
Calculates the Lennard-Jones potential for N particles after a single particle moved
r_old: old positions
r_new: new positions
VLJ_old: old LJ potential
part_index: index of the particle that moved
N: number of particles
*/
double LJ_potential(double *r_old, double *r_new, double VH_old, int part_index, int N);

/*
Calculates the potential energy for N particles (LJ + harmonic) by brute force
r: positions
var_param: variational parameters
N: number of particles
*/
double potential_bruteforce(double *r, double *var_param, int N);

#endif