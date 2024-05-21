#ifndef FUNCTIONS_H
#define FUNCTIONS_H

/*
Returns the scalar product of two vectors a and b (pointers)
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
Returns the Lennard-Jones potential V(r) for two particles at distance r
[eps, sigma] parameters
*/
double lennard_jones(double r, double *param);

/*
Calculates the kinetic energy for N particles
*/
double kinetic_energy(double *r, double *param, int N);

#endif