#ifndef METROPOLIS_H
#define METROPOLIS_H

/*
Returns the acceptance value for the M(RT)^2 step
E: energy
E_proposed: energy of the proposed configuration
*/
double acceptance(double E, double E_proposed);

#endif