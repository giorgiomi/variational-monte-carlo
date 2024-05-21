# Variational Monte Carlo
Computational physics simulation using the Variational Monte Carlo method to calculate the ground state energy for a system of $^4\text{He}$ at temperature $T = 0$ trapped in an external harmonic potential. 

## Theory

The hamiltonian used is
$$H = -\frac{\hbar^2}{2m} \sum_{i = 1}^N \nabla_i^2 + \frac{1}{2} \omega^2 \sum_{i = 1}^N r_i^2 + \sum_{i < j} V(r_{ij}),$$
where $V(r)$ is the classical Lennard-Jones potential:
$$V(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right].$$

The variational wave function used is
$$\Psi(\mathbf{r}_1,\dots,\mathbf{r}_N) = \exp{\left(-\frac{1}{2\alpha} \sum_{i = 1}^N r_i^2 - \frac{1}{2}\sum_{i < j} u_\beta(r_{ij})\right)},$$
with $u_\beta(r) = (\beta_1/r)^{\beta_2}$.

## Compilation and execution
Move to the `src` folder and compile using:
```
    make
```
An executable named `run` will be generated. Use
```
    ./run
```
to execute it.