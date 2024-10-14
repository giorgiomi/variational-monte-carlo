# Variational Monte Carlo
Computational physics simulation using the Variational Monte Carlo method to calculate the ground state energy for a system of $^4\text{He}$ at temperature $T = 0$ trapped in an external harmonic potential. 

## Theory

The hamiltonian used is
```math
H = -\frac{\hbar^2}{2m} \sum_{i = 1}^N \nabla_i^2 + \frac{1}{2} m\omega^2 \sum_{i = 1}^N r_i^2 + \sum_{i < j} V(r_{ij})\ ,
```

where $V(r)$ is the classical Lennard-Jones potential:
```math
V(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]\ .
```

The variational wave function used is
```math
\Psi(r_1,\dots,r_N)=\exp{\left(-\frac{1}{2\alpha}\sum_{i = 1}^Nr_i^2-\frac{1}{2}\sum_{i < j}u_\beta(r_{ij})\right)}\ ,
```

with $u_\beta(r) = (\beta_1/r)^{\beta_2}$.

## Compilation and execution
Compile in the repository directory using:
```
    make
```
An executable named `run` will be generated. Use
```
    ./run N n_steps [alpha_saved]
```
to execute it, where `N` is the number of particles, `n_steps` is the number of Metropolis steps. `alpha_saved` is optional but it's needed to look at the observables values at a certain value of $\alpha$. For example, executing
```
    ./run 2 100000 25.0
```
will save the observable values of the simulation for $\alpha = 25\ \text{Å}^2$. Then, to plot the results one can use
```
    python3 plots/energy_plot.py N n_steps INT show
```
where `INT` and `show` are needed for plotting the interacting simulation results and to show the graph. If the simulation is run with multiple values of the parameter $\alpha$, one can look at the energy as a function of that, by running `plots/variational_plot.py`. At the moment, this can be achieved by modifying `alpha_start` and `alpha_end` in `main.c`.