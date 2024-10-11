import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Set default font size
plt.rcParams.update({'font.size': 14})

# Parameters
N = sys.argv[1]
n_steps = sys.argv[2]

# Load the data
data = pd.read_csv(f'data/NOINT_{N}_{n_steps}/variational.csv')

alpha = data['alpha']
T = data['T']
T_lap = data['T_lap']
T_grad = data['T_grad']
V = data['V']
E = data['E']

T_std = data['T_std']
T_lap_std = data['T_lap_std']
T_grad_std = data['T_grad_std']
V_std = data['V_std']
E_std = data['E_std']

# Print energy
print(f'Energy estimate for N={N} and {n_steps} steps: E = {data["E"][data["alpha"] == 25.0].values[0]:.4f} +/- {data["E_std"][data["alpha"] == 25.0].values[0]:.4f} K')

# Create subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))

# First subplot: Plot the data
axs[0].errorbar(alpha, T, yerr=T_std, capsize=2, label=r'$\langle T \rangle$', fmt='-')
axs[0].errorbar(alpha, T_lap, yerr=T_lap_std, capsize=2, label=r'$\langle T \rangle_{lap}$', fmt='-')
axs[0].errorbar(alpha, T_grad, yerr=T_grad_std, capsize=2, label=r'$\langle T \rangle_{grad}$', fmt='-')
axs[0].errorbar(alpha, V, yerr=V_std, capsize=2, label=r'$\langle V \rangle$', fmt='-')
axs[0].errorbar(alpha, E, yerr=E_std, capsize=2, label=r'$\langle E \rangle$', fmt='-')

# Add labels and title to the first subplot
axs[0].set_xlabel(r'$\alpha$ $[Å^2]$')
axs[0].set_ylabel('Energy [K]')
axs[0].set_title(f'T, V, E for N={N} and {n_steps} steps')
axs[0].legend(loc='upper right')

# Second subplot: Plot the standard deviations
axs[1].plot(alpha, T_std, label=r'$\sigma(T)$')
axs[1].plot(alpha, np.zeros(len(T_std)), label=r'$\sigma(T_{lap})$')
axs[1].plot(alpha, T_grad_std, label=r'$\sigma(T_{grad})$')
axs[1].plot(alpha, V_std, label=r'$\sigma(V)$')
axs[1].plot(alpha, E_std, label=r'$\sigma(E)$')

# Add labels and title to the second subplot
axs[1].set_xlabel(r'$\alpha$ $[Å^2]$')
axs[1].set_ylabel(r'$\sigma$ [K]')
axs[1].set_title(f'StDev for N={N} and {n_steps} steps')
axs[1].legend(loc='upper right')

# Adjust layout and show the plot
plt.tight_layout()
plt.savefig(f'report/figures/NOINT/variational_{N}_{n_steps}.png', dpi=500)
# plt.show()