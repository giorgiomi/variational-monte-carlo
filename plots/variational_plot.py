import pandas as pd

import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv('data/variational.csv')

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

# Load parameters
parameters = pd.read_csv('data/param.csv')
N = parameters['N'][0]
n_steps = parameters['n_steps'][0]

# Create subplots
fig, axs = plt.subplots(2, 1, figsize=(8, 8))

# First subplot: Plot the data
axs[0].errorbar(alpha, T, yerr=T_std, capsize=2, label='T', fmt='-')
axs[0].errorbar(alpha, T_lap, yerr=T_lap_std, capsize=2, label='T_lap', fmt='-')
axs[0].errorbar(alpha, T_grad, yerr=T_grad_std, capsize=2, label='T_grad', fmt='-')
axs[0].errorbar(alpha, V, yerr=V_std, capsize=2, label='V', fmt='-')
axs[0].errorbar(alpha, E, yerr=E_std, capsize=2, label='E', fmt='-')

# Add labels and title to the first subplot
axs[0].set_xlabel(r'$\alpha$ $[Å^2]$')
axs[0].set_ylabel('Energy [K]')
axs[0].set_title(f'T, V, E for N={N} and {n_steps} steps')
axs[0].legend()

# Second subplot: Plot the standard deviations
axs[1].plot(alpha, T_std, label='T_std')
axs[1].plot(alpha, T_lap_std, label='T_lap_std')
axs[1].plot(alpha, T_grad_std, label='T_grad_std')
axs[1].plot(alpha, V_std, label='V_std')
axs[1].plot(alpha, E_std, label='E_std')

# Add labels and title to the second subplot
axs[1].set_xlabel(r'$\alpha$ $[Å^2]$')
axs[1].set_ylabel('Standard Deviation [K]')
axs[1].set_title(f'Standard Deviations for N={N} and {n_steps} steps')
axs[1].legend()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()