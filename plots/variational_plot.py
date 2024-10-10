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

# Plot the data
plt.figure(figsize=(10, 6))
plt.errorbar(alpha, T, yerr=T_std, capsize=2, label='T', fmt='-')
plt.errorbar(alpha, T_lap, yerr=T_lap_std, capsize=2, label='T_lap', fmt='-')
plt.errorbar(alpha, T_grad, yerr=T_grad_std, capsize=2, label='T_grad', fmt='-')
# plt.errorbar(alpha, V, yerr=V_std, capsize=2, label='V', fmt='-o')
# plt.errorbar(alpha, E, yerr=E_std, capsize=2, label='E', fmt='-o')

# Add labels and title
plt.xlabel(r'Alpha $[Å^2]$')
plt.ylabel('Energy [K]')
plt.title(f'T, V, E for N={N} and {n_steps} steps')
plt.legend()

# Show the plot
plt.show()