import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D

# Set default font size
plt.rcParams.update({'font.size': 14})

# Parameters
N = sys.argv[1]
n_steps = sys.argv[2]
case = sys.argv[3]

# Load the data
data = pd.read_csv(f'data/{case}_{N}_{n_steps}/variational.csv')

alpha = data['alpha']
beta1 = data['beta1']  # Assuming beta1 is another column in your CSV
T = data['T']
T_lap = data['T_lap']
T_grad = data['T_grad']
VHO = data['VHO']
VLJ = data['VLJ']
E = data['E']

T_std = data['T_std']
T_lap_std = data['T_lap_std']
T_grad_std = data['T_grad_std']
VHO_std = data['VHO_std']
VLJ_std = data['VLJ_std']
E_std = data['E_std']

# Find minimum of E
min_E = np.min(E)
alpha_min = alpha[np.argmin(E)]
beta1_min = beta1[np.argmin(E)]
E_std_min = E_std[np.argmin(E)]

# Create a 3D plot
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_proj_type('ortho')

# Plot the data
# ax.plot_trisurf(alpha, beta1, T, label=r'$\langle T \rangle$', alpha=0.7)
# ax.plot_trisurf(alpha, beta1, T_lap, label=r'$\langle T \rangle_{lap}$', alpha=0.7)
# ax.plot_trisurf(alpha, beta1, T_grad, label=r'$\langle T \rangle_{grad}$', alpha=0.7)
# ax.plot_trisurf(alpha, beta1, VHO, label=r'$\langle VHO \rangle$', alpha=0.7)
# ax.plot_trisurf(alpha, beta1, VLJ, label=r'$\langle VLJ \rangle$', alpha=0.7)
ax.plot_trisurf(alpha, beta1, E, label=r'$\langle E \rangle$', alpha=0.7)

# Draw a point at the minimum energy
ax.scatter(alpha_min, beta1_min, min_E, color='red', s=100, label='Min E')
# Add a text box with alpha_min and beta1_min
textstr = f'α_min = {alpha_min:.2f} Å^2\nβ1_min = {beta1_min:.2f} Å\nE_min = {min_E:.3f} ± {E_std_min:.3f} K'
props = dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.5)
ax.text2D(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
          verticalalignment='top', bbox=props)

# Add labels and title
ax.set_xlabel(r'$\alpha$ $[Å^2]$')
ax.set_ylabel(r'$\beta_1$ $[Å]$')
ax.set_zlabel('Energy [K]')
ax.set_title(f'3D Plot for N={N} and {n_steps} steps')

# Show the plot
plt.legend()
plt.tight_layout()

# Plot E_std
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_proj_type('ortho')

# Plot the data
ax.plot_trisurf(alpha, beta1, E_std, label=r'$\langle E \rangle$', alpha=0.7)

# Add labels and title
ax.set_xlabel(r'$\alpha$ $[Å^2]$')
ax.set_ylabel(r'$\beta_1$ $[Å]$')
ax.set_zlabel('Energy std [K]')
ax.set_title(f'3D Plot for N={N} and {n_steps} steps')

# Show the plot
plt.legend()
plt.tight_layout()


plt.show()
