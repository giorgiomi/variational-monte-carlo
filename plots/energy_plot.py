import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import glob
import sys

# Parameters
N = sys.argv[1]
n_steps = sys.argv[2]

# Load data from CSV file
path_pattern = f"data/NOINT_{N}_{n_steps}/energy_*.csv"
all_files = glob.glob(path_pattern)
all_files = glob.glob(path_pattern)
if not all_files:
	raise FileNotFoundError(f"No files found for pattern: {path_pattern}")
data = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
i = data['i']
T = data['T']
V = data['V']
E = data['E']

# Extract alpha from the filenames
alpha_values = []
for file in all_files:
	alpha_match = re.search(r'_(\d+\.\d+)', file)
	if alpha_match:
		alpha_values.append(float(alpha_match.group(1)))

alpha = np.mean(alpha_values) if alpha_values else 0.0

# compute averages and standard deviations
T_avg = np.mean(T)
T_std = np.std(T)
V_avg = np.mean(V)
V_std = np.std(V)
E_avg = np.mean(E)
E_std = np.std(E)

# Plot the data
plt.figure(figsize=(10, 6))

plt.plot(i, T, linestyle='-', label=r'$\langle T \rangle$ = {:.3f} $\pm$ {:.3f}'.format(T_avg, T_std/np.sqrt(len(T))))
plt.plot(i, V, linestyle='-', label=r'$\langle V \rangle$ = {:.3f} $\pm$ {:.3f}'.format(V_avg, V_std/np.sqrt(len(V))))
plt.plot(i, E, linestyle='-', label=r'$\langle E \rangle$ = {:.3f} $\pm$ {:.3f}'.format(E_avg, E_std/np.sqrt(len(E))))

plt.title(f'Energy plot, N = {N}, steps = {n_steps}, alpha = {alpha}')
plt.xlabel('steps')
plt.ylabel('energy')

# plt.yscale('log')

plt.legend()
plt.grid(True)
plt.show()