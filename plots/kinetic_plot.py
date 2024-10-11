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
path_pattern = f"data/NOINT_{N}_{n_steps}/kinetic_avg_*.csv"
all_files = glob.glob(path_pattern)
if not all_files:
	raise FileNotFoundError(f"No files found for pattern: {path_pattern}")
data = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
i = data['i']
T_lap = data['T_lap']
T_grad = data['T_grad']

# Extract alpha from the filenames
alpha_values = []
for file in all_files:
	alpha_match = re.search(r'_(\d+\.\d+)', file)
	if alpha_match:
		alpha_values.append(float(alpha_match.group(1)))

alpha = np.mean(alpha_values) if alpha_values else 0.0

# compute averages and standard deviations
T_lap_avg = np.mean(T_lap)
T_lap_std = np.std(T_lap)

T_grad_avg = np.mean(T_grad)
T_grad_std = np.std(T_grad)

# Add average and standard deviation to the legend
T_lap_label = f'T_lap = {T_lap_avg:.2f} ± {T_lap_std/np.sqrt(len(T_lap)):.2f} K'
T_grad_label = f'T_grad = {T_grad_avg:.2f} ± {T_grad_std/np.sqrt(len(T_grad)):.2f} K'

# Plot the data
plt.figure(figsize=(10, 6))

plt.plot(i, T_lap, linestyle='-', color='b', label=T_lap_label)
plt.plot(i, T_grad, linestyle='-', color='g', label=T_grad_label, alpha=0.7)

plt.title(f'Kinetic energy averages, N = {N}, steps = {n_steps}, alpha = {alpha}')
plt.xlabel('steps')
plt.ylabel('T [K]')

plt.legend()
plt.grid(True)
plt.show()