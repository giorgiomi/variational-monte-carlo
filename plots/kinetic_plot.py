import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import glob
import sys

# Parameters
N = sys.argv[1]
n_steps = sys.argv[2]
case = sys.argv[3]
show = sys.argv[4]

# Load data from CSV file
path_pattern = f"data/{case}_{N}_{n_steps}/kinetic_avg_*.csv"
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
T_lap_label =r'$\langle T \rangle_{lap}$' + f' = {T_lap_avg:.3f} $\pm$ {T_lap_std/np.sqrt(len(T_lap)):.3f} K'
T_grad_label =r'$\langle T \rangle_{grad}$' + f' = {T_grad_avg:.3f} $\pm$ {T_grad_std/np.sqrt(len(T_grad)):.3f} K'

# Plot the data
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(8, 6))

plt.plot(i, T_grad, linestyle='-', label=T_grad_label)
plt.plot(i, T_lap, linestyle='-', label=T_lap_label)

plt.title(f'Kinetic energy estimators with N = {N}, n_steps = {n_steps}, alpha = {alpha}')
plt.xlabel('steps')
plt.ylabel('T [K]')

plt.legend()
plt.tight_layout()
#plt.grid(True)
plt.savefig(f'report/figures/{case}/kinetic_{N}_{n_steps}.png', dpi=500)
if (show == 'show'):
    plt.show()