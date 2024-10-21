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
path_pattern = f"data/{case}_{N}_{n_steps}/energy_*.csv"
all_files = glob.glob(path_pattern)
all_files = glob.glob(path_pattern)
if not all_files:
	raise FileNotFoundError(f"No files found for pattern: {path_pattern}")
data = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
i = data['i']
T = data['T']
VHO = data['VHO']
VLJ = data['VLJ']
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
VHO_avg = np.mean(VHO)
VHO_std = np.std(VHO)
VLJ_avg = np.mean(VLJ)
VLJ_std = np.std(VLJ)
E_avg = np.mean(E)
E_std = np.std(E)

# Plot VHO and VLJ
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(8, 6))

plt.plot(i, VHO, linestyle='-', alpha=0.5, label=r'$\langle V_{{HO}} \rangle$ = {:.3f} $\pm$ {:.3f} K'.format(VHO_avg, VHO_std/np.sqrt(len(VHO))))
plt.plot(i, VLJ, linestyle='-', alpha=0.5, label=r'$\langle V_{{LJ}} \rangle$ = {:.3f} $\pm$ {:.3f} K'.format(VLJ_avg, VLJ_std/np.sqrt(len(VLJ))))

plt.title(f'VHO and VLJ plot with N = {N}, steps = {n_steps}, alpha = {alpha}')
plt.xlabel('steps')
plt.ylabel('energy [K]')

plt.legend()
plt.tight_layout()
plt.savefig(f'report/figures/{case}/VHO_VLJ_{N}_{n_steps}.png', dpi=500)
if (show == 'show'): 
	plt.show()

# Plot T and E
plt.figure(figsize=(8, 6))

plt.plot(i, T, linestyle='-', alpha=0.5, label=r'$\langle T \rangle$ = {:.3f} $\pm$ {:.3f} K'.format(T_avg, T_std/np.sqrt(len(T))))
plt.plot(i, E, linestyle='-', label=r'$\langle E \rangle$ = {:.3f} $\pm$ {:.3f} K'.format(E_avg, E_std/np.sqrt(len(E))))

plt.title(f'T and E plot with N = {N}, steps = {n_steps}, alpha = {alpha}')
plt.xlabel('steps')
plt.ylabel('energy [K]')

plt.legend()
plt.tight_layout()
plt.savefig(f'report/figures/{case}/T_E_{N}_{n_steps}.png', dpi=500)
if (show == 'show'): 
	plt.show()