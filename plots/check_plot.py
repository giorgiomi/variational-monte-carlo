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
path_pattern = f"data/{case}_{N}_{n_steps}/check_*.csv"
all_files = glob.glob(path_pattern)
all_files = glob.glob(path_pattern)
if not all_files:
	raise FileNotFoundError(f"No files found for pattern: {path_pattern}")
data = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
i = data['i']
a = data['a']
T_old = data['T_old']
T_new = data['T_new']
VHO_old = data['VHO_old']
VHO_new = data['VHO_new']
VLJ_old = data['VLJ_old']
VLJ_new = data['VLJ_new']

# Extract alpha from the filenames
alpha_values = []
for file in all_files:
	alpha_match = re.search(r'_(\d+\.\d+)', file)
	if alpha_match:
		alpha_values.append(float(alpha_match.group(1)))

alpha = np.mean(alpha_values) if alpha_values else 0.0

# # Plot VLJ_new/VLJ_old and a
# plt.rcParams.update({'font.size': 14})
# plt.figure(figsize=(8, 6))

# plt.plot(i, np.abs(VLJ_new/VLJ_old), linestyle='-', alpha=0.5, label=r'$\langle V_{{LJ,new}}/V_{{LJ,old}} \rangle$')
# plt.plot(i, a, linestyle='-', alpha=0.5, label=r'$\langle a \rangle$')

# plt.title(f'VLJ_new/VLJ_old and a plot with N = {N}, steps = {n_steps}, alpha = {alpha}')
# plt.xlabel('steps')
# plt.ylabel('check')
# plt.legend()
# plt.tight_layout()
# plt.yscale('log')
	
# # Plot VHO_new/VHO_old and a
# plt.figure(figsize=(8, 6))

# plt.plot(i, np.abs(VHO_new/VHO_old), linestyle='-', alpha=0.5, label=r'$\langle V_{{HO,new}}/V_{{HO,old}} \rangle$')
# plt.plot(i, a, linestyle='-', alpha=0.5, label=r'$\langle a \rangle$')

# plt.title(f'VHO_new/VHO_old and a plot with N = {N}, steps = {n_steps}, alpha = {alpha}')
# plt.xlabel('steps')
# plt.ylabel('check')
# plt.legend()
# plt.tight_layout()
# #plt.yscale('log')

# Plot VHO_new/VHO_old and VLJ_new/VLJ_old
plt.figure(figsize=(8, 6))

plt.plot(i, np.abs(VHO_new/VHO_old), linestyle='-', alpha=0.5, label=r'$\langle V_{{HO,new}}/V_{{HO,old}} \rangle$')
plt.plot(i, np.abs(VLJ_new/VLJ_old), linestyle='-', alpha=0.5, label=r'$\langle V_{{LJ,new}}/V_{{LJ,old}} \rangle$')

plt.title(f'VHO_new/VHO_old and VLJ_new/VLJ_old plot with N = {N}, steps = {n_steps}, alpha = {alpha}')
plt.xlabel('steps')
plt.ylabel('check')
plt.legend()
plt.tight_layout()
plt.yscale('log')

if (show == 'show'): 
    plt.show()