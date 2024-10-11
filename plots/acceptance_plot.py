import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import re
import glob

# Parameters
N = sys.argv[1]
n_steps = sys.argv[2]
case = sys.argv[3]

# Load data from CSV file
path_pattern = f"data/{case}_{N}_{n_steps}/acceptance_*.csv"
all_files = glob.glob(path_pattern)
all_files = glob.glob(path_pattern)
if not all_files:
	raise FileNotFoundError(f"No files found for pattern: {path_pattern}")
data = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
i = data['i']
a = data['a']

# Extract alpha from the filenames
alpha_values = []
for file in all_files:
	alpha_match = re.search(r'_(\d+\.\d+)', file)
	if alpha_match:
		alpha_values.append(float(alpha_match.group(1)))

alpha = np.mean(alpha_values) if alpha_values else 0.0

# Plot the data
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(8, 6))

plt.plot(i, a, linestyle='-', label='a')
plt.axhline(y=0.5, linestyle='--', color='tab:orange', label='50%')

plt.title(f'Acceptance rate plot with N = {N}, steps = {n_steps}, alpha = {alpha}')
plt.xlabel('steps')
plt.ylabel('acceptance rate')

plt.legend()
plt.yticks(np.arange(0, 1.1, 0.1))
plt.tight_layout()
#plt.grid(True)
plt.savefig(f'report/figures/{case}/acceptance_{N}_{n_steps}.png', dpi=500)
# plt.show()