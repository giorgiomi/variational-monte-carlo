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
show = sys.argv[4]

# Load data from CSV file
path_pattern = f"data/{case}_{N}_{n_steps}/pos_*.csv"
all_files = glob.glob(path_pattern)
if not all_files:
    raise FileNotFoundError(f"No files found for pattern: {path_pattern}")
data = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)

# Extract alpha from the filenames
alpha_values = []
for file in all_files:
    alpha_match = re.search(r'_(\d+\.\d+)', file)
    if alpha_match:
        alpha_values.append(float(alpha_match.group(1)))

alpha = np.mean(alpha_values) if alpha_values else 0.0

# Extract x, y, z positions
x_positions = data[[f'x{i}' for i in range(int(N))]].values
y_positions = data[[f'y{i}' for i in range(int(N))]].values
z_positions = data[[f'z{i}' for i in range(int(N))]].values

# Plot x1, x2
time_steps = np.arange(int(n_steps) + 1)
fig, axs = plt.subplots(3, 1, figsize=(12, 8))
for i in range(int(N)):  
    axs[0].plot(time_steps, x_positions[:, i], label=f'x{i}')
    axs[1].plot(time_steps, y_positions[:, i], label=f'y{i}')
    axs[2].plot(time_steps, z_positions[:, i], label=f'z{i}')

# Set labels and titles for each subplot
axs[0].set_xlabel('Time Step')
axs[0].set_ylabel('X Position')
axs[0].set_title(f'X Position vs Time Step for {N} Particles (Alpha = {alpha})')
axs[0].legend()
axs[0].grid(True)

axs[1].set_xlabel('Time Step')
axs[1].set_ylabel('Y Position')
axs[1].set_title(f'Y Position vs Time Step for {N} Particles (Alpha = {alpha})')
axs[1].legend()
axs[1].grid(True)

axs[2].set_xlabel('Time Step')
axs[2].set_ylabel('Z Position')
axs[2].set_title(f'Z Position vs Time Step for {N} Particles (Alpha = {alpha})')
axs[2].legend()
axs[2].grid(True)

plt.tight_layout()

# Calculate the distance between the first two particles
distances = np.sqrt((x_positions[:, 0] - x_positions[:, 1])**2 +
                    (y_positions[:, 0] - y_positions[:, 1])**2 +
                    (z_positions[:, 0] - z_positions[:, 1])**2)

# Create a new figure for the distance plot
fig_dist, ax_dist = plt.subplots(figsize=(10, 5))
ax_dist.plot(time_steps, distances, label='Distance between Particle 1 and 2')
ax_dist.axhline(y=2.556, linestyle='--', label=r'$\sigma$')
ax_dist.set_xlabel('Time Step')
ax_dist.set_ylabel('Distance')
ax_dist.set_title(f'Distance between Particle 1 and 2 vs Time Step (Alpha = {alpha})')
ax_dist.legend()
ax_dist.grid(True)

plt.tight_layout()

# Calculate LJ potential
eps = 10.22
sigma = 2.556
def lj_potential(r):
    return 4 * eps * ((sigma/r)**12 - (sigma/r)**6)

VLJ = lj_potential(distances)
VLJ_avg = np.mean(VLJ)

# Calculate harmonic potential


# Create a new figure for the LJ potential plot
fig_lj, ax_lj = plt.subplots(figsize=(10, 5))
ax_lj.plot(time_steps, VLJ, label='LJ Potential between Particle 1 and 2')
ax_lj.axhline(y=VLJ_avg, linestyle='--', label='LJ_avg = {:.2f}'.format(VLJ_avg))
ax_lj.set_xlabel('Time Step')
ax_lj.set_ylabel('LJ Potential')
ax_lj.set_title(f'LJ Potential between Particle 1 and 2 vs Time Step (Alpha = {alpha})')
ax_lj.legend()
ax_lj.grid(True)
#ax_lj.set_yscale('log')

plt.tight_layout()

if show.lower() == 'show':
    plt.show()
else:
    plt.savefig(f"motion_plot_{case}_{N}_{n_steps}.png")
