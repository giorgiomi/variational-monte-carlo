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
path_pattern = f"data/{case}_{N}_{n_steps}/positions_*.csv"
all_files = glob.glob(path_pattern)
if not all_files:
    raise FileNotFoundError(f"No files found for pattern: {path_pattern}")
data = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
pos = data['position']

# Extract alpha from the filenames
alpha_values = []
for file in all_files:
    alpha_match = re.search(r'_(\d+\.\d+)', file)
    if alpha_match:
        alpha_values.append(float(alpha_match.group(1)))

alpha = np.mean(alpha_values) if alpha_values else 0.0

# Plot the data
plt.rcParams.update({'font.size': 14})
fig, axs = plt.subplots(3, 1, figsize=(15, 8))

# Reshape the position array to (n_steps, 3, N)
positions_reshaped = pos.values.reshape((int(n_steps), 3, int(N)))

# Extract x, y, z positions for each particle
x_positions = positions_reshaped[:, 0, :]
y_positions = positions_reshaped[:, 1, :]
z_positions = positions_reshaped[:, 2, :]

# Plot x1, x2
time_steps = np.arange(int(n_steps))
for i in range(min(2, int(N))):  # Plot only the first two particles
    axs[0].plot(time_steps, x_positions[:, i], label=f'x{i+1}')
    axs[1].plot(time_steps, y_positions[:, i], label=f'y{i+1}')
    axs[2].plot(time_steps, z_positions[:, i], label=f'z{i+1}')

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
ax_dist.set_xlabel('Time Step')
ax_dist.set_ylabel('Distance')
ax_dist.set_title(f'Distance between Particle 1 and 2 vs Time Step (Alpha = {alpha})')
ax_dist.legend()
ax_dist.grid(True)

plt.tight_layout()

if show.lower() == 'show':
    plt.show()
else:
    plt.savefig(f"motion_plot_{case}_{N}_{n_steps}.png")
