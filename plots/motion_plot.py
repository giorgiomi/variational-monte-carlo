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
all_files = glob.glob(path_pattern)
if not all_files:
	raise FileNotFoundError(f"No files found for pattern: {path_pattern}")
data = pd.concat((pd.read_csv(f) for f in all_files), ignore_index=True)
#i = data['i']
pos = data['position']

# Extract alpha from the filenames
alpha_values = []
for file in all_files:
	alpha_match = re.search(r'_(\d+\.\d+)', file)
	if alpha_match:
		alpha_values.append(float(alpha_match.group(1)))

alpha = np.mean(alpha_values) if alpha_values else 0.0

import matplotlib.animation as animation

# Plot the data
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(8, 6))

# Reshape the position array to (n_steps, 3, N)
positions_reshaped = pos.values.reshape((int(n_steps), 3, int(N)))

# Extract x, y, z positions for each particle
x_positions = positions_reshaped[:, 0, :]
y_positions = positions_reshaped[:, 1, :]
z_positions = positions_reshaped[:, 2, :]

# Initialize the plot
lines = [ax.plot([], [], label=f'Particle {i+1}')[0] for i in range(int(N))]

ax.set_xlim(np.min(x_positions), np.max(x_positions))
ax.set_ylim(np.min(y_positions), np.max(y_positions))
ax.set_xlabel('X Position')
ax.set_ylabel('Y Position')
ax.set_title(f'Motion Plot for {N} Particles over {n_steps} Steps (Alpha = {alpha})')
ax.legend()
ax.grid(True)

# Animation update function
def update(frame):
    for i, line in enumerate(lines):
        line.set_data(x_positions[:frame, i], y_positions[:frame, i])
    return lines

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=int(n_steps), blit=True)

if show.lower() == 'show':
    plt.show()
else:
    ani.save(f"motion_plot_{case}_{N}_{n_steps}.mp4", writer='ffmpeg')

if show.lower() == 'show':
    plt.show()
else:
    plt.savefig(f"motion_plot_{case}_{N}_{n_steps}.png")
