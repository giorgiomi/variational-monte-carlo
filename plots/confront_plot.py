import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Set default font size
plt.rcParams.update({'font.size': 14})

# Load the data
data2 = pd.read_csv(f'data/NOINT_2_1000000/variational.csv')
data4 = pd.read_csv(f'data/NOINT_4_1000000/variational.csv')
data6 = pd.read_csv(f'data/NOINT_6_1000000/variational.csv')
data8 = pd.read_csv(f'data/NOINT_8_1000000/variational.csv')

alpha = data2['alpha']

E2 = data2['E']
E2_std = data2['E_std']

E4 = data4['E']
E4_std = data4['E_std']

E6 = data6['E']
E6_std = data6['E_std']

E8 = data8['E']
E8_std = data8['E_std']

# Find minimum energy
min_E2 = np.min(E2/2)
min_E4 = np.min(E4/4)
min_E6 = np.min(E6/6)
min_E8 = np.min(E8/8)

pos_min_E2 = np.argmin(E2)
pos_min_E4 = np.argmin(E4)
pos_min_E6 = np.argmin(E6)
pos_min_E8 = np.argmin(E8)

alpha_min_E2 = alpha[pos_min_E2]
alpha_min_E4 = alpha[pos_min_E4]
alpha_min_E6 = alpha[pos_min_E6]
alpha_min_E8 = alpha[pos_min_E8]

print(f"Minimum E2: {min_E2} at alpha {alpha_min_E2}")
print(f"Minimum E4: {min_E4} at alpha {alpha_min_E4}")
print(f"Minimum E6: {min_E6} at alpha {alpha_min_E6}")
print(f"Minimum E8: {min_E8} at alpha {alpha_min_E8}")

# Create subplots
fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))

# First subplot: Plot the data
axs[0].errorbar(alpha, E2/2, yerr=E2_std, capsize=2, label=r'$N=2$', fmt='-')
axs[0].errorbar(alpha, E4/4, yerr=E4_std, capsize=2, label=r'$N=4$', fmt='-')
axs[0].errorbar(alpha, E6/6, yerr=E6_std, capsize=2, label=r'$N=6$', fmt='-')
axs[0].errorbar(alpha, E8/8, yerr=E8_std, capsize=2, label=r'$N=8$', fmt='-')

# Plot the minimum energy values
axs[0].scatter(alpha_min_E2, min_E2, zorder=5, color='tab:gray')
axs[0].hlines(min_E2, alpha_min_E2 - 10, alpha_min_E2 + 10, colors='tab:gray', linestyles='dashed', zorder=4)
axs[0].vlines(alpha_min_E2, min_E2 - 0.1, min_E2, colors='tab:gray', linestyles='dashed', zorder=4)

# Add a separate box for the minimum energy label
bbox_props = dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.7)
axs[0].text(alpha_min_E2 + 9, min_E2 - 0.05, f'$(E/N)_{{min}} \simeq {min_E2:.3f}$ K', ha="center", va="bottom", bbox=bbox_props)

# Add another box for the alpha value
bbox_props = dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white", alpha=0.7)
axs[0].text(alpha_min_E2 + 9, min_E2 - 0.105, f'$\\alpha_{{min}} \simeq {alpha_min_E2:.3f}$ $Å^2$', ha="center", va="bottom", bbox=bbox_props)



# Add labels and title to the first subplot
axs[0].set_xlabel(r'$\alpha$ $[Å^2]$')
axs[0].set_ylabel('E/N [K]')
axs[0].set_title(f'E/N for N=2,4,6,8 and 1000000 steps')
axs[0].legend(loc='upper right')

# Second subplot: Plot the standard deviations
axs[1].plot(alpha, E2_std/np.sqrt(2), label=r'$\sigma(E)_{N=2}$')
axs[1].plot(alpha, E4_std/np.sqrt(4), label=r'$\sigma(E)_{N=4}$')
axs[1].plot(alpha, E6_std/np.sqrt(6), label=r'$\sigma(E)_{N=6}$')
axs[1].plot(alpha, E8_std/np.sqrt(8), label=r'$\sigma(E)_{N=8}$')

# Add labels and title to the second subplot
axs[1].set_xlabel(r'$\alpha$ $[Å^2]$')
axs[1].set_ylabel(r'$\sigma/N$ [K]')
axs[1].set_title(f'StDev for N=2,4,6,8 and 1000000 steps')
axs[1].legend(loc='upper right')

# Adjust layout and show the plot
plt.tight_layout()
plt.savefig('report/figures/NOINT/confront_1000000.png', dpi=500)
plt.show()