import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data from CSV file
data = pd.read_csv('data/energy.csv')
i = data['i']
T = data['T']
V = data['V']
E = data['E']

# compute averages and standard deviations
T_avg = np.mean(T)
T_std = np.std(T)

# Plot the data
plt.figure(figsize=(10, 6))

plt.plot(i, T, linestyle='-', color='b', label='T')
# plt.axhline(y=T_avg, color='r', linestyle='-', label=r'$\langle T \rangle$')
# plt.axhline(y=T_avg + T_std, color='r', linestyle='--', label=r'$\langle T \rangle + \sigma$')
# plt.axhline(y=T_avg - T_std, color='r', linestyle='--', label=r'$\langle T \rangle - \sigma$')

# plt.plot(i, V, linestyle='-', color='g', label='V')
# plt.plot(i, E, linestyle='-', color='r', label='E')

plt.title('Energy plot')
plt.xlabel('steps')
plt.ylabel('energy')

# plt.yscale('log')

plt.legend()
plt.grid(True)
plt.show()