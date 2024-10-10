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
V_avg = np.mean(V)
V_std = np.std(V)
E_avg = np.mean(E)
E_std = np.std(E)

# Plot the data
plt.figure(figsize=(10, 6))

plt.plot(i, T, linestyle='-', color='b', label=r'$\langle T \rangle$ = {:.3f} $\pm$ {:.3f}'.format(T_avg, T_std/np.sqrt(len(T))))
plt.plot(i, V, linestyle='-', color='g', label=r'$\langle V \rangle$ = {:.3f} $\pm$ {:.3f}'.format(V_avg, V_std/np.sqrt(len(V))))
plt.plot(i, E, linestyle='-', color='r', label=r'$\langle E \rangle$ = {:.3f} $\pm$ {:.3f}'.format(E_avg, E_std/np.sqrt(len(E))))

plt.title('Energy plot')
plt.xlabel('steps')
plt.ylabel('energy')

# plt.yscale('log')

plt.legend()
plt.grid(True)
plt.show()