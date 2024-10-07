import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data from CSV file
data = pd.read_csv('data/kinetic_avg.csv')
i = data['i']
T_lap = data['T_lap']
T_grad = data['T_grad']

# compute averages and standard deviations
T_lap_avg = np.mean(T_lap)
T_lap_std = np.std(T_lap)

T_grad_avg = np.mean(T_grad)
T_grad_std = np.std(T_grad)

# Add average and standard deviation to the legend
T_lap_label = f'T_lap = {T_lap_avg:.2f} ± {T_lap_std/np.sqrt(len(T_lap)):.2f} K'
T_grad_label = f'T_grad = {T_grad_avg:.2f} ± {T_grad_std/np.sqrt(len(T_grad)):.2f} K'

# Plot the data
plt.figure(figsize=(10, 6))

plt.plot(i, T_lap, linestyle='-', color='b', label=T_lap_label)
plt.plot(i, T_grad, linestyle='-', color='g', label=T_grad_label, alpha=0.7)

plt.title('Kinetic energy averages')
plt.xlabel('steps')
plt.ylabel('T [K]')

plt.legend()
plt.grid(True)
plt.show()