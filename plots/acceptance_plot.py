import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load data from CSV file
data = pd.read_csv('data/acceptance.csv')
i = data['i']
a = data['a']

# Plot the data
plt.figure(figsize=(10, 6))

plt.plot(i, a, linestyle='-', color='b', label='a')

plt.title('Acceptance rate')
plt.xlabel('steps')
plt.ylabel('a')

plt.legend()
plt.grid(True)
plt.show()