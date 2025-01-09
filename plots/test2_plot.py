import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("data/test/MC.csv")
m = data['m']
E = data['E']
print(E)

plt.plot(m,E)
plt.show()