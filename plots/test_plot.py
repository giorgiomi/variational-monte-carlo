import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

data = pd.read_csv(f'data/test.csv')
AlphaValues = np.arange(21.0, 36.0, 1.0)
BetaValues = np.arange(1.95, 2.66, 0.05)
Energies = np.array(data['energies'])
Energies = Energies.reshape((15, 15))
print(len(Energies))
print(len(AlphaValues))
print(len(BetaValues))


# Prepare for plots
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Plot the surface.
X, Y = np.meshgrid(AlphaValues, BetaValues)
surf = ax.plot_surface(X, Y, Energies,cmap=cm.coolwarm,linewidth=0, antialiased=False)
# Customize the z axis.
zmin = np.matrix(Energies).min()
zmax = np.matrix(Energies).max()
ax.set_zlim(zmin, zmax)
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'$\beta$')
ax.set_zlabel(r'$\langle E \rangle$')
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()