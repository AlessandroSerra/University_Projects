import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

L = 1
Nx, Ny = 1000, 1000
x, y = np.linspace(-L, L, Nx), np.linspace(-L, L, Ny)
x, y = np.meshgrid(x, y)        #crea una griglia con tutte le righe uguali a partire da una lista

z1 = np.sin(2 * np.pi * (x + y))
z2 = np.sin(2 * np.pi * x * y)

fig1 = plt.figure()
ax3D1 = fig1.add_subplot(111, projection = '3d')
ax3D1.plot_wireframe(x, y, z1)

fig2 = plt.figure()
ax3D2 = fig2.add_subplot(111, projection = '3d')
ax3D2.plot_wireframe(x, y, z2)

plt.show()