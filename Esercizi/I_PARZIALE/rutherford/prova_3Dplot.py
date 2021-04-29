import numpy as np
import matplotlib.pyplot as plt

z = np.linspace(0, 15, 1000)
x = np.sin(z)
y = np.cos(z)

fig = plt.subplots()
ax = plt.axes(projection = '3d')
ax.plot3D(x, y, z)
ax.set_zlabel('prova asse z')

plt.show()