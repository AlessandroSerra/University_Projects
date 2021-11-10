import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0,10000)
y = x

fig, ax = plt.subplots()
ax.plot(x,y, label = 'dati')
ax.set_xscale('log')
ax.set_yscale('log')
plt.show()