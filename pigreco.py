#programma per calcolare pigreco con una buona approssimazione

import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

pi_greco = np.pi
x_list = []
y_list = []
n = 0               #contatore dei cicli
N = 10000000           #cicli da fare

np.random.seed()    #decide numeri casuali in base allo stato del computer rendendoli casuali

for i in range(N):
    x = np.random.rand()
    y = np.random.rand()
    x_list.append(x)
    y_list.append(y)
    r = np.sqrt(x**2 + y**2)

    if r <= 1:
        n += 1

print(4*(n/N))

x_circle = np.linspace(0,1,1000)
y_circle = np.sqrt(1-x_circle**2)

fig, ax = plt.subplots()
ax.axis('equal')
ax.set_xlabel('x (A.U.)')
ax.set_ylabel('y (A.U.)')
#ax.add_patch(patches.Rectangle(0,0), 1,1, linewidth=1, edgecolor='b')
ax.plot(x_list,y_list, linewidth = 0, marker = '.', markersize = 1)
ax.plot(x_circle,y_circle, color = 'r')


plt.show()
