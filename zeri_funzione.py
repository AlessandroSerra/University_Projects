import numpy as np
import matplotlib.pyplot as plt
import time

def f(x):
    return np.log(x)

inter = [1e-3, 1.5]     #modo compatto per indicare le potenze di 10
soglia = 1e-3

if f(inter[0]*inter[1]) > 0:
    print('metodo non utilizzabile')

else:
    fm = 1
    while np.abs(fm) > soglia:
        m = np.mean(inter)
        fm = f(m)

        if f(inter[0]) * fm < 0:
            inter = [inter[0],m]

        else:
            inter = [m,inter[1]]

print('zero at m:' + str(m))

x = np.linspace(1e-3, 1, 1000)
y = f(x)

plt.plot(x,y)
plt.show()
