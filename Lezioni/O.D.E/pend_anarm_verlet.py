import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import g

## NOTE: definiamo la ODE del pendolo
def f(theta):
    l = .1           #lunghezza della corda
    return -g/l * np.sin(theta)

"""
NOTE: definiamo l'algoritmo di verlet, è necessario svolgere 
il primo passo con eulero o un altro algoritmo
"""
def eulero(xold, vold, tau):
    xnew = xold + vold * tau
    vnew = vold + f(xold) * tau
    return xnew, vnew

def verlet(x1, xold, vold, tau):
    xnew = 2 * xold + x1 + tau**2 * f(xnew)
    vnew = (xnew - xold) / (2 * tau)
    return xnew, vnew

## NOTE: definiamo l'algritmo di velocity-verlet per la ODE
def velocity_verlet(xold, vold, tau):
    xnew = xold + tau * vold + f(xold) * tau**2 / 2
    vnew = vold + (f(xold) + f(xnew)) * tau / 2
    return xnew, vnew

# NOTE: definiamo tau, x0, v0, N cicli e tempi n*tau
x0 = np.pi / 2 + np.pi / 6
v0 = 0
tau = 1e-3                          #time step
N = 2000                            #passi di integrazione
t = [tau * i for i in range(N+1)]   #modo veloce per riempire lista

## NOTE: definiamo le liste dei valori di x e v
x = [x0]
v = [v0]

##NOTE: ciclo per velocity verlet

for i in range(N):
    xnew, vnew = velocity_verlet(x[i], v[i], tau)
    x.append(xnew)
    v.append(vnew)

"""
##NOTE: ciclo per verlet
x1 = x0 + v0 * tau
x2, v2 = verlet(x1, x0, v0, tau)
x.append(x1, x2)
v.append(v0, v2)

for i in range(N-1):
    xnew = verlet(x1, )

"""

## NOTE: rappresentiamo quello che abbiamo ottenuto con eulero
fig, ax = plt.subplots()
ax.plot(t, x, label = 'posizione')
ax.plot(t, v, label = 'velocità')
plt.legend()

plt.show()