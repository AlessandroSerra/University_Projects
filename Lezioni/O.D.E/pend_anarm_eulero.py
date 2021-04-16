import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import g       #accelerazione di gravità

## NOTE: definiamo la ODE del pendolo
def f(theta):
    l = .1         #lunghezza della corda
    return - g/l * np.sin(theta)

## NOTE: definiamo l'algoritmo di eulero la ODE
def eulero(xold, vold, tau):
    xnew = xold + vold * tau
    vnew = vold + f(xold) * tau
    return xnew, vnew

##NOTE: definiamo eulero cromer per la ODE
def eulero_cromer(xold, vold, tau):
    vnew = vold + tau * f(xold)
    xnew = xold + tau * vnew
    return xnew, vnew

## NOTE: definiamo tau, x0, v0, N cicli e tempi n*tau
x0 = np.pi / 2 + np.pi / 6
v0 = 0
tau = 1e-3                          #time step
N = 2000                            #passi di integrazione
t = [tau * i for i in range(N+1)]   #modo veloce per riempire lista

## NOTE: definiamo le liste dei valori di x e v
x = [x0]
v = [v0]

"""
## NOTE: ciclo per eulero
for i in range(N):
    xnew, vnew = eulero(x[i], v[i], tau)
    x.append(xnew)
    v.append(vnew)

"""

##NOTE: ciclo per eulero cromer 
for i in range(N):
    xnew, vnew = eulero_cromer(x[i], v[i], tau)
    x.append(xnew)
    v.append(vnew)

## NOTE: rappresentiamo quello che abbiamo ottenuto con eulero
fig, ax = plt.subplots()
ax.plot(t, x, label = 'posizione')
ax.plot(t, v, label = 'velocità')
plt.legend()

plt.show()