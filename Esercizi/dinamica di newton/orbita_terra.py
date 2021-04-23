import numpy as np
import matplotlib.pyplot as plt

## NOTE: definiamo l'equazione per ottenere il raggio dell'orbita
def r(x, y):
    return np.sqrt(x**2 + y**2)

## NOTE: definiamo la ODE dell'orbita della terra 
def fx(x, y):
    return GMs * x / r(x,y)**3

def fy(x, y):
    return GMs * y / r(x,y)**3

## NOTE: useremo l'algoritmo velocity-verlet per calcolare le grandezze 
def velocity_verlet(xold, yold, vxold, vyold, tau):
    xnew = xold + vxold * tau + fx(xold, yold) * tau**2 / 2
    ynew = yold + vyold * tau + fy(xold, yold) * tau**2 / 2
    vxnew = vxold + tau / 2 * (fx(xold, yold) + fx(xnew, ynew))
    vynew = vyold + tau / 2 * (fy(xold, yold) + fy(xnew, ynew))
    return xnew, ynew, vxnew, vynew


## NOTE: costanti del moto conosciute
GMs = 4 * np.pi**2
Mterra = 5.9722e24

## NOTE: posizioni iniziali della terra (arbitrarie), le vel iniziali e lo step e tau
x0 = 0.98           #perielio
y0 = 0
vx0 = 0
vy0 = 2*np.pi
N = 1000
tau = 1e-3

x_pos = [x0]
y_pos = [y0]
x_vel = [vx0]
y_vel = [vy0]

## NOTE: ciclo per velocity-verlet
for i in range(N):
    x, y, vx, vy = velocity_verlet(x_pos[i], y_pos[i], x_vel[i], y_vel[i], tau)
    x_pos.append(x)
    y_pos.append(y)
    x_vel.append(vx)
    y_vel.append(vy)

## NOTE: plottiamo il risultato
fig, ax = plt.subplots()
ax.plot(x_pos, y_pos, label = 'orbita terra')
ax.set_xlabel('x (AU)')
ax.set_ylabel('y (AU)')
ax.axis('equal')
plt.legend()

plt.show()