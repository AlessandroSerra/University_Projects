import numpy as np
import matplotlib.pyplot as plt

## NOTE: definiamo la ODE del pendolo anarmonico
def f(theta):
    g = 9.81
    l = 0.1
    return - g/l * np.sin(theta)

## NOTE: algoritmo di eulero per la ODE
def eulero(xold, vold, fold, tau):
    xnew = xold + vold * tau
    vnew = vold + fold * tau
    return xnew, vnew

## NOTE: algoritmo di velocity_verlet per le ODE
def velocity_verlet(xold, vold, tau):
    xnew = xold + vold * tau + f(xold) * tau**2 / 2
    vnew = vold + tau/2 * (f(xold) + f(xnew))
    return xnew, vnew

def K(v):
    l = 0.1
    return 1/2 * v**2 * l**2

def U(x):
    l = 0.1
    g = 9.81
    return - g * l * np.cos(x)

tau = 1e-3
N = 2000
t = [tau * i for i in range(N+1)]

x0 = np.pi / 2 + np.pi / 6
v0 = 0

x = [x0]
v = [v0]

## NOTE: ciclo per eulero
#for i in range(N):
#    xnew, vnew = eulero(x[i], v[i], f(x[i]), tau)
#    x.append(xnew)
#    v.append(vnew)

## NOTE: ciclo per velocity_verlet
for i in range(N):
    xnew, vnew = velocity_verlet(x[i], v[i], tau)
    x.append(xnew)
    v.append(vnew)

## NOTE: valori di energia meccanica del sistema
cinetica = [K(v[i]) for i in range(len(x))]
potenziale = [U(x[i]) for i in range(len(x))]
Etot = [cinetica[i] + potenziale[i] for i in range(len(cinetica))]

fig, ax = plt.subplots()
ax.plot(t, x, label = 'posizione')
ax.plot(t, v, label = 'velocità')
plt.legend()

fig1, ax1 = plt.subplots()
ax1.plot(t, cinetica, label = 'E cinetica')
ax1.plot(t, potenziale, label = 'E potenziale')
ax1.plot(t, Etot, label = 'E meccanica')
plt.legend()

## NOTE: per due figure distinte va usato un solo plt.show, sennò non apre la seconda
plt.show()
