import numpy as np
import matplotlib.pyplot as plt

## NOTE: equazione differenziale del pendolo anarmonico, x = theta

def f(x):
    g = 9.81
    l = 0.1         #note: lunghezza del braccio del pendolo
    return - g/l * np.sin(x)

## NOTE: algoritmo di eulero per una ODE, tau rappresenta l'incremento
# della x (tempo), definiamo f_n+1 
def eulero(xold, vold, fold, tau):
    xnew = xold + vold * tau
    vnew = vold + f(xold) * tau
    return xnew, vnew

## NOTE: defiiamo un intervallo temporale sufficientemente piccolo
tau = 1e-3

## NOTE: definiamo il numero di cicli da fare
N = 2000

## NOTE: definiamo una lista che contenga tutti i valori della x (t)
t = [tau * i for i in range(N+1)]

## NOTE: definiamo velocità e posizione inizieli
x0 = np.pi / 2 + np.pi / 6      #random
v0 = 0

## NOTE: definiamo una lista che contenga i vaori di velocità e posizione
x = [x0]
v = [v0]

## NOTE: ciclo per il metodo di eulero
for i in range(N):
    xnew, vnew = eulero(x[i], v[i], f(x[i]), tau)
    x.append(xnew)
    v.append(vnew)

fig, ax = plt.subplots()
ax.plot(t,x, label='posizione')
ax.plot(t,v, label='velocità')
plt.legend()
plt.show()