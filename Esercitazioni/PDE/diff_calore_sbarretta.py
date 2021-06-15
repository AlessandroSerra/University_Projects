'''
simuliamo la diffusione del calore in una sbarretta di lunghezza L = h * N_x = 1m con condizioni al contorno
di Dirichlet o di Newmann 
'''

import numpy as np
import matplotlib.pyplot as plt

##NOTE: definiamo una funzione che implementi le consizioni al contorno a seconda della condizione richiesta
def cond_cont(old_T, new_T, condition = 'dir'):     #di defoult restituisce le condizioni di dirichlet

    if condition == 'dir':      
        return old_T[0], old_T[-1]      #temp costante agli estremi

    elif condition == 'neu':
        return new_T[1], new_T[-2]      #temp ai punti 1 e penultimo peer fissare il flusso di calore

##NOTE: definiamo la funzione che esegua uno step di integrazione
def step(old_T, k, h, dt, condizione = 'dir'):

    N_x = len(old_T)
    new_T = np.zeros(N_x)

    for i in range(1, N_x - 1):         #primo e ultimo punto voncolati dalle cond al contorno

        new_T[i] = old_T[i] + dt * k / h**2 * (old_T[i+1] + old_T[i-1] - 2 * old_T[i])

    new_T[0], new_T[-1] = cond_cont(old_T, new_T, condizione)

    return new_T

##NOTE: definiamo la funzione che esegua tutti gli step necessari
def solve(N_t, T, k, h, tau, condizione = 'dir'):

    for n in range(1, N_t):         #la prima posizione contienenla temp all'istante iniziale t0

        T[n] = step(T[n-1], k, h, tau, condizione)

    return T


N_x = 100       #divisioni della sbarretta
N_t = 10000      #divisioni del tempo complessivo
k = 5000           #m**2/s
h = 0.01        #m
dt = 0.5 * h**2 / (2 * k)      #s, 0.5 per il termine sigma(t)

T = np.zeros([N_t, N_x])        #array delle temperature

##NOTE: inizializziamo l'array delle temperature con le cond al contorno di Dirichlet o Neumann
T = T + 300         #tutti gli elementi di T sono a 300K temperatura ambiente
T[0, 0] = 330       #K, temperatura a t=0 e a x=0
T[0, -1] = 300      #K, temperatura a t=0 e all'ultima posizione

T = solve(N_t, T, k, h, dt, 'dir')     #risulviamo l'equazione differenziale del sistema

fig, ax = plt.subplots()
x = [h * i for i in range(N_x)]             #lista dei valori della x

for n in range(0, N_t, int(N_t / 5)):       #loop per visualizzare solo 5 profili di temperatura al fronte dei 1000

    ax.plot(x, T[n], label = 't = %1.2f s' % (n * dt))   #passo un numero float con 1 cifra non decimale e 2 decimali

ax.set_xlabel('x (m)')
ax.set_ylabel('T (K)')

ax.legend(frameon = False)
plt.show()