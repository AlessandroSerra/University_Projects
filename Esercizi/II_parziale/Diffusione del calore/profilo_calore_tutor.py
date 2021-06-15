'''
simuliamo la diffusione del calore in una sbarretta di lunghezza L = h * N_x = 1m con condizioni al contorno
di Dirichlet o di Newmann 
'''

import numpy as np
import matplotlib.pyplot as plt

##NOTE: funzione del calore al variare del tempo
def init_profile(x, N_x, N_t, init_T, init_A, wavelen):

    T_0 = np.zeros((N_t, N_x))

    T_list = [init_T for i in range(len(x))]
    T_0[0, :] = T_list - init_A * np.cos(2 * np.pi * x / wavelen)

    return T_0

##NOTE: definiamo una funzione che implementi le consizioni al contorno a seconda della condizione richiesta
def cond_cont(old_T, new_T, condition = 'dir'):     #di defoult restituisce le condizioni di dirichlet

    if condition == 'dir':      
        return old_T[0], old_T[-1]      #temp costante agli estremi

    elif condition == 'neu':
        return new_T[1], new_T[-2]      #temp ai punti 1 e penultimo peer fissare il flusso di calore

##NOTE: definiamo la funzione che esegua uno step di integrazione
def step(old_T, k, h, tau, condizione = 'dir'):

    N_x = len(old_T)
    new_T = np.zeros(N_x)

    for i in range(1, N_x - 1):         #primo e ultimo punto voncolati dalle cond al contorno

        new_T[i] = old_T[i] + tau * k / h**2 * (old_T[i+1] + old_T[i-1] - 2 * old_T[i])

    new_T[0], new_T[-1] = cond_cont(old_T, new_T, condizione)

    return new_T

##NOTE: definiamo la funzione che esegua tutti gli step necessari
def solve(N_t, T, k, h, tau, condizione = 'dir'):

    for n in range(1, N_t):         #la prima posizione contienenla temp all'istante iniziale t0

        T[n] = step(T[n-1], k, h, tau, condizione)

    return T


k = 1           #m**2/s
h = 1        #m
L = 100
wavelen = L
tau = 0.5 * h**2 / (2 * k)      #s, 0.5 per il termine sigma(t)
init_T = 300
init_A = 10
N_x = int(L / h)       #divisioni della sbarretta
N_t = 5000      #divisioni del tempo complessivo

T = np.zeros([N_t, N_x])        #array delle temperature
x = np.array([h * i for i in range(N_x)])             #lista dei valori della x

##NOTE: inizializziamo l'array delle temperature con le cond al contorno di Dirichlet o Neumann
T = init_profile(x, N_x, N_t, init_T, init_A, wavelen)

T = solve(N_t, T, k, h, tau, 'neu')     #risulviamo l'equazione differenziale del sistema

fig, ax = plt.subplots()

for n in range(0, N_t, int(N_t / 5)):       #loop per visualizzare solo 5 profili di temperatura al fronte dei 1000

    ax.plot(x, T[n], label = 't = %1.2f s' % (n * tau))   #passo un numero float con 1 cifra non decimale e 2 decimali

ax.plot(x, T[-1])

ax.set_xlabel('x (m)')
ax.set_ylabel('T (K)')

ax.legend(frameon = False)
plt.show()