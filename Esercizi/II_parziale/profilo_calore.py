import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time as t

t_s = t.time()


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
def step(old_T, N_x, k, h, tau, condizione = 'dir'):

    new_T = np.zeros(N_x)

    for i in range(1, N_x - 1):         #primo e ultimo punto voncolati dalle cond al contorno

        new_T[i] = old_T[i] + tau * k / h**2 * (old_T[i+1] + old_T[i-1] - 2 * old_T[i])

    new_T[0], new_T[-1] = cond_cont(old_T, new_T, condizione)

    return new_T


##NOTE: definiamo la funzione che esegua tutti gli step necessari
def solve(N_t, N_x, T, k, h, tau, condizione = 'dir'):

    for n in range(1, N_t):         #la prima posizione contienenla temp all'istante iniziale t0

        T[n] = step(T[n-1], N_x, k, h, tau, condizione)

    return T


init_T = 300        #K
init_A = 10         #K
L = 100             #m
wavelen = L         #m
k = .5              #m**2/s
h = 1               #m
tau = 0.9 * h**2 / (2 * k)      #s
N_x = 1000                      #divisioni della sbarretta
N_t = 10000                     #divisioni temporali

X = np.linspace(0, L, N_x)
T_0 = init_profile(X, N_x, N_t, init_T, init_A, wavelen)
T = solve(N_t, N_x, T_0, k, h, tau, 'neu')

fig, ax = plt.subplots()

'''
for n in range(0, N_t, int(N_t / 5)):       #loop per visualizzare solo 5 profili di temperatura al fronte dei 1000

    ax.plot(X, T[n], label = 't = %1.2f s' % (n * tau))   #passo un numero float con 1 cifra non decimale e 2 decimali
'''


fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
gridx , gridy = np.meshgrid(range(N_t),range(N_x)) 
ax.plot_surface(gridx,gridy,T,cmap=plt.cm.coolwarm, vmax=250,linewidth=0,rstride=2, cstride=100)

ax.set_xlabel('x (m)')
ax.set_ylabel('T (K)')

ax.legend(frameon = False)
plt.show()