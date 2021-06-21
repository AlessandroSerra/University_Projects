import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import time as t

t_s = t.time()      #tempo di inizio simulazione


##NOTE: funzione di inizializzazione del profilo del calore
def init_profile(x, N_x, N_t, init_T, init_A, wavelen):

    T_0 = np.zeros((N_t, N_x))                  #array delle temperature

    T_list = [init_T for i in range(len(x))]    #lista delle temperature iniziali dei punti
    T_0[0, :] = T_list - init_A * np.cos(2 * np.pi * x / wavelen)      #inizializzazione del profilo

    return T_0


##NOTE: definiamo una funzione che implementi le consizioni al contorno a seconda della condizione richiesta
def cond_cont(old_T, new_T, condition = 'dir'):     #di defoult restituisce le condizioni di dirichlet

    if condition == 'dir':      
        return old_T[0], old_T[-1]      #temp costante agli estremi

    elif condition == 'neu':
        return new_T[1], new_T[-2]      #temp ai punti 1 e penultimo per fissare il flusso di calore


##NOTE: definiamo la funzione che esegua uno step di integrazione
def step(old_T, N_x, k, h, dt, condizione = 'dir'):

    new_T = np.zeros(N_x)               #array delle nuove temperature

    for i in range(1, N_x - 1):         #primo e ultimo punto voncolati dalle cond al contorno

        new_T[i] = old_T[i] + dt * k / h**2 * (old_T[i+1] + old_T[i-1] - 2 * old_T[i])      #FTCS

    new_T[0], new_T[-1] = cond_cont(old_T, new_T, condizione)      #implementazione delle condizioni al contorno

    return new_T


##NOTE: definiamo la funzione che esegua tutti gli step necessari
def solve(N_t, N_x, T, k, h, dt, condizione = 'dir'):

    for n in range(1, N_t):         #la prima posizione contienenla temp all'istante iniziale t0

        T[n] = step(T[n-1], N_x, k, h, dt, condizione)

    return T

##NOTE: funzione per fittare l'ampiezza
def profile_width_fit(T, init_A, L, dt, N_t, wavelen):

    def width_func(t, A0, tau, c):  #funzione da fittare

        return A0 * np.exp(-t / tau) + c

    def diff_termica(tau):          #espressione di k (diffusività termica)

        return wavelen**2 / (4 * np.pi**2 * tau)

    A = np.amax(T, axis = 1)        #restituisce il valore max per ogni colonna
    time = np.array([dt * i for i in range(N_t)])   #array dei tempi

    p, cov = spo.curve_fit(width_func, time, A, p0 = [init_A, L**2 / time[-1], 300])
    y_fit = width_func(time, p[0], p[1], p[2])

    print('Tempo di rilassamento:\t', p[1], 's')
    print('Diffusività termica sperimentale:\t', diff_termica(p[1]), 'm^2/s')

    return y_fit, time, A


##NOTE: funzione di run del programma
def run_calor_profile(N_t, N_x, init_T, init_A, wavelen, dt, k, h):

    X = np.array([h * i for i in range(N_x)])             #lista dei valori della x
    T_0 = init_profile(X, N_x, N_t, init_T, init_A, wavelen)
    T = solve(N_t, N_x, T_0, k, h, dt, 'neu')

    y_fit, time, A = profile_width_fit(T, init_A, L, dt, N_t, wavelen)

    ##NOTE: rappresentazione grafica dei profili ti temperatura al variare del tempo
    fig, ax = plt.subplots()

    for n in range(0, N_t, int(N_t / 5)):       #loop per visualizzare solo 5 profili di temperatura al fronte dei 1000

        ax.plot(X, T[n], label = 't = %1.2f s' % (n * dt))   

    ax.plot(X, T[-1], label = 't = %1.2f s' % (N_t * dt))   #profilo di temperatura al tempo finale
    ax.legend(frameon = False)

    ##NOTE: rappresentazione grafica delle ampiezze in funzione del tempo e del loro fit
    fig1, ax1 = plt.subplots()
    ax1.plot(time, A, marker = '.', markersize = .1, linewidth = .1, label = 'Ampiezze in funzione del tempo')
    ax1.plot(time, y_fit, label = 'Fit delle ampiezze in funzione del tempo')
    ax1.legend(frameon = False)

    ##NOTE: rappresentazione 3D dell'andamento del profilo in funzione del tempo
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection = '3d')
    gridx , gridy = np.meshgrid(range(0, N_x), range(0, N_t))
    img = ax2.plot_surface(gridx, gridy, T, cmap = plt.cm.hot, linewidth = 0)
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('y (t)')
    ax2.set_zlabel('T (K)')
    ax2.set_title('Temperature 3D plot')
    ax2.legend(frameon = False)
    fig2.colorbar(img, orientation = 'vertical')

    ##NOTE: il tempo di simulazione non include l'intervallo necessario a visualizzare il grafico
    print('Time taken by the simulation:\t', t.time() - t_s, 'seconds')

    plt.show()


'''
Main del programma
'''

##NOTE: parametri di simulazione

init_T = 300        #K
init_A = 10         #K
L = 100             #m
wavelen = L     #m
k = 0.5             #m**2/s
h = 1               #m
dt = 0.5 * h**2 / (2 * k)      #s
N_x = int(L / h)                      #divisioni della sbarretta
N_t = L * 3 * 4**(4 - (wavelen + 25) / 25)

if wavelen == L:
    N_t = 3600

elif wavelen == L / 2:
    N_t = 1200

else:
    N_t = 300

##NOTE: funzione di run del programma
run_calor_profile(N_t, N_x, init_T, init_A, wavelen, dt, k, h)