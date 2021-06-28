import numpy as np
import matplotlib.pyplot as plt

N_x = 100
wavelen1 = 100
wavelen2 = 50
wavelen3 = 25

X = np.array([1 * i for i in range(N_x)])

T_0 = np.zeros((N_x))                  #array delle temperature

T_list = [300 for i in range(len(X))]    #lista delle temperature iniziali dei punti
T_1 = T_list - 10 * np.cos(2 * np.pi * X / wavelen1)      #inizializzazione del profilo
T_2 = T_list - 10 * np.cos(2 * np.pi * X / wavelen2)      #inizializzazione del profilo
T_3 = T_list - 10 * np.cos(2 * np.pi * X / wavelen3)      #inizializzazione del profilo

fig, ax = plt.subplots()
ax. plot([k for k in range(N_x)], [300 for i in range(N_x)], linewidth = .5, color = 'black')
ax.plot(X, T_1, linewidth = .7, label = 'Profilo iniziale $\\lambda = L$')
ax.plot(X, T_2, linewidth = .7, label = 'Profilo iniziale $\\lambda = L/2$')
ax.plot(X, T_3, linewidth = .7, label = 'Profilo iniziale $\\lambda = L/4$')
ax.set_title('Profili iniziali di temperatura al variare di $\\lambda$')
ax.set_xlabel('x (m)')
ax.set_ylabel('T(K)')
ax.legend(frameon = False)

plt.show()