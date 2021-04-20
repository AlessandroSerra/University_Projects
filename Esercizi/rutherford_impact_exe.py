'''
programma che riproduce in simulazione l'esperimento di Rutherford sul bombardamento di un foglio d'oro da
parte di particelle alpha. si studia l'angolo di riflessione delle particelle alpha dopo l'eventuale impatto
con i nuclei d'oro.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import scipy.optimize as spo
import vec2d as v2d                #importiamo la classe vec2d creata in vec2d.py (deve essere nella stessa cartella)
import time 

t_s = time.time()

Nparticles = 1

## NOTE: numero atomico rispettivamente dell'elio e dell'oro
Ze, Zo = 2, 79
alpha_mass = 2 * spc.proton_mass + 2 * spc.neutron_mass               #massa di protone e neutrone

## NOTE: forza di cuolomb in forma scalare
Fc = Ze * Zo * spc.e**2 / (4 * np.pi * spc.epsilon_0 * alpha_mass)    #carica elettrone, costante dielettrica nel vuoto
E = 5e5 * spc.electron_volt                                           #energia di 5MeV convertita in Joule
dis = Fc / E * alpha_mass                                              #distanza minima tra particella sparata e nucleo
vel = np.sqrt(2 * E / alpha_mass)
tau = dis / vel

def F(r):
    mod = r.mod()
    return Fc / mod**3 * r      #F(r) restituisce la forza di coulomb vettoriale

def step(r, v, tau):                            #velocity-verlet
    rnew = r + v * tau + F(r) * tau**2 / 2
    vnew = v + (F(r) + F(rnew)) * tau / 2
    return rnew, vnew

## NOTE: funzione che esegue i vari step del velocity verlet
def solve(r0, v0, tau, Nsteps):
    t, r, v = [0], [r0], [v0]

    for i in range(Nsteps - 1):         #dato che abbiamo già inserito il primo elemento nelle liste
        rnew, vnew = step(r[i], v[i], tau)
        t.append(t[i] + tau)
        r.append(rnew)
        v.append(vnew)

    return t, r, v

r_list = []
v_list = []
imp_list = np.linspace(-100,100, Nparticles)

for i in range(Nparticles):
    r0 = v2d.vec2d(-100 * dis, imp_list[i])
    v0 = v2d.vec2d(vel, 0)
    r_list.append(r0)
    v_list.append(v0)                                              #all'inizio vogliamo moto solo lungo x
    Nsteps = 2 * int(r_list[i].mod()  / dis)

    t, r, v = solve(r_list[i], v_list[i], tau, Nsteps)


fig, ax = plt.subplots()
ax.plot([pos.x / dis for pos in r], [pos.y / dis for pos in r])

ax.plot(0, 0, marker = '.', color = 'k')                            #disegnamo l'atomo d'oro nell'origine
ax.set_xlabel('x (d)')
ax.set_ylabel('y (d)')

t_f = time.time() - t_s
print(t_f, 'secondi')

plt.legend()
plt.show()
