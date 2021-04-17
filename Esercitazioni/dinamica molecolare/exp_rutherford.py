'''
programma che riproduce in simulazione l'esperimento di Rutherford sul bombardamento di un foglio d'oro da
parte di particelle alpha. si studia l'angolo di riflessione delle particelle alpha dopo l'eventuale impatto
con i nuclei d'oro.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import vec2d as v2d                 #importiamo la classe vec2d

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

## NOTE: funzione che esegue vari step del velocity verlet
def solve(r0, v0, tau, Nsteps):
    t, r, v = [0], [r0], [v0]

    for i in range(Nsteps - 1):         #dato che abbiamo già onserito il primo elemento nelle liste
        rnew, vnew = step(r[i], v[i], tau)
        t.append(t[i] + tau)
        r.append(rnew)
        v.append(vnew)

    return t, r, v

theta = []
Nparticles = 1000

for i in range(Nparticles):
    r0 = v2d.vec2d(-100 * dis, 100 *(2 * np.random.rand() - 1) * dis) #rand() restituisce un numero tra 0 ed 1 ma noi lo vogliamo tra -100 e 100 d
    v0 = v2d.vec2d(vel, 0)        #all'inizio vogliamo moto solo lungo x

    Nsteps = 2 * int(r0.mod()  / dis)

    t, r, v = solve(r0, v0, tau, Nsteps)
    theta.append(v0.get_angle(v[-1]))

## NOTE: vogliamo vedere come sono distribuiti gli angoli theta delle N particelle alpha
fig, ax = plt.subplots()
ax.hist(theta, histtype = 'step')       #creiamo un istogramma
ax.set_yscale('log')                    #usiamo scala log per vedere più suddivisioni
ax.set_xlabel('$\\theta$')
ax.set_ylabel('counts')

'''
fig, ax = plt.subplots()
ax.plot([pos.x / dis for pos in r], [pos.y / dis for pos in r])
ax.plot(0, 0, marker = '.', color = 'k')        #disegnamo l'atomo d'oro nell'origine
ax.set_xlabel('x (d)')
ax.set_ylabel('y (d)')
'''
#plt.legend()
plt.show()
