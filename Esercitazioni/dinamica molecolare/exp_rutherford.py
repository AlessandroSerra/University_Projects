'''
programma che riproduce in simulazione l'esperimento di Rutherford sul bombardamento di un foglio d'oro da 
parte di particelle alpha. si studia l'angolo di riflessione delle particelle alpha dopo l'eventuale impatto
con i nuclei d'oro.
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import vec2d as v                   #importiamo la classe vec2d

## NOTE: numero atomico rispettivamente dell'elio e dell'oro
Ze, Zo = 2, 79
alpha_mass = 2 * spc.proton_mass + 2 * spc.neutron_mass               #massa di protone e neutrone

## NOTE: forza di cuolomb in forma scalare 
Fc = Ze * Zo * spc.e**2 / (4 * np.pi * spc.epsilon_0 * alpha_mass)    #carica elettrone, costante dielettrica nel vuoto
E = 5e5 * spc.electron_volt                                           #energia di 5MeV convertita in Joule
dis = F0 / E * alpha_mass                                              #distanza minima tra particella sparata e nucleo 
vel = np.sqrt(2 * E / alpha_mass)
tau = dis / vel

def F(r):
    mod = r.mod()
    return F0 / mod**3 * r      #F(r) restituisce la forza di coulomb vettoriale

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

r0 = v.vec2d(-100 * dis, 100 *(2 * np.random.rand() - 1) * dis) #rand() restituisce un numero tra 0 ed 1 ma noi lo vogliamo tra -100 e 100 d
v0 = v.vec2d(vel, 0)        #all'inizio vogliamo moto solo lungo x

Nsteps = 2 * int(r0  / dis)