import numpy as np
import matplotlib.pyplot as plt
import scipy.constants  as spc
#from sklearn.metrics import mean_squared_error
import random
import math
import sympy as sp
from sympy.core.evalf import N
from sympy.series.limits import limit

mass_of_argon = 39.948 # Massa dell'argon in unità di massa atomiche
Nsteps=200
Nparts=3 #numero di atomi di Argon
L=20 #amstrong dimensioni box
tau = 0.1 #passo temporale
T=300 #K kelvin
sigma=3.46 #Angstrong 10 alla - 10 metri Distanza in cui la potenziale energia è zero
epsilon = 0.0103 #eV è l'energia potenziale alla separazione di equilibrio
cutoff=20 #se due particelle sono molto distanti tra loro, non è necessario calcolare lj


def lj_force(r, epsilon, sigma):
    """Definisco una funzione per calcolare la forza dal potenziale di Lennard-Jones. Questo dipende
    dalla distanza tra le due particelle. Restituisce il valore della forza derivante dal potenziale
    ovvero la forza di interazione tra due particelle di argon,
    Con il cut-off se due particelle sono molto distanti tra loro, non è necessario 
    calcolare l'energia / forza. Invece il valore viene semplicemente preso come 0"""
    
    if r < cutoff:
        return 48 * epsilon * np.power(
            sigma / r, 13) - 24 * epsilon * np.power(
            sigma / r, 7)
    else:
        return 0        
    

def init_velocity(T, number_of_particles):
    """La seguente funzione inizializza le velocità, prende come argomenti la temperatura T ,
    il Numero di particelle e gli step e restituisce un array di velocità pari al numero di particelle
    si calcola la velocità media delle particelle in equilibrio termico a temperatura T. Poichè l'energia
    media a temperatura T è Kb*T calola la radice quadrata dell'energia divisa la massa ottenendo così
    una stima delle velocità medie delle particelle inoltre effettua la caonversione da J a eV """      
    R = np.random.rand(number_of_particles) - 0.5
    return R * np.sqrt(spc.Boltzmann * T / (
        mass_of_argon * 1.602e-19))

def get_accelerations(positions):
    """
    Calcola l'accelerazione su ogni particella
     come risultato l'una dell'altra particella."""
   
    accel_x = np.zeros((positions.size, positions.size))
    for i in range(0, positions.size - 1):
        for j in range(i + 1, positions.size):
            r_x = positions[j] - positions[i]
            rmag = np.sqrt(r_x * r_x)
            force_scalar = lj_force(rmag, epsilon, sigma)
            force_x = force_scalar * r_x / rmag
            accel_x[i, j] = force_x / mass_of_argon
            accel_x[j, i] = - force_x / mass_of_argon
    return np.sum(accel_x, axis=0)

def update_pos(x, v, a, dt):
    """Queste due funzioni aggiornano le posizione e le velocità in
    ingresso prende le vecchie posizioni, le vecchie velocità, le forze
    e l'intervallo di tempo"""
    
    return x + v * dt + 0.5 * a * dt * dt

def update_velo(v, a, a1, dt):

    return v + 0.5 * (a + a1) * dt

def get_mean_quad_dev(positions):

    mean_quad_dev = []
    quad_dev = 0

    for i in range(Nparts):
            for k in range(Nsteps):

                quad_dev += (positions[k, i]- positions[k, 0])**2
                
            mean_quad_dev.append(quad_dev)

    return mean_quad_dev

def get_selfdiff_par(mean_quad_dev):

    x = sp.Symbol('x')
    a = []

    for i in range(Nparts):

        a.append(sp.limit((mean_quad_dev / (2*x)), x, sp.oo))

    return a

def run_md(dt, Nparts, number_of_steps, initial_temp, x):

    """In questa parte del programma si esegue la simulazione molecolare e infine
    restituisce gli array delle posizioni delle paticelle in input prende il passo
    temporale tau, il numero di steps, una temperatura iniziale """
    positions = np.zeros((number_of_steps, Nparts))
    v = init_velocity(initial_temp, Nparts)
    a = get_accelerations(x)
    for i in range(number_of_steps):
        x = update_pos(x, v, a, dt)
        a1 = get_accelerations(x)
        v = update_velo(v, a, a1, dt)
        a = np.array(a1)
        positions[i, :] = x

    mean_quad_dev = get_mean_quad_dev(positions)
    selfdiff_par = get_selfdiff_par(Nparts)
    print(mean_quad_dev)
    print(selfdiff_par)

    return positions

x = np.linspace(0, L, Nparts)
sim_pos = run_md(tau, Nparts, Nsteps, T, x)

atoms=[Nparts]

'''
def distribuzione_Maxwell_Boltzman(mass_of_argon,Nparts):
        """Applica la distribuzione di Boltzmann alle velocità atomiche"""
        normDist = []
        scaling_factor = np.sqrt(spc.Boltzmann*T/mass_of_argon)

        # Stabilisce una distribuzione normale
        for i in range(0, 3*Nparts):
            normDist.append(random.gauss(0,1))
      
        # Applica il fattore di scala alla distribuzione
        for number in range(0, 3*Nparts):
            normDist[number] = normDist[number]*scaling_factor
            
        # Distribuisce le velocità
        for atom in range(0, Nparts):
            vx = normDist[atom*3]
            #atoms[atom].vy = normDist[atom*3+1]
        plt.plot(normDist,linestyle='', marker=',',markerfacecolor="blue" )

#max=distribuzione_Maxwell_Boltzman(mass_of_argon,Nparts)
'''


for i in range(sim_pos.shape[1]):
    plt.plot(sim_pos[:, i], '.', label='atom {}'.format(i))
plt.xlabel(r'Step')
plt.ylabel(r'$x$-Position/Å')
plt.legend(frameon=False)
plt.show()