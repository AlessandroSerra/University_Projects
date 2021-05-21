import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as spc 


"COSTANTI"
mass_argon = 39.984 #uma #qui scriveremo in termini di uma, angstrog e eV
epsilon = 0.0103 #eV #costante di accoppiamento per l'argon
sigma = 3.46 #distanza entro cui l'interazione è più efficace

dt = 0.01 #passo temporale necessario per farlo runnare
N_steps = 100000
T = 300 #K
x = np.array([1, 5, 10])
L = 15 * sigma #angstrom

def lj_force(d, L): #definizione del potenziale di Lennar-Jones
    d = d % (np.sqrt(2) * L / 2)
    #controllare se le particelle sono a distanza max 3 o 5 sigma
    if d <= 5 * sigma:
        return 4 * epsilon * (12 * sigma**12 / d**13 - 6 * sigma**6 / d**7)
    return 0 #posso scrivere così anche senza else, perché python legge fintato che le condizioni non sono soddisfatte!

def init_vel(T, N_partcl):
    R = 2 * np.random.rand() - 1 #in questo modo definiamo un numero random compreso fra -1 e 1
    return R * np.sqrt(spc.Boltzmann * T / (mass_argon * spc.electron_volt))

def compute_forces(positions, L):
    N_partcl = len(positions)
    forces = np.zeros([N_partcl, N_partcl]) #matrice quadrata
    for i in range(N_partcl - 1): #in questo modo evito che l'ultima partcl sia conteggiata due volte
        for j in range(i + 1, N_partcl):
            d = positions[j] - positions[i]
            force = lj_force(d, L) / mass_argon
            forces[i, j] = force #i=riche; j=colone !!!!!!!!!!
            forces[j, i] = - force
    return np.sum(forces, axis = 0) #axis = 0 vuol dire che deve sommare tutti gli elementi che hanno gli stessi indici di colonna (con axis=2 vuol dire sommare gli stessi sulla riga)

def update_positions(positions, velocities, forces, dt, L):
    new_pos = positions + velocities * dt + forces * dt**2 / 2
    return new_pos % L

def update_velocities(velocities, forces, new_forces, dt):
    return velocities + (forces + new_forces) * dt / 2

def run_md(dt, N_steps, T, initial_pos, L): #questa sarà la funzione che girerà. dt è il passo temporale. T è temperatura, initial_pos sono le posizioni iniziali.
    N_partcl = len(initial_pos)
    positions = np.zeros([N_steps, N_partcl])
    positions[0, :] = initial_pos #questa dicitura indica di prendere la prima colonna e considerarne tutti gli elementi. 
    x = initial_pos
    vel = init_vel(T, N_partcl)
    forces = compute_forces(initial_pos, L)
    for i in range(N_steps):
        #aggiornare posizioni
        x = update_positions(x, vel, forces, dt, L)
        #aggiornare forze con nuove posizioni
        new_forces = compute_forces(x, L) #lo chiamo newforces per evitare sovrascrizioni
        #calcolare velocità
        vel = update_velocities(vel, forces, new_forces, dt)
        #aggiorno array di posizione
        positions[i, :] = x
    return positions 
    #semplicemente usando la stessa sintassi si può appendere e scrivere anche la velocità 


positions = run_md(dt, N_steps, T, x, L)

fig, ax = plt.subplots()
for i in range(len(x)): #per un array, len restituisce il numero di array (n° righe)
    ax.plot(positions[:, i], label = "Atomo %g" % i)
ax.set_xlabel("Step")
ax.set_ylabel("Positions (Angstrom)")
ax.legend()

plt.show()