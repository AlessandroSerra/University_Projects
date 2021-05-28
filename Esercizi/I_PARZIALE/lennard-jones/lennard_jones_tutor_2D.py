import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as spc 
import vec2d as v2d


"COSTANTI"
mass_argon = 39.984 #uma qui scriveremo in termini di uma, angstrog e eV
epsilon = 0.0103 #eV costante di accoppiamento per l'argon
sigma = 3.46 #distanza entro cui l'interazione è più efficace

dt = 0.01 #passo temporale necessario per farlo runnare
N_steps = 10000
T = 300 #K
L = 25 * sigma #angstrom

def lj_force(r, L): #definizione del potenziale di Lennar-Jones
    d = r.mod()
    #controllare se le particelle sono a distanza max 3 o 5 sigma
    return 4 * epsilon * (12 * sigma**12 / d**13 - 6 * sigma**6 / d**7) * r.unitary()
    
def init_vel(T, N_partcl):
    np.random.seed()
    R1 = 2 * (np.random.rand(N_partcl) - 1) #in questo modo definiamo un numero random compreso fra -1 e 1
    R2 = 2 * (np.random.rand(N_partcl) - 1)

    const = np.sqrt(spc.Boltzmann * T / (mass_argon * spc.electron_volt))
    vel_list = [const * v2d.vec2d(R1[i], R2[i]) for i in range(N_partcl)]
    return np.array(vel_list)


def compute_forces(positions, L):
    N_partcl = len(positions)
    forces = np.full([N_partcl, N_partcl], fill_value = v2d.vec2d(0, 0)) #matrice quadrata
    for i in range(N_partcl - 1): #in questo modo evito che l'ultima partcl sia conteggiata due volte
        for j in range(i + 1, N_partcl):
            d = positions[j] - positions[i]
            force = lj_force(d, L) / mass_argon
            forces[i, j] = force #i=righe; j=colone !!!!!!!!!!
            forces[j, i] = - force
    return np.sum(forces, axis = 0) #axis = 0 vuol dire che deve sommare tutti gli elementi che hanno gli stessi indici di colonna (con axis=2 vuol dire sommare gli stessi sulla riga)

def update_positions(positions, velocities, forces, dt, L):
    new_pos = (positions + velocities * dt + forces * dt**2 / 2) % L
    return new_pos

def update_velocities(velocities, forces, new_forces, dt):
    return velocities + (forces + new_forces) * dt / 2

def run_md(dt, N_steps, T, L): #questa sarà la funzione che girerà. dt è il passo temporale. T è temperatura, initial_pos sono le posizioni iniziali.
    N_partcl = 25
    initial_pos = np.array([v2d.vec2d(i + 10, j + 10) * sigma * 1.1 for i in range(int(np.sqrt(N_partcl))) for j in range(int(np.sqrt(N_partcl)))])
    positions = np.full([N_steps, N_partcl], fill_value = v2d.vec2d(0, 0))
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


positions = run_md(dt, N_steps, T, L)

fig, ax = plt.subplots()
for i in range(positions.shape[1]): #per un array, len restituisce il numero di array (n° righe)
    pos_x, pos_y = [positions[j, i].x for j in range(N_steps)], [positions[j, i].y for j in range(N_steps)]
    ax.plot(pos_x, pos_y)
ax.set_xlabel("Step")
ax.set_ylabel("Positions (Angstrom)")
ax.legend()

plt.show()


'''
le velocità sono in angstrom al nanosecondo
quando fittiamo maxwell boltzman dobbiamo mettere la velocita finale meno quella 
iniziale perche senno abbiamo una distribuzione che non parte da 0
'''