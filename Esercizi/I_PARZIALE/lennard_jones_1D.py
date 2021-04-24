import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc

## NOTE: costanti della simulazione
mass_of_argon = 39.948      #amu
epsilon = .0103             
sigma = 3.4                 #Angstron
L_box = 20                  #Angstron

## NOTE: variabili della simulazione
#N_particles = int(input('Inserire il numero di particelle da far interagire\n'))
N_particles = 30
N_steps = 100
tau = .1                    #femtosecondi
initial_T = 300             #Kelvin

## NOTE: definiamo la forza di interazione come la derivata parziale opposta del potenziale
def lj_force(r, epsilon, sigma):

    return 48 * epsilon * np.power(
        sigma, 12) / np.power(
        r, 13) - 24 * epsilon * np.power(
        sigma, 6) / np.power(r, 7)

## NOTE: inizializziamo le velocità in funzione della temperatura del sistema
def init_velocity(initial_T, N_particles):

    P = np.random.rand(N_particles) - 0.5
    return P * np.sqrt(spc.Boltzmann * initial_T
        / (mass_of_argon * spc.e))

## NOTE: ricaviamo l'acelerazione dalla forza di interazione
def get_accelerations(positions):

    accel_x = np.zeros((N_particles, N_particles))

    for i in range(0, N_particles - 1):
        for k in range(i + 1, N_particles):
            r_x = positions[k] - positions[i]
            rmag = np.sqrt(r_x**2)
            force_scalar = lj_force(r_x, epsilon, sigma)
            force_vec = force_scalar * r_x / rmag
            accel_x[i, k] = force_vec / mass_of_argon
            accel_x[k, i] = - accel_x[i, k]

    return np.sum(accel_x, axis = 0)

## NOTE: usiamo velocity-verlet per integrare posizione e velocità delle particelle
def update_pos(x, v, a, tau):

    return x + v * tau + 0.5 * a * tau**2

def update_vel(v, a, anew, tau):

    return v + 0.5 * (a + anew) * tau

## NOTE: funzione che fa girare velocity-verlet
def run_md(tau, N_steps, initial_T, x):

    positions = np.zeros((N_steps, N_particles))
    v = init_velocity(initial_T, N_particles)
    a = get_accelerations(x)
    for i in range(N_steps):
        x = update_pos(x, v, a, tau)
        anew = get_accelerations(x)
        v = update_vel(v, a, anew, tau)
        a = np.array(anew)
        positions[i, :] = x

    return positions

x = np.array([5 * i for i in range(N_particles)])
sim_pos = run_md(tau, N_steps, initial_T, x)

fig, ax = plt.subplots()

## NOTE: ciclo che plotta le posizioni delle particelle
for i in range(sim_pos.shape[1]):
    plt.plot(sim_pos[:, i], '.', label='atom {}'.format(i + 1))

ax.set_xlabel(r'Step')
ax.set_ylabel(r'$x$-Position/Å')
plt.legend(frameon=False)
plt.show()