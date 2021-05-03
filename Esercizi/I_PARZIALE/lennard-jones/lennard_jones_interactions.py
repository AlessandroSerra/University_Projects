'''
codice che simula la traiettoria di tre atomi di argon in funzione del tempo di simulazione e del time step scelto
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import *
mass_of_argon = 39.948 # amu

def lj_force(r, box_L, epsilon, sigma):

    r = r % (np.sqrt(2) * box_L / 2)

    if r <= 5 * sigma:

        return 48 * epsilon * np.power(
            sigma, 12) / np.power(
            r, 13) - 24 * epsilon * np.power(
            sigma, 6) / np.power(r, 7)

    return 0

def init_velocity(T, number_of_particles):

    R = np.random.rand(number_of_particles) - 0.5
    return R * np.sqrt(Boltzmann * T / (
        mass_of_argon * 1.602e-19))

def get_accelerations(positions, box_L):

    accel_x = np.zeros((positions.size, positions.size))
    for i in range(0, positions.size - 1):
        for j in range(i + 1, positions.size):
            r_x = positions[j] - positions[i]
            rmag = np.sqrt(r_x * r_x)
            force_scalar = lj_force(r_x, box_L, 0.0103, 3.4)
            force_x = force_scalar * r_x / rmag
            accel_x[i, j] = force_x / mass_of_argon
            accel_x[j, i] = - force_x / mass_of_argon

#np.sum resituisce la somma di ogni riga 
    return np.sum(accel_x, axis=0)

def update_pos(x, v, a, box_L, dt):

    new_pos = x + v * dt + 0.5 * a * dt * dt
    return new_pos % box_L

def update_velo(v, a, a1, dt):

    return v + 0.5 * (a + a1) * dt

def run_md(dt, number_of_steps, initial_temp, box_L, x):

    positions = np.zeros((number_of_steps, 4))
    v = init_velocity(initial_temp, 4)
    print(v, 'm/s')
    a = get_accelerations(x, box_L)
    for i in range(number_of_steps):
        x = update_pos(x, v, a, box_L, dt)
        a1 = get_accelerations(x, box_L)
        v = update_velo(v, a, a1, dt)
        a = np.array(a1)
        positions[i, :] = x

    return positions

box_L = 20
x = np.array([1, 5, 10, 15])
sim_pos = run_md(0.1, 10000, 300, box_L, x)

    
for i in range(sim_pos.shape[1]):
    plt.plot(sim_pos[:, i], '.', label='atom {}'.format(i))
plt.xlabel(r'Step')
plt.ylabel(r'$x$-Position/Å')
plt.legend(frameon=False)
plt.show()