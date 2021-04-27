import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import vec2d as v2d
import time 

t_s = time.time()

def acc_coulomb(r):
    
    return acc_scalar / r.mod()**3 * r

def update_pos(r, v, tau):

    r_new = r + v * tau + acc_coulomb(r) * tau**2 / 2
    return r_new

def update_vel(r, r_new, v, tau):

    v_new = v + (acc_coulomb(r) + acc_coulomb(r_new)) * tau / 2
    return v_new

def run_rtfrd(r0, v0, N_steps, tau, k):

    '''
    t_list = [tau]
    r_list = [r0]
    v_list = [v0]
    '''

    r_new = update_pos(r0, v0, tau)
    v_new = update_vel(r0, r_new, v0, tau)

    return r_new, v_new

##NOTE: costanti di simulazione
Ze, Zo = 2, 79  #numero atomico elio e oro
energy = 5e5 * spc.electron_volt    #5e5 eV convertiti in joule
alpha_mass = 2 * spc.proton_mass + 2 * spc.neutron_mass

acc_scalar = Ze * Zo * spc.e**2 / (4 * np.pi * spc.epsilon_0 * alpha_mass)
dis = acc_scalar / energy * alpha_mass
vel = np.sqrt(2 * energy / alpha_mass)
tau = dis / vel

r0 = v2d.vec2d(-100 * dis, 100 * (2 * np.random.rand() - 1) * dis)
v0 = v2d.vec2d(vel, 0)
N_steps = 2 * int(r0.mod()/dis)
N_particles = 2
particles_array = np.zeros((N_particles, N_steps))
vel_array 

for i in range(N_particles):
    
    for k in range(1, N_steps):

        pos_array[i][k], vel_array[i][k] = run_rtfrd(r0, v0, N_steps, tau, k)
        r0 = pos_array[i][k]
        v0 = vel_array[i][k]

print(particles_array)

