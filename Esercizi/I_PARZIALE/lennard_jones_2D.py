import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import vec2d as v2d
import time 

t_s = time.time()

## NOTE: funzione che inizializza le posizioni degli atomi
def init_pos(N_particles, N_steps):

    pos_array = np.array([N_particles], [N_steps])

    ## NOTE: loop che assegna la posizione iniziale a tutte le particelle
    for i in range(N_particles):
        
        pos_array[0, i] = v2d.vec2d(1, 20 * (2 * np.random.rand() - 1))
    
    return pos_array

def init_vel(N_particles, N_steps, initial_T):

    vel_array = np.zeros((N_steps, N_particles))

    ## NOTE: loop che assegna la vlocità iniziale a tutte le particelle    
    for i in range(N_particles):

        P = np,random.rand(N_particles) - 0.5
        R = np.random.rand(N_particles) - 0.5
        vel_array[0, i] = (v2d.vec2d(P * np.sqrt(spc.Boltzmann * initial_T / 
                (mass_of_argon / spc.e)), R * np.sqrt(spc.Boltzmann * initial_T / 
                (mass_of_argon / spc.e))))
    
    return vel_array

def get_radius(pos_array, i):

    dist_array = np.zeros((N_particles, N_particles))

    for j in range(0, N_particles - 1):
            for k in range(j + 1, N_particles):

                dist_array[j, k] = pos_array[i, j] - pos_array[i, k]

    return dist_array

# NOTE: definiamo la forza di interazione come la derivata parziale opposta del potenziale
def lj_acceleration(dist_array, epsilon, sigma, i):

    acc_array = np.zeros(N_particles, N_particles)
    force_array = np.zeros(N_particles, N_particles)

    for j in range(0, N_particles - 1):     #dist_array contiene il vettore differenza tra i vettori posizione
        for k in range(j + 1, N_particles):
            acc_array[j, k] = (48 * epsilon * np.power(
                sigma, 12) / np.power(
                dist_array[j, k].mod(), 13) - 24 * epsilon * np.power(
                sigma, 6) / np.power(dist_array[j, k].mod(), 7)
                ) / mass_of_argon * dist_array[j, k] / dist_array[j, k].mod()

            acc_array[k, j] = - acc_array[j, k]     #terza legge di newton

    return acc_array   #accelerazione lj vettoriale

def update_pos(dist_array, vel_array, accel, tau, i):

    temp_list = []

    for k in range(pos_array, N_particles):

        temp_list.append(dist_array[i, k] + vel_array[i, k] * tau + accel * tau**2 / 2)
    
    return temp_list

def update_vel(vel_array, accel, accel_new, tau, i):

    temp_list = []

    for k in range(N_particles):

        temp_list.append(vel_array[i, k] + (accel + accel_new) * tau**2 / 2)

    return temp_list


## NOTE: costanti della simulazione
mass_of_argon = 39.948      #amu
epsilon = .0103             
sigma = 3.4                 #Angstron
L_box = 20                  #Angstron

## NOTE: variabili della simulazione
#N_particles = int(input('Inserire il numero di particelle da far interagire\n'))
N_particles = 2
N_steps = 5
tau = .1                    #femtosecondi
initial_T = 300             #Kelvin

## NOTE: inizializziamo le posizioni, le velocità, le distanze e le accelerazioni
pos_array = init_pos(N_particles, N_steps)
vel_array = init_vel(N_particles, N_steps, initial_T)
dist_array = get_radius(pos_array, 0)
accel = lj_acceleration(dist_array, epsilon, sigma, 0)

## NOTE: eseguiamo il programma con un loop for, il primo punto è già inserito
for i in range(1, N_steps):

    dist_array = get_radius(pos_array, i)
    accel_new = lj_acceleration(dist_array, epsilon, sigma, i)
    pos_list = update_pos(dist_array, vel_array, accel, tau, i)
    vel_list = update_vel(vel_array, accel, accel_new, tau, i)

    accel = accel_new
    pos_array[:, i] = pos_list
    vel_array[:, i] = vel_list

print(pos_array)
print(vel_array)