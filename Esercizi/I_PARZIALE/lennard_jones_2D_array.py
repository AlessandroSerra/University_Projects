import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import vec2d as v2d
import time 

t_s = time.time()

## NOTE: funzione che inizializza le posizioni degli atomi
def init_pos(N_particles, N_steps):

    pos_array_x = np.zeros((N_particles, N_steps))
    pos_array_y = np.zeros((N_particles, N_steps))

    ## NOTE: loop che assegna la posizione iniziale a tutte le particelle
    for i in range(N_particles):
        
        pos_array_x[0, i] = 20 * (2 * np.random.rand() - 1)
        pos_array_y[0, i] = 1
    
    return pos_array_x, pos_array_y

def init_vel(N_particles, N_steps, initial_T):

    vel_array_x= np.zeros((N_steps, N_particles))
    vel_array_y= np.zeros((N_steps, N_particles))

    ## NOTE: loop che assegna la vlocità iniziale a tutte le particelle    
    for i in range(N_particles):

        P = np.random.rand(N_particles) - 0.5
        R = np.random.rand(N_particles) - 0.5
        vel_array_x[0, i] = 1 #P * np.sqrt(spc.Boltzmann * initial_T / (mass_of_argon / spc.e))
        vel_array_y[0, i] = 1 #R * np.sqrt(spc.Boltzmann * initial_T / (mass_of_argon / spc.e))

    return vel_array_x, vel_array_y

#RIVEDI QUESTA FUNZIONE PER LA DISTANZA FRA DUE PARTICELLE
def get_radius(pos_array_x, pos_array_y, i):

    dist_array_x = np.zeros((N_particles, N_particles))
    dist_array_y = np.zeros((N_particles, N_particles))

    for j in range(0, N_particles - 1):
            for k in range(j + 1, N_particles):

                dist_array_x[j, k] = pos_array_x[j, i] - pos_array_x[k, i]
                dist_array_y[j, k] = pos_array_y[j, i] - pos_array_y[k, i]

    return dist_array_x, pos_array_y

# NOTE: definiamo la forza di interazione come la derivata parziale opposta del potenziale
def lj_acceleration(dist_array_x, dist_array_y, epsilon, sigma, i):

    acc_array_x = np.zeros((N_particles, N_steps))
    acc_array_y = np.zeros((N_particles, N_steps))
    force_array_x = np.zeros((N_particles, N_steps))
    force_array_y = np.zeros((N_particles, N_steps))

    for j in range(0, N_particles - 1):     #dist_array contiene il vettore differenza tra i vettori posizione
        for k in range(j + 1, N_particles):
            acc_array_x[j, k] = (48 * epsilon * np.power(
                sigma, 12) / np.power(
                dist_array_x[j, k], 13) - 24 * epsilon * np.power(
                sigma, 6) / np.power(dist_array_x[j, k], 7)
                ) / mass_of_argon

            acc_array_y[j, k] = (48 * epsilon * np.power(
                sigma, 12) / np.power(
                dist_array_y[j, k], 13) - 24 * epsilon * np.power(
                sigma, 6) / np.power(dist_array_y[j, k], 7)
                ) / mass_of_argon

            acc_array_x[k, j] = - acc_array_x[j, k]     #terza legge di newton
            acc_array_y[k, j] = - acc_array_y[j, k]

    return acc_array_x, acc_array_y   #accelerazione lj vettoriale

def update_pos(dist_array_x, dist_array_y, vel_array_x, vel_array_y, accel_x, accel_y, tau, i):

    temp_list_x = []
    temp_list_y = []

    for k in range(N_particles):

        temp_list_x.append(dist_array_x[i, k] + vel_array_x[i, k] * tau + accel_x * tau**2 / 2)
        temp_list_y.append(dist_array_y[i, k] + vel_array_y[i, k] * tau + accel_y * tau**2 / 2)
    

    return temp_list_x, temp_list_y

def update_vel(vel_array_x, vel_array_y, accel_x, accel_y, accel_new_x, accel_new_y, tau, i):

    temp_list_x = []
    temp_list_y = []

    for k in range(N_particles):

        temp_list_x.append(vel_array_x[i, k] + (accel_x + accel_new_x) * tau**2 / 2)
        temp_list_y.append(vel_array_y[i, k] + (accel_y + accel_new_y) * tau**2 / 2)

    return temp_list_x, temp_list_y


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
pos_array_x, pos_array_y = init_pos(N_particles, N_steps)
vel_array_x, vel_array_y = init_vel(N_particles, N_steps, initial_T)
dist_array_x, dist_array_y = get_radius(pos_array_x, pos_array_y, 0)
accel_x, accel_y = lj_acceleration(dist_array_x, dist_array_y, epsilon, sigma, 0)

## NOTE: eseguiamo il programma con un loop for, il primo punto è già inserito
for i in range(1, N_steps):

    dist_array_x, dist_array_y = get_radius(pos_array_x, pos_array_y, i)
    accel_new_x, accel_new_y = lj_acceleration(dist_array_x, dist_array_y, epsilon, sigma, i)
    pos_list_x, pos_list_y = update_pos(dist_array_x, dist_array_y, vel_array_x, vel_array_y, accel_x, accel_y, tau, i)
    vel_list_x, vel_list_y = update_vel(vel_array_x, vel_array_y, accel_x, accel_y, accel_new_x, accel_new_y, tau, i)

    accel_x, accel_y = accel_new_x, accel_new_y

'''
    for m in range(N_particles):
        pos_array_x[m, i] = pos_list_x[m]
        pos_array_y[m, i] = pos_list_y[m]
        vel_array_x[m, i] = vel_list_x[m]
        vel_array_y[m, i] = vel_list_y[m]

print(pos_array_x)
print(pos_array_y)
print(vel_array_x)
print(vel_array_y)
'''