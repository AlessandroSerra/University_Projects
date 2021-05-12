import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import scipy.optimize as spo
import vec2d as v2d
import time

t_s = time.time()       #tempo di inizio simulazione

'''Funzioni necessarie ad eseguire il programma'''

## NOTE: funzione che inizializza le posizioni degli atomi
def init_values(N_particles, N_steps, initial_T, mass_of_argon, L_box, sigma, status):

    pos_array = np.full((N_particles, N_steps), fill_value = v2d.vec2d(0, 0))
    vel_array = np.full((N_particles, N_steps), fill_value = v2d.vec2d(0, 0))
    
    np.random.seed()
    pos_y = np.linspace(0, L_box, N_particles)      #posizioni lungo y egualmente distanziate
    v_max = np.sqrt(spc.Boltzmann * initial_T / (mass_of_argon * spc.e))

    for i in range(N_particles):

        pos_array[i, 0] = v2d.vec2d(0, pos_y[i])
        
        if status == 'rand':    #caso di inizializzazione con velocità random ma dipendenti dalla temperature del sistema

            P = np.random.rand()
            R = 2 * np.random.rand() - 1
            vel_x = P * v_max
            vel_y = R * v_max

        else:   #caso di inizializzazione con velocità mono-modali dipendenti dalla temperatura del sistema
            
            vel_x = v_max
            vel_y = v_max

        vel_array[i, 0] = v2d.vec2d(vel_x, vel_y)

    tau = (sigma / v_max) * 1e-2

    return pos_array, vel_array, tau


##NOTE: funzione che recupera le distanze relative delle particelle step per step
def get_distance(pos_array, N_particles, i):

    dist_array = np.full((N_particles, N_particles), fill_value = v2d.vec2d(0, 0))

    ##NOTE: loop che restituiscono i vettori distanza tra le particelle
    for j in range(N_particles - 1):
        for k in range(j + 1, N_particles):

            dist_array[j, k] = (pos_array[j, i] - pos_array[k, i]).abs_value()  #funzione valore assoluto implementata in vec2d
            dist_array[k, j] = dist_array[j, k]

    return dist_array


##NOTE: funzione che recupera le accelerazioni relative delle particelle step per step
def get_acceleration(dist_array, epsilon, sigma, i):

    acc_rel_array = np.full((N_particles, N_particles), fill_value = v2d.vec2d(0, 0))
    acc_list = []

    ##NOTE: loop che restiruiscono i vettori accelerazione delle particelle
    for j in range(N_particles - 1):
        for k in range(j + 1, N_particles):

            acc_rel_array[j, k] = ((48 * epsilon * np.power(
                sigma, 12) / np.power(
                dist_array[j, k].mod(), 13)) - (24 * epsilon * np.power(
                sigma, 6) / np.power(dist_array[j, k].mod(), 7)
                )) / mass_of_argon * dist_array[j, k].unitary()

            acc_rel_array[k, j] = -1 * acc_rel_array[j, k]

    ##NOTE: loop che restituisce la somma vettoriale delle accelerazioni agenti su una particella
    for i in range(N_particles):

        acc_list.append(np.sum(acc_rel_array[i, :], axis = 0))

    return acc_list


##NOTE: funzione che aggiorna le posizioni con velocity-verlet\
def update_pos(pos_array, vel_array, acc, N_particles, L_box, tau, i):

    pos_list = []

    for k in range(N_particles):

        pos_list.append(pos_array[k, i] + vel_array[k, i] * tau + 0.5 * acc[k] * tau**2)

    return pos_list


##NOTE: funzione che aggiorna le velocità con velocity-verlet
def update_vel(vel_array, acc, acc_new, N_particles, tau, i):

    vel_list = []

    for k in range(N_particles):

        vel_list.append(vel_array[k, i] + 0.5 * (acc[k] + acc_new[k]) * tau)

    return vel_list


##NOTE: funzione di run del programma
def run_lennard_jones2D(N_particles, N_steps, initial_T, mass_of_argon, L_box, status):

    ##NOTE: primo ciclo di inizialzzazione
    pos_array, vel_array, tau = init_values(N_particles, N_steps, initial_T, mass_of_argon, L_box, sigma, status)
    dist_array = get_distance(pos_array, N_particles, 0)
    acc = get_acceleration(dist_array, epsilon, sigma, 0)

    ##NOTE: ciclo che gestisce il run del programma per N-1 passi
    for i in range(N_steps - 1):

        pos_array[:, i + 1] = update_pos(pos_array, vel_array, acc, N_particles, L_box, tau, i)
        dist_array = get_distance(pos_array, N_particles, i + 1)
        acc_new = get_acceleration(dist_array, epsilon, sigma, i + 1)
        vel_array[:, i + 1] = update_vel(vel_array, acc, acc_new, N_particles, tau, i)

        acc = acc_new

    return pos_array, vel_array


'''Main del programma'''

## NOTE: parametri di simulazione
mass_of_argon = 39.948      #amu
epsilon = .0103             
sigma = 3.4                 #Angstron
L_box = 50                  #Angstron
#N_particles = int(input('Inserire il numero di particelle da far interagire\n'))
N_particles = 3
N_steps = 1000
initial_T = 300              #Kelvin

##NOTE: parametro per inizializzare le velocità random o mono-modali
status = 'rand' #input('Select initial random velocities "r" or mono-modal initial velocities "m":\t')

##NOTE: funzione di run del programma
pos_array, vel_array = run_lennard_jones2D(N_particles, N_steps, initial_T, mass_of_argon, L_box, status)


'''parte del programma dedicata alla rappresentazione grafica'''

fig, ax = plt.subplots()

for i in range(N_steps):
    for k in range(N_particles):

        ax.plot(pos_array[k, i].x, pos_array[k, i].y, marker = '.', color = 'b', markersize = 1.5)

ax.set_title('Lennard-Jones Interactions for Argon')
ax.set_xlabel('x (Å)')
ax.set_ylabel('y (Å)')

##NOTE: il tempo di simulazione non include l'intervallo necessario a visualizzare il grafico
print('Time taken by the simulation:', time.time() - t_s, 'seconds')

plt.show()
