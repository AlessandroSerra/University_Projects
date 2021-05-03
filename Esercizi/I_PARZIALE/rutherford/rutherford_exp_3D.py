import numpy as np
import matplotlib.pyplot as plt 
import scipy.constants as spc
import scipy.optimize as spo
import vec3d as v3d
import time

t_s = time.time()       #tempo di inizio simulazione

'''
Funzioni necessarie ad eseguire il programma
'''

def init_values(N_particles, N_steps, Ze, Zo, energy, alpha_mass):

    pos_array = np.full((N_particles, N_steps), fill_value = v3d.vec3d(0, 0, 0))
    vel_array = np.full((N_particles, N_steps), fill_value = v3d.vec3d(0, 0, 0))
    inter_distance = Ze * Zo * spc.e**2 / (4 * np.pi * spc.epsilon_0 * alpha_mass) * alpha_mass / energy
    angles_list = []

    emittent_radius = 5
    pos_x = -100 * inter_distance
    vel_x = np.sqrt(2 * energy / alpha_mass)
    tau = inter_distance / vel_x
    np.random.seed()                            #diverso seed per la funzione random.rand() ad ogni run del programma

    for k in range(N_particles):
        
        pos_y = emittent_radius * (2 * np.random.rand() - 1) * inter_distance
        pos_z = emittent_radius * (2 * np.random.rand() - 1) * inter_distance
        pos_array[k, 0] = v3d.vec3d(pos_x, pos_y, pos_z)
        vel_array[k, 0] = v3d.vec3d(vel_x, 0, 0)

    return pos_array, vel_array, angles_list, inter_distance, tau


def get_acceleration(pos_array, Ze, Zo, alpha_mass, i):

    acc_list = []

    for k in range(N_particles):

        acc_list.append(((Ze * Zo * spc.e**2) / (4 * np.pi * spc.epsilon_0 * alpha_mass)) * (pos_array[k, i] / pos_array[k, i].mod()**3))

    return acc_list


def update_pos(pos_array, vel_array, acc, N_particles, tau, i):

    pos_temp_list = []

    for k in range(N_particles):

        pos_temp_list.append(pos_array[k, i] + vel_array[k, i] * tau + 0.5 * acc[k] * tau**2)

    return pos_temp_list


def update_vel(vel_array, acc, acc_new, N_particles, tau, i):

    vel_temp_list = []

    for k in range(N_particles):

        vel_temp_list.append(vel_array[k, i] + 0.5 * (acc[k] + acc_new[k]) * tau)

    return vel_temp_list


def scattering_angles_curve_fit(angles_list, N_particles):

    def theor_distribution(angles_list, N_particles, alpha):

        return N_particles / np.power(np.sin(angles_list / 2), alpha)

    counts, bins = np.histogram(angles_list, bins = 30)   
    bins = bins[1:] - (bins[1] - bins[0]) / 2
    p, cov = spo.curve_fit(theor_distribution, bins, counts, p0 = [1, 4], sigma = counts)
    x  = np.linspace(bins[0], bins[-1], 1000)
    y = theor_distribution(x, p[0], p[1])

    return x, y


def run_rutherford_exp3D(N_particles, N_steps, Ze, Zo, energy, alpha_mass):

    pos_array, vel_array, angles_list, inter_distance, tau = init_values(N_particles, N_steps, Ze, Zo, energy, alpha_mass)
    acc = get_acceleration(pos_array, Ze, Zo, alpha_mass, 0)

    for i in range(N_steps - 1):

        pos_array[:, i + 1] = update_pos(pos_array, vel_array, acc, N_particles, tau, i)
        acc_new = get_acceleration(pos_array, Ze, Zo, alpha_mass, i + 1)
        vel_array[:, i + 1] = update_vel(vel_array, acc, acc_new, N_particles, tau, i)
        acc = acc_new

    for i in range(N_particles):

        angles_list.append(vel_array[i, 0].get_angle(vel_array[i, -1], 'rad'))

    return pos_array, vel_array, angles_list, inter_distance

'''
Corpo principale del programma
'''

## NOTE: parametri di simulazione
N_particles = 100000
N_steps = 200
energy = 5e5 * spc.electron_volt    #energia di 5MeV convertita in Joule                 #secondi
Ze, Zo = 2, 79                                              #numero atomico di elio (Ze) ed oro (Zo)
alpha_mass = 2 * spc.proton_mass + 2 * spc.neutron_mass     #massa atomica elio in uma

what_to_plot = input('Insert "t" to see the alpha particles trajectories or "a" to see the scattering angles:\t')

pos_array, vel_array, angles_list, inter_distance = run_rutherford_exp3D(N_particles, N_steps, Ze, Zo, energy, alpha_mass)
x_fit, y_fit = scattering_angles_curve_fit(angles_list, N_particles)

'''
Parte del codice dedicata esclusivamente alla rappresentazione grafica
'''

## NOTE: segmento dedicato alla rappresentazione delle traiettorie delle particelle
if what_to_plot == 't':

    fig = plt.subplots()
    ax = plt.axes(projection = '3d')
    ax.plot3D(0, 0, 0, marker = 'o', color = '#FFD700')         #atomo d'oro
    
    for i in range(N_steps):
        for k in range(N_particles):

            ax.plot3D([pos_array[k, i].x / inter_distance],[pos_array[k, i].y / inter_distance], [pos_array[k, i].z / inter_distance], marker = '.', markersize = .3, color = '#39FF14')

    ax.set_title('3D Rutherford Experiment')
    ax.set_xlabel('x / impact parameter')
    ax.set_ylabel('y / impact parameter')
    ax.set_zlabel('z / impact parameter')

##NOTE: segmento dedicato ala rappresentazione degli angoli di scattering
elif what_to_plot == 'a': 

        fig, ax = plt.subplots()
        ax.hist(angles_list, histtype = 'step', bins = 30, label = 'Data')
        ax.set_yscale('log')                #scala log per visualizzare meglio alti numeri di conteggi
        ax.set_xlabel('$\\theta$')
        ax.set_ylabel('Counts')
        ax.set_title('Distribution of the Scattered Angles in 3 dimentions')
        ax.plot(x_fit, y_fit, label = 'Scattered angles fit', linewidth = 1)
        plt.legend()

##NOTE: messaggio di inserimento errato del codice what_to_do
else:
    print('WARNING: wrong code, use "t" to see the trajectories or "a" to see the angles distribution\n')
    exit()

##NOTE: il tempo di simulazione non include l'intervallo necessario a visualizzare il grafico
print('Time taken by the simulation:\t', time.time() - t_s, 'seconds\n')

plt.show()