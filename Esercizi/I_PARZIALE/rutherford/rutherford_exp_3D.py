from operator import pos
import numpy as np
import matplotlib.pyplot as plt 
import scipy.constants as spc
import scipy.optimize as spo
import vec3d as v3d
import time

t_s = time.time()       #tempo di inizio simulazione

'''Funzioni necessarie ad eseguire il programma'''

## NOTE: funzione che inizializza le posizioni delle particelle alpha
def init_values(N_particles, N_steps, Ze, Zo, energy, alpha_mass):

    pos_array = np.full((N_particles, N_steps), fill_value = v3d.vec3d(0, 0, 0))
    vel_array = np.full((N_particles, N_steps), fill_value = v3d.vec3d(0, 0, 0))
    inter_distance = Ze * Zo * spc.e**2 / (4 * np.pi * spc.epsilon_0) / energy
    angles_list = []

    emittent_radius = 100 * inter_distance
    pos_x = -100 * inter_distance
    vel_x = np.sqrt(2 * energy / alpha_mass)
    tau = inter_distance / vel_x
    k = 0                                       #inizializo l'indice del loop while per riampire le posizioni delle particelle
    np.random.seed()                            #diverso seed per la funzione random.rand() ad ogni run del programma

    while k < N_particles:
        
        pos_y = emittent_radius * (2 * np.random.rand() - 1)
        pos_z = emittent_radius * (2 * np.random.rand() - 1)

        if np.sqrt(pos_y**2 + pos_z**2) <= (emittent_radius):      #condizione di raggio cilindrico-

            pos_array[k, 0] = v3d.vec3d(pos_x, pos_y, pos_z)
            vel_array[k, 0] = v3d.vec3d(vel_x, 0, 0)
            k += 1                                                                  #incrementa SOLO se la posizione generata random sta dentro il cerchio

    return pos_array, vel_array, angles_list, inter_distance, tau


##NOTE: funzione che recupera le accelerazioni relative delle particelle step per step
def get_acceleration(pos_array, Ze, Zo, alpha_mass, i):

    acc_list = []

    for k in range(N_particles):

        acc_list.append(((Ze * Zo * spc.e**2) / (4 * np.pi * spc.epsilon_0 * alpha_mass)) * (pos_array[k, i] / pos_array[k, i].mod()**3))

    return acc_list


##NOTE: funzione che aggiorna le posizioni con velocity-verlet
def update_pos(pos_array, vel_array, acc, N_particles, tau, i):

    pos_temp_list = []

    for k in range(N_particles):

        pos_temp_list.append(pos_array[k, i] + vel_array[k, i] * tau + 0.5 * acc[k] * tau**2)

    return pos_temp_list


##NOTE: funzione che aggiorna le velocità con velocity-verlet
def update_vel(vel_array, acc, acc_new, N_particles, tau, i):

    vel_temp_list = []

    for k in range(N_particles):

        vel_temp_list.append(vel_array[k, i] + 0.5 * (acc[k] + acc_new[k]) * tau)

    return vel_temp_list


##NOTE: funzione che fitta la distribuzione attesa degli angoli di scattering
def scattering_angles_curve_fit(angles_list):

    def theor_distribution(angles_list, alpha, N0):

        return N0 / np.power(np.sin(angles_list / 2), alpha)        #N0 parametro del fin NON numero di particelle inziale

    counts, bins = np.histogram(angles_list, bins = 50)   
    bins = bins[1:] - (bins[1] - bins[0]) / 2
    p, cov = spo.curve_fit(theor_distribution, bins, counts, p0 = [4, 1e-3])
    x_fit  = np.linspace(bins[0], bins[-1], 1000)
    y_fit = theor_distribution(x_fit, p[0], p[1])

    return x_fit, y_fit, p[0], p[1], cov


##NOTE: funzione di run del programma
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


    what_to_plot = 'a' #input('Insert "t" to see the alpha particles trajectories or "a" to see the scattering angles:\t')

    ##NOTE: funzione di fit degli angoli
    x_fit, y_fit, alpha, N0, cov = scattering_angles_curve_fit(angles_list) 


    ## NOTE: segmento dedicato alla rappresentazione delle traiettorie delle particelle
    if what_to_plot == 't':

        fig = plt.subplots()
        ax = plt.axes(projection = '3d')
        ax.plot3D(0, 0, 0, marker = 'o', color = '#FFD700')         #atomo d'oro
    
        for k in range(pos_array.shape[0]):

            ax.plot3D([(pos_array[k, i].x / inter_distance) for i in range(pos_array.shape[1])],
            [(pos_array[k, i].y / inter_distance) for i in range(pos_array.shape[1])], 
            [pos_array[k, i].z / inter_distance for i in range(pos_array.shape[1])], color = 'purple')

        ax.set_title('3D Rutherford Experiment')
        ax.set_xlabel('x / interaction distance')
        ax.set_ylabel('y / interaction distance')
        ax.set_zlabel('z / interaction distance')

    ##NOTE: segmento dedicato ala rappresentazione degli angoli di scattering
    elif what_to_plot == 'a': 

        print('N0:\t', N0)
        print('Esponente del seno:\t', alpha)
        print('Covarianza:\n', cov)

        fig, ax = plt.subplots()
        ax.hist(angles_list, histtype = 'step', bins = 50, label = 'Angles')
        ax.plot(x_fit, y_fit, label = 'Scattered angles fit', linewidth = 1)
        ax.set_yscale('log')                #scala log per visualizzare meglio alti numeri di conteggi
        ax.set_xlabel('$\\theta$')
        ax.set_ylabel('Counts')
        ax.set_title('Distribution of the Scattered Angles in 3 dimentions')

        plt.legend(frameon = False)
        
    ##NOTE: messaggio di inserimento errato del codice what_to_do
    else:
        print('WARNING: wrong code, use "t" to see the trajectories or "a" to see the angles distribution\n')
        exit()

    print('Number of particels:\t', pos_array.shape[0])

    ##NOTE: il tempo di simulazione non include l'intervallo necessario a visualizzare il grafico
    print('Time taken by the simulation:\t', time.time() - t_s, 'seconds')

    plt.show()


'''Main del programma'''

## NOTE: parametri di simulazione
N_particles = 100000 
N_steps = 200
energy = 5e5 * spc.electron_volt    #energia di 5MeV convertita in Joule                 
Ze, Zo = 2, 79                                          #numero atomico di elio (Ze) ed oro (Zo)
alpha_mass = 2 * spc.proton_mass + 2 * spc.neutron_mass     #massa atomica elio in uma


##NOTE: funzione di run del programma
run_rutherford_exp3D(N_particles, N_steps, Ze, Zo, energy, alpha_mass)