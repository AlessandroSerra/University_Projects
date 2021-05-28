import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import scipy.optimize as spo
import vec2d as v2d
import time as t

t_s = t.time()       #tempo di inizio simulazione

'''Funzioni necessarie ad eseguire il programma'''

## NOTE: funzione che inizializza le posizioni degli atomi
def init_values(N_particles, N_steps, initial_T, mass_of_argon, L_box, sigma, status):

    pos_array = np.full((N_particles, N_steps), fill_value = v2d.vec2d(0, 0))
    vel_array = np.full((N_particles, N_steps), fill_value = v2d.vec2d(0, 0))
    v_max = np.sqrt(2 * spc.Boltzmann * initial_T / (mass_of_argon * spc.electron_volt))    #tolto il 2
    
    np.random.seed()
    M = 2 * (np.random.rand(N_particles) - 1)
    N = 2 * (np.random.rand(N_particles) - 1)

    pos_array[:, 0] = np.array([v2d.vec2d(i + 10, j + 10) * sigma * 1.1 for i in range(int(np.sqrt(N_particles))) for j in range(int(np.sqrt(N_particles)))])

    if status == 'r':
        for i in range(N_particles):

            vel_x = M[i] * v_max
            vel_y = N[i] * v_max
            vel_array[i, 0] = v2d.vec2d(vel_x, vel_y)   #angstrom / femtosecondi

    elif status == 'm':
        for i in range(N_particles):

            if M[i] < 0:            #queato parametro cambia le direzioni di alcuni atomi senza variarne la velocita'
                M[0] = - M[0]

            vel_x = M[0] * v_max
            vel_y = M[1] * v_max
            vel_array[i, 0] = v2d.vec2d(vel_x, vel_y)

    tau = (sigma / v_max) * 1e-4        #femtosecondi

    return pos_array, vel_array, tau


##NOTE: funzione che recupera le accelerazioni relative delle particelle step per step
def get_acceleration(pos_array, N_particles, epsilon, sigma, i):

    acc_rel_array = np.full((N_particles, N_particles), fill_value = v2d.vec2d(0, 0))

    ##NOTE: loop che restiruiscono i vettori accelerazione delle particelle
    for j in range(N_particles - 1):
        for k in range(j + 1, N_particles):

            dist = pos_array[j, i] - pos_array[k, i]  #distanza relativa tra particella j-esima e k-esima all'istante i-esimo

            acc_rel_array[j, k] = 4 * epsilon * ((12 * sigma**12 / dist.mod()**13) - (6 * sigma**6 / dist.mod()**7)) / mass_of_argon * dist.unitary()
            acc_rel_array[k, j] = - acc_rel_array[j, k]     #3a legge di Newton

    ##NOTE: loop che restituisce la somma vettoriale delle accelerazioni agenti su una particella
    acc_array = np.sum(acc_rel_array, axis = 1)     #axis = 1 somma gli elementi sulle righe dell'array

    return acc_array


##NOTE: funzione che aggiorna le posizioni con velocity-verlet\
def update_pos(pos_array, vel_array, acc, N_particles, L_box, tau, i):

    positions = np.full((N_particles), fill_value = v2d.vec2d(0, 0))

    for k in range(N_particles):

        positions[k] = (pos_array[k, i] + vel_array[k, i] * tau + 0.5 * acc[k] * tau**2) % L_box

    return positions


##NOTE: funzione che aggiorna le velocità con velocity-verlet
def update_vel(vel_array, acc, acc_new, N_particles, tau, i):

    velocities = np.full((N_particles), fill_value = v2d.vec2d(0, 0))

    for k in range(N_particles):

        velocities[k] = vel_array[k, i] + 0.5 * (acc[k] + acc_new[k]) * tau

    return velocities


##NOTE: funzione per il fit delle velocità secondo Maxwell-Boltzmann
def velocity_curve_fit(vel_array, status):

    def M_B_distribution(v, C, A):

        return C * v**2 * np.exp(- A * v**2)        #funzione da fittare

    velocity_list = []

    for k in range(N_particles):

        if status == 'r':
            velocity_list.append((vel_array[k, -1]).mod())

        else:
            velocity_list.append((vel_array[k, -1] - vel_array[k, 0]).mod())
    
    velocity_list = np.array(velocity_list)

    ##NOTE: segmento di codice per inizializzare curve-fit
    counts, bins = np.histogram(velocity_list, bins = 10)
    bins = bins[1:] - (bins[1] - bins[0]) / 2
    p, cov = spo.curve_fit(M_B_distribution, bins, counts, p0 = [1e-2, 1], maxfev = 5000)
    x_fit = np.linspace(bins[0], bins[-1], 1000)
    y_fit = M_B_distribution(x_fit, p[0], p[1])

    return velocity_list, p[0], p[1], x_fit, y_fit


##NOTE: funzione che calcola lo spostamento quadratico medio
def get_mean_quadratic_deviation(pos_array, N_particles, N_steps):

    quadratic_deviation = 0
    mean_quadratic_deviation = np.zeros((N_steps))

    for i in range(N_steps):
        for k in range(N_particles):

            quadratic_deviation += ((pos_array[k, i] - pos_array[k, 0]).mod())**2

        mean_quadratic_deviation[i] = quadratic_deviation / N_particles

    return mean_quadratic_deviation


##NOTE: funzione che fitta la curva dello spostamento quadratico medio in un intervallo pseudo-lineare
def get_diffusion_par(mean_quadratic_deviation, tau):

    x_fit = np.array([tau * i for i in range(int(2 * N_steps / 3), N_steps)])
    y_fit = np.array(mean_quadratic_deviation[int(2 * N_steps / 3): N_steps])
    p = np.polyfit(x_fit, y_fit, 1)

    x_fit1 = np.array([tau * i for i in range(0, int(N_steps / 3))])
    y_fit1 = np.array(mean_quadratic_deviation[0: int(N_steps / 3)])
    p1 = np.polyfit(x_fit1, y_fit1, 1)

    x_fit2 = np.array([tau * i for i in range(int(N_steps / 3), int(2 * N_steps / 3))])
    y_fit2 = np.array(mean_quadratic_deviation[int(N_steps / 3): int(2 * N_steps / 3)])
    p2 = np.polyfit(x_fit2, y_fit2, 1)

    return x_fit, x_fit1, x_fit2, p[0], p[1], p1[0], p1[1], p2[0], p2[1]


##NOTE: funzione di run del programma
def run_lennard_jones2D(N_particles, N_steps, initial_T, mass_of_argon, L_box):

    status =  input('Insert "r" for random initial velocities or "m" for mono-modal velocities:\t')

    ##NOTE: primo ciclo di inizialzzazione
    time_list = [0]
    pos_array, vel_array, tau = init_values(N_particles, N_steps, initial_T, mass_of_argon, L_box, sigma, status)
    acc = get_acceleration(pos_array, N_particles, epsilon, sigma, 0)

    ##NOTE: ciclo che gestisce il run del programma per N-1 passi
    for i in range(N_steps - 1):

        pos_array[:, 1] = update_pos(pos_array, vel_array, acc, N_particles, L_box, tau, 0)
        pos_array[:, i + 1] = update_pos(pos_array, vel_array, acc, N_particles, L_box, tau, i)
        acc_new = get_acceleration(pos_array, N_particles, epsilon, sigma, i + 1)
        vel_array[:, i + 1] = update_vel(vel_array, acc, acc_new, N_particles, tau, i)

        acc = acc_new
        time_list.append(tau * (i + 1))     #lista dei passi temporali

    ##NOTE: funzioni di fit e per il RMD
    velocity_list, C, A, x_fit, y_fit = velocity_curve_fit(vel_array, status)
    mean_quadratic_deviation = get_mean_quadratic_deviation(pos_array, N_particles, N_steps)
    x_diff_fit, x_diff_fit1, x_diff_fit2, coeff_ang, coeff_retta, coeff_ang1, coeff_retta1, coeff_ang2, coeff_retta2 = get_diffusion_par(mean_quadratic_deviation, tau)


    ##NOTE: segmento per la creazioe delle figure
    fig, ax = plt.subplots()
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    ##NOTE: rappresentazione delle  traiettorie delle particelle
    for i in range(N_particles):

        ax.plot([pos_array[i, step].x for step in range(N_steps)], [pos_array[i, step].y for step in range(N_steps)], linewidth = 2)

    ax.set_title('Lennard-Jones Interactions for Argon')
    ax.set_xlabel('x (Å)')
    ax.set_ylabel('y (Å)')
    
    ##NOTE: rappresntazione grafica dello spostamento quadratico medio e del suo fit
    coeff_retta = np.array([coeff_retta for i in range(len(x_diff_fit))])
    coeff_retta1 = np.array([coeff_retta1 for i in range(len(x_diff_fit1))])
    coeff_retta2 = np.array([coeff_retta2 for i in range(len(x_diff_fit2))])
    y_diff_fit = coeff_retta + coeff_ang * x_diff_fit
    y_diff_fit1 = coeff_retta1 + coeff_ang1 * x_diff_fit1
    y_diff_fit2 = coeff_retta2 + coeff_ang2 * x_diff_fit2

    ax1.plot(time_list, mean_quadratic_deviation, label = 'Root Mean Displacement')
    ax1.plot(x_diff_fit, y_diff_fit, label = 'Fit 2a regione pseudo-lineare')
    ax1.plot(x_diff_fit1, y_diff_fit1, label = 'Fit 1a regione pseudo-lineare')
    ax1.plot(x_diff_fit2, y_diff_fit2, label = 'Fit 3a regione pseudo-lineare')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('RMD')
    
    print('Parametro di autodiffusione 1a regione pseudo-lineare:\t', coeff_ang1 / 4)
    print('Parametro di autodiffusione 2a regione lineare\t', coeff_ang2 / 4)
    print('Parametro di autodiffusione 3a regione pseudo-lineare:\t', coeff_ang / 4)

    ##NOTE: rappresentazione grafica della distribuzione delle velocità
    ax2.hist(velocity_list, histtype = 'step', bins = 10, label = 'Particles velocities')
    ax2.plot(x_fit, y_fit, label = 'Maxwell-Boltzmann velocities distribution fit', linewidth = 1)
    ax2.set_xlabel('Velocity')
    ax2.set_ylabel('Counts')
    ax2.set_title('Distribution of the Argon Atoms Velocities')


    ##NOTE: il tempo di simulazione non include l'intervallo necessario a visualizzare il grafico
    print('Time taken by the simulation:', t.time() - t_s, 'seconds')

    plt.legend(frameon = False)     #nessun bordo nella legenda
    plt.show()



'''Main del programma'''

## NOTE: parametri di simulazione
mass_of_argon = 39.948                      #amu
epsilon = .0103                             #eV
sigma = 3.46                                #Angstrom
N_particles = 64
L_box = N_particles * sigma                 #Angstrom
N_steps = 10000              
initial_T = 300                             #Kelvin

##NOTE: funzione di run del programma
run_lennard_jones2D(N_particles, N_steps, initial_T, mass_of_argon, L_box)