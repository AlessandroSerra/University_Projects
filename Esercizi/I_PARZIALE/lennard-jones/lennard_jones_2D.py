import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import scipy.optimize as spo
import vec2d as v2d
import time as t

t_s = t.time()       #tempo di inizio simulazione

'''Funzioni necessarie ad eseguire il programma'''

## NOTE: funzione che inizializza le posizioni degli atomi
def init_values(N_particles, N_steps, initial_T, mass_of_argon, L_box, sigma):

    pos_array = np.full((N_particles, N_steps), fill_value = v2d.vec2d(0, 0))
    vel_array = np.full((N_particles, N_steps), fill_value = v2d.vec2d(0, 0))
    
    np.random.seed()
    M = np.random.rand()
    N = 2 - (np.random.rand() - 1)

    pos_x = [1 for i in range(N_particles)]
    pos_y = np.linspace(2, L_box - 2, N_particles)      #posizioni lungo y egualmente distanziate
    v_max = np.sqrt(spc.Boltzmann * initial_T / (mass_of_argon * spc.e))

    for i in range(N_particles):

        pos_array[i, 0] = v2d.vec2d(pos_x[i], pos_y[i])
            
        vel_x = M * v_max
        vel_y = N * v_max

        vel_array[i, 0] = v2d.vec2d(vel_x, vel_y)

    tau = (sigma / v_max) * 1e-3
    print('tau = ', tau)

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

        pos_list.append((pos_array[k, i] + vel_array[k, i] * tau + 0.5 * acc[k] * tau**2) % L_box)

    return pos_list


##NOTE: funzione che aggiorna le velocità con velocity-verlet
def update_vel(vel_array, acc, acc_new, N_particles, tau, i):

    vel_list = []

    for k in range(N_particles):

        vel_list.append(vel_array[k, i] + 0.5 * (acc[k] + acc_new[k]) * tau)

    return vel_list


##NOTE: funzione per il fit delle velocità secondo Maxwell-Boltzmann
def velocity_curve_fit(initial_T, vel_array, mass_of_argon):

    def M_B_distribution(velocity_list, initial_T, mass_of_argon, C):

        T_Kb = initial_T * spc.Boltzmann
        return C / T_Kb * velocity_list * np.exp(- mass_of_argon / (2 * T_Kb) * velocity_list)

    velocity_list = []

    for k in range(N_particles):

        velocity_list.append(vel_array[k ,-1]*vel_array[k, -1])        #fa la differenza dell'ultima velocità dalla prima

    counts, bins = np.histogram(velocity_list, bins = 30)
    bins = bins[1:] - (bins[1] - bins[0]) / 2
    p, cov = spo.curve_fit(M_B_distribution, bins, counts, p0 = [30])
    x_fit = np.linspace(bins[0], bins[-1], 1000)
    y_fit = M_B_distribution(x_fit, initial_T, mass_of_argon, p[0])

    return x_fit, y_fit, p[0], velocity_list


##NOTE: funzione che calcola lo spostamento quadratico medio
def get_mean_quadratic_deviation(pos_array, N_particles, N_steps):

    quadratic_deviation = 0
    mean_quadratic_deviation = []

    for i in range(N_steps):
        for k in range(N_particles):

            quadratic_deviation = quadratic_deviation + ((pos_array[k, i] - pos_array[k, 0]).mod())**2

        mean_quadratic_deviation.append(quadratic_deviation / N_particles)

    return mean_quadratic_deviation


##NOTE: funzione di run del programma
def run_lennard_jones2D(N_particles, N_steps, initial_T, mass_of_argon, L_box):

    ##NOTE: primo ciclo di inizialzzazione
    time_list = [0]
    pos_array, vel_array, tau = init_values(N_particles, N_steps, initial_T, mass_of_argon, L_box, sigma)
    dist_array = get_distance(pos_array, N_particles, 0)
    acc = get_acceleration(dist_array, epsilon, sigma, 0)

    ##NOTE: ciclo che gestisce il run del programma per N-1 passi
    for i in range(N_steps - 1):

        pos_array[:, i + 1] = update_pos(pos_array, vel_array, acc, N_particles, L_box, tau, i)
        dist_array = get_distance(pos_array, N_particles, i + 1)
        acc_new = get_acceleration(dist_array, epsilon, sigma, i + 1)
        vel_array[:, i + 1] = update_vel(vel_array, acc, acc_new, N_particles, tau, i)

        acc = acc_new
        time_list.append(tau * (i + 1))

    return pos_array, vel_array, time_list


'''Main del programma'''

## NOTE: parametri di simulazione
mass_of_argon = 39.948      #amu
epsilon = .0103             
sigma = 3.4                 #Angstron
L_box = 10                  #Angstron
#N_particles = int(input('Inserire il numero di particelle da far interagire\n'))
N_particles = 2
N_steps = 1000
initial_T = 300              #Kelvin

what_to_plot = 't' #input('Insert "t" to see the alpha particles trajectories or "a" to see the scattering angles:\t')


##NOTE: funzione di run del programma
pos_array, vel_array, time_list = run_lennard_jones2D(N_particles, N_steps, initial_T, mass_of_argon, L_box)

#x_fit, y_fit, C_fit, velocity_list = velocity_curve_fit(initial_T, vel_array, mass_of_argon)

quadratic_mean_deviation = get_mean_quadratic_deviation(pos_array, N_particles, N_steps)

print(quadratic_mean_deviation[999])

'''parte del programma dedicata alla rappresentazione grafica'''

fig, ax = plt.subplots()

if what_to_plot == 't':

    for i in range(N_particles):

        ax.plot([pos_array[i, step].x for step in range(N_steps)], [pos_array[i, step].y for step in range(N_steps)], label = 'atom {}'.format(i + 1), linewidth = 1.5)

    ax.set_title('Lennard-Jones Interactions for Argon')
    ax.set_xlabel('x (Å)')
    ax.set_ylabel('y (Å)')

elif what_to_plot == 'v':

    print('Costante C della distribuzione:', C_fit)

    ax.hist(velocity_list, histtype = 'step', bins = 30, label = 'Data')
    ax.set_yscale('log')                #scala log per visualizzare meglio alti numeri di conteggi
    ax.set_xlabel('Velocity')
    ax.set_ylabel('Counts')
    ax.set_title('Distribution of the Argon Atoms Velocities')
    ax.plot(x_fit, y_fit, label = 'Maxwell-Boltzmann Velocity fit', linewidth = 1)

elif what_to_plot == 'rmd':

    for i in range(N_steps):

        ax.plot(time_list[i], quadratic_mean_deviation[i])

##NOTE: il tempo di simulazione non include l'intervallo necessario a visualizzare il grafico
print('Time taken by the simulation:', t.time() - t_s, 'seconds')

plt.legend(frameon = False)
plt.show()

