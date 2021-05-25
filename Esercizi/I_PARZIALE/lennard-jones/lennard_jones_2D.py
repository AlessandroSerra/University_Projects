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
    v_max = np.sqrt(2 * spc.Boltzmann * initial_T / (mass_of_argon * spc.electron_volt))
    
    np.random.seed()
    M = 2 * (np.random.rand(N_particles) - 1)
    N = 2 * (np.random.rand(N_particles) - 1)

    pos_array[:, 0] = np.array([v2d.vec2d(i + 10, j + 10) * sigma * 1.1 for i in range(int(np.sqrt(N_particles))) for j in range(int(np.sqrt(N_particles)))])

    if status == 'r':
        for i in range(N_particles):

            vel_x = M[i] * v_max
            vel_y = N[i] * v_max
            vel_array[i, 0] = v2d.vec2d(vel_x, vel_y)   #angstrom / femtosecondi

    else:
        for i in range(N_particles):

            vel_x = M[0] * v_max
            vel_y = M[0] * v_max
            vel_array[i, 0] = v2d.vec2d(vel_x, vel_y)

    tau = (sigma / v_max) * 1e-4        #femtosecondi

    return pos_array, vel_array, tau


##NOTE: funzione che recupera le accelerazioni relative delle particelle step per step
def get_acceleration(pos_array, N_particles, epsilon, sigma, i):

    dist_array = np.full((N_particles, N_particles), fill_value = v2d.vec2d(0, 0))
    acc_rel_array = np.full((N_particles, N_particles), fill_value = v2d.vec2d(0, 0))
    acc_list = []

    ##NOTE: loop che restiruiscono i vettori accelerazione delle particelle
    for j in range(N_particles - 1):
        for k in range(j + 1, N_particles):

            dist_array[j, k] = (pos_array[j, i] - pos_array[k, i]).abs_value()  #funzione valore assoluto implementata in vec2d
            dist_array[k, j] = dist_array[j, k]

        if dist_array[j, k].mod() <= 5 * sigma:

            acc_rel_array[j, k] = ((48 * epsilon * np.power(
                sigma, 12) / np.power(
                dist_array[j, k].mod(), 13)) - (24 * epsilon * np.power(
                sigma, 6) / np.power(dist_array[j, k].mod(), 7)
                )) / mass_of_argon * dist_array[j, k].unitary()

            acc_rel_array[k, j] = - acc_rel_array[j, k]

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
def velocity_curve_fit(vel_array, status):

    def M_B_distribution(v, C, A):

        return C * v**2 * np.exp(- A * v**2)

    velocity_list = []

    for k in range(N_particles):

        if status == 'r':
            velocity_list.append((vel_array[k, -1]).mod())

        else:
            velocity_list.append((vel_array[k, -1] - vel_array[k, 0]).mod())
    
    velocity_list = np.array(velocity_list)

    counts, bins = np.histogram(velocity_list, bins = 10)
    bins = bins[1:] - (bins[1] - bins[0]) / 2
    p, cov = spo.curve_fit(M_B_distribution, bins, counts, p0 = [1, 1], maxfev = 5000)
    x_fit = np.linspace(bins[0], bins[-1], 1000)
    y_fit = M_B_distribution(x_fit, p[0], p[1])

    return velocity_list, p[0], p[1], x_fit, y_fit


##NOTE: funzione che calcola lo spostamento quadratico medio
def get_mean_quadratic_deviation(pos_array, N_particles, N_steps):

    quadratic_deviation = 0
    mean_quadratic_deviation = []

    for i in range(N_steps):
        for k in range(N_particles):

            quadratic_deviation = quadratic_deviation + ((pos_array[k, i] - pos_array[k, 0]).mod())**2

        mean_quadratic_deviation.append(quadratic_deviation / N_particles)

    return mean_quadratic_deviation


##NOTE:
def get_diffusion_coeff(mean_quadratic_deviation, tau):

    def fit_func(B, x):

        return B * x

    x_fit = np.array([tau * i for i in range(int(N_steps / 2), N_steps)])
    y_fit = np.array(mean_quadratic_deviation[int(N_steps/2): N_steps])
    p = spo.curve_fit(fit_func, x_fit, y_fit)
    p = np.polyfit(x_fit, y_fit, 1)

    return x_fit, p[0], p[1]



##NOTE: funzione di run del programma
def run_lennard_jones2D(N_particles, N_steps, initial_T, mass_of_argon, L_box):

    what_to_plot = 'rmd' #input('Insert "t" to see the alpha particles trajectories, "v"to see the velocity distributiion or "rmd" to see the root mean deviation:\t')
    status = 'r' #input('Insert "r" for random initial velocities or "m" for mono-modal velocities:\t')

    ##NOTE: primo ciclo di inizialzzazione
    time_list = [0]
    pos_array, vel_array, tau = init_values(N_particles, N_steps, initial_T, mass_of_argon, L_box, sigma, status)
    acc = get_acceleration(pos_array, N_particles, epsilon, sigma, 0)

    ##NOTE: ciclo che gestisce il run del programma per N-1 passi
    for i in range(N_steps - 1):

        pos_array[:, i + 1] = update_pos(pos_array, vel_array, acc, N_particles, L_box, tau, i)
        acc_new = get_acceleration(pos_array, N_particles, epsilon, sigma, i + 1)
        vel_array[:, i + 1] = update_vel(vel_array, acc, acc_new, N_particles, tau, i)

        acc = acc_new
        time_list.append(tau * (i + 1))

    velocity_list, C, A, x_fit, y_fit = velocity_curve_fit(vel_array, status)

    mean_quadratic_deviation = get_mean_quadratic_deviation(pos_array, N_particles, N_steps)
    x_diff_fit, self_diffusion_par, coeff_retta = get_diffusion_coeff(mean_quadratic_deviation, tau)

    fig, ax = plt.subplots()

    if what_to_plot == 't':

        for i in range(N_particles):

            ax.plot([pos_array[i, step].x for step in range(N_steps)], [pos_array[i, step].y for step in range(N_steps)], linewidth = 2)

        ax.set_title('Lennard-Jones Interactions for Argon')
        ax.set_xlabel('x (Å)')
        ax.set_ylabel('y (Å)')

    elif what_to_plot == 'v':
        
        print('Coefficiente moltiplicativo:\t', C)
        print('Coefficiente esponenziale:\t', A)

        ax.hist(velocity_list, histtype = 'step', bins = 10, label = 'Data')
        ax.plot(x_fit, y_fit, label = 'Maxwell-Boltzmann Velocity fit', linewidth = 1)
        ax.set_xlabel('Velocity')
        ax.set_ylabel('Counts')
        ax.set_title('Distribution of the Argon Atoms Velocities')

    elif what_to_plot == 'rmd':

        print('coeff diff:', self_diffusion_par)
        coeff_retta = np.array([coeff_retta for i in range(len(x_diff_fit))])
        y_diff_fit = coeff_retta + self_diffusion_par * x_diff_fit

        ax.plot(time_list, mean_quadratic_deviation)
        ax.plot(x_diff_fit, y_diff_fit)
        ax.set_xlabel('Time')
        ax.set_ylabel('RMD')

        print('parametro di autodiffusione:\t', self_diffusion_par / 4)

    ##NOTE: il tempo di simulazione non include l'intervallo necessario a visualizzare il grafico
    print('Time taken by the simulation:', t.time() - t_s, 'seconds')

    plt.legend(frameon = False)
    plt.show()



'''Main del programma'''

## NOTE: parametri di simulazione
mass_of_argon = 39.948      #amu
epsilon = .0103             #eV
sigma = 3.46                #Angstrom
N_particles = 64
L_box = 150                 #angstrom
N_steps = 1000
initial_T = 300              #Kelvin

##NOTE: funzione di run del programma
run_lennard_jones2D(N_particles, N_steps, initial_T, mass_of_argon, L_box)


#0.423 \pm 0.003 cm^2 / s coeff diff argon in argon