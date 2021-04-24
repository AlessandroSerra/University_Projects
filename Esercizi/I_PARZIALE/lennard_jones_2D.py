import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import vec2d as v2d
import time 

t_s = time.time()

## NOTE: costanti della simulazione
mass_of_argon = 39.948      #amu
epsilon = .0103             
sigma = 3.4                 #Angstron
L_box = 20                  #Angstron

## NOTE: variabili della simulazione
#N_particles = int(input('Inserire il numero di particelle da far interagire\n'))
N_particles = 5
N_steps = 10000
tau = .1                    #femtosecondi
initial_T = 300             #Kelvin

## NOTE: funzione che inizializza le posizioni degli atomi
def init_pos(N_particles):

    r_init[i] = v2d.vec2d(1, 20 * (2 * np.random.rand() - 1))
    
    return r_init

def init_vel(N_particles, initial_T):
    
    P = np,random.rand(N_particles) - 0.5
    R = np.random.rand(N_particles) - 0.5
    
    v_init[i] = v2d.vec2d(P * np.sqrt(spc.Boltzmann * initial_T / 
            (mass_of_argon / spc.e)), R * np.sqrt(spc.Boltzmann * initial_T / 
            (mass_of_argon / spc.e)))
    
    return v_init

# NOTE: definiamo la forza di interazione come la derivata parziale opposta del potenziale
def lj_force(r, epsilon, sigma):

    f_scalar = 48 * epsilon * np.power(
        sigma, 12) / np.power(
        r.mod(), 13) - 24 * epsilon * np.power(
        sigma, 6) / np.power(r.mod(), 7)

    return f_scalar * r

def update_pos(r, v, tau):

    new_r = r + v * tau + lj_force(r, epsilon, sigma) * tau**2 / 2
    return new_r

def update_vel(v, a_part, anew_part, tau):
    new_v = v + (a_part + anew_part) * tau**2 / 2

def run_md(tau, N_steps, initial_T, ):
    dnf