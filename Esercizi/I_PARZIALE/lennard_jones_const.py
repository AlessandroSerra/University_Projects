'''
programma per fittare il potenziale doi lennard-jones al fine di ottenere le
migliori stime possibili delle costanti epsilon e sigma
'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo

'''
the following cod e should calculate the energy values of the 
lj potential in order to get some valuable points to later
implement the lj fit to get the constants
'''

radius_list = np.arange(3.5, 7, 0.5)
energy_list = np.array()
print(radius_list)

def lj_potential(r, epsilon, sigma):
    '''
    Implementation of the Lennard-Jones potential 
    to calculate the energy of the interaction.
    
    Parameters
    ----------
    r: float
        Distance between two particles (Å)
    epsilon: float 
        Potential energy at the equilibrium bond 
        length (eV)
    sigma: float 
        Distance at which the potential energy is 
        zero (Å)
    
    Returns
    -------
    float
        Energy of the van der Waals interaction (eV)
    '''

    return 4 * epsilon * np.power(sigma / r, 12) - 4 * epsilon * np.power(sigma / r, 6)

'''
the scipy function curve_fit should return the best value that fits the
lj_potential curve
'''

p, cov = spo.curve_fit(lj_potential, r, epsilon, sigma)