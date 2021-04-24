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
energy_list = np.array([0.1374, -0.0195, -0.0218, -0.0133, -0.0076, -0.0043, -0.0025])
energy_err_list = energy_list * 0.1

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

p, cov = spo.curve_fit(lj_potential, radius_list, energy_list, sigma = energy_err_list)
print('best value for sigma: ', p[0])
print('best value for epsilon: ', p[1])

fig, ax = plt.subplots()
ax.errorbar(radius_list, energy_list, yerr = energy_err_list, marker = 'o', ls = '', label = 'exp values lj potential')
ax.set_xlabel('r/A')
ax.set_ylabel('E/eV')

x_list = np.linspace(3.5, 7, 1000)
ax.plot(x_list, lj_potential(x_list, p[0], p[1]), label = 'fit lj potential')

plt.legend()
plt.show()