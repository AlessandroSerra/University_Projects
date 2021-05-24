import numpy as np
import matplotlib.pyplot as plt

def attractive_energy(r, epsilon, sigma):
    """
    Attractive component of the Lennard-Jones 
    interactionenergy.
    
    Parameters
    ----------
    r: float
        Distance between two particles (Å)
    epsilon: float 
        Negative of the potential energy at the 
        equilibrium bond length (eV)
    sigma: float 
        Distance at which the potential energy is 
        zero (Å)
    
    Returns
    -------
    float
        Energy of attractive component of 
        Lennard-Jones interaction (eV)
    """
    return -4 * epsilon * np.power(sigma / r, 6)

def repulsive_energy(r, epsilon, sigma):
    """
    Repulsive component of the Lennard-Jones 
    interactionenergy.
    
    Parameters
    ----------
    r: float
        Distance between two particles (Å)
    epsilon: float 
        Negative of the potential energy at the 
        equilibrium bond length (eV)
    sigma: float 
        Distance at which the potential energy is 
        zero (Å)
    
    Returns
    -------
    float
        Energy of repulsive component of 
        Lennard-Jones interaction (eV)
    """
    return 4 * epsilon * np.power(sigma / r, 12)

def lj_energy(r, epsilon, sigma):
    """
    Implementation of the Lennard-Jones potential 
    to calculate the energy of the interaction.
    
    Parameters
    ----------
    r: float
        Distance between two particles (Å)
    epsilon: float 
        Negative of the potential energy at the 
        equilibrium bond length (eV)
    sigma: float 
        Distance at which the potential energy is 
        zero (Å)
    
    Returns
    -------
    float
        Energy of the Lennard-Jones potential 
        model (eV)
    """
    return repulsive_energy(
        r, epsilon, sigma) + attractive_energy(
        r, epsilon, sigma)

x = np.linspace(2.9, 5, 100)
y = np.linspace(-.10, .20, 100)

r = np.linspace(3, 5, 100)
plt.plot(r, attractive_energy(r, 0.0103, 3.4),
         label='Attractive')
plt.plot(r, repulsive_energy(r, 0.0103, 3.4), 
         label='Repulsive')
plt.plot(r, lj_energy(r, 0.0103, 3.4), 
         label='Lennard-Jones')
plt.plot(x, [0 for i in range(len(x))], color = 'black', linewidth = .6)
plt.plot([3 for i in range(len(y))], y, color = 'black', linewidth = .6)
plt.xlabel(r'$r$/Å')
plt.ylabel(r'$V_{LJ}$/eV')
plt.legend(frameon=False)
plt.show()