import numpy as np
import matplotlib.pyplot as plt

mass_of_argon = 39.948 # amu

def get_accelerations(positions):
    """
    Calculate the acceleration on each 
    particle as a  result of each other 
    particle. 
    N.B. We use the Python convention of
    numbering from 0.
    
    Parameters
    ----------
    positions: ndarray of floats
        The positions, in a single dimension, 
        for all of the particles (Å)
        
    Returns
    -------
    ndarray of floats
        The acceleration on each particle (eV/Åamu)
    """
    accel_x = np.zeros((positions.size, positions.size))
    for i in range(0, positions.size - 1):
        for j in range(i + 1, positions.size):
            r_x = positions[j] - positions[i]
            rmag = np.sqrt(r_x * r_x)
            force_scalar = lj_force(rmag, 0.0103, 3.4)
            force_x = force_scalar * r_x / rmag
            accel_x[i, j] = force_x / mass_of_argon #eV Å-1 amu-1
            # appling Newton's third law
            accel_x[j, i] = - force_x / mass_of_argon
    return np.sum(accel_x, axis=0)

accel = get_accelerations(np.array([1, 5, 10]))
print('Acceleration on particle 0 = {:.3e} eV/Åamu'.format(
    accel[0]))
print('Acceleration on particle 1 = {:.3e} eV/Åamu'.format(
    accel[1]))
print('Acceleration on particle 2 = {:.3e} eV/Åamu'.format(
    accel[2]))