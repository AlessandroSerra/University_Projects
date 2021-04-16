cutoff = 15 

def lj_force_cutoff(r, epsilon, sigma):
    """
    Implementation of the Lennard-Jones potential 
    to calculate the force of the interaction which 
    is considerate of the cut-off.
    
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
        Force of the van der Waals interaction (eV/Å)
    """
    if r < cutoff:
        return 48 * epsilon * np.power(
            sigma / r, 13) - 24 * epsilon * np.power(
                sigma / r, 7)
        
    else:
        return 0