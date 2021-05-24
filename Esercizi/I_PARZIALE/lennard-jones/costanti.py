import numpy as np
import scipy.constants as spc

massa_argon = 39.984
temperatura = 300
particelle = 10 

A = massa_argon / (2 * spc.Boltzmann * temperatura)
C = 4 * np.pi * particelle * np.power(massa_argon / (2 * np.pi * spc.Boltzmann * temperatura), 1.5)

print('A=\t', A)
print('C=\t', C)