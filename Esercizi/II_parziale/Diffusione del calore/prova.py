import numpy as np


x = np.random.randint(0,10, size = 5)
y = np.random.randint(0,10, size = 5)
z = np.random.randint(0,10, size = 5)
T = np.array((x, y, z))
print(T)

A = np.amax(T, axis = 1)
print(A)