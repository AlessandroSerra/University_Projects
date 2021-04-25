import numpy as np
import vec2d as v

N_particles = 3
N_steps = 2
a = np.zeros((N_particles, N_steps))
b = [1, 2, 3]

for i in range(N_steps):
    for l in range(N_particles):
        a[l, i] = b[l]

print(a)