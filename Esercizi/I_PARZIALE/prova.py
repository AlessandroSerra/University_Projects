import numpy as np
import vec2d as v

N_particles = 3
N_steps = 3
temp_list = []

a = np.full((N_particles, N_steps), fill_value = v.vec2d)

for i in range(N_particles):
    for j in range(N_steps):

        a[i, j] = v.vec2d(i, j)

temp_list.append(a[1, 1])


print(len(temp_list))

print(temp_list)