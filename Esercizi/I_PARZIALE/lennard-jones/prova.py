import numpy as np
import vec2d as v

a = np.array([[v.vec2d(1, 1), v.vec2d(2, 2)], [v.vec2d(3, 3), v.vec2d(4, 4)]])
print(a)

b = []

for i in range(2):

    b.append(np.sum(a[i, :], axis = 0))

print(b[0], b[1])