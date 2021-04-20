import vec3d as v3d
import vec2d as v2d
import numpy as np

c = v2d.vec2d(1, 2)
d = v2d.vec2d(3, 4)

a = v3d.vec3d(1, 2, 3)
b = v3d.vec3d(4, 5, 6)

s = 2

angle = a.get_angle(b, 'deg')
print(angle)