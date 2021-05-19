import numpy as np
import vec2d as v

a0 = v.vec2d(1,2)
b0 = v.vec2d(1,18)

a = v.vec2d(1.0003634725679749,2.0087991324596537)
b = v.vec2d(1.0006843416079425,18.00758183098231)

quad_dev1 = ((a - a0).mod())**2
quad_dev2 = ((b - b0).mod())**2

quad_dev = (quad_dev1 + quad_dev1) / 2

print(quad_dev)