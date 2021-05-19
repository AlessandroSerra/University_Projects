from re import X
import numpy as np
import sympy as sp
import vec2d as v

w = v.vec2d(1, 2)
x = sp.Symbol('x')
a = sp.limit((1 + 1/x)**x, x, sp.oo)

print(w*w)