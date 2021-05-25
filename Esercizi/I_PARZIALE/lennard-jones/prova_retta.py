import numpy as np
import matplotlib.pyplot as plt

x = np.array([0, 1, 2, 3, 4, 5])
a = 2
q = np.array([1 for i in range(len(x))])

y = q + a * x

print(x + q)