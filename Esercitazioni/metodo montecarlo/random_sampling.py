import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return np.exp(-x**2 / 8)        #gaussiana con ampiezza 1 e sigma 2

def rand_range(min, max, n = 1):

    return (max - min) * np.random.rand(n) + min

def random_sampling(f, n, min, max):

    x = rand_range(min, max, n)
    y = f(x)
    mean = np.mean(y)

    return x, y, mean * (max - min)

min, max, n = -10, 10, 1000
analytic = np.sqrt(8 * np.pi)
x, y, numerical = random_sampling(func, n, min, max)
x_ord = np.linspace(min, max, n)
y_ord = func(x_ord)

print('analitico:\t', analytic)
print('numerico:\t', numerical)
print('diff percentuale:\t', (numerical - analytic) / analytic)

fig, ax = plt.subplots()
ax.plot(x, [0 for i in range(len(x))])
ax.plot(x, func(x), marker = '.', linewidth = 0, markersize = 1)
ax.plot(x_ord, y_ord, linewidth = 1)

plt.show()