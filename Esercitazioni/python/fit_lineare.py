import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo

## NOTE: fit lineare di dati che si adattano ad una retta

def retta(x, A, B):
    return A * x + B

A = 2
B = 1
sigma = .1
x = np.linspace(-1, 1, 20)
y = retta(x, A, B) + np.random.normal(loc = 0, scale = sigma, size = 20)    #genera numeri random normalmente distribuiti
weights = 1 / sigma * np.ones(20)

fig, ax = plt.subplots()
ax.errorbar(x, y, yerr = sigma, elinewidth = 1, linewidth = 0, marker = '.', label = 'data')
ax.set_xlabel('x')
ax.set_ylabel('y')

p, cov = np.polyfit(x, y, 1, w = weights, cov = True)     #la funzione fitta i dati x ed y con un polinomio di grado n (lineare = 1), w sta per pesi
ax.plot(x, retta(x, p[0], p[1]), label = 'fit')
print(p[0], p[1])
print(cov)

plt.legend()
plt.show()
