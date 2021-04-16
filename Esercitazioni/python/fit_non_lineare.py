import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo

## NOTE: definiamo una funzione che si adatti ad un fit NON lineare
def f(x, A, tau, omega):
    return A * np.exp(-x / tau) * np.cos(omega * x)

x = np.linspace(0, 10, 100)
x_new = np.linspace(0, 10, 1000)
tau = 5
A = 1
omega =2 * np.pi / 2
sigma = .1

## NOTE: loc indica l'eventuale decentramento della distribuzione normale
y = f(x, A, tau, omega) + np.random.normal(loc = 0, scale = sigma, size = 100)

fig, ax = plt.subplots()
ax.errorbar(x, y, yerr = sigma * np.ones(100), elinewidth = 1, linewidth = 0, marker = '.', label = 'data')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.legend()

## NOTE:passiamo alla funzinìone quelli che crediamo essere parametri possibilie gli errori

p, cov = spo.curve_fit(f, x, y,  p0 = [A, tau, omega], sigma = sigma * np.ones(100))
ax.plot(x_new, f(x_new, p[0], p[1], p[2]), label = 'fit')

plt.show()
