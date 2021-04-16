import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo

## NOTE: definiamo una funzione che si adatti ad un fit NON lineare
def f(x, A, tau, omega):
    return A * np.exp(-x / tau) * np.cos(omega * x)

x = np.linspace(0, 10, 100)
tau = 5A = 1
omega =2 * np.pi

y = f(x, A, tau, omega)

plt.sublìplots()