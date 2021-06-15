import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

def f(t, a, b):

    return a * np.exp(-t / b) * np.sin(30 * 2 * np.pi * t) / 10

def Fourier(y, T, t):

    N = len(y)
    freq = np.array([i / T for i in range(N)])     #array delle frequenze della funzione in ingresso

    FT = np.array(
        [np.sum(y * np.exp(2j * np.pi * freq[i] * t)) for i in range(N)]       
)                                                                               #FT di una funzione y in input
    return freq, FT


N = 1000
t = np.linspace(0, 1, N)
y = f(t, 1, 0.1)
T = t[-1]
tau = t[1] - t[0]

omega = [t[i] / (2 * T * tau) for i in range(N)]
F = np.abs(fft(y))[:N//2] * tau
freq, my_F = np.abs(Fourier(y, T, t))

fig, ax = plt.subplots()
ax.plot(t, y)
ax.set_xlabel('t (s)')
ax.set_ylabel('y(t)')

fig1, ax1 = plt.subplots()
ax1.plot(omega[:N//2], F, color = 'orange')
ax1.set_xlabel('$\\nu$')
ax1.set_ylabel('$\mathcal{F})(\\nu)$')

fig2, ax2 = plt.subplots()
ax2.plot(freq[:N//2], my_F[:N//2], color = 'blue')
ax2.set_xlabel('$\\nu$')
ax2.set_ylabel('$\mathcal{F})(\\nu)$')

plt.show()