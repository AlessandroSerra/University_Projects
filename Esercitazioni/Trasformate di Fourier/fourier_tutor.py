import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

def Fourier(y, T, t):

    N = len(y)
    freq = np.array([(i - N//2) / T for i in range(N)])     #array delle frequenze della funzione in ingresso
    tau = t[1] - t[0]

    FT = np.array(
        [np.sum(y * np.exp(2j * np.pi * freq[i] * t)) for i in range(N)]   #FT di una funzione y in input
    ) * tau

    fast_FT = fft(y) * tau

    return freq, FT, fast_FT

def InvFourier(y, T, nu):

    N = len(y)
    nu_max = N / T
    t = np.array([i / nu_max for i in range(N)])
    delta_nu = nu[1] - nu[0]

    inv_FT = np.array(
        [np.sum(y * np.exp(- 2j * np.pi * t[i] * nu)) for i in range(N)]
    ) * delta_nu

    fast_IFFT = ifft(y) * delta_nu

    return t, inv_FT, fast_IFFT


prices = np.loadtxt('/Users/aleserra/Università/Fisica Computazionale/University_Projects/Esercizi/II_parziale/Trasformate di Fourier/dow_jones.txt')
y = prices
N = len(prices)
t = np.array([i for i in range(N)])
T = t[-1]

nu, FT, fast_FFT = Fourier(y, T, t)
new_t, inv_FT, fast_IFFT = InvFourier(FT, T, nu)

fig, ax = plt.subplots()
ax.plot(t, y, label = 'Original data')
ax.plot(new_t, inv_FT.real, label = 'Inverse Fourier transform')
#ax.plot(new_t, fast_IFFT.real, label = 'Fast Inverse Fourier transform')
ax.set_xlabel('t(s)')
ax.set_ylabel('f(t)')
ax.legend()

fig1, ax1 = plt.subplots()
ax1.plot(nu, np.abs(FT), label = 'Fourier transform')
#ax1.plot(nu, np.abs(fast_FFT), label = 'Fast Fourier transform')
ax1.set_xlabel('$\\nu (Hz)$')
ax1.set_ylabel('$\mathcal{F}(\\nu)$')
ax1.legend()


plt.show()