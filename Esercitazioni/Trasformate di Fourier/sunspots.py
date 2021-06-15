import numpy as np
import matplotlib.pyplot as plt

def Fourier_trans(y, T, t):

    N = len(y)
    nu = np.array([i / T for i in range(N)])
    tau = t[1] - t[0]

    FT = np.array(
        [np.sum(y * np.exp(2j * np.pi * nu[i] * t)) for i in range(N)]
    ) * tau

    fast_FT = np.fft.rfft(y) * tau

    return FT, fast_FT, nu


def Inv_Fourier_trans(y, T, nu):

    N = len(y)
    nu_max = N / T
    new_t = np.array([i / nu_max for i in range (N)])
    delta_nu = nu[1] - nu[0]

    inv_FT = np.array(
        [np.sum(y * np.exp(2j * np.pi * new_t * nu[i])) for i in range(N)]
    ) * delta_nu

    fast_IFT = np.fft.irfft(y) * delta_nu

    return inv_FT, fast_IFT, new_t


t = []
y = []

with open('/Users/aleserra/Università/Fisica Computazionale/University_Projects/Esercitazioni/Trasformate di Fourier/sunspots.txt', 'r') as file:
    data = file.readlines()

for line in data:

    current_t, current_y = line.split()         #divide gli elementi di una riga
    t.append(float(current_t))
    y.append(float(current_y))

t, y = np.array(t), np.array(y)
T = t[-1]
nu, F = np.fft.rfftfreq(len(t)), np.fft.rfft(y)
new_y = np.fft.irfft(F)

#FT, fast_FT, my_nu = Fourier_trans(y, T, t)
#inv_FT, fast_IFT, new_t = Inv_Fourier_trans(y, T, nu)

fig, ax = plt.subplots()
ax.plot(t, y, label = 'Original data')
ax.plot(t[:-1], new_y, label = 'Fast inverse Fourier transform')
ax.set_xlabel('t (A.U.)')                   #Arbitrary Units
ax.set_ylabel('y (A.U.)')
ax.legend(frameon = False)

fig1, ax1 = plt.subplots()
ax1.plot(nu, np.abs(F), label = 'Fast Fourier transform')
ax1.set_xlabel('$\\nu$')
ax1.set_ylabel('$\mathcal{F}(\\nu)$')
ax1.legend(frameon = False)

plt.show()