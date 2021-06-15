import numpy as np
import matplotlib.pyplot as plt

def init_data():

    price = np.loadtxt('/Users/aleserra/Università/Fisica Computazionale/University_Projects/Esercizi/II_parziale/Trasformate di Fourier/dow_jones.txt')
    time = np.array([i for i in range(price.shape[0])])

    return time, price


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


time, price = init_data()

T = time[-1]                #non è periodica quindi prendiamo l'ultimo valore del tempo come periodo

F_trans, fast_F_trans, nu = Fourier_trans(price, T, time)
inv_F_trans, gast_inv_F_trans, new_t = Inv_Fourier_trans(F_trans, T, nu)

fig, ax = plt.subplots()
ax.plot(time, price, label = 'Original dow_jones data')
ax.plot(time, inv_F_trans, label = 'Inverse Fourier transform')
ax.set_xlabel('Time')
ax.set_ylabel('Price')
ax.legend(frameon = False)

fig1, ax1 = plt.subplots()
ax1.plot(nu, np.abs(F_trans), label = 'Fourier transform')
ax1.set_xlabel('$\\nu$')
ax1.set_ylabel('$\mathcal{F}(\\nu)$')
ax1.legend(frameon = False)

plt.show()