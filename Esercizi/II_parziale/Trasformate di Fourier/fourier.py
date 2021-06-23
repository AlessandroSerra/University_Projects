import numpy as np
import matplotlib.pyplot as plt

def init_data():

    prices = np.loadtxt('/Users/aleserra/Università/Fisica Computazionale/University_Projects/Esercizi/II_parziale/Trasformate di Fourier/dow_jones.txt')
    time = np.array([i for i in range(len(prices))])

    return time, prices


def Fourier_trans(y, T, t):

    N = len(y)
    nu = np.array([(i - N//2) / T for i in range(N)])     #array delle frequenze della funzione in ingresso
    tau = t[1] - t[0]

    FT = np.array(
        [np.sum(y * np.exp(2j * np.pi * nu[i] * t)) for i in range(N)]   #FT di una funzione y in input
    ) * tau

    return FT, nu


def Inv_Fourier_trans(y, T, nu):

    N = len(y)
    nu_max = N / T
    t = np.array([i / nu_max for i in range(N)])
    delta_nu = nu[1] - nu[0]

    inv_FT = np.array(
        [np.sum(y * np.exp(- 2j * np.pi * t[i] * nu)) for i in range(N)]
    ) * delta_nu

    y2 = y.copy()
    for i in range(N):
        if np.abs(nu[i]) > 0.02 * nu_max:
            y2[i] = 0

    inv_FT2 = np.array(
        [np.sum(y2 * np.exp(- 2j * np.pi * t[i] * nu)) for i in range(N)]
    ) * delta_nu

    return inv_FT, inv_FT2, t


##NOTE: funzione di run del programma
def run_Dow_Jones_Analisys(time, prices):

    N = len(prices)
    T = time[-1]                #non è periodica quindi prendiamo l'ultimo valore del tempo come periodo

    windowed_prices = np.hanning(N) * prices
    F_trans, nu = Fourier_trans(prices, T, time)
    inv_F_trans, inv_F_trans2, new_t = Inv_Fourier_trans(F_trans, T, nu)

    fast_nu = np.fft.rfftfreq(N)
    fast_F_trans = np.fft.rfft(prices)
    fast_inv_F_trans = np.fft.irfft(fast_F_trans)

    fig, ax = plt.subplots()
    ax.plot(time, prices, label = 'Original data')
    ax.plot(new_t[1:], np.abs(inv_F_trans)[1:], label = 'Inverse Fourier transform')
    #ax.plot(time, fast_inv_F_trans, label = 'Fast inverse Fourier transform')
    ax.set_xlabel('Time')
    ax.set_ylabel('Prices')
    ax.set_ylim(9000, 14500)
    ax.legend(frameon = False)

    fig1, ax1 = plt.subplots()
    ax1.plot(nu, np.abs(F_trans), label = 'Fourier transform')
    #ax1.plot(fast_nu, np.abs(fast_F_trans), label = 'Fast Fourier transform')
    ax1.set_xlabel('$\\nu$')
    ax1.set_ylabel('$\mathcal{F}(\\nu)$')
    ax1.set_yscale('log')
    ax1.set_xlim(0, nu[-1])
    ax1.legend(frameon = False)

    fig2, ax2 = plt.subplots()
    ax2.plot(time, prices, label = 'Original data')
    ax2.plot(new_t[10:-5], np.abs(inv_F_trans2)[10:-5], label = 'Inverse Fourier transform, only first 2%')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Prices')
    ax2.set_ylim(9000, 14500)
    ax2.legend(frameon = False)

    plt.show()


'''
Main del programma
'''

time, prices = init_data()
run_Dow_Jones_Analisys(time, prices)