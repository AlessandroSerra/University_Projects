import numpy as np
import matplotlib.pyplot as plt
import time as t

t_s = t.time()


##NOTE: funzione per leggere ed inizializzare i dati
def init_data():

    prices = np.loadtxt('/Users/aleserra/Università/Fisica Computazionale/University_Projects/Esercizi/II_parziale/Trasformate di Fourier/dow_jones.txt')
    time = np.array([(i + 1) for i in range(len(prices))])

    return time, prices


#NOTE: funzione che definisce la FT discreta
def Fourier_trans(y, T, t):

    N = len(y)
    nu = np.array([(i - N//2) / T for i in range(N)])     #array delle frequenze della funzione in ingresso
    tau = t[1] - t[0]                                     #intervallino temporale (1 giorno)

    FT = np.array(
        [np.sum(y * np.exp(2j * np.pi * nu[i] * t)) for i in range(N)]   #FT di una funzione y in input
    ) * tau

    return FT, nu


##NOTE: funzione che definisce la anti-FT discreta
def Inv_Fourier_trans(y, T, nu):

    N = len(y)
    nu_max = N / T                                  #frequenza di campionamento
    t = np.array([i / nu_max for i in range(N)])
    delta_nu = nu[1] - nu[0]                        #intervallino di frequenza

    inv_FT = np.array(
        [np.sum(y * np.exp(- 2j * np.pi * t[i] * nu)) for i in range(N)]    #anti_FT di una funzione y in input
    ) * delta_nu

    y2 = y.copy()               #funzione per copiare l'array y in y2

    ##NOTE: segmento di codice per mantenere diversi da 0 solo i valori della FT relativi al primo 2% delle frequenze
    for i in range(N):
        if np.abs(nu[i]) > 0.01 * nu_max:       #primo 2% delle frequenze sia positive che negative
            y2[i] = 0

    inv_FT2 = np.array(
        [np.sum(y2 * np.exp(- 2j * np.pi * t[i] * nu)) for i in range(N)] #anti-FT del nuovo array y2
    ) * delta_nu

    return inv_FT, inv_FT2, t, y2


##NOTE: funzione di run del programma
def run_Dow_Jones_Analisys(time, prices):

    N = len(prices)
    T = time[-1]                #il periodo non è noto quindi si prende l'ultimo valore tel tempo come periodo

    F_trans, nu = Fourier_trans(prices, T, time)
    inv_F_trans, inv_F_trans2, new_t, y2 = Inv_Fourier_trans(F_trans, T, nu)

    #NOTE: segmento di codice relativo all'mplementazione delle fft di numpy
    fast_nu = np.fft.rfftfreq(N)
    fast_F_trans = np.fft.rfft(prices)
    fast_inv_F_trans = np.fft.irfft(fast_F_trans)

    #NOTE: rappresentazione grafica dei dati originali, della anti-FT e della fast anti-FT
    fig, ax = plt.subplots()
    ax.plot(time, prices, label = 'Original data')
    ax.plot(new_t[1:], np.abs(inv_F_trans)[1:], linewidth = .3, label = 'Inverse Fourier transform')
    ax.plot(time, np.abs(fast_inv_F_trans), linewidth = .3, label = 'Fast inverse Fourier transform')
    ax.set_title('Anti Trasformate di Fourier delle Trasformate complete')
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Prices ($)')
    ax.set_ylim(9000, 14500)
    ax.legend(frameon = False)

    #NOTE: rappresentazione grafica della FT e della fast FT
    fig1, ax1 = plt.subplots()
    ax1.plot(nu, np.abs(F_trans), label = 'Fourier transform')
    ax1.plot(fast_nu, np.abs(fast_F_trans), linewidth = .8, label = 'Fast Fourier transform')
    ax1.set_title('Trasformate di Fourier dei dati in input')
    ax1.set_xlabel('$\\nu (days^{-1})$')
    ax1.set_ylabel('$log(\mathcal{F}(\\nu))$')
    ax1.set_yscale('log')                   #scala logaritmica per visualizzare meglio i valori
    ax1.set_xlim(0, nu[-1])
    ax1.legend(frameon = False)

    #NOTE: rappresentazione grafica dei dati originali, della anti-FT e della fast anti-FT usando solo il 2% dei valori
    fig2, ax2 = plt.subplots()
    ax2.plot(time, prices, label = 'Original data')
    ax2.plot(new_t, np.abs(inv_F_trans2), label = 'Inverse Fourier transform, only first 2%')
    ax2.set_title('Anti Trasformate di Fourier delle Trasformate ridotte')
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Prices ($)')
    ax2.set_ylim(9000, 14500)
    ax2.legend(frameon = False)

    print('TIme taken by the simulation:\t', t.time() - t_s)

    plt.show()


'''
Main del programma
'''

time, prices = init_data()
run_Dow_Jones_Analisys(time, prices)