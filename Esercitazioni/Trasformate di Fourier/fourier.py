import numpy as np
import matplotlib.pyplot as plt

def Fourier(y, T, t):

    N = len(y)
    freq = np.array([i / T for i in range(N)])     #array delle frequenze della funzione in ingresso

    FT = np.array(
        [np.sum(y * np.exp(2j * np.pi * freq[i] * t)) for i in range(N)]       
)                                                                               #FT di una funzione y in input
    return freq, FT


T, N = 1, 100                       #periodo di 1 secondo, 100 divisioni dell'intervallo di tempo
t = np.linspace(0, 1, 1000)         #secondi, freq massima = 100 Hz
y = np.sin(2 * np.pi * 50 * t)      #usiamo una frequenza di 50 Hz

freq, FT = Fourier(y, T, t)

fig, ax = plt.subplots()
ax.plot(t, y)
ax.set_xlabel('t (s)')
ax.set_ylabel('y(t)')

fig1, ax1 = plt.subplots()
ax1.plot(freq[:N//2], np.abs(FT)[:N//2])
ax1.set_xlabel('$\\nu$')
ax1.set_ylabel('$\mathcal{F})(\\nu)$')

plt.show()