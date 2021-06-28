import matplotlib
import numpy as np 
import matplotlib.pyplot as plt

prices = np.loadtxt('/Users/aleserra/Università/Fisica Computazionale/University_Projects/Esercizi/II_parziale/Trasformate di Fourier/dow_jones.txt')
time = np.array([(i + 1) for i in range(len(prices))])

fig, ax = plt.subplots()
ax.plot(time, prices)
ax.set_title('Valore indice Dow Jones (2004 - 2008)')
ax.set_ylabel('Punti (A.U.)')
ax.set_xlabel('Time (days)')

plt.show()