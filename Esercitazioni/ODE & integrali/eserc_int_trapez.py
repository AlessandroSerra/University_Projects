import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

## NOTE: definiamo la funzione da integrare
def f(x):
    return np.exp(-x**2)

def trap(y, h):
    return h/2 * (y[0]+y[-1]) + h * np.sum(y[1:-1])

def simp(y, h):
    return h/3 * (y[0]+y[-1] + np.sum(4 * y[1:-1:2]) + np.sum(2 * y[2:-1:2]))

## NOTE: definiamo l'intervallo di integrazione
int = [0,1]

## NOTE: modo veloce di definire una lista
N_list = [100 + i * 100 for i in range(20)]
exact = sp.erf(int[-1]) * np.sqrt(np.pi)/2

diff_trap = []
diff_simp = []

for N in N_list:
    x = np.linspace(int[0],int[1],N+1)
    y = f(x)
    h = x[1]-x[0]
    int_trap = trap(y,h)
    int_simp = simp(y,h)
    diff_trap.append(exact - int_trap)
    diff_simp.append(exact - int_simp)

print('trapezi', int_trap)
print('simpson', int_simp)
print('esatto', exact)

## NOTE: la funzione subplot restituisce una figura ed un asse
fig, ax = plt.subplots()
ax.plot(N_list, diff_trap, label = 'trapezi')
ax.plot(N_list, diff_simp, label = 'simpson')
ax.legend()
ax.set_yscale('log')
ax.set_xscale('log')

## NOTE: potevamo anche direttamente fare ax.loglog() per avere il grafico in scala logar.

plt.show()
