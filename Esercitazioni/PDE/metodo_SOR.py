import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

##NOTE: funzione che verifica la consizione al bordo
def fuori_bordo(i, j, N, bordo):

    i, j = int(N / 2), int(N / 2)
    
    return ((i - c_x)**2 + (j - c_y)**2) > bordo

##NOTE: funzione che calcola un singolo punti della griglia conoscendo quelli adiacenti, \ manda a capo la riga
def SOR_alg(V, rho, h, w, N, i, j, bordo):

    if fuori_bordo(i, j, N, bordo):
        return 0

    new_V = w * (V[i+1, j] + V[i-1, j] + V[i, j+1] + V[i, j-1] + h**2 * rho[i, j]) / 4\
                    + (1 - w) * V[i, j]

    return new_V

##NOTE: funzine che esegue uno step dell'algoritmo SOR
def step(V, rho, h, w, N, bordo):

    for i in range(1, N - 1):
        for j in range(1, N - 1):

            V[i, j] = SOR_alg(V, rho, h, w, N, i, j, bordo)

    return V

##NOTE: funzione che verifichi quando si è arrivati a convergenza
def conver_cond(diff, soglia):

    if diff < soglia:
        return True
    
    return False

##NOTE: funzione per iterare il SOR sino ad N step, il parametro soglia indica il livello di convergenza che vogliamo
def run_SOR(V, rho, h, w, soglia, bordo):

    diff = 10
    N = len(V)

    while not conver_cond(diff, soglia):

        old_V = V.copy()        #copia di un array in un altro
        V = step(V, rho, h, w, N, bordo)
        diff = np.max(np.abs(V - old_V)) / (np.max(rho) * h**2)     #differenza percentuale

    return V


##NOTE: costanti
N = 100              #punti su x e y (griglia quadrata)
L = 1               #m, lunghezza griglia
h = L / N           #m, lunghezza passo di integrazione
w = 1.8             #peso del metodo SOR (parametro di rilassamento)
soglia = 1e-3       #condizione di convergenza del SOR
bordo = N / 2

V = np.zeros([N, N])        #array del potenziale
rho = np.zeros([N, N])      #array della densità di carica

##NOTE: condizioni iniziali
c_x, c_y = int(N / 2), int(N / 2)
rho[c_x, c_y] = 1e-3 / h**2                #densità di carica iniziale

V = run_SOR(V, rho, h, w, soglia, bordo)
x, y = np.linspace(-L / 2, L / 2, N), np.linspace(-L / 2, L / 2, N)
x, y = np.meshgrid(x, y)

##NOTE: plottiamo in 3d il potenziale
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot_wireframe(x, y, V)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('V')

fig1, ax1 = plt.subplots()
ax1.contour(x, y, V, levels = 20)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.axis('equal')

plt.show()