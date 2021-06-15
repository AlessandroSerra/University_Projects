#ESERCIZIO 2: TRASFORMAZIONI DI FOURIER
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as spc 

#FUNZIONI
def Fourier(data, T, t):
    N = len(data)
    nu = np.array([(i-N//2)/T for i in range(N) ]) #frequenze, create con un ciclo for, array
    F = np.array( 
        [np.sum( data * np.exp(2j * np.pi *nu[i] * t)) for i in range(N) ]  # j è l'unità immaginaria
        )
    delta_t = t[1] - t[0]
    return nu, F* delta_t

def invFourier(data, nu_max, nu):
    N = len(data)
    t = np.array([i/nu_max for i in range(N)])

    invF = np.array(
        [np.sum( data * np.exp(- 2j * np.pi * nu * t[i])) for i in range (N) ]
        )
    delta_nu = nu[1]- nu[0]
    return t, invF * delta_nu



#TRASFORMATA
with open("/Users/aleserra/Università/Fisica Computazionale/University_Projects/Esercizi/II_parziale/Trasformate di Fourier/dow_jones.txt", "r") as f:
    data = list(map(float, f.readlines()))

t=np.arange(0,1024,1)
T= 1024

nu, F = Fourier(data, T , t)   
nu1, F1 =  np.fft.rfftfreq(len(t)), np.fft.rfft(data)

#INVERSA TRASFORMATA
j = int(0.02 * len(data))
Finv = F.copy()
for i in range(j, len(Finv)):
    Finv[i] = 0

tinv, InvF = invFourier(Finv, len(t)/T , nu)
FFinv = np.fft.irfft(Finv) #trasformata inversa reale 

t1=np.arange(0,1023,0.5)

#GRAFICI

fig1,ax1 = plt.subplots()
ax1.plot(nu, np.abs(F), label= "Trasformata", color = "red")
ax1.plot(nu1, np.abs(F1) ,label ="Fast Fourier", color = "blue")
ax1.set_xlabel("$\\nu$")
ax1.set_ylabel("$\mathcal{F}(\\nu)$")
plt.title("Trasformate")
plt.legend()

fi2,ax2, = plt.subplots()  
ax2.plot(tinv, abs(InvF), label = "Inversa Trasformata", color = "red" )
ax2.plot( t1 ,abs(FFinv) , label = "Fast Fourier inversa", color="orange")
ax2.plot(t , data , label = "dati", color="blue")
plt.title("Trasformate Inverse e Dati")
ax2.set_xlabel("x (s)")
ax2.set_ylabel("y (t)")
plt.legend()

plt.show()
