import matplotlib.pyplot as plt #questo in modo da non avere funzioni chiamate uguali

#abbiamo una certa serie numerica e voglis conoscerne la somma, approssimiamo la somma come
#la somma parziale sino ad un termine n sufficientemente grande

def geom(x,a,n):        #funzione che restituiscce il termine n-esimo della serie
    return a * x**n

x = 0.9                 #valore nella x nella serie gemoetrica (scelto da noi)
a = 1                   #vediamo la somma dei primi 100 elementi
exact = a / (1-x)
N_values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200]
errors = []             #lista vuota per gli errori

for N in N_values:
    sum = 0
    for n in range(N+1):
        sum += geom(x, a, n)
    errors.append(exact-sum)

fig, ax = plt.subplots()        #la funzione subplots restituisce una finestra con i dati e gli assi x ed y
ax.plot(N_values, errors)       #ax.plot indica quali sono i valori da plottare sugli assi
ax.set_yscale('log')            #indica la scala scelta per l'asse y
ax.set_xlabel('N')
ax.set_ylabel('$\>Delta S$')    #si può utilizzare la sintassi di LaTeX per usare simboli ed altro

plt.show()             #funzione null per visualizzare il grafico
