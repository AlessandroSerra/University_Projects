import numpy as np
import matplotlib.pyplot as plt

## NOTE: creiamo un array che contenga la funzione f(x)=x^4-2x+1 per x in [0,2]
## NOTE: maggiore N, minore l'errore tra esatto e i due metodi
N = 1000
x = np.linspace(0,2,N+1)
y = x**4 - 2*x + 1

## NOTE: metodo TRAPEZIOIDALE
h = np.abs(x[1]-x[0])
## NOTE: l'indice -1 indica l'ultimo elemento dell'array
trap = 0.5*h*(y[0] + y[-1])
## NOTE: sommiamo tutti gli elementi di y dal secondo all'ultimo
trap += np.sum(y[1:-1])*h

## NOTE: metodo di SIMPSON
## NOTE: h rimane invariata
simp = y[0] + y[-1]
## NOTE: usiamo la slicing indicando questa volta anche il passo, solo i numeri pari vengono contati
simp += 4*np.sum(y[1:-1:2]) + 2*np.sum(y[2:-1:2])
simp *= h/3

## NOTE: valore esatto dell'integrale
exact = 4.4

print('esatto, trapez, simpson:', exact, trap, simp)
