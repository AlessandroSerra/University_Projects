'''
l'algoritmo per la creazione di numeri random consiste nel moltiplicare u valore a per l'ultimo numero random generato 
(nel caso del primo si passa il seed) sommato ad un altro valore e poi si prende il modulo di tutto secondo un 
numero molto grande
'''

import numpy as np
import matplotlib.pyplot as plt

def rand_gen(a, c, M, r0, n = 1):       #dando un valore ad un argomento si crea un default in caso non definissimo n

    rand_list = []
    r = r0              #seed del random

    if n == 1:
        return ((a * r + c) % M ) / float(M)         #otteniamo un numero compeso tra 0 ed 1

    for i in range(n):

        r = (a * r + c) % M     #vecchio r (numero random precedente)
        rand_list.append(r / float(M))     #dividendo per M otteniamo un numero random compreso fra 0 ed 1

    return rand_list


a, c, M, r0, n = 5, 8, 104729, 9, 1000
rand_list = rand_gen(a, c, M, r0, n)

fig, ax = plt.subplots()
ax.plot(rand_list, marker = '.', linewidth = 0)

fig2, ax2 = plt.subplots()
ax2.hist(rand_list, histtype = 'step', bins = 30)

plt.show()