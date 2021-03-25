from pylab import *
ioff()

x = linspace(0,2*pi,10)       #linspace prende anche l'elemento n mentre arange sino ad n-1
y1 = sin(x)

plot(x,y1)
show()
savefig('grafic_prova.png',dpi=500)
