from pylab import *

esatto = 4.4
Nmax=10000
trap_err=zeros(int((Nmax-1)/2))
simp_err=zeros(int((Nmax-1)/2))
i=0
for N in range(2,Nmax,2):
    x = linspace(0,2,N+1)
    y = x**4 - 2*x + 1
    h = abs( x[1] - x[0] )
    trap = h* ( 0.5*y[0] + 0.5*y[-1] + sum(y[1:-1]) )
    simp = h/3 *( y[0] + y[-1] + 4*sum(y[1:-1:2])+ 2*sum(y[2:-1:2]) )
    trap_err[i] = abs(trap - esatto)
    simp_err[i] = abs(simp - esatto)
    i=i+1

N_arr = range(2,Nmax,2)
loglog(N_arr, trap_err)
loglog(N_arr, simp_err)

deltaY1 = log(trap_err[350])-log(trap_err[300])
deltaX1 = log(N_arr[350])-log(N_arr[300])
slope_trap = deltaY1 / deltaX1
deltaY2 = log(simp_err[350])-log(simp_err[300])
deltaX2 = log(N_arr[350])-log(N_arr[300])
slope_simp = deltaY2 / deltaX2
print (slope_simp, slope_trap)

plot(deltaX1,deltaY1)
plot(deltaX2,deltaY2)
show()

## NOTE: facendo così siamo in grado di confrontare gli errori tra il metodo TRAPEZIOIDALE e di SIMPSON all'aumentare dei cicli
