from pylab import *

#Costanti moto
GM = 4*pi**2 # G*M in unita astronomiche (UA)
m=5.9722E+24 # massa terra (kg)

#Posizioni e velocita' iniziali in UA, al perielio
x=0.98
y=0.0
r=sqrt(x**2+y**2)
vx=0
vy=2*pi

# Nstep*tau e' il tempo totale  in anni
Nstep=1000 # Numero di passi temporali
tau=0.001 # lunghezza singolo passo temporale

#arrays posizione e velocit
x_arr=zeros(Nstep+1)
y_arr=zeros(Nstep+1)
vx_arr=zeros(Nstep+1)
vy_arr=zeros(Nstep+1)

#Salvataggio delle condizioni iniziali
x_arr[0]=x
y_arr[0]=y
vx_arr[0]=vx
vy_arr[0]=vy

for istep in range(Nstep): 
    Fx=-GM*x/r**3 
    Fy=-GM*y/r**3
  #Eulero-Cromer
  #vx+=tau*Fx
  #vy+=tau*Fy
  #x+=tau*vx
  #y+=tau*vy
  #r=sqrt(x**2+y**2)
  #velocity-Verlet
    x+=vx*tau+0.5*Fx*tau**2
    y+=vy*tau+0.5*Fy*tau**2
    r=sqrt(x**2+y**2)
    F1x=-GM*x/r**3
    F1y=-GM*y/r**3
    vx+=0.5*tau*(Fx+F1x)
    vy+=0.5*tau*(Fy+F1y)
    x_arr[istep+1]=x
    y_arr[istep+1]=y
    vx_arr[istep+1]=vx
    vy_arr[istep+1]=vy

xlabel("x (AU)") # label dell'asse X
ylabel("y (AU)") # label dell'asse Y   
plot(x_arr,y_arr)
axis("equal")

show()