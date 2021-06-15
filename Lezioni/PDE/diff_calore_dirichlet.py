import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N=100
tstep=1000
T=np.zeros((N,tstep))
#Temperatura iniziale
T[50,0]=10.
ts=1
tau=0.9*ts

for time in range(1,tstep):
    for i in range(1,N-1): 
        T[i,time]=T[i,time-1]+0.5*tau/ts*(T[i-1,time-1]+T[i+1,time-1]-2*T[i,time-1])


fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
gridx , gridy = np.meshgrid(range(0,200),range(0,N,2)) 
#ax.plot_wireframe(gridx,gridy,T[::2,0:200])
ax.plot_surface(gridx,gridy,T[::2, 0:200],cmap=plt.cm.coolwarm, vmax=25,linewidth=0,rstride=2, cstride=100)


plt.show()
