import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N=100
tstep1=3000
Temp=np.zeros((N,tstep1))
Temp[35:65,0]=500.

#Temperatura iniziale
ts=1
tau=0.9*ts

#ciclo for su i    

for time in range(1,tstep1):
    for i in range(1,N-1): 
        Temp[i,time]=Temp[i,time-1]+0.5*tau/ts*(Temp[i-1,time-1]+Temp[i+1,time-1]-2*Temp[i,time-1])
        
        Temp[0,time]=Temp[1,time] 
        Temp[99,time]=Temp[98,time]

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
gridx , gridy = np.meshgrid(range(tstep1),range(N)) 
ax.plot_surface(gridx,gridy,Temp,cmap=plt.cm.coolwarm, vmax=250,linewidth=0,rstride=2, cstride=100)

plt.show()
