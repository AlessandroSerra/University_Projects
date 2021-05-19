import numpy as np
import matplotlib.pyplot as plt
import scipy.constants  as spc

"""Inizialmente nel programmma ho definito le costanti utili alla simulazione"""

mass_of_argon = 39.948 # Massa dell'argon in unità di massa atomiche
Nsteps=100
Nparts=64 #numero di atomi di Argon 
L=20 #amstrong dimensioni box
dt = 0.1 #passo temporale 
T=300 #K kelvin
sigma=3.46 #Angstrong 10 alla - 10 metri Distanza in cui la potenziale energia è zero 
epsilon = 0.0103 #eV è l'energia potenziale alla separazione di equilibrio


"""Definisco una funzione che calcola il potenziale di Lennard-Jones per calcolare la forza 
dell'interazione. """
def lj_force(r, epsilon, sigma):
    
    if r.all()<=5*sigma:
        
        return 48 * epsilon * np.power(
            
            sigma, 12) / np.power(
                
            r, 13) - 24 * epsilon * np.power(
                
            sigma, 6) / np.power(r, 7)
                
    return 0 



""" init_position è una funzione che inizializza le posizioni delle particelle
di argon """
def init_positions(L,Nparts):
    
    pos_array_x=np.zeros((Nsteps,Nparts))
    
    pos_array_y=np.zeros((Nsteps,Nparts))
    
    pos_x=0
    
    pos_y=np.linspace(0,L,Nparts)
    
    for i in range(Nparts):
        
        pos_array_y[0][i]=pos_y[i]
        
    return pos_array_x,pos_array_y    



"""La seguente funzione inizializza le velocità """      
def init_velocity(T, number_of_particles,Nsteps):
    
    v=np.sqrt((spc.Boltzmann / spc.electron_volt)*T *2 /mass_of_argon)
    
    velox_x=np.zeros([Nsteps,Nparts])
    
    velox_y=np.zeros([Nsteps,Nparts])
    
    np.random.seed()
    
    for i in range(Nparts):
        
       theta=2 *np.pi * np.random.rand()
       
       vx=v*np.cos(theta)
       
       vy=v*np.sin(theta)
       
       velox_x[0,i]=vx
       
       velox_y[0,i]=vy
       
    return velox_x,velox_y


       
"""La seguente funzione calcola le forze tra particelle"""
def calcolo_forze(positions_x,positions_y,Nparts):
    
    forza_x = np.zeros((Nparts, Nparts))
    
    acel_listx=[]
    
    acel_listy=[]
    
    forza_y = np.zeros((Nparts, Nparts))
    
    for i in range(0, Nparts - 1):
        
        for j in range(i + 1,Nparts):
            
            r_x =positions_x[j] - positions_x[i]
            
            r_y = positions_y[j] - positions_y[i]
            
            modulo = np.sqrt(r_x **2+r_y**2)
            
            force_scalar = lj_force(modulo,epsilon, sigma)
            force_x = force_scalar * r_x / modulo
            
            force_y = force_scalar * r_y / modulo
            
            accel_x[i, j] = force_x / mass_of_argon
            
            accel_x[j, i] = - force_x / mass_of_argon
            
            accel_y[i, j] = force_y / mass_of_argon
            
            accel_y[j, i] = - force_y / mass_of_argon
            
    for i in range(Nparts):
        
         acel_listx.append(np.sum(force_x[i,:]))
         
         acel_listy.append(np.sum(force_y[i,:]))
         
    return acel_listx, acel_listy




"""Qeste due funzioni aggiornano le posizione e le velocità in
ingresso prende le vecchie posizioni, le vecchie velocità, le forze
e l'intervallo di tempo"""
def update_pos(x,y, vx,vy, ax,ay, dt,i):
    
    for k in range(Nparts):
        
        posi_x= x[i,k] + vx[i,k] * dt + 0.5 * ax[k] * dt * dt
        
        posi_y= y[i,k] + vy[i,k] * dt + 0.5 * ay[k] * dt * dt
        
    return posi_x,posi_y

def update_velo(vx,vy,ax,ay, a1x,a1y,dt):
    
    for h in range(Nparts):
        
        velo_x= vx[i,h] + 0.5 * (ax[h] + a1x[h]) * dt
        
        velo_y= vy[i,h] + 0.5 * (ay[h] + a1y[h]) * dt
        
    return velo_x, velo_y

            

"""In questa parte del programma si esegue la simulazione molecolare e infine 
restituisce gli array delle posizioni delle paticelle in input prende il passo
 temporale tau, il numero di steps, una temperatura iniziale """
def simula_md(dt, number_of_steps, T, Nparts):
    
    positions_x = np.zeros((number_of_steps, Nparts))
    
    positions_y = np.zeros((number_of_steps, Nparts))
    
    v_x, v_y = init_velocity(T, Nparts,Nsteps)
    
    position_x,positions_y=init_positions(L,Nparts)
    
    a_x,a_y = calcolo_forze(positions_x,positions_y,Nparts)
    
    for i in range(1,number_of_steps):
        
        #fare stessa cosa per vx #sto aggiornando la colonna iesima con la lista delle accele delle particell 
        
        position_x[i,:],positions_y[i,:] = update_pos(positions_x,positions_y,v_x,v_y,a_x,a_y,dt,i-1)
        
        a1x,a1y = calcolo_forze(positions_x,positions_y,Nparts)
       
        vx,vy = update_velo(v_x,v_y,a_x,a_y, a1x,a1y,dt)
        
        a_x = np.array(a1x)
        
        a_y = np.array(a1y)
        
        positions_x[i, :] = x
        
        positions_y[i, :] = y
   
    return positions_x,positions_y


x,y = simula_md(dt,Nsteps,T,Nparts)



for i in range(simula_md.shape[1]):
    plt.plot(simula_md[:, i], '.', label='atom {}'.format(i))
plt.xlabel("Step")
plt.ylabel("Position (Angstrom)")
plt.legend()
plt.show()