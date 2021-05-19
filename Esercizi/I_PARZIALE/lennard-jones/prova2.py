import numpy as np
import matplotlib.pyplot as plt
import scipy.constants  as spc
#from sklearn.metrics import mean_squared_error
import random
import sympy as sp


mass_of_argon = 39.948 # Massa dell'argon in unità di massa atomiche
Nsteps=100
Nparts=64 #numero di atomi di Argon
L=35 #amstrong dimensioni box
dt = 0.1 #passo temporale
T=300 #K kelvin
sigma=3.46 #Angstrong 10 alla - 10 metri Distanza in cui la potenziale energia è zero
epsilon = 0.0103 #eV è l'energia potenziale alla separazione di equilibrio


"""Definisco una funzione per calcolare la forza dal potenziale di Lennard-Jones. Questo dipende
dalla distanza tra le due particelle. Restituisce il valore della forza derivante dal potenziale
ovvero la forza di interazione tra due particelle di argon """
def lj_force(r, epsilon, sigma):

    return 48 * epsilon * np.power(
           
        sigma, 12) / np.power(
               
        r, 13) - 24 * epsilon * np.power(
               
        sigma, 6) / np.power(r, 7)


""" init_position è una funzione che inizializza le posizioni delle particelle
di argon """
def init_positions(L,Nparts):
   
    pos_array_x=np.zeros((Nsteps,Nparts))
    pos_array_y=np.zeros((Nsteps,Nparts))
   
    pos_array_y[0,:]=np.linspace(0,L,Nparts)
       
    return pos_array_x,pos_array_y
    #return pos_x,pos_y    


"""La seguente funzione inizializza le velocità, prende come argomenti la temperatura T ,
il Numero di particelle e gli step e restituisce array di velocità pari al numero di particelle
si calcola la velocità media delle particelle in equilibrio termico a temperatura T. Poichè l'energia
media a temperatura T è Kb*T calola la radice quadrata dell'energia divisa la massa ottenendo così
una stima delle velocità medie delle particelle inoltre effettua la caonversione da J a eV """      
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
    #DEBUGG:  
    return velox_x,velox_y
    #return vx,vy

     
"""La seguente funzione calcola tutte le forze di interazione tra le varie particelle,
si sono create 2 matrici che contengono tutte le forze che le particelle esercitano tra di loro
e si è inizializzata a zero."""
def calcolo_forze(positions_x,positions_y,Nparts,i):
   
    force_x = np.zeros((Nparts, Nparts))
    #in questa matrice si salvano tutte le forze lungo la componente x
   
    force_y = np.zeros((Nparts, Nparts))
    #in questa matrice si salvano tutte le forze lungo la componente y
   
    acel_listx=[]
   
    acel_listy=[]
   
    r_x=np.zeros((Nparts,Nparts))
   
    r_y=np.zeros((Nparts,Nparts))
   
    accel_x=np.zeros((Nparts,Nparts))
   
    accel_y=np.zeros((Nparts,Nparts))

    """Attraverso il seguente ciclo for in k si calcola la forza agente su ogni singola particella
    dalla prima alla n-unesima. il secondo in j Per evitare che l'ultima particella sia conteggiata due volte,
    successivamente calcola su ciascuna particella le forze dovute a tutte le altre particelle
    evitando le forze già calcolate uso solo le coppie uniche (i+1) . """


    for k in range(0, Nparts):
       
        for j in range(k+1 ,Nparts):

            """Calcolo le componenti delle distanze tra le particelle e successivamente ne faccio il modulo """
           
            r_x[k,j] =np.abs(positions_x[i, j] - positions_x[i, k])
            r_y[k,j] = np.abs(positions_y[i,j] - positions_y[i,k])
           
            modulo = np.sqrt(r_x[k,j] **2 + r_y[k,j]**2)
           
            """ Ora calcolo le forze tramite la funzione relativa al potenziale di lennard Jones e determino
            le accelerazioni utili per l'algoritmo della velocity di verlet successivamente salvo le forze
            all'interno della matrice """
           
            force_scalar = lj_force(modulo,epsilon, sigma)

            force_x[k, j] = force_scalar * r_x[k,j] / modulo
            force_y[k, j] = force_scalar * r_y[k,j] / modulo  
            force_x[j, k] = -force_scalar * r_x[k,j] / modulo
            force_y[j, k] = -force_scalar * r_y[k,j] / modulo  
           
            accel_x[k, j] = force_x[k,j] / mass_of_argon
            accel_x[j, k] = - force_x[k,j] / mass_of_argon
           
            accel_y[k, j] = force_y[k,j] / mass_of_argon
            accel_y[j, k] = - force_y[k,j] / mass_of_argon
           
    for i in range(Nparts):
       
         acel_listx.append(np.sum(force_x[i,:]))
         
         acel_listy.append(np.sum(force_y[i,:]))
         
    return acel_listx, acel_listy


"""Qeste due funzioni aggiornano le posizione e le velocità in
ingresso prende le vecchie posizioni, le vecchie velocità, le forze
e l'intervallo di tempo"""
def update_pos(x,y, vx,vy, ax,ay, dt,i):
   
    for k in range(Nparts):
       
        posi_x= x[k] + vx[k] * dt + 0.5 * ax[k] * dt * dt
       
        posi_y= y[k] + vy[k] * dt + 0.5 * ay[k] * dt * dt
       
    return posi_x,posi_y

def update_velo(vx,vy,ax,ay, a1x,a1y,dt, i):

    velo_x, velo_y = np.zeros(Nparts), np.zeros(Nparts)
   
    for h in range(Nparts):
       
        velo_x[h]= vx[i,h] + 0.5 * (ax[h] + a1x[h]) * dt
       
        velo_y[h]= vy[i,h] + 0.5 * (ay[h] + a1y[h]) * dt
       
    return velo_x, velo_y

#i=np.linspace(0,Nparts,1)

"""In questa parte del programma si esegue la simulazione molecolare e infine
restituisce gli array delle posizioni delle paticelle in input prende il passo
 temporale tau, il numero di steps, una temperatura iniziale """
def simula_md(dt, number_of_steps, T, Nparts):
   
    v_x, v_y = init_velocity(T, Nparts,Nsteps)
   
    positions_x,positions_y=init_positions(L,Nparts)
    print(v_x ,":Vx",v_y,":Vy")
    print('eccomi')
    print(positions_x,positions_y)
    a_x,a_y = calcolo_forze(positions_x,positions_y,Nparts, 0)
   
    """Dopo aver inizializzato le posizioni, le velocità , e le forze , il seguente
    ciclo for riaggiorna volta per volta posizioni e velocità ."""
   
    for i in range(1,number_of_steps):
       
        positions_x[i,:], positions_y[i,:]= update_pos(positions_x[i-1,:], positions_y[i-1,:], v_x[i-1,:],v_y[i-1,:],a_x,a_y,dt,i-1)

        a1x,a1y = calcolo_forze(positions_x,positions_y,Nparts,i)
       
        v_x[i,:], v_y[i,:]= update_velo(v_x,v_y,a_x,a_y,a1x, a1y, dt, i)

        a_x[i]=a1x
       
        a_y[i]=a1y
       
    return positions_x,positions_y

x,y = simula_md(dt,Nsteps,T,Nparts)

atoms=[]

def distribuzione_Maxwell_Boltzman(mass_of_argon,Nparts):
        """Applica la distribuzione di Boltzmann alle velocità atomiche"""
        normDist = []
        scaling_factor = np.sqrt(spc.Boltzmann*T/mass_of_argon)

        # Stabilisce una distribuzione normale
        for i in range(0, 3*Nparts):
            normDist.append(random.gauss(0,1))
            
        # Applica il fattore di scala alla distribuzione
        for number in range(0, 3*Nparts):
            normDist[number] = normDist[number]*scaling_factor
            
        # Distribuisce le velocità
        for atom in range(0, Nparts):
            vx = normDist[atom*3]
            vy = normDist[atom*3+1]
        
            
maxwell=distribuzione_Maxwell_Boltzman(mass_of_argon,Nparts)            

"""Questa funzione calcola lo spostamento quadratico medio in funzione del tempo"""            
def spotamento_quadratico_medio(particle_iesima, particle_0):
    return np.sqrt(((particle_iesima - particle_0) ** 2).mean())

delta_R = spotamento_quadratico_medio(np.array(d), np.array(p))
print("Lo spostamento quadratico medio è: " + str(delta_R))



def coeff_autodiff_Einstein(Delta_R,dt):
    x = sp.Symbol('x')
    y = Delta_R/2*dt
    limite=sp.limit(y, x, sp.oo)
    print("Il coefficiente di diffusione di Einstein è :\n",limite)

einstein= coeff_autodiff_Einstein(Delta_R,dt)



for i in range(simula_md.shape[1]):
    plt.plot(simula_md[:, i], '.', label='atom {}'.format(i))
plt.xlabel("Step")
plt.ylabel("Position (Angstrom)")
plt.legend()
plt.show()
