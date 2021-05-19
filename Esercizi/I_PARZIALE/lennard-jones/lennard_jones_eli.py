import numpy as np #per gli array
import matplotlib.pyplot as plt #per plottare
import scipy.constants as spc #libreria delle costanti
import scipy.optimize as spo

#Atomi di Argon
mass_argon = 39.984 #uma
epsilon = 0.0103 #eV
sigma = 3.46 #A

def lj_force(d, L, epsilon, sigma):
                                                                            #print("d: " + str(d) + "; L: "+str(L))
    d= d % (np.sqrt(2) *L/2)
    if d<= 5*sigma:
        return 4*epsilon * (12 * sigma**12 / d**13 - 6*sigma**6/ d**7)
    return 0

def init_vel(T, N_part):
    R = 2* np.random.rand(N_part) - 1 #numero random -1<r<1
    return R * np.sqrt(spc.Boltzmann * T /(mass_argon* spc.electron_volt))

def compute_forces(positions, L, mass_argon, epsilon, sigma):
    N_part=len(positions)
    forces = np.zeros([N_part, N_part])
    for i in range(N_part - 1):
        for j in range(i+1, N_part):
            d= positions[j] - positions[i]
            force = lj_force(d, L, epsilon, sigma)/mass_argon
            forces[i,j] = force
            forces[j,i] = - force
    return np.sum(forces, axis = 0) #sommare tutti gli elementi aventi stesso indice di colonna

def update_positions(positions, velocities, forces, dt, L):
    new_pos= positions + velocities *dt + forces*dt**2/2
    return new_pos % L

def update_velocities(velocities, forces, new_forces, dt):
    return velocities + (forces + new_forces)*dt/2


def run_md(dt, N_steps, T, initial_pos, L, mass_argon, epsilon, sigma):
    N_part = len(initial_pos)
    #n steps riga e n° di particelle colonna
    positions = np.zeros([N_steps, N_part])
    positions[0, :] = initial_pos
    velocities = np.zeros([N_steps, N_part])
    velocities[0,:] = init_vel(T, N_part)
    x= initial_pos
    vel = init_vel(T, N_part)
    forces = compute_forces(initial_pos, L, mass_argon, epsilon, sigma)
    for i in range(1, N_steps):
        #Aggiornare posizioni
        x= update_positions(x, vel, forces, dt, L)
        #Aggiornare le forze con le nuove posizioni
        new_forces = compute_forces(x, L, mass_argon, epsilon, sigma) #??
        #Calcolare le velocità
        vel= update_velocities(vel, forces, new_forces, dt)
        #Aggiorno array posizione
        positions[i, :] = x
        velocities[1,:]=vel #'i' per 1
    return positions, velocities


dt = 0.01
N_steps = 5
N_part = 3
T=300 #K
L=11 #A

x = np.array([1, 5, 10]) #array di posizioni
#print("dt: " + str(dt) + "; N_steps: " + str(N_steps)+ "; T: "+str(T)+"; x: " + str(x) +"; L: "+str(L))
positions, velocities = run_md(dt, N_steps, T, x, L, mass_argon, epsilon, sigma)

fig, ax = plt.subplots()
for i in range(len(x)):
    ax.plot(positions[:,i], label="Atomo %g" % i)
ax.set_xlabel("Step") #t
ax.set_ylabel("Positions (Angstrom)")
ax.legend()

for i in range(N_part):
    for k in range(N_steps):
        print('velocità atomi:', velocities[k, i])



'''
#Distribuzione di Maxwell Boltzmann
def fit(v,A,alpha):
    return A*v**2*np.exp(-alpha*v**2)

#istogramma
fig1, ax1 = plt.subplots()
ax1.hist(np.abs(velocities[-1,:]), histtype="step")
ax1.set_xlabel("v")
ax1.set_ylabel("Counts")

counts, bins= np.histogram(np.abs(velocities[-1,:]))
b_width = (bins[1]-bins[0])/2
bins = bins[1:]-b_width
p, cov = spo.curve_fit(fit,bins,counts,p0=[1, 100]) #A e alpha all'interno di p0
x = np.linspace(bins[0]-b_width, bins[-1], 1000)
ax1.plot(x, fit(x, p[0], p[1])) #p0=A, p1=aplha
print(p[1]/(mass_argon/(spc.Boltzmann*4*np.pi)))
'''



plt.show()