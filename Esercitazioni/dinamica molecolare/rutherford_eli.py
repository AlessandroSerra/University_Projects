import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import vec2d as v2d

Z1, Z2 = 2 , 79
alpha_mass = 2 * spc.proton_mass + 2 * spc.neutron_mass
F0 = Z1 + Z2 *spc.e **2 / (4 * np.pi * spc.epsilon_0 * alpha_mass)              #N * m**2/kg
E = 5 *1e6 * spc.electron_volt                                                  #N * m
d = F0 / E *alpha_mass
vel= np.sqrt(2 * E / alpha_mass)
tau = d/ vel

def F(r):
    mod = r.mod()
    return F0/ mod**3 * r

def step(r, v, tau):
    new_r = r + v * tau + F(r) * tau **2 /2
    new_v = v + (F(r) + F(new_r)) * tau /2
    return new_r, new_v

def solve(r0, v0, tau, Nsteps):
    t, r, v = [0], [r0], [v0]
    for i in range(Nsteps-1):       #mettenedo il meno 1 ho il numero di posiszioni uguali a Nsteps
                                    #se non lo metto qavrò Nsteps +1 perchè il punto iniziale lo avevo già inserito
        new_r, new_v = step(r[i], v[i], tau)
        t.append(t[i]+tau)
        r.append(new_r)
        v.append(new_v)
    return t, r, v

theta = []
Nparticels = 1000

for i in range(Nparticels):
    r0 = v2d.vec2d(-100 * d, 100*(2 * np.random.rand() -1)*d)
    v0 = v2d.vec2d(vel, 0)

    Nsteps = 2 * int(r0.mod()/d)

    t, r, v = solve(r0, v0, tau, Nsteps)
    theta.append(v0.get_angle(v[-1]))

#lo uso se ho una sola particella
#fig, ax = plt.subplots()
#ax.plot([pos.x /d for pos in r], [pos.y / d for pos in r])
#ax.plot(0, 0, marker=".", color= "y")
#plt.show()

fig, ax = plt.subplots()
ax.hist(theta, histtype = "step")
ax.set_yscale("log")
ax.set_xlabel("$\\theta$")
ax.set_ylabel("Counts")

plt.show()