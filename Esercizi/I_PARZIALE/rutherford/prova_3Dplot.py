import numpy as np
import matplotlib.pyplot as plt

'''
N_part = 5
i = 0
num_list = []
np.random.seed()

while i < N_part:

    a = np.random.rand()

    if a <= 0.5:
        print('numero random: ', a)
        num_list.append(1)
        i += 1

print('lista:', num_list)
'''

z = np.linspace(0, 2*np.pi, 100)
x = np.sin(z)
y = np.cos(z)
a = [0 for i in range(100)]

fig = plt.subplots()
ax = plt.axes(projection = '3d')
ax.plot3D(x,y,[0 for i in range(100)])
plt.show()