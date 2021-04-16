#import numpy as np
#import matplotlib.pyplot as plt
from pylab import *

ion()

x = linspace(0,2*np.pi,10)
y = sin(x)

plot(x, y)
