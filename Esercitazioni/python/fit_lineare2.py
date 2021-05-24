import scipy.optimize as spo
import numpy as np

def retta(x, A, B):
    return A + B * x

x = np.linspace(0, 10, 100)
y = retta
p, cov = np.polyfit(x, retta, 1)