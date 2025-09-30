import tde
import numpy as np


x = np.ones(3)
y = np.ones(3)
z = np.ones(3)

x[0], y[0], z[0] = 0, 0, 1
x[1], y[1], z[1] = 0, 5, 3
x[2], y[2], z[2] = 5, 5, 3

e = np.linspace(0,5,10)
n = np.linspace(0,5,10)
z = np.zeros(len(e))
u=tde.calc_tri_displacements(e, n, z, x, y, z, 0.25, 0, 0, 1)
