import bvls
import numpy as np
G = np.random.random((10,2))
m = np.array([1.0,2.0])
d = G.dot(m)
s = []
x = []
lower_bounds = np.array([0.0,0.0])
upper_bounds = np.array([1.5,1.5])
bounds = [lower_bounds,upper_bounds]
a=bvls.bvls(0, G,d,lower_bounds,upper_bounds)
print(a)
