#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA. Dec. 23, 2018.
import numpy as np
from numpy import arctan

class strikeslip2d(object):
    '''
    Finite dislocation in elastic half-space.
                 x2 (up)
                 |
                 |
                 |
                 O-----------+--------------+ x1 
                /           /              /|
               /           /F             / |
              /           /A             /  |
             /           /U             /   |
            /           /L             /   /|
           /           /T             /   / |
          /           /              /   /  |
       x3+-----------+--------------+   /   |
         |                          |  /   /|
         |          locked          | /   / |
         |                          |/   /  |
         +--------------------------+ d2/   |
         |                          |  /    |
         |           slip           | /    /
         |                          |/    /
         +--------------------------+ d1 /
         |                          |   /
         |            locked        |  /
         |                          | /
         +--------------------------+

    if d2 == 0  : coseismic rupture
    if d1 == inf: interseismic phase
    if 0<d2<d1  : afterslip 
    '''

    def __init__(self, s, d1, d2, x2, constant):
        '''
        Constractor.

        Input:
            s    = slip/slip rate
            d1   = depth 1 in km
            d2   = depth 2 in km
            x3   = observation depth km
            constant = constant
        '''
        if d2 == 0:
            raise ValueError('d2 == 0 is the coseismic case and is not '
                             'handled here; for the interseismic case use '
                             'd1=inf, d2=D')
        self.d2       = d2
        self.s        = s
        self.d1       = d1
        self.x2       = x2
        self.constant = constant

    def forward(self, x1):
        '''
        Compute dislocation deformation along fault-normal distance x1

        Input:
            x1   = a arrary of distance from fault trace
        '''
        u3 = -self.s*(arctan(x1/(self.x2+self.d1))
                    - arctan(x1/(self.x2-self.d1))
                    - arctan(x1/(self.x2+self.d2))
                    + arctan(x1/(self.x2-self.d2)))/(2*np.pi)+self.constant
        return u3

    def plotProfile(self, x1, u3):
        '''
        Plot fault-normal displacements
        '''
        import matplotlib.pyplot as plt
        plt.figure(figsize=(15,5))
        plt.scatter(x1, u3, s=1)
        plt.plot(x1, V)
        plt.xlabel('Distance from the fault (km)')
        plt.ylabel('Fault-normal motion (mm/yr)')
        plt.show()
