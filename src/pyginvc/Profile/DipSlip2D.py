#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA. Dec. 23, 2018.
import numpy as np
from numpy import arctan

class dipslip2d(object):
    '''
    Finite dislocation in elastic half-space -- dip-slip (thrust / normal).

    Cross section perpendicular to strike.  x1 is the horizontal distance
    from the fault trace; u is the horizontal (fault-normal) and w the
    vertical displacement.

                 surface
        ---------+-------------------------------- x1
                 |                          /
                 |  footwall               /
                 |                        /  dip
                 |                       /
                 |                      /
                 |       locked        /
                 |                    /
                 +-------------------+ D
                 |                  /
                 |       slip      /
                 |                /
                 |               /
                 |              /
                 |             /
                 +------------+

    the fault is locked above the depth D and slips freely below it.
    Surface displacement (see the note 震间断层滑动速率_文献总览 3.3,
    profile method - dip-slip fault):

        u(x) = s/pi [ cos(dip) atan(x/D) + sin(dip) x D/(x^2+D^2) ]
        w(x) = s/pi [ sin(dip) atan(x/D) - cos(dip) x D/(x^2+D^2) ]

    Far field (x -> +-inf):  u -> +- s cos(dip)/2 ,  w -> +- s sin(dip)/2 ,
    i.e. the across-fault velocity difference is s cos(dip) horizontally and
    s sin(dip) vertically.  The sign of s selects the sense of slip
    (thrust <-> normal).
    '''

    def __init__(self, s, dip, D, constant=0.0):
        '''
        Constractor.

        Input:
            s        = slip/slip rate along the dip
            dip      = fault dip in degree
            D        = locking depth in km
            constant = constant
        '''
        self.s        = s
        self.dip      = dip
        self.D        = D
        self.constant = constant

    @staticmethod
    def kernel(x1, D):
        '''
        The two shape functions of one locking depth D:
            atan(x/D)      and      x D / (x^2 + D^2)

        D == inf -> fault slips from the surface down (no locking)
        D == 0   -> fault fully locked
        '''
        if D == np.inf:
            return np.zeros_like(x1), np.zeros_like(x1)
        if D == 0.0:
            return np.sign(x1)*np.pi/2, np.zeros_like(x1)
        return arctan(x1/D), x1*D/(x1**2 + D**2)

    def forward(self, x1):
        '''
        Compute the surface deformation along fault-normal distance x1

        Input:
            x1   = a array of distance from fault trace in km
        Output:
            u    = fault-normal horizontal displacement in km or mm/yr
            w    = vertical displacement in km or mm/yr
        '''
        x1 = np.atleast_1d(np.asarray(x1, dtype=float))
        cd, sd = np.cos(np.deg2rad(self.dip)), np.sin(np.deg2rad(self.dip))
        S1, S2 = self.kernel(x1, self.D)

        u = self.s*(cd*S1 + sd*S2)/np.pi + self.constant
        w = self.s*(sd*S1 - cd*S2)/np.pi
        return u, w

    def plotProfile(self, x1, u, w):
        '''
        Plot fault-normal and vertical displacements
        '''
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        ax[0].plot(x1, u, 'k')
        ax[0].set_ylabel('Fault-normal motion')
        ax[1].plot(x1, w, 'k')
        ax[1].set_ylabel('Vertical motion')
        ax[1].set_xlabel('Distance from the fault (km)')
        plt.tight_layout()
        plt.show()


if __name__ == '__main__':
    x = np.arange(-150, 150, 0.1)
    ds = dipslip2d(10, 30, 15, 0.0)
    u, w = ds.forward(x)
    ds.plotProfile(x, u, w)
