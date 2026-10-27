#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
PyGINVC's rectangular dislocation with the numpy-2 fix.

`pyginvc.Greens.Okada.Okada.disloc` is correct apart from one line,
`x = np.mat(x)`, which numpy 2 removed.  This subclass changes only that
line; the geometry convention is unchanged:

    disgeom = [len, wid, dep, dip, strike, delE, delN, ss, ds, op]

where `dep` is the depth of the lower edge and `(delE, delN)` is the map
midpoint of the lower edge.  `x` is [East, North] in km (a single point or
an (m, 2) array); the result is an (m, 3) [E, N, U] array.

The module is imported lazily by models.Okada3D, so a missing `okada85`
never prevents the 2D models from being imported.
'''
import numpy as np

from pyginvc.Greens import okada85 as okada
from pyginvc.Greens.Okada import Okada as _Okada


class OkadaKernel(_Okada):
    '''
    Drop-in Okada.disloc that runs under numpy 2.
    '''

    @staticmethod
    def disloc(nu, disgeom, x):
        '''
        Surface displacement of one rectangular dislocation.

        Input:
            nu       = Poisson's ratio
            disgeom  = [len, wid, dep, dip, strike, delE, delN, ss, ds, op]
                       lengths km, angles degree, slip in any unit
            x        = [East, North] km, shape (2,) or (m, 2)
        Output:
            u        = (m, 3) [E, N, U] in the unit of ss/ds
        '''
        alp = 1 - 2*nu
        dip = disgeom[3]
        if dip == 90 or dip == 270:
            cs = 0.0
            sn = 1.0
        else:
            dipr = np.deg2rad(dip)
            cs = np.cos(dipr)
            sn = np.sin(dipr)

        strkr = np.deg2rad(90 - disgeom[4])

        # the only change from PyGINVC: np.mat -> np.atleast_2d
        x = np.atleast_2d(np.asarray(x, dtype=float))
        m = x.shape[0]

        delE = disgeom[5]
        delN = disgeom[6]

        u  = np.zeros((m, 3))
        ut = np.zeros(9)
        for k in range(m):
            zlocal = (complex(x[k, 0], x[k, 1]) - complex(delE, delN)) \
                     * np.exp(complex(0, -strkr)) + 0.5*disgeom[0]
            xrec = np.array([zlocal.real, zlocal.imag])

            okada.okadaio(alp, [xrec[0], xrec[1]],
                          [disgeom[2], 0, disgeom[0], 0, disgeom[1],
                           sn, cs, disgeom[7], disgeom[8], disgeom[9]], ut)

            zglob = complex(ut[0], ut[1])*np.exp(complex(0, strkr))
            u[k, 0] = zglob.real
            u[k, 1] = zglob.imag
            u[k, 2] = ut[2]

        return u
