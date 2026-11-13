# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 09:34:53 2019

@author: zhao
"""
import numpy as np
import logging
import h5py
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")


class RectPlaneLap(object):
    '''
    RectPlaneLap is a class representing Lapacian for rectangualr fault
    '''
    def __init__(self, nsegs, ndeps, bcs, method='1'):
        '''
        Constructor.
        
        Input:
            nsegs = number of patchs along strike
            ndeps = number of patchs along depth
            bcs   = boundary conditions
        '''
        self.nsegs = nsegs
        self.ndeps = ndeps
        self.bcs   = bcs
        method     = str(method).strip()
        if method not in ('1', '4'):
            raise ValueError(
                "Unsupported Laplacian method %r; use '1' (Laplacian) or "
                "'4' (Tikhonov)." % (method,))
        self.method = method

    def _check_bcs(self, bcs):
        '''Normalise and validate the boundary-condition list.'''
        bcs = np.asarray([] if bcs is None else bcs, dtype=int).ravel()
        if bcs.size not in (0, 4):
            raise ValueError(
                'bcs must be empty or a four-element list (top, bottom, '
                'end, start), got %r.' % (bcs.tolist(),))
        if bcs.size and not np.isin(bcs, [0, 1]).all():
            raise ValueError(
                'bcs entries must be 0 or 1, got %r.' % (bcs.tolist(),))
        return bcs

    def build_laps(self):
        if self.method == '1':
            self.GenLap(self.nsegs, self.ndeps, self.bcs)
            logging.info('smoothing method 1 is used.')
        elif self.method == '4':
            self.tikhonov()
            logging.info('smoothing method 4 is used.')


    def GenLap(self, nsegsa, ndepsa, bcs, nplanes=1):
        '''
        Generate Green functions for GPS, InSAR and Level data
        Input:
            nplanes   :  number of planes, I'm not sure what does this mean
            nsegsa    :  number of patches along strike direction
            ndepsa    :  number of patches along dip dircetion
            bcs       :  python list [0,0,0,0]
                
        Output:
            G_lap     :  array, shape(m,m) Lap 
                
    
        MOD by Zhao Bin, Jan 1, 2017. if bcs = []
        '''
        
        if nsegsa < 2 or ndepsa < 2:
            raise ValueError(
                'nsegs and ndeps must both be >= 2 for Laplacian smoothing; '
                'got %d and %d. Use Tikhonov (method 4) instead.'
                % (nsegsa, ndepsa))

        bcs = self._check_bcs(bcs)

        # init parameters
        G_lapn = np.array([])
        G_lap  = np.array([])

        # for each plane
        for i in range(nplanes):
            if bcs.size == 4:
                nlap   = self.MakeLapF_edge(nsegsa, ndepsa, bcs)
            else:
                nlap   = self.MakeLapF(nsegsa, ndepsa)
            G_lapn = np.array([])
            
            if i > 0:
                for k in range(0,i-1):
                    zarray   = np.zeros((3*nsegsa[i]*ndepsa[i], 3*nsegsa[k]*ndepsa[k]))
                    G_lapn   = np.hstack(G_lapn, zarray)
            if len(G_lapn) > 0:
                G_lapn  = np.hstack(G_lapn, nlap)
            else:
                G_lapn  = nlap
                    
            
            if i+2 <= nplanes:
                for k in range(0,nplanes):
                    zarray   = np.zeros((3*nsegsa[i]*ndepsa[i], 3*nsegsa[k]*ndepsa[k]))
                    G_lapn   = np.hstack(G_lapn, zarray)    
            
            if len(G_lapn) & len(G_lap) > 0:
                G_lap  = np.vstack(G_lap, G_lapn)
            else:
                G_lap  = G_lapn
            
        # print the status
        logging.info('Laplacian matrix is computed.')
        self.G_lap =  G_lap
    
    
    @staticmethod
    def MakeLapF(nsegs, ndeps):
        '''
        Make a Laplacian smoothing filter for 3 slip components.
        Note: The model vector must be presorted in depth and along strike to use this filter.
        So, the model vector has elements m[d0, as0], m[d0, as1], m[d0, as2], m[d0, as3]
    
        Written by Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley.
        Dec 9, 2016
        '''
        
        # number of model vectors
        nmod = nsegs * ndeps
        n3   = 3*nmod
        LAP  = np.zeros((n3, n3))

        def put(row0, col0, value):
            '''Set the three component-diagonal cells of one stencil entry.

            Columns outside the model vector are skipped, mirroring the
            original loop, which only visited m in range(3*nmod).
            '''
            if 0 <= col0 < n3:
                LAP[row0,     col0]     = value
                LAP[row0 + 1, col0 + 1] = value
                LAP[row0 + 2, col0 + 2] = value

        # for each depth
        for i in range(ndeps):
            # for each along-strike
            for j in range(nsegs):
                indm = 3*nsegs*i     + 3*j
                inda = 3*nsegs*(i-1) + 3*j
                indb = 3*nsegs*(i+1) + 3*j
                indl = 3*nsegs*i + 3*j - 1
                indr = 3*nsegs*i + 3*j + 3

                # centre
                put(indm, indm, -1)

                # up neighbour
                if 0 <= inda < n3:
                    if i == 0 and j == 0:
                        value = 0.0
                    elif i == 0 and j > 0 and j < nsegs-1:
                        value = 0.0
                    elif i == 0 and j == nsegs-1:
                        value = 0.0
                    elif j == 0 and i > 0 and i < ndeps-1:
                        value = 1.0/3
                    elif j == 0 and i == ndeps-1:
                        value = 0.5
                    elif i == ndeps-1 and j > 0 and j < nsegs-1:
                        value = 1./3
                    elif i == ndeps-1 and j == nsegs-1:
                        value = 0.5
                    elif j == nsegs-1 and i > 0 and i < ndeps-1:
                        value = 1./3
                    else:
                        value = 0.25
                    put(indm, inda, value)

                # down neighbour
                if 0 <= indb < n3:
                    if i == 0 and j == 0:
                        value = 0.5
                    elif i == 0 and j > 0 and j < nsegs-1:
                        value = 1./3
                    elif i == 0 and j == nsegs-1:
                        value = 0.5
                    elif j == 0 and i > 0 and i < ndeps-1:
                        value = 1./3
                    elif j == 0 and i == ndeps-1:
                        value = 0.0
                    elif i == ndeps-1 and j > 0 and j < nsegs-1:
                        value = 0.0
                    elif i == ndeps-1 and j == nsegs-1:
                        value = 0.0
                    elif j == nsegs-1 and i > 0 and i < ndeps-1:
                        value = 1./3
                    else:
                        value = 0.25
                    put(indm, indb, value)

                # left neighbour
                if 0 <= indl < n3:
                    if i == 0 and j == 0:
                        value = 0.0
                    elif i == 0 and j > 0 and i < nsegs-1:
                        value = 1./3
                    elif i == 0 and j == nsegs-1:
                        value = 0.5
                    elif j == 0 and i > 0 and i < ndeps-1:
                        value = 0.0
                    elif j == 0 and i == ndeps:
                        value = 0.0
                    elif i == ndeps-1 and j > 0 and j < nsegs-1:
                        value = 1./3
                    elif i == ndeps-1 and j == nsegs-1:
                        value = 0.5
                    elif j == nsegs-1 and i > 0 and i < ndeps-1:
                        value = 1./3
                    else:
                        value = 0.25
                    put(indm, indl, value)

                # right neighbour
                if 0 <= indr < n3:
                    if i == 0 and j == 0:
                        value = 0.5
                    elif i == 0 and j > 0 and j < nsegs-1:
                        value = 1./3
                    elif i == 0 and j == nsegs-1:
                        value = 0.0
                    elif j == 0 and i > 0 and i < ndeps-1:
                        value = 1./3
                    elif j == 0 and i == ndeps-1:
                        value = 0.5
                    elif i == ndeps-1 and j > 0 and j < nsegs-1:
                        value = 1./3
                    elif i == ndeps-1 and j == nsegs-1:
                        value = 0.0
                    elif j == nsegs-1 and i > 0 and i < ndeps-1:
                        value = 0.0
                    else:
                        value = 0.25
                    put(indm, indr, value)
        return LAP
    
    
    @staticmethod
    def MakeLapF_edge(nsegs, ndeps, bcs):
        '''
        Make Lap
        Written by Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley April 16 2016
        MOD by Zhao Bin, Jan 1, 2017. Correct a bug when nsegs != ndeps the result is wrong.
        MOD by Zhao Bin, June 26, 2018. Qi Yujie found the bug and fixed it.
            Input:
                nsegs   : number of segments along strike
                ndeps   : number of segments along dip
                bcs     : boundary conditions for Lap matrix. If bs is set to 1, sysmmetric bcs are used,
                        : whereas if zero, zero slip is assumed at fault termination.
                        : bcs must be a four element binary vector, with each element representing
                        : differenent sides of the mesh.
                        : The order of bcs is: top, bottom, left, right.
        '''
        
        nx     = nsegs
        nz     = ndeps
        
        # make smoothing matrix in real coordinates:
        LAP    = np.zeros((nx*nz*3, nx*nz*3))
        
        k      = range(nx*nz)
    #   k      = np.reshape(k, (nx, nz))
        k      = np.reshape(k, (nz, nx))
        k_zero = k
        
        counter= 0
        bcs    = np.array(bcs)*2-1 
        for i in range(0,len(k_zero[:,0])):
            for j in range(0,len(k_zero[0,:])):
    #           edgei  = ((i==0) | (i==len(k_zero[:,0])))
    #           edgej  = ((j==0) | (j==len(k_zero[0,:])))
                
                totalnodes = 4
                
                if i == 0:
                    im1  = k_zero[i+1,j]*bcs[0]
                else:
                    im1  = k_zero[i-1,j]
                    
                if i == len(k_zero[:,0])-1:
                    ip1  = k_zero[i-1,j]*bcs[1]
                else:
                    ip1  = k_zero[i+1,j]
                    
                pos = k_zero[i,j]
                
                if j == 0:
                    jm1  = k_zero[i,j+1]*bcs[2]
                else:
                    jm1  = k_zero[i,j-1]
    
                if j == len(k_zero[0,:])-1:
                    jp1  = k_zero[i,j-1]*bcs[3]
                else:
                    jp1  = k_zero[i,j+1]
                    
                # Now, this smoothing has to be put into the three components
    #           if im1 != 0:
                if im1 >= 0:
                    LAP[counter, 3*im1+0] = LAP[counter, 3*im1+0] + 1./totalnodes
                    LAP[counter+1, 3*im1+1] = LAP[counter+1, 3*im1+1] + 1./totalnodes
                    LAP[counter+2, 3*im1+2] = LAP[counter+2, 3*im1+2] + 1./totalnodes
        
    #           if ip1 != 0:
                if ip1 >= 0:
                    LAP[counter, 3*ip1+0] = LAP[counter, 3*ip1+0] + 1./totalnodes
                    LAP[counter+1, 3*ip1+1] = LAP[counter+1, 3*ip1+1] + 1./totalnodes
                    LAP[counter+2, 3*ip1+2] = LAP[counter+2, 3*ip1+2] + 1./totalnodes
                    
    #           if jp1 != 0:
                if jp1 >= 0:
                    LAP[counter, 3*jp1+0] = LAP[counter, 3*jp1+0] + 1./totalnodes
                    LAP[counter+1, 3*jp1+1] = LAP[counter+1, 3*jp1+1] + 1./totalnodes
                    LAP[counter+2, 3*jp1+2] = LAP[counter+2, 3*jp1+2] + 1./totalnodes
                    
    #           if jm1 != 0:
                if jm1 >= 0:
                    LAP[counter, 3*jm1+0] = LAP[counter, 3*jm1+0] + 1./totalnodes
                    LAP[counter+1, 3*jm1+1] = LAP[counter+1, 3*jm1+1] + 1./totalnodes
                    LAP[counter+2, 3*jm1+2] = LAP[counter+2, 3*jm1+2] + 1./totalnodes
                    
                    
                LAP[counter,3*pos+0] = -1
                LAP[counter+1,3*pos+1] = -1
                LAP[counter+2,3*pos+2] = -1
                
                counter = counter+3
        
        return LAP
    
    

    def tikhonov(self):
        '''
        Tikhonov (minimum norm) smoothing.
        '''
    
        nelem = self.nsegs*self.ndeps
        Lap   = np.eye(3*nelem)
        self.G_lap = Lap
        return
    

    def DumpLap(self):
        with h5py.File('lap.h5', 'w') as h5:
            h5.create_dataset('lap', data=self.G_lap, compression='gzip')
