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
        
        if method == '1':
            self.GenLap(self.nsegs, self.ndeps, self.bcs)
            logging.info('smoothing method 1 is used.')
        elif method == '4':
            self.tikhonov()
            logging.info('smoothing method 4 is used.')
        
        return

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
        
        # init parameters
        G_lapn = np.array([])
        G_lap  = np.array([])
    
        # for each plane
        for i in range(nplanes):
            if np.size(bcs) == 4:
                nlap   = self.MakeLapF_edge(nsegsa, ndepsa, bcs)
            elif np.size(bcs) == 0:
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
        LAP  = np.zeros((3*nmod, 3*nmod))
    
        # for each depth
        for i in range(ndeps):
            # for each along-strike
            for j in range(nsegs):
                indm = 3*nsegs*i     + 3*j
                inda = 3*nsegs*(i-1) + 3*j
                indb = 3*nsegs*(i+1) + 3*j
                indl = 3*nsegs*i + 3*j - 1
                indr = 3*nsegs*i + 3*j + 3
                for m in range(3*nmod):
                    if m == indm:
                        LAP[indm, indm]     = -1
                        LAP[indm+1, indm+1] = -1
                        LAP[indm+2, indm+2] = -1
                    elif m == inda:
                        if i == 0 and j == 0:
                            LAP[indm, inda]     = 0.0
                            LAP[indm+1, inda+1] = 0.0
                            LAP[indm+2, inda+2] = 0.0
                        elif i == 0 and j > 0 and j < nsegs-1:
                            LAP[indm, inda]     = 0.0
                            LAP[indm+1, inda+1] = 0.0
                            LAP[indm+2, inda+2] = 0.0
                        elif i == 0 and j == nsegs-1:
                            LAP[indm, inda]     = 0.0
                            LAP[indm+1, inda+1] = 0.0
                            LAP[indm+2, inda+2] = 0.0
                        elif j == 0 and i > 0 and i < ndeps-1:
                            LAP[indm, inda]     = 1.0/3
                            LAP[indm+1, inda+1] = 1.0/3
                            LAP[indm+2, inda+2] = 1.0/3
                        elif j == 0 and i == ndeps-1:
                            LAP[indm, inda]     = 0.5
                            LAP[indm+1, inda+1] = 0.5
                            LAP[indm+2, inda+2] = 0.5
                        elif i == ndeps-1 and j > 0 and j < nsegs-1:
                            LAP[indm, inda]     = 1./3
                            LAP[indm+1, inda+1] = 1./3
                            LAP[indm+2, inda+2] = 1./3
                        elif i == ndeps-1 and j == nsegs-1:
                            LAP[indm, inda]     = 0.5
                            LAP[indm+1, inda+1] = 0.5
                            LAP[indm+2, inda+2] = 0.5
                        elif j == nsegs-1 and i > 0 and i < ndeps-1:
                            LAP[indm, inda]     = 1./3
                            LAP[indm+1, inda+1] = 1./3
                            LAP[indm+2, inda+2] = 1./3
                        else:
                            LAP[indm, inda]     = 0.25
                            LAP[indm+1, inda+1] = 0.25
                            LAP[indm+2, inda+2] = 0.25
                    elif m == indb:
                        if i == 0 and j == 0:
                            LAP[indm, indb]     = 0.5
                            LAP[indm+1, indb+1] = 0.5
                            LAP[indm+2, indb+2] = 0.5
                        elif i == 0 and j > 0 and j < nsegs-1:
                            LAP[indm, indb]     = 1./3
                            LAP[indm+1, indb+1] = 1./3
                            LAP[indm+2, indb+2] = 1./3
                        elif i == 0 and j == nsegs-1:
                            LAP[indm, indb]     = 0.5
                            LAP[indm+1, indb+1] = 0.5
                            LAP[indm+2, indb+2] = 0.5
                        elif j == 0 and i > 0 and i < ndeps-1:
                            LAP[indm, indb]     = 1./3
                            LAP[indm+1, indb+1] = 1./3
                            LAP[indm+2, indb+2] = 1./3
                        elif j == 0 and i == ndeps-1:
                            LAP[indm, indb]     = 0.0
                            LAP[indm+1, indb+1] = 0.0
                            LAP[indm+2, indb+2] = 0.0
                        elif i == ndeps-1 and j > 0 and j < nsegs-1:
                            LAP[indm, indb]     = 0.0
                            LAP[indm+1, indb+1] = 0.0
                            LAP[indm+2, indb+2] = 0.0
                        elif i == ndeps-1 and j == nsegs-1:
                            LAP[indm, indb]     = 0.0
                            LAP[indm+1, indb+1] = 0.0
                            LAP[indm+2, indb+2] = 0.0
                        elif j == nsegs-1 and i > 0 and i < ndeps-1:
                            LAP[indm, indb]     = 1./3
                            LAP[indm+1, indb+1] = 1./3
                            LAP[indm+2, indb+2] = 1./3
                        else:
                            LAP[indm, indb]     = 0.25
                            LAP[indm+1, indb+1] = 0.25
                            LAP[indm+2, indb+2] = 0.25
                    elif m == indl:
                        if i == 0 and j == 0:
                            LAP[indm, indl]     = 0.0
                            LAP[indm+1, indl+1] = 0.0
                            LAP[indm+2, indl+2] = 0.0
                        elif i == 0 and j > 0 and i < nsegs-1:
                            LAP[indm, indl]     = 1./3
                            LAP[indm+1, indl+1] = 1./3
                            LAP[indm+2, indl+2] = 1./3
                        elif i == 0 and j == nsegs-1:
                            LAP[indm, indl]     = 0.5
                            LAP[indm+1, indl+1] = 0.5
                            LAP[indm+2, indl+2] = 0.5
                        elif j == 0 and i > 0 and i < ndeps-1:
                            LAP[indm, indl]     = 0.0
                            LAP[indm+1, indl+1] = 0.0
                            LAP[indm+2, indl+2] = 0.0
                        elif j == 0 and i == ndeps:
                            LAP[indm, indl]     = 0.0
                            LAP[indm+1, indl+1] = 0.0
                            LAP[indm+2, indl+2] = 0.0
                        elif i == ndeps-1 and j > 0 and j < nsegs-1:
                            LAP[indm, indl]     = 1./3
                            LAP[indm+1, indl+1] = 1./3
                            LAP[indm+2, indl+2] = 1./3
                        elif i == ndeps-1 and j == nsegs-1:
                            LAP[indm, indl]     = 0.5
                            LAP[indm+1, indl+1] = 0.5
                            LAP[indm+2, indl+2] = 0.5
                        elif j == nsegs-1 and i > 0 and i < ndeps-1:
                            LAP[indm, indl]     = 1./3
                            LAP[indm+1, indl+1] = 1./3
                            LAP[indm+2, indl+2] = 1./3
                        else:
                            LAP[indm, indl]     = 0.25
                            LAP[indm+1, indl+1] = 0.25
                            LAP[indm+2, indl+2] = 0.25
                    elif m == indr:
                        if i == 0 and j == 0:
                            LAP[indm, indr]     = 0.5
                            LAP[indm+1, indr+1] = 0.5
                            LAP[indm+2, indr+2] = 0.5
                        elif i == 0 and j > 0 and j < nsegs-1:
                            LAP[indm, indr]     = 1./3
                            LAP[indm+1, indr+1] = 1./3
                            LAP[indm+2, indr+2] = 1./3
                        elif i == 0 and j == nsegs-1:
                            LAP[indm, indr]     = 0.0
                            LAP[indm+1, indr+1] = 0.0
                            LAP[indm+2, indr+2] = 0.0
                        elif j == 0 and i > 0 and i < ndeps-1:
                            LAP[indm, indr]     = 1./3
                            LAP[indm+1, indr+1] = 1./3
                            LAP[indm+2, indr+2] = 1./3
                        elif j == 0 and i == ndeps-1:
                            LAP[indm, indr]     = 0.5
                            LAP[indm+1, indr+1] = 0.5
                            LAP[indm+2, indr+2] = 0.5
                        elif i == ndeps-1 and j > 0 and j < nsegs-1:
                            LAP[indm, indr]     = 1./3
                            LAP[indm+1, indr+1] = 1./3
                            LAP[indm+2, indr+2] = 1./3
                        elif i == ndeps-1 and j == nsegs-1:
                            LAP[indm, indr]     = 0.0
                            LAP[indm+1, indr+1] = 0.0
                            LAP[indm+2, indr+2] = 0.0
                        elif j == nsegs-1 and i > 0 and i < ndeps-1:
                            LAP[indm, indr]     = 0.0
                            LAP[indm+1, indr+1] = 0.0
                            LAP[indm+2, indr+2] = 0.0
                        else:
                            LAP[indm, indr]     = 0.25
                            LAP[indm+1, indr+1] = 0.25
                            LAP[indm+2, indr+2] = 0.25
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
