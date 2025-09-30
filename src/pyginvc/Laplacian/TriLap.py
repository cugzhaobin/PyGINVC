# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:38:43 2019

@author: zhao
"""
import logging, h5py
import numpy as np
from numpy import logical_and

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")


class TriLap(object):
    '''
    '''
    
    def __init__(self, vertex, element, method='1'):
        '''
        '''
        self.vertex  = vertex
        self.element = element
        if method == '1':
            self.tri_smooth()
            logging.info('smoothing method 1 is used.')
        elif method == '2':
            # the smoothing factor should be large in this case
            self.tri_laplacian(allow=True)
            logging.info('smoothing method 2 is used.')
        elif method == '3':
            self.build_laplacian()
            logging.info('smoothing method 3 is used.')
        elif method == '4':
            self.tikhonov()
            logging.info('smoothing method 4 is used.')
        else:
            self.build_laplacian()
    
    
    @staticmethod    
    def neighbors(i, tri):
        '''
        Find adjacent element for given element in the triangular mesh list
        Revised from neighbors.m in NIF program.
        
        Modified by Zhao Bin, June 28, 2021

        Parameters
        ----------
        i : int
            index of tri element
        tri : array
            triangle element

        Returns
        -------
        k : list
            index of neighors
        '''
        
        
        if i > np.size(tri,0):
            logging.warning("error")
        elif np.size(tri,1) < 3:
            logging.warning("error")
        elif np.size(tri,0) == 3 and np.size(tri,1) > 3:
            logging.warning("error")
        elif np.size(tri,1) > 3:
            logging.warning("error")
    
        v1 = tri[i,0]
        v2 = tri[i,1]
        v3 = tri[i,2]
        
    
        # find element adjacemnt to side 3, across vertice 3
        m1 = np.where(logical_and(tri[:,0]==v1, tri[:,1]==v2))[0]
        m2 = np.where(logical_and(tri[:,1]==v1, tri[:,2]==v2))[0]
        m3 = np.where(logical_and(tri[:,2]==v1, tri[:,0]==v2))[0]
        m4 = np.where(logical_and(tri[:,1]==v1, tri[:,0]==v2))[0]
        m5 = np.where(logical_and(tri[:,2]==v1, tri[:,1]==v2))[0]
        m6 = np.where(logical_and(tri[:,0]==v1, tri[:,2]==v2))[0]
        
        
        n3 = np.unique(np.hstack((m1, m2, m3, m4, m5, m6)))
        n3 = np.setdiff1d(n3, i)
        ed3 = 3*np.ones((1,len(n3)),dtype=int)
        ed3 = np.unique(ed3)
    
        # find element adjacemnt to side 2, across vertice 2
        m1 = np.where(logical_and(tri[:,0]==v3, tri[:,1]==v1))[0]
        m2 = np.where(logical_and(tri[:,1]==v3, tri[:,2]==v1))[0]
        m3 = np.where(logical_and(tri[:,2]==v3, tri[:,0]==v1))[0]
        m4 = np.where(logical_and(tri[:,1]==v3, tri[:,0]==v1))[0]
        m5 = np.where(logical_and(tri[:,2]==v3, tri[:,1]==v1))[0]
        m6 = np.where(logical_and(tri[:,0]==v3, tri[:,2]==v1))[0]
        n2 = np.unique(np.hstack((m1, m2, m3, m4, m5, m6)))
        n2 = np.setdiff1d(n2, i)
        ed2 = 2*np.ones((1,len(n2)),dtype=int)
        ed2 = np.unique(ed2)
    
        # find element adjacemnt to side 1, across vertice 1
        m1 = np.where(logical_and(tri[:,0]==v2, tri[:,1]==v3))[0]
        m2 = np.where(logical_and(tri[:,1]==v2, tri[:,2]==v3))[0]
        m3 = np.where(logical_and(tri[:,2]==v2, tri[:,0]==v3))[0]
        m4 = np.where(logical_and(tri[:,1]==v2, tri[:,0]==v3))[0]
        m5 = np.where(logical_and(tri[:,2]==v2, tri[:,1]==v3))[0]
        m6 = np.where(logical_and(tri[:,0]==v2, tri[:,2]==v3))[0]
        n1 = np.unique(np.hstack((m1, m2, m3, m4, m5, m6)))
        n1 = np.setdiff1d(n1, i)
        ed1 = 1*np.ones((1, len(n1)),dtype=int)
        ed1 = np.unique(ed1)
    
        edge = np.hstack((ed1, ed2, ed3))-1
        k    = np.hstack((n1,n2,n3))
        return k, edge
    
    
    
    def tri_smooth(self):
        '''
        Calculate Laplcaian for triangular dislocation
        Written by Zhao Bin, Institute of Seismology, CEA. May 12 2017
        Ref: Barnhat & Lohman, 2010. Faultresampler code.
        Input:
        element  --- 
        '''
    
        element = self.element
        smooth = np.eye(3*len(element))
        for i in range(len(element)):
            common = []
            a      = element[i,:]
            for j in range(len(element)):
                b  = element[j,:]
                tmp= np.intersect1d(a, b)
                if len(tmp) == 2:
                    common.append(j)
            nneighbor = len(common)
            common = np.array(common)
    
            if nneighbor == 1:
                ind = common
                smooth[3*i, 3*ind] = -1
                smooth[3*i+1, 3*ind+1] = -1
                smooth[3*i+2, 3*ind+2] = -1
            elif nneighbor == 2:
                smooth[3*i, 3*common] = -1.0/2
                smooth[3*i+1, 3*common+1] = -1.0/2
                smooth[3*i+2, 3*common+2] = -1.0/2
            elif nneighbor == 3:
                smooth[3*i, 3*common] = -1.0/3
                smooth[3*i+1, 3*common+1] = -1.0/3
                smooth[3*i+2, 3*common+2] = -1.0/3
            smooth[3*i, 3*common] = -1.0/3
            smooth[3*i+1, 3*common+1] = -1.0/3
            smooth[3*i+2, 3*common+2] = -1.0/3
    
        self.G_lap = smooth


    def tri_laplacian(self, allow=True):
        '''
        Make Laplcaian for triangular dislocation, method from NIF by Liu Zhen
        Ref. Maerten et al., 2005, BSSA
        
        Modified by Zhao Bin, June 28, 2021. We discard using matplotlib.tri.Triangulation
        to obtain neighbors, since it is limited up to three edges.

        Modified by Zhao Bin, Aug. 3, 2021. adding code when allow is zero
        
    
        Input:
        vertex: array, n*3, vertice number
        tri   : array, n*3, element number
        allow : True/False, allow slip at the bounary, False, does not allow
    
        Output:
        Lap     :  m*m
        Lap_inv :  m*m
        '''
    
        vertex    = self.vertex
        tri       = self.element
        nelem     = len(tri)
    
        # calculate center coordinates of elements
        tri  = tri-1
        cxyz = np.zeros((len(tri),3))
        cxyz[:,0] = (vertex[tri[:,0],0]+vertex[tri[:,1],0]+vertex[tri[:,2],0])/3
        cxyz[:,1] = (vertex[tri[:,0],1]+vertex[tri[:,1],1]+vertex[tri[:,2],1])/3
        cxyz[:,2] = (vertex[tri[:,0],2]+vertex[tri[:,1],2]+vertex[tri[:,2],2])/3

        x    = vertex[tri[:,0:3],0].reshape((nelem,3))
        y    = vertex[tri[:,0:3],1].reshape((nelem,3))
        z    = vertex[tri[:,0:3],2].reshape((nelem,3))
    
        # trian = Triangulation(vertex[:,0], vertex[:,1], tri)
        # neig  = trian.neighbors
    
        Lap  = np.zeros((3*nelem, 3*nelem))
        for i in range(nelem):
            k = []
            neig,edge = self.neighbors(i, tri)
            for j in neig:
                if j!=-1: k.append(j)
        
            # calculate distance between centers
            diff = np.tile(cxyz[i,:], (len(k),1))-cxyz[k,:]
            hh   = np.sqrt(np.sum(diff**2,1))
       
            
            # Laplacian term
            Li               = np.sum(hh)
            Lap[3*i,3*i]     = -2.0/Li*np.sum(1.0/hh)
            Lap[3*i+1,3*i+1] = -2.0/Li*np.sum(1.0/hh)
            Lap[3*i+2,3*i+2] = -2.0/Li*np.sum(1.0/hh)
            
            for j in range(len(k)):
                Lap[3*i+0, 3*k[j]+0] = 2/Li/hh[j]
                Lap[3*i+1, 3*k[j]+1] = 2/Li/hh[j]
                Lap[3*i+2, 3*k[j]+2] = 2/Li/hh[j]
            if allow != True:
                nb = len(k)
                if nb < 3:
                    r0 = np.vstack((x[i], y[i], z[i]))
                    r1 = r0[:,[1,2,0]]
                    redg = (r0+r1)/2.0
                    redg = redg.T
                    redg = redg[[1,2,0]]
                    diffe = np.tile(cxyz[i], (3,1))-redg
                    dist = np.sqrt(np.sum(diffe**2,1))
                    cc   = np.setdiff1d([0,1,2], edge)
                    ii   = np.empty(0, dtype=int)
                    for c in cc:
                        ii = np.hstack((ii, np.where(np.array([0,1,2],dtype=int)==c)[0]))
                    h1   = np.zeros(3)
                    h1[ii]   = dist[cc]
                    h1[edge] = hh
                    Li       = np.sum(h1)
                    Lap[3*i+0, 3*i+0] = -2.0/Li*np.sum(1.0/h1)
                    Lap[3*i+1, 3*i+1] = -2.0/Li*np.sum(1.0/h1)
                    Lap[3*i+2, 3*i+2] = -2.0/Li*np.sum(1.0/h1)
                    for j in range(len(k)):
                        Lap[3*i+0, 3*k[j]+0] = 2/Li/hh[j]
                        Lap[3*i+1, 3*k[j]+1] = 2/Li/hh[j]
                        Lap[3*i+2, 3*k[j]+2] = 2/Li/hh[j]

        self.G_lap = Lap
        return

    def build_laplacian(self):
        '''
        Make Laplcaian for triangular dislocation
        
        Modified by Zhao Bin, June 28, 2021. We discard using matplotlib.tri.Triangulation
        to obtain neighbors, since it is limited up to three edges.
        Adopted from CSI code (https://github.com/jolivetr/csi)
    
        Input:
        vertex: array, n*3, vertice number
        tri   : array, n*3, element number
        allow : 1, allow slip at the bounary, ~=1, does not allow
    
        Output:
        Lap     :  m*m
        Lap_inv :  m*m
        '''
    
        vertex    = self.vertex
        tri       = self.element
        nelem     = len(tri)
    
        # calculate center coordinates of elements
        tri  = tri-1
        cxyz = np.zeros((len(tri),3))
        cxyz[:,0] = (vertex[tri[:,0],0]+vertex[tri[:,1],0]+vertex[tri[:,2],0])/3
        cxyz[:,1] = (vertex[tri[:,0],1]+vertex[tri[:,1],1]+vertex[tri[:,2],1])/3
        cxyz[:,2] = (vertex[tri[:,0],2]+vertex[tri[:,1],2]+vertex[tri[:,2],2])/3
    
#       trian = Triangulation(vertex[:,0], vertex[:,1], tri)
#       neig  = trian.neighbors
        sumP = 0 
        Lap  = np.zeros((3*nelem, 3*nelem))
        for i in range(nelem):
            neig,edge = self.neighbors(i, tri)
            k = []
            for j in neig:
                if j!=-1: k.append(j)
        
            # calculate distance between centers
            diff = np.tile(cxyz[i,:], (len(k),1))-cxyz[k,:]
            dist = np.sqrt(np.sum(diff**2,1))
       
            
            if len(dist) == 3:
                h12, h13, h14 = dist
                Lap[3*i,3*k[0]]     = -h13*h14
                Lap[3*i+1,3*k[0]+1] = -h13*h14
                Lap[3*i+2,3*k[0]+2] = -h13*h14
                Lap[3*i,3*k[1]]     = -h12*h14
                Lap[3*i+1,3*k[1]+1] = -h12*h14
                Lap[3*i+2,3*k[1]+2] = -h12*h14
                Lap[3*i,3*k[2]]     = -h12*h13
                Lap[3*i+1,3*k[2]+1] = -h12*h13
                Lap[3*i+2,3*k[2]+2] = -h12*h13
        
                sumP = h13*h14 + h12*h14 + h12*h13
            if len(dist) == 2:
                h12, h13 = dist
                h14 = max(h12, h13)
                Lap[3*i,3*k[0]]     = -h13*h14
                Lap[3*i+1,3*k[0]+1] = -h13*h14
                Lap[3*i+2,3*k[0]+2] = -h13*h14
                Lap[3*i,3*k[1]]     = -h12*h14
                Lap[3*i+1,3*k[1]+1] = -h12*h14
                Lap[3*i+2,3*k[1]+2] = -h12*h14
                sumP = h13*h14+h12*h14+h12*h13
            Lap[3*i,3*i]     = sumP
            Lap[3*i+1,3*i+1] = sumP
            Lap[3*i+2,3*i+2] = sumP
        Lap = Lap/np.max(np.abs(np.diag(Lap)))
        self.G_lap = Lap
        return

    def tikhonov(self):
        '''
        Tikhonov (minimum norm) smoothing.
        '''
        tri   = self.element
        nelem = len(tri)
        print('********************')
        print(nelem)
        Lap   = np.eye(3*nelem)
        self.G_lap = Lap
        return
    
    def DumpLap(self):

        with h5py.File('lap.h5', 'w') as h5:
            h5.create_dataset('lap', data=self.G_lap, compression='gzip')
