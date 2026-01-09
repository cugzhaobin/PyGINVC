#!/usr/bin/env python
# Written by Zhao Bin, Aug. 26, 2020.
import numpy as np
import logging, os, h5py
import pyginvc.libs.geotools as gt

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")

class BaseGreen(object):
    '''
    BaseGreen is a class representing general Green's function 
    '''    
    
    def __init__(self, flt, data, dict_green):
        '''
        Constructor.
        
        Input:
            flt        = an instance of class Fault
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
        '''
        
        greenfile = dict_green['greenfile']
        beta      = 0.0
        if os.path.isfile(greenfile) == True:
            self.LoadGreens(greenfile)
        else:
            self.GenGreens(flt, data, dict_green)
            self.SaveGreens()

        if 'rake_beta' in dict_green.keys():
            beta = dict_green['rake_beta']
            self.RotateGreens(beta)

        self.G_gps_ramp = None
        self.G_sar_ramp = None
        self.rake_beta = beta
        self.modulus   = float(dict_green['modulus'])
        return

    def LoadGreens(self, greenfile):
        '''
        Load Green's function.

        Input:
            grnfile    = file name of Green function
        '''
        if os.path.isfile(greenfile) == False:
            logging.info('Green file {} does not exist.'.format(greenfile))
        ext = greenfile.split('.')[-1]

        if ext == 'h5':
            with h5py.File(greenfile, 'r') as h5:
                self.G     = h5['G'][()]
                self.G_sar = h5['G_sar'][()]
                logging.info('Load Greens function from {}'.format(greenfile))
        elif ext == 'npy' or ext == 'npz':
            dat    = np.load('greenfile')
            self.G     = dat['G']
            self.G_sar = dat['G_sar']
            logging.info('Load Greens function from {}'.format(greenfile))
        return

    def SaveGreens(self, greenfile='green.h5'):
        '''
        '''
        if os.path.isfile('green.h5') == True:
            os.remove('green.h5')
        with h5py.File('green.h5', 'w') as h5:
            h5.create_dataset('G', data = self.G, compression='gzip')
            h5.create_dataset('G_sar', data = self.G_sar, compression='gzip')
            logging.info('Green function is saved to green.h5')

    def RotateGreens(self, beta, greenfile='green_rot.h5'):
        '''
        Rotate the green function.
           ds'     ds      ss'
            \      |      /
              \    |    /
                \  |  /  deta
        __________\|/____________ss
        '''

        if beta == 0: return
        rbeta = np.deg2rad(beta)
        R     = np.array([[np.cos(rbeta), np.sin(rbeta)],
                         [-np.sin(rbeta), np.cos(rbeta)]])

        Gnew    = np.empty(0)
        Gsarnew = np.empty(0)
        if len(self.G) > 0:
            G    = self.G
            Gnew = np.copy(G)
            for i in range(int(G.shape[1]/3)):
                Gnew[:, 3*i+0] = G[:,3*i+0]*R[0,0]+G[:,3*i+1]*R[0,1]
                Gnew[:, 3*i+1] = G[:,3*i+0]*R[1,0]+G[:,3*i+1]*R[1,1]
        if len(self.G_sar) > 0:
            G_sar     = self.G_sar
            Gsarnew   = np.copy(G_sar)
            for i in range(int(G_sar.shape[1]/3)):
                Gsarnew[:, 3*i+0] = G_sar[:,3*i+0]*R[0,0]+G_sar[:,3*i+1]*R[0,1]
                Gsarnew[:, 3*i+1] = G_sar[:,3*i+0]*R[1,0]+G_sar[:,3*i+1]*R[1,1]
        self.G     = Gnew
        self.G_sar = Gsarnew
        with h5py.File('./green_new.h5', 'w') as h5:
            h5.create_dataset('G',     data=Gnew,    compression='gzip')
            h5.create_dataset('G_sar', data=Gsarnew, compression='gzip')
        logging.info('Rotate Greens function finished.')

    def MakeGSARRamp(self, xy, dim):
        G = np.ones((len(xy),3))
        for i in range(len(xy)):
            G[i,1] = xy[i,0]
            G[i,2] = xy[i,1]
            
    def MakeGGPSRamp(self, xy, dim):
        G = np.ones((dim*len(xy), 3))
        for i in range(len(xy)):
            if dim == 2:
                G[2*i+0] = np.array([1, 0, -xy[i,0]])
                G[2*i+1] = np.array([0, 1,  xy[i,1]])
            elif dim == 3:
                G[3*i+0] = np.array([1, 0, -xy[i,0]])
                G[3*i+1] = np.array([0, 1,  xy[i,1]])
        return G
    
    def _convert_to_local(self, llh, origin):
        xy = np.zeros((len(llh),2))
        for i in range(len(llh)):
            xy[i,:] = gt.llh2localxy(llh[i], origin)
        return xy
        
