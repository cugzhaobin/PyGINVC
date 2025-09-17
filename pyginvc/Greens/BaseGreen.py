#!/usr/bin/env python
# Written by Zhao Bin, Aug. 26, 2020.
import logging, os, h5py
import numpy as np
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


#   def GenGreens(self, flt, data, dict_green):
#       '''
#       Generate Green's Function.

#       Input:
#           flt        = an instance of class Fault
#           data       = an instance of class GeoData
#           dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
#       Output:
#           self.G     = Green's function for GPS
#           self.G_sar = Green's function for InSAR
#       '''
#       import numpy as np
#       import geotools as gt

#       llh_gps       = data.llh_gps
#       llh_lev       = data.llh_lev
#       llh_sar       = data.llh_sar
#       ndim          = data.ndim
#       unit          = data.unit
#       origin        = flt.origin
#       dis_geom_grid = flt.dis_geom_grid
        
#       greentype     = dict_green['greentype']
#       nu            = dict_green['nu']
#       if 'verbose' not in dict_green.keys():
#           verbose   = True
#       else:
#           verbose   = dict_green['verbose']
    
        # convert the LLH to local coordinate
#       xy_gps        = np.zeros((len(llh_gps),2))
#       xy_lev        = np.zeros((len(llh_lev),2))
#       xy_sar        = np.zeros((len(llh_sar),2))
#       if len(llh_gps) > 0:
#           for i in range(len(llh_gps)):
#               xy_gps[i,:] = gt.llh2localxy(llh_gps[i], origin)
#               xy_gps[i,:] = gt.llh2utm(llh_gps[i], origin)
#           
#       if len(llh_lev) > 0:
#           for i in range(len(llh_lev)):
#               xy_lev[i,:] = gt.llh2localxy(llh_lev[i], origin)
#               xy_lev[i,:] = gt.llh2utm(llh_lev[i], origin)
            
#       if len(llh_sar) > 0:
#           for i in range(len(llh_sar)):
#               xy_sar[i,:] = gt.llh2localxy(llh_sar[i], origin)
#               xy_sar[i,:] = gt.llh2utm(llh_sar[i], origin)
    
        # if we have GPS data
#       G_dis = np.array([])
#       if len(xy_gps) > 0 and len(dis_geom_grid)>0:
#           G_dis = self.MakeGGPS(dis_geom_grid, xy_gps, nu, greentype[0], greentype[1], greentype[2], ndim)
    
            # print the status
#           if verbose:
#               logging.info('Green function for %d GPS stations are computed.' %(len(xy_gps)))
#       G = G_dis
            
        # if we have level data
#       G_dis = np.array([])    
#       if len(xy_lev) > 0:
#           G_dis = self.MakeGLEV(dis_geom_grid, xy_lev, nu, greentype[0], greentype[1], greentype[2])
#           G = np.hstack((G, G_dis))
    
            # print the status
#           if verbose:
#               logging.info('Green function for %d Level stations are computed.' %(len(xy_lev)))
        
        # if we have SAR data
#       G_sar   = np.array([])
#       if len(xy_sar) > 0:
#           G_sar = self.MakeGSAR(unit, dis_geom_grid, xy_sar, nu, greentype[0], greentype[1], greentype[2])
#   
            # print the status
#           if verbose:
#               logging.info('Green function for %d InSAR points are computed.' %(len(xy_sar)))
        
#       self.G     = G
#       self.G_sar = G_sar

        # print the status
#       if verbose:
#           logging.info('Green function for geodetic data are computed using Okada model.')
        
#       return


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
