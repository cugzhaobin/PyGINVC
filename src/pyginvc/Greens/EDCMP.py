#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Mar 7 2016
# Revised by Zhao Bin, Apr. 26, 2019. to Object Oriented Code.

import numpy as np
import subprocess, logging
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")


class EDCMP(object):
    '''
    Okada is a class representing rectangualr dislocation
    '''    
    
    def __init__(self, flt, data, dict_green):
        '''
        Constructor.
        
        Input:
            flt        = an instance of class Fault
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
        '''
        
        import h5py, os
        greenfile = dict_green['greenfile']
        if greenfile == "":
            self.GenGreens(flt, data, dict_green)
        elif greenfile == "SAVE":
            self.GenGreens(flt, data, dict_green)
            with h5py.File('greenfunc.h5', 'w') as h5:
                h5.create_dataset('G', data = self.G, compression='gzip')
                h5.create_dataset('G_sar', data = self.G_sar, compression='gzip')
        elif os.path.isfile(greenfile):
            with h5py.File(greenfile, 'r') as h5:
                self.G     = h5['G'][()]
                self.G_sar = h5['G_sar'][()]
        else:
            self.GenGreens(flt, data, dict_green)
            
        return


    def GenGreens(self, flt, data, dict_green):
        '''
        Generate Green's Function.

        Input:
            flt        = an instance of class Fault
            data       = an instance of class GeoData
            dict_green = a dict containing 'greentype', 'nu', 'bcs', 'greenfile'
        Output:
            self.G     = Green's function for GPS
            self.G_sar = Green's function for InSAR
        '''
        import numpy as np
        from pyginvc.libs import geotools as gt

        llh_gps       = data.llh_gps
        llh_lev       = data.llh_lev
        llh_sar       = data.llh_sar
        ndim          = data.ndim
        unit          = data.unit
        origin        = flt.origin
        edcmp_flt     = flt.FaultGeom2EDCMP()
        greentype     = dict_green['greentype']
        edgrndir      = dict_green['grndir']
        edgrn_ssfile  = dict_green['grn_ssfile']
        edgrn_dsfile  = dict_green['grn_dsfile']
        edgrn_clfile  = dict_green['grn_clfile']

    
        # convert the LLH to local coordinate
        xy_gps        = np.zeros((len(llh_gps),2))
        xy_lev        = np.zeros((len(llh_lev),2))
        xy_sar        = np.zeros((len(llh_sar),2))
        for i in range(len(llh_gps)):
            xy_gps[i,:] = gt.llh2localxy(llh_gps[i], origin)
        for i in range(len(llh_lev)):
            xy_lev[i,:] = gt.llh2localxy(llh_lev[i], origin)
        for i in range(len(llh_sar)):
            xy_sar[i,:] = gt.llh2localxy(llh_sar[i], origin)

        # if we have GPS data
        G_dis = np.array([])
        if len(xy_gps) > 0 and len(edcmp_flt)>0:
            G_dis = self.MakeGGPS(edcmp_flt, xy_gps, greentype, ndim, edgrndir, edgrn_ssfile, edgrn_dsfile, edgrn_clfile)
    
            # print the status
            logging.info('Green function for %d GPS stations are computed.' %(len(xy_gps)))
        G = G_dis
            
        # if we have level data
        G_dis = np.array([])    
        if len(xy_lev) > 0:
            G_dis = self.MakeGLEV(edcmp_flt, xy_lev, greentype, edgrndir, edgrn_ssfile, edgrn_dsfile, edgrn_clfile)
            G = np.hstack((G, G_dis))
    
            # print the status
            logging.info('Green function for %d Level stations are computed.' %(len(xy_lev)))
        
        # if we have SAR data
        G_sar   = np.array([])
        if len(xy_sar) > 0:
            G_sar = self.MakeGSAR(unit, edcmp_flt, xy_sar, greentype, edgrndir, edgrn_ssfile, edgrn_dsfile, edgrn_clfile)
    
            # print the status
            logging.info('Green function for %d InSAR points are computed.' %(len(xy_sar)))
        
        self.G     = G
        self.G_sar = G_sar

        # print the status
        logging.info('Green function for geodetic data are computed using Okada model.')
        
        return
        
    def MakeGGPS(self, edcmp_flt, xy, greentype, gdim, grndir, grn_ssfile, grn_dsfile, grn_clfile):
        '''
        Make green's function for GPS data
        Written by Zhao Bin when his is in UC Berkeley, Jan 28 2016
        Subroutine to make design matrix for slip
    
        Input:
            dis_geom  = Dislocation Geometry
            xy        = xy coordinates of GPS stations
            nu        = Poisson's ratio
            ss        = 0/1
            ds        = 0/1
            op        = 0/1
            gdim      = 2/3
        Output:
            G         = Green's Function for GPS
     
        '''
        ss = greentype[0]
        ds = greentype[1]
        op = greentype[2]
        if edcmp_flt.ndim == 1:
            edcmp_flt = edcmp_flt.reshape(len(edcmp_flt),10)
        if xy.ndim == 1:
            xy = xy.reshape(len(xy),2)

        nf   = edcmp_flt.shape[0]
        nsta = xy.shape[0]
    
        # whether we used 3D data
        if gdim == 3:
            ndata = nsta*3
        else:
            ndata = nsta*2
    
        # init for G
        G = np.zeros((ndata,3*nf))
    
        # for each subfaults
        for i in range(nf):
            # if we want to constrain strike slip 
            if bool(ss) is True:
                edcmp_flt[i,1] = 1.0
                edcmp_flt[i,9] = 0
                gen_edcmp_inp(edcmp_flt[i], xy[:,1], xy[:,0], np.zeros(len(xy)), grndir, grn_ssfile, grn_dsfile, grn_clfile)
                subprocess.call("echo tmp.inp | edcmp", shell=True)
                if gdim == 3:
                    dat = np.genfromtxt('tmp.disp', usecols=[2,3,4], comments='#')
                    dat[:,2] = -dat[:,2]
                if gdim == 2:
                    dat = np.genfromtxt('tmp.disp', usecols=[2,3])
                dat = dat.flatten()
                G[:,3*i+0] = dat

         
            # if we want to constrain dip slip 
            if bool(ds) is True:
                edcmp_flt[i,1] = 1.0
                edcmp_flt[i,9] = 90
                gen_edcmp_inp(edcmp_flt[i], xy[:,1], xy[:,0], np.zeros(len(xy)), grndir, grn_ssfile, grn_dsfile, grn_clfile)
                subprocess.call("echo tmp.inp | edcmp", shell=True)
                if gdim == 3:
                    dat = np.genfromtxt('tmp.disp', usecols=[2,3,4])
                    dat[:,2] = -dat[:,2]
                if gdim == 2:
                    dat = np.genfromtxt('tmp.disp', use_cols=[2,3])
                dat = dat.flatten()
                G[:,3*i+1] = dat
    
            # if we want to constrain openning slip 
            if bool(op) is True:
                edcmp_flt[i,1] = 0.0
                edcmp_flt[i,9] = 0
                gen_edcmp_inp(edcmp_flt[i], xy[:,1], xy[:,0], np.zeros(len(xy)), grndir, grn_ssfile, grn_dsfile, grn_clfile)
                subprocess.call("echo tmp.inp | edcmp", shell=True)
                if gdim == 3:
                    dat = np.genfromtxt('tmp.disp', usecols=[2,3,4])
                    dat[:,2] = -dat[:,2]
                if gdim == 2:
                    dat = np.genfromtxt('tmp.disp', usecols=[2,3])
                dat = dat.flatten()
                G[:,3*i+2] = dat
    
        return G
        
        
    def MakeGSAR(self, unit, edcmp_flt, xy, greentype, edgrndir, edgrn_ssfile, edgrn_dsfile, edgrn_clfile):
        '''
        Make Green's Function for InSAR data
        This treats range change data as absolute with respect to a reference
        level of zero deformation.  In practice the reference level is uncertain
        and the inversion should allow for a best-fitting reference level.
        Hence the appearance of 3*nf+1 parameters below.
    
        Input:
            unit         -  unit vector of InSAR data
            dis_geom     -  numpy array, Dislocation Geometry
            xy           -  numpy array, xy coordinates of stations
            nu           -  Poisson's ratio
            ss           -  if we will estimate strike   slip, ss=1, else ss=0
            ds           -  if we will estimate dip      slip, ds=1, else ds=0
            op           -  if we will estimate openning slip, op=1, else op=0
    
        Output:
            G            -  numpy array, desigin matrix
        
        '''
    
        import numpy as np
        from numpy import cos, sin, deg2rad, hstack
    
        # 
        ss = greentype[0]
        ds = greentype[1]
        op = greentype[2]
        if edcmp_flt.ndim == 1:
            edcmp_flt = edcmp_flt.reshape(len(edcmp_flt),10)
        if xy.ndim == 1:
            xy = xy.reshape(len(xy),2)

        nf     = edcmp_flt.shape[0]
        ndata  = xy.shape[0]

        # Descending track unit vector
        if len(unit) == 0:
            lookangle = 23
            track     = 13.9
            unit      = np.array([-cos(deg2rad(track))*sin(deg2rad(lookangle)), 
                                  sin(deg2rad(track))*sin(deg2rad(lookangle)), 
                                 -cos(deg2rad(lookangle))])
            unit      = -1.0*unit
        G = np.zeros((ndata, 3*nf+3))
    
        # for each subfaults
        for i in range(nf):   
            # SS motion.
            if ss != 0:
                edcmp_flt[i,1] = 1.0
                edcmp_flt[i,9] = 0
                gen_edcmp_inp(edcmp_flt[i], xy[:,1], xy[:,0], np.zeros(len(xy)), edgrndir, ss, ds, cl)
                subprocess.call("echo tmp.inp | edcmp", shell=True)
                dat = np.genfromtxt('tmp.disp', use_cols=[2,3,4])
                dat[:,2] = -dat[:,2]
                G[:,3*i] = unit.dot(dat).T
    
    	    # DS motion.
            if ds != 0: 
                edcmp_flt[i,1] = 1.0
                edcmp_flt[i,9] = 90
                gen_edcmp_inp(edcmp_flt[i], xy[:,1], xy[:,0], np.zeros(len(xy)), edgrndir, ss, ds, cl)
                subprocess.call("echo tmp.inp | edcmp", shell=True)
                dat = np.genfromtxt('tmp.disp', use_cols=[2,3,4])
                dat[:,2] = -dat[:,2]
                G[:,3*i] = unit.dot(dat).T
    
    	    # Opening
            if op != 0: 
                edcmp_flt[i,1] = 0.0
                edcmp_flt[i,9] = 0
                gen_edcmp_inp(edcmp_flt[i], xy[:,1], xy[:,0], np.zeros(len(xy)), edgrndir, ss, ds, cl)
                subprocess.call("echo tmp.inp | edcmp", shell=True)
                dat = np.genfromtxt('tmp.disp', use_cols=[2,3,4])
                dat[:,2] = -dat[:,2]
                G[:,3*i] = unit.dot(dat).T
    
        # The following is for adding a constant and slope to the SAR data
#       for j in range(ndata):
#           G[j,3*nf+0] = 1.0
#           G[j,3*nf+1] = xy[i,0]
#           G[j,3*nf+2] = xy[i,1]
     
        return G


    def MakeGLEV(self, edcmp_flt, xy_lev, greentype, edgrndir, edgrn_ssfile, edgrn_dsfile, edgrn_clfile):
        pass

def gen_edcmp_inp(edcmp_flt, n, e, u, edgrndir="", ss="", ds="", cl=""):
    '''
    make input file for EDCMP program when we want to calculate the stress changes at a single point
    By Zhao Bin, Institute of Seismology, CEA. jan  9 2018

    Input:
        sfault:   source fault rupture in EDCMP format, output from faultgeom2edcmp.py
        n     :   north component in km   
        e     :   east component in km   
        u     :   up component in km. downward is positive
    Output:
    '''

    if edcmp_flt.ndim == 1:
        edcmp_flt = edcmp_flt.reshape(1,10)
    
    # input file name for EDCMP
    infile = 'tmp.inp'
    fid = open(infile, 'w')

    # write the file
    # 1. receiver pms
    fid.write("#receiver points\n")
    fid.write('0\n')
    fid.write("%d\n" %(len(n)))
    for i in range(len(n)):
        line = "(%e, %e)\n" %(n[i]*1000.0, e[i]*1000.0)
        fid.write(line)

    # 2. output file
    fid.write('# output\n')
    fid.write("'./'\n")
    fid.write('1               0              0                0\n')
    line =  "'tmp.disp'    'tmp.strn'  'tmp.strss'  'tmp.tilt'\n"
    fid.write(line)

    # 3. source fault
    fid.write('# source fault\n')
    fid.write(str(len(edcmp_flt))+"\n")
    sfault = edcmp_flt
    for j in range(len(sfault)):
        line =  "%d %f %e %e %e %e %e %f %f %f\n" %(sfault[j,0], sfault[j,1], sfault[j,2], sfault[j,3], sfault[j,4], sfault[j,5], sfault[j,6], sfault[j,7], sfault[j,8], sfault[j,9])
        fid.write(line)

    # 4. earth model
    fid.write('# earth model\n')
    fid.write('1\n')
    line = "'{}' '{}' '{}' '{}'".format(edgrndir, ss, ds, cl)
    fid.write(line)

    # close the file
    fid.close()
