#!/usr/bin/env python

import numpy as np
import pandas as pd
import os, sys, logging
from numpy import mat, math, dot 
from numpy import zeros, ones, deg2rad, rad2deg, sin, cos, sqrt
from pyginvc.libs import geotools as gt

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")


class Fault(object):
    '''
    Fault represents rectangular fault plane 
    '''
    def __init__(self, faultfile, nsegs, ndeps, doSubFault, origin=[]):
        '''
        Constructor.
        
        Input:
            faultfile   = Fault Geometry file name
            nsegs       = number of segments along the strike direction
            ndeps       = number of segments along the dip direction
            doSubFault  = True or False
            origin      = [lat, lon]
        '''
        self.faultfile     = faultfile
        self.nf            = None
        self.doSubFault    = doSubFault
        self.origin        = np.asarray(origin)
        self.nsegs         = nsegs
        self.ndeps         = ndeps
        self.dis_geom_grid = np.array([])
        self.geom_grid     = np.array([])
        self.load_fault()    
    

    def load_fault(self):
        # load the original fault geometry
        geom = self.LoadFaultGeom()

        
        # convert the fault into DISMODEL format
        if len(self.origin) != 2:
            self.origin = np.array([np.mean(geom[:,[3,5]]), np.mean(geom[:,[4,6]])])
        logging.info(f'origin is set to {self.origin[0]:10.4f} {self.origin[1]:10.4f}')
        dis_geom = self.SetOrigin(geom)
        
        # if we need subdivide the fault
        if not self.doSubFault:
            self.dis_geom_grid = dis_geom
            self.geom_grid     = geom
            self.nf            = len(geom)
            logging.info(f'Fault geometry is loaded. There are {self.nsegs} * {self.ndeps} patches in it.')
        else:
            pm = self.PatchFault(dis_geom)
            self.nf = len(pm)
            self.dis_geom_grid = np.column_stack((pm, zeros((self.nf,3))))            
            self.geom_grid     = self.disgeom2geom(self.dis_geom_grid)
            logging.info(f'Fault geometry is loaded and subdivided into {self.nsegs} * {self.ndeps} patches.')
        return
    
    def LoadFaultGeom(self):
        '''
        Load fault geometry file, the format is as following
    
        Input:
           filename  = filename of fault geometry
           The format of the file is:
           Width  Depth  Dip  SLat SLon ELat Elon Strike-Slip Dip-slip Openning
           ( km    km    deg  deg  deg  deg  deg  mm          mm       mm)
                                                  +rightslip  +updip
                                                  
        Output:
           geom      = numpy array, shape(m,10). The structure is same as input
        '''
    
        # check the file is exit
        isfile = os.path.isfile(self.faultfile)
        if isfile is False:
           logging.fatal('fault file does not exist! Please input another file!')
           sys.exit()
    
        # load the data
        geom = np.genfromtxt(self.faultfile, comments='#')
        if geom.ndim == 1: geom = geom.reshape((1, len(geom)))
        logging.info('Fault Geometry is loaded.')
        return geom

    def SetOrigin(self, geom):
        '''
        Written by Zhao Bin when his is in UC Berkeley, Feb 20 2016
    
        Input:
            geom     = Fault Geometry 
                      Width  Depth  Dip  SLon SLat ELon Elat Strike-Slip Dip-slip Openning
        Output:
            dis_geom = numpy array, containing Length Width Depth Dip Strike dE dN Strike-Slip Dip-slip Openning
        '''
        
        nf, ncol = geom.shape
        if ncol != 10:
            raise ValueError("geom must have 10 columns")
            
        dis_geom   = zeros((nf,10))
       
        for i in range(nf):
            # convert latitude and longitude to XY with respect to origin
            xy1  = gt.llh2utm([geom[i,3], geom[i,4]], self.origin)
            xy2  = gt.llh2utm([geom[i,5], geom[i,6]], self.origin)
            
            # Length of the fault
            length = np.linalg.norm(xy1-xy2)
            
            # strike angle of the fault
            strik = rad2deg(np.arctan2((xy2[1]-xy1[1]), (xy2[0]-xy1[0])))
            strik = 90-strik
            
            # local coordinates with respect to origin
            dE = 0.5*(xy1[0] + xy2[0])
            dN = 0.5*(xy1[1] + xy2[1])
            
            dis_geom[i,:] = np.array([length,
                                   geom[i,0], 
                                   geom[i,1], 
                                   geom[i,2], 
                                   strik, 
                                   dE, 
                                   dN, 
                                   geom[i,7], 
                                   geom[i,8], 
                                   geom[i,9]])    
        return dis_geom
        
    def PatchFault(self, m):
        '''
        This functions discretizes a fault model into i*j distinct patches.
    
        Input:
            m = numpy array, 1x10 vector defined as follows
                m[0] = fault length along the strike direction (km)
                m[1] = fault width in dip direction (km)
                m[2] = depth of lower edge of fault (km)
                m[3] = dip angle, from the horizontal (degrees) 
                m[4] = strike, clockwise from N (degrees)
                m[5] = East offset of midpoint of lower edge from origin (km)
                m[6] = North offset of midpoint of lower edge from origin (km) 
            i = number of patches along fault length 
            j = number of patches along fault width
    
        Output:
            pm = mx7 matrix of patch models, where m=i*j
        '''
        
        i, j    = self.nsegs, self.ndeps
        n       = i*j
        dip     = deg2rad(m[0,3])
        strike  = deg2rad(-m[0,4])
        sin_dip = sin(dip)
        cos_dip = cos(dip)
        iw      = m[0,0]/i
        jw      = m[0,1]/j
        iis     = np.arange(1, i+1)
        jjs     = np.transpose(np.arange(1, j+1))
        
        
        # Calculate midpoints, depths of each patch    
        p         = np.mat(cos_dip*(jw*jjs - m[0,1])).T * ones((1, i))   
        q         = np.mat(ones((j,1)) * ((iw*iis) - 0.5 * (m[0,0] + iw))) 
        r         = np.mat((m[0,2] - jw * sin_dip * (j - jjs))).T * ones((1, i))
        p         = p.T.flatten().T
        q         = q.T.flatten().T
        r         = r.T.flatten().T
    
        mp        = np.array([p, q, r]).reshape(3, p.size).T
        
        # Adjust midpoints for strike
        R         = np.array([[cos(strike), -sin(strike), 0],
                           [sin(strike), cos(strike), 0],
                           [ 0, 0, 1]])
    
        mp        = np.dot(mp, R.T)
    
        # Adjust midpoints for offset from origin
        mp[:,0]   = mp[:,0] + m[0,5]
        mp[:,1]   = mp[:,1] + m[0,6]
    
        # Form patch-models
        pm        = zeros((n,7))
        pm[:,0]   = ones(n) * iw
        pm[:,1]   = ones(n) * jw
        pm[:,2]   = mp[:,2]
        pm[:,3:5] = ones((n,1))*m[:,3:5]
        pm[:,5:7] = mp[:,0:2]
        
    
        # sort the pm
        idx = np.lexsort((pm[:,5], pm[:,2]))
        pm  = pm[idx]    
        return pm
        
    def disgeom2geom(self, dis_geom):
        '''
        Convert fault geometry from xy into llh, the output format is the same as FaultGeom
        Input:
            dis_geom    = numpy array, 1x10 vector defined as follows
            dis_geom[0] = fault length along the strike direction (km)
            dis_geom[1] = fault width in dip direction (km)
            dis_geom[2] = depth of lower edge of fault (km)
            dis_geom[3] = dip angle, from the horizontal (degrees) 
            dis_geom[4] = strike, clockwise from N (degrees)
            dis_geom[5] = East offset of midpoint of lower edge from origin (km)
            dis_geom[6] = North offset of midpoint of lower edge from origin (km)        
        Output:
            geom        = numpy array, nf*7 vector defined as following
                            Width  Depth  Dip  SLon SLat ELon Elat
                            ( km    km    deg  deg  deg  deg  deg)
            
        '''
        geom = zeros((self.nf,7))
        
        if self.nf > 0:
            for i in range(0, self.nf):
                # convert middle of fault plane to upper corner (lat1, lon1)
                len_x     = dis_geom[i,0]*sin(deg2rad(dis_geom[i,4]))
                len_y     = dis_geom[i,0]*cos(deg2rad(dis_geom[i,4]))
                tmpx      = dis_geom[i,5] - 0.5*len_x
                tmpy      = dis_geom[i,6] - 0.5*len_y
                
                fltllh    = gt.xy2llh([tmpy, tmpx], self.origin)
                geom[i,3] = fltllh[0,0]
                geom[i,4] = fltllh[0,1]
                
                tmpx      = dis_geom[i,5] + 0.5*len_x
                tmpy      = dis_geom[i,6] + 0.5*len_y
                fltllh    = gt.xy2llh([tmpy, tmpx], self.origin)
                geom[i,5] = fltllh[0,0]
                geom[i,6] = fltllh[0,1]
                geom[i,0] = dis_geom[i,1]
                geom[i,1] = dis_geom[i,2]
                geom[i,2] = dis_geom[i,3]
        return geom
        
    def FaultGeom2AllVertex(self):
        '''
        convert FaultGeom of DISMODEL/GINVPY to get all vertex of rectangular
        Written by Zhao Bin, Institute of Seismology, CEA. July 8, 2017
        Written by Zhao Bin, Institute of Seismology, CEA. July 31, 2017. Output neu vertex
        Mod by Zhao Bin, Feb. 21, 2019. Fix bug for calculating rake angle
        Mod by Zhao Bin, Sep. 26, 2022. adding parameter origin
    
        Input:
            geom    = numpy array, shape(m,10). The structure is same as FaultGeom
            disgeom = numpy array, shape(m,10). Leng, wid, dep, dip, strk, dE, dN, ss, ds, op
    
        Output:
            flt_elem_all = numpy array, shape(m,33). Leng, wid, dip, strk, ss, ds, op, ts, rake, top_left_llh, top_right_llh, bottom_right_llh, bottom_left_llh, top_left_neu, top_right_neu, bottom_right_neu, bottom_left_neu
        '''

        N            = 6378206.4
        e            = 1.0/298.257
        felems       = zeros((self.nf, 33))
    
        for i in range(self.nf):
            W     = sqrt(1-e**2*sin(deg2rad(self.geom_grid[i,3])))
            Rp    = (N/W)*cos(deg2rad(self.geom_grid[i,3]))
            Rm    = N*(1-e**2)/W**3
    
            dip   = self.dis_geom_grid[i,3]
            strk  = self.dis_geom_grid[i,4]
            dipr  = deg2rad(dip)
            strkr = deg2rad(strk)
            dipdir= np.pi - strkr
            ss    =-self.dis_geom_grid[i,7]
            ds    = self.dis_geom_grid[i,8]
            op    = self.dis_geom_grid[i,9]
            ts    = sqrt(ss**2 + ds**2)
            rake  = 0.0 if ts == 0 else rad2deg(np.arctan2(ds, ss))

            length   = self.dis_geom_grid[i,0]
            width    = self.dis_geom_grid[i,1]
            bottom   =-self.dis_geom_grid[i,2]
            top      = bottom + width*sin(dipr)
            
            gll      = complex(self.geom_grid[i,3], self.geom_grid[i,4])
            glr      = complex(self.geom_grid[i,5], self.geom_grid[i,6])
    
            # factor of 1e3 converts km to meters
            horproj  = self.dis_geom_grid[i,1]*cos(dipr)*1e3
            offset   = rad2deg(horproj)*complex( sin(dipdir)/Rm, cos(dipdir)/Rp )
            
            # convert offset from radians to degrees
            gul      = gll + offset
            gur      = glr + offset

            def llh_to_enu(c, z):
                e, n = gt.llh2utm([c.real, c.imag], self.origin)
                return e, n, z
    
            # convert llh to neu
            gll_e, gll_n, gll_u = llh_to_enu(gll, bottom)
            glr_e, glr_n, glr_u = llh_to_enu(glr, bottom)
            gur_e, gur_n, gur_u = llh_to_enu(gur, top)
            gul_e, gul_n, gul_u = llh_to_enu(gul, top)
            
            
            # assign to flt_elem_all
            felems[i,:] = np.array([
                length, width, strk, dip, rake,
                ss, ds, op, ts,
                gll.imag, gll.real, bottom, 
                glr.imag, glr.real, bottom, 
                gur.imag, gur.real, top, 
                gul.imag, gul.real, top, 
                gll_e, gll_n, gll_u, 
                glr_e, glr_n, glr_u, 
                gur_e, gur_n, gur_u, 
                gul_e, gul_n, gul_u])
        
        columns = ["length", "width", "strike", "dip", "rake",
                   "strike_slip", "dip_slip", "tensile_slip", "total_slip",
                   "lt_lon", "lt_lat", "lt_dep",
                   "rt_lon", "rt_lat", "rt_dep",
                   "rb_lon", "rb_lat", "rb_dep",
                   "lb_lon", "lb_lat", "lb_dep",
                   "lt_e", "lt_n", "lt_u",
                   "rt_e", "rt_n", "rt_u",
                   "rb_e", "rb_n", "rb_u",
                   "tl_e", "tr_n", "br_u"]
        df_felems = pd.DataFrame(felems, columns=columns)
        return df_felems

    @staticmethod
    def geom2x(geom, origin):
        '''
        Written by Zhao Bin when his is in UC Berkeley, Feb 20 2016

        Input:
            geom   = Fault Geometry 
                    Width  Depth  Dip  SLon SLat ELon Elat Strike-Slip Dip-slip Openning  
            origin = [lat0, lon0]
        Output:
            dis_geom - numpy array, containing Length Width Depth Dip Strike midLat midLon Strike-Slip Dip-slip Openning

    '''
        # convert array to matrix
        geom     = np.asarray(geom)
        nf, ncol = geom.shape
        x        = zeros((nf,7))
   
	    # for each subfaults
        for i in range(nf):
            xy1    = gt.llh2utm([geom[i,3], geom[i,4]], origin)
            xy2    = gt.llh2utm([geom[i,5], geom[i,6]], origin)
            length = np.linalg.norm(xy1-xy2)
	        
	        
	        # strike angle of the fault
            strike = rad2deg(np.atan2((xy2[1]-xy1[1]),
                                     (xy2[0]-xy1[0])))
            strike = 90-strike
	        
	        # midpoint
            midlat = 0.5*(geom[i,3] + geom[i,5])
            midlon = 0.5*(geom[i,4] + geom[i,6])
            x[i,:] = np.array([length, geom[i,0], geom[i,1], geom[i,2], strike, midlat, midlon])    
	    
        return x
        
    @staticmethod
    def x2geom(x, origin):
        '''
        Convert model vector x to fault geometry and dis_geom arrays
        Modified by Zhao Bin, Aug. 29, 2018.
        
        Input:
            x      = number of faults * 7 array
                     length, width, depth, dip, strike, lat, lon
            origin = [lat0, lon0]
        
        '''

        fmt = 7*"%10.4f"
        print(fmt %(x[0,0], x[0,1], x[0,2], x[0,3], x[0,4], x[0,5], x[0,6]))
        nf       = len(x)
        geom     = zeros((nf, 7))
        dis_geom = zeros((nf, 7))
    
        N        = 6378206.4
        e        = 1./298.257
    
        for i in range(nf):
            length, width, depth, dip, strike, lat, lon = x[i]
            theta  = deg2rad(90 - strike)
            dE     = 0.5*length*cos(theta)
            dN     = 0.5*length*sin(theta)
    
            W        = sqrt(1-e**2*sin(deg2rad(lat)))
            Rp       = (N/W)*cos(deg2rad(lat))
            Rm       = N*(1-e**2)/W**3
    
            dlat     = rad2deg(1.0e3*dN/Rm)
            dlon     = rad2deg(1.0e3*dE/Rp)
    
            lat1, lon1 = lat-dlat, lon-dlon
            lat2, lon2 = lat+dlat, lon+dlon
            geom[i,:]  = [width, depth, dip, lat1, lon1, lat2, lon2]
            
            xy1 = gt.llh2utm([geom[i,3], geom[i,4]], origin)
            xy2 = gt.llh2utm([geom[i,5], geom[i,6]], origin)
    
            midE    = 0.5*(xy1[0] + xy2[0])
            midN    = 0.5*(xy1[1] + xy2[1])
    
            dis_geom[i,:] = [length, width, depth, dip, strike, midE, midN] 
        return geom, dis_geom


    @staticmethod
    def GetFaultBounds(filename, nf):
	    '''
	    Read in FalutBound file to get FaultGeom bound for GeomEst.py
	    By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley. Dec. 1, 2016
        Modified by Zhao Bin, write to 'coulomb.inp'

	    Input:
	    faultbound    -  file name of FaultBound
	                  -   Len, Wid, depth, dip, strike, latitude, longitude
	    nf            -  number of fault
	    Output:
	    xlb           -  lower bound
	    xub           -  upper bound
	    
	
	    # check the file is exit
        if not os.path.isfile(filename):
        if not os.path.isfile(filename):
	       logging.warning(f'{filename} is not exist! Please input another file!')
           sys.exit()
	
	    # load the FaultBound
	    fltbnds = np.genfromtxt(filename, comments='#')
	    if len(fltbnds)%2 != 0:
	        logging.warning('{} has uneven lines for fault bounds!'.format(filename))
	        sys.exit()
	    
	    # get the size of the matrix
	    [nf2, ncol] = fltbnds.shape
	
	    # 
	    xlb = np.zeros((nf2//2,ncol))
	    xub = np.zeros((nf2//2,ncol))
	    if nf2 == 2*nf:
	        for i in range(nf):
	            for j in range(0, 7):
	                xlb[i,j] = fltbnds[2*i,j]
	                xub[i,j] = fltbnds[2*i+1,j]
	#               if xlb[i,6] > 180:
	#                   xlb[i,6]=xlb[i,6]-180
	#               if xub[i,6] > 180:
	#                   xub[i,6]=xub[i,6]-180
	    
	    return xlb, xub'''

    def FaultGeom2Coulomb(self):
        '''
        Convert FaultGeom file into the format Coulomb 
        By Zhao Bin, Mar 1st, 2016
        
        Notice that the in Coulomb x(east) and y(north)
        '''
    
        # compute some values for ouput file
        geom     = self.geom_grid
        depth    = np.mean(geom[:,1])
        nf       = len(geom)
        
        # output using the format of Coulomb program
        with open('coulomb.inp', 'w') as fid:
            fid.write("This is an example for Coulomb 3.2.\n")
            fid.write("You can use first two lines (rows) for some comments or notes.\n")
            fid.write("#reg1=  0  #reg2=  0  #fixed=  %d  sym=  1" %(nf))
            fid.write(" PR1=       0.250     PR2=       0.250   DEPTH=     %.3f\n" %(depth))
            fid.write("  E1=      8.000e+05   E2=      8.000e+05\n")
            fid.write("XSYM=       .000     YSYM=       .000\n")
            fid.write("FRIC=          0.400\n")
            fid.write("S1DR=         19.000 S1DP=         -0.010 S1IN=        100.000 S1GD=          0.000\n")
            fid.write("S2DR=         89.990 S2DP=         89.990 S2IN=         30.000 S2GD=          0.000\n")
            fid.write("S3DR=        109.000 S3DP=         -0.010 S3IN=          0.000 S3GD=          0.000\n")
            fid.write("\n")
            fid.write("  #   X-start    Y-start     X-fin      Y-fin   Kode  rt.lat    reverse   dip angle     top      bot\n")
            fid.write("xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx xxxxxxxxxx\n")
    
        # for each subfault
        with open('coulomb.inp', 'a') as fid:
            for i in range (nf):            
                if geom[i,2] > 180: rdip = deg2rad(geom[i,2]-180)
                
                ll_1end  = np.array([geom[i,3], geom[i,4]])
                ll_2end  = np.array([geom[i,5], geom[i,6]])
                xy_1end  = gt.llh2localxy(ll_1end, self.origin)
                xy_2end  = gt.llh2localxy(ll_2end, self.origin)
                fid.write("%3d %10.4f %10.4f %10.4f %10.4f %3d %10.4f %10.4f %10.4f %10.4f %10.4f\n" 
                    %(1, xy_1end[0], xy_1end[1], xy_2end[0], xy_2end[1], 100, 
                      geom[i,7]/1e3, geom[i,8]/1e3, geom[i,2]-180, geom[i,1],
                      geom[i,1]+geom[i,0]*sin(rdip)))
    
        # get lbllh and rtllh
        lbllh = gt.xy2llh([-200, -200], self.origin)
        rtllh = gt.xy2llh([200, 200], self.origin)
    
        # output other information
        with open('coulomb.inp', 'a') as fid:
            fid.write("\n")
            fid.write("     Grid Parameters\n")
            fid.write("  1  ----------------------------  Start-x =%17.4f\n"  %-500)
            fid.write("  2  ----------------------------  Start-y =%17.4f\n"  %-500)
            fid.write("  3  --------------------------   Finish-x =%17.4f\n"  %500)
            fid.write("  4  --------------------------   Finish-y =%17.4f\n"  %500)
            fid.write("  5  -----------------------   x-increment =%17.4f\n"  %50 )
            fid.write("  6  -----------------------   y-increment =%17.4f\n"  %50 )
            fid.write("\n")
            fid.write("     Size Parameters\n")
            fid.write("  1  --------------------------  Plot size =%16.7f\n"   %6.0)
            fid.write("  2  --------------  Shade/Color increment =%16.7f\n"   %1.0)
            fid.write("  3  ------  Exaggeration for disp.& dist. =%16.7f\n"   %10000.0)
            fid.write("\n")
            fid.write("     Cross section default\n")
            fid.write("  1  ----------------------------  Start-x =%16.7f\n"   %-200)
            fid.write("  2  ----------------------------  Start-y =%16.7f\n"   %-200)
            fid.write("  3  ---------------------------  Finish-x =%16.7f\n"   %200)
            fid.write("  4  ---------------------------  Finish-y =%16.7f\n"   %200)
            fid.write("  5  ------------------  Distant-increment =%16.7f\n"   %20)
            fid.write("  6  ----------------------------- Z-depth =%16.7f\n"   %0.0)
            fid.write("  7  ------------------------  Z-increment =%16.7f\n"   %1.0)
            fid.write("   Map info\n")
            fid.write("  1  ---------------------------  min. lon =%16.7f\n"   %lbllh[:,1])
            fid.write("  2  ---------------------------  max. lon =%16.7f\n"   %rtllh[:,1])
            fid.write("  3  ---------------------------  zero lon =%16.7f\n"   %self.origin[1])
            fid.write("  4  ---------------------------  min. lat =%16.7f\n"   %lbllh[:,0])
            fid.write("  5  ---------------------------  max. lat =%16.7f\n"   %rtllh[:,0])
            fid.write("  6  ---------------------------  zero lat =%16.7f\n"   %self.origin[0])
            fid.write("  7  ------------------------  Z-increment =%16.7f\n"   %1.0)
        return
 
    def FaultGeom2EDCMP(self):
        '''
        Convert FaultGeom of DISMODEL program into the format of EDCMP progam
        Written By Zhao Bin, Mar 1, 2016
        Mod by Zhao Bin, Mar 3, 2016 : output to file
    
        Notice:
             In EDCMP x = north  and y = east
             The unit is meter
    
        Input file:
            file name of FaultGeom
        Output file:
            file in the format of EDCMP program
        '''
    
        felem = self.FaultGeom2AllVertex()
    
        # open output file
        with open('edcmp.fault', 'w') as fid:
            fid.write("# origin: {:10.4f} {:10.4f} ".format(self.origin[0], self.origin[1]))
            for i, row in felem.iterrows():
                strike = row['strike'] % 360
                dip    = row['dip'] % 360
      
                # total slip, convert mm in FautlGeom to m in EDCMP
                slip   = row['total_slip']/1e3
                
                length = row['length']*1e3
                width  = row['width']*1e3
                north  = row['lt_n']*1e3
                east   = row['lt_e']*1e3
                dep    =-row['lt_u']*1e3
            
                fid.write("%4d %10.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3f %10.3f %10.3f\n" 
                    %(i+1, slip, north, east, dep, length, width, strike, dip, row['rake']))


    def FaultGeom2VISCO1D(self):
        '''
        Convert FaultGeom of DISMODEL program into the format of VISCO1D progam
        Written By Zhao Bin, Apr 8, 2016
        MOD by Zhao Bin, Nov 2, 2016. correct for Rp 
        MOD by Zhao Bin, July 8, 2017. call FaultGeom2AllVertex function
    
        Input file:
            file name of FaultGeom
    
        Output file:
            file in the format of VISCO1D program, the following is a sample of the file:
    
            15.00 0.0 30.0                      % lower and upper depth of fault in km, and dip angle in deg
            1975. 1975. 2080.82 1.0             % inital time, start and end time in decimal year, and ISRATE
            1                                   % number of fault plane
            0.899361 0.233661 200. 0. 90. 100.  % latitude(deg) and longitude(deg) of the lowermost corner of the fault closest to the strike direction, length(km), strke angle(deg), rake angle(deg), slip(cm)
        '''
    
        ################################# convert them to VISCO1D format ################################
        felem = self.FaultGeom2AllVertex()
        for _, r  in felem.iterrows():
            print("%10.3f %10.3f %10.3f %8.2f %8.2f %10.3f %10.3f %10.3f %8.2f" 
                %(r['rb_lat'], r['rb_lon'], r['width'], r['strike'], r['rake'],
                  r['total_slip']/10, -r['lt_dep'], -r['lb_dep'], r['dip']))
        return


    def FaultGeom2VISCO2PT5D(self):
        '''
        Convert FaultGeom of DISMODEL program into the format of VISCO2.5D program
        Written By Zhao Bin, June 10, 2016
        MOD by Zhao Bin, Nov 2, 2016. correct for Rp
    
        Input file:
            file name of FaultGeom
    
        Output file:
            file in the format of VISCO2PT5D program, the following is a sample of the file:
    
            1                                   % number of fault plane
            15.00 0.0 30.0                      % lower and upper depth of fault in km, and dip angle in deg
            0.899361 0.233661 200. 0. 90. 100.  % latitude(deg) and longitude(deg) of the lowermost corner of the fault closest to the strike direction, length(km), strke angle(deg), rake angle(deg), slip(cm)
        '''
        felem  = self.FaultGeom2AllVertex()
        active = felem[felem['total_slip']>0] 
    
        print("{}".format(len(active)))
        for _, row in active:
            if row['dip']<180:
                dip  = row['dip']
                rake = row['rake']
            else:
                dip  = row['dip']-180
                rake = row[rake]-90
            top_dep = -row['lt_dep']
            bot_dep = -row['lb_dep']
            print(f"{top_dep:10.3f} {bot_dep:10.3f} {dip:5.2f}")
            print(f"{row['lb_lat']:10.3f} {row['lb_lon']:10.3f} {row['length']:10.3f} "
              f"{row['strike']:8.2f} {rake:8.2f} {row['total_slip']/10.:10.3f}")
        return


    def FaultGeom2RELAX(self):
        '''
        Convert FaultGeom of DISMODEL program into the format of RELAX program
        Written By Zhao Bin, Apr 8, 2016
        MOD by Zhao Bin, Jan 14, 2017. rake angle
    
        Input file:
            file name of FaultGeom -- contain slip
    
        Output file:
            fault slip model for RELAX software
            No.  Slip  X1  X2  X3  Length  Width  Strike  Dip  Rake 
            #    m     km  km  km  km      km     deg     deg  deg 
        '''
        felem = self.FaultGeom2AllVertex()
        # open output file
        with open('relax.fault', 'w') as fid:
            fid.write("# origin: {:10.4f} {:10.4f}\n".format(self.origin[0], self.origin[1]))
            # for each patches
            for i, row in felem.iterrows():
                fid.write("%4d %10.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3f %10.3f %10.3f\n" 
                    %(i+1, row['total_slip']/1e3,
                      row['lt_n']*1e3,
                      row['lt_e']*1e3,
                     -row['lt_u']*1e3,
                      row['length']*1e3,
                      row['width']*1e3,
                      row['strike'],
                      row['dip'] % 180,
                      row['rake']))
        return
                      
                         
    def FaultGeom2PSCMP(self):
        '''
        Convert FaultGeom of DISMODEL program into the format of PSCMP program
        Written By Zhao Bin, Apr 10, 2016
    
        Input file:
            file name of FaultGeom -- contain slip
    
        Output file:
            fault slip model for PSCMP software
            n   O_lat   O_lon   O_depth length  width strike dip   np_st np_di start_time
            [-] [deg]   [deg]   [km]    [km]     [km] [deg]  [deg] [-]   [-]   [day]
            pos_s   pos_d   slp_stk slp_ddip open
            [km]    [km]    [m]     [m]      [m]
        '''
        
    
        np_st      = 1
        np_di      = 1
        start_time = 0
        
        felem = self.FaultGeom2AllVertex()
        for i, r in felem.iterrows():
            strike = r['strike'] % 360
            dip    = r['dip'] % 180
            print("%5d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %5d %5d %5d"
                %(i+1, r['lt_lat'], r['lt_lon'], -r['lt_dep'],
                     r['length'], r['width'],
                     strike, dip,
                     np_st, np_di, start_time))
            print("     %10.3f %10.3f %10.3f %10.3f %10.3f" 
                %(r['length']/2.0, r['width']/2,
                  r['strike_slip']/1e3, -r['dip_slip']/1e3, r['tensile_slip']/1e3))
        return
            
    def FaultGeom2GMT(self):
        '''
        Convert fault slip model in standard FaultGeom to GMT plotting format
        Written by Zhao Bin, @ UC Berkeley April 22 2016
        '''
        def write_gmtfile(filename, slip_col, df):
            """Write a GMT fault patch file for a given slip column"""
            with open(filename, 'w') as fid:
                for _, row in df.iterrows():
                    fid.write(f"> -Z{row[slip_col]:.6f}\n")
                    # write corners: top-left, top-right, bottom-right, bottom-left
                    corners = [
                        ("lt_lon", "lt_lat", "lt_dep"),
                        ("rt_lon", "rt_lat", "rt_dep"),
                        ("rb_lon", "rb_lat", "rb_dep"),
                        ("lb_lon", "lb_lat", "lb_dep")
                    ]
                    for lon_col, lat_col, dep_col in corners:
                        fid.write(f"{row[lon_col]:10.5f} {row[lat_col]:10.5f} {row[dep_col]:10.4f}\n")

        felem = self.FaultGeom2AllVertex()

        # Define output files and corresponding slip columns
        slip_mapping = {
                "ss_slip_llh.gmtlin": "strike_slip",
                "ds_slip_llh.gmtlin": "dip_slip",
                "ts_slip_llh.gmtlin": "total_slip"}

        # Write each slip type to its file
        for fname, slip_col in slip_mapping.items():
            write_gmtfile(fname, slip_col, felem)
        return

    def Moment(self, shear_modulus=3e10):
        '''
        '''

        Mo   = np.zeros(self.nf)
        geom = self.geom_grid
        slip = np.sqrt(geom[:,7]**2+geom[:,8]**2)
        
        for i in range(self.nf):
            Mo[i] = 1e3*self.dis_geom_grid[i,0] * \
                    1e3*self.dis_geom_grid[i,1] * \
                    slip[i] * shear_modulus
        
        Mo_total = np.sum(Mo)/1e3
        Mw_total = 2.0/3.0*np.log10(Mo_total) - 6.067

        return Mo_total, Mw_total

    def FaultGeom2VTK(self, scale=1.0):
        '''
        Write the fault geometry and slip distribution to VTK file for Paraview.
        '''
        felem = self.FaultGeom2AllVertex()
        with open('faultgeom.vtp', 'w') as fid:
            fid.writelines(["<?xml version=\"1.0\"?>\n",
                       "<VTKFile type=\"PolyData\" version=\"0.1\">\n",
                       "  <PolyData>\n"])

            for i, r in felem.iterrows():
                fid.write("    <Piece NumberOfPoints=\"4\" NumberOfPolys=\"1\">\n")
                fid.write("      <Points>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"Fault Patch\" NumberOfComponents=\"3\" format=\"ascii\">\n")
                fid.write("         {} {} {}\n".format(r['lt_e']*scale, r['lt_n']*scale, r['lt_u']*scale))
                fid.write("         {} {} {}\n".format(r['re_e']*scale, r['rt_n']*scale, r['rt_u']*scale))
                fid.write("         {} {} {}\n".format(r['rb_e']*scale, r['rb_n']*scale, r['rb_u']*scale))
                fid.write("         {} {} {}\n".format(r['lb_e']*scale, r['lb_n']*scale, r['lb_u']*scale))
                fid.write("         </DataArray>\n")
                fid.write("      </Points>\n")
                fid.write("      <Polys>\n")
                fid.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"{}\">\n".format(self.nf-1))
                fid.write("0 1 2 3\n")
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"4\" RangeMax=\"4\">\n")
                fid.write("          4\n")
                fid.write("        </DataArray>\n")
                fid.write("      </Polys>\n")
                fid.write("      <CellData Scalar=\"geometry\">\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"strike\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(r['strike']))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"dip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(r['dip']))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"strike slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(r['strike_slip']))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"dip slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(r['dip_slip']))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(r['total_slip']))
                fid.write("        </DataArray>\n")
                fid.write("      </CellData>\n")
                fid.write("    </Piece>\n")
            fid.write("  </PolyData>\n")
            fid.write("</VTKFile>\n")

    def plot_faultgeom(self, afsfile="", coordtype="llh", show_idx=False, azimuth=-40, elevation=40, sarfile="", gpsfile=""):
        """"""
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        felem = self.FaultGeom2AllVertex()
        slip  = felem['total_slip'] / felem["total_slip"].max()
    
        if coordtype == 'llh':
            x = felem[['lt_lon', 'rt_lon', 'rb_lon', 'lb_lon']]
            y = felem[['lt_lat', 'rt_lat', 'rb_lat', 'lb_lat']]
            z = felem[['lt_dep', 'rt_dep', 'rb_dep', 'lb_dep']]
        elif coordtype == 'enu':
            x = felem[['lt_e', 'rt_e', 'rb_e', 'lb_e']]
            y = felem[['lt_n', 'rt_n', 'rb_n', 'lb_n']]
            z = felem[['lt_u', 'rt_u', 'rb_u', 'lb_u']]

        # Create 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.view_init(elev=elevation, azim=azimuth)
        ax.scatter(x, y, z, s=0.0)

        def plot_aftershocks(ax, afsfile):
            """Plot aftershocks if the file exists."""
            if os.path.exists(afsfile):
                afs = np.genfromtxt(afsfile)
                ax.scatter(afs[:, 0], afs[:, 1], afs[:, 2], s=4, marker='s')
            else:
                print('Input aftershock catalog does not exist!')
                
        def plot_sar(ax, sarfile):
            """Plot aftershocks if the file exists."""
            if os.path.exists(sarfile):
                sar = np.genfromtxt(sarfile)
                s   = ax.scatter(sar[:,0], sar[:,1], np.zeros(len(sar)), c=sar[:,2], cmap=plt.cm.rainbow)
                plt.colorbar(s, shrink=0.8, pad=0.1, label='LOS (mm)')
            else:
                print('Input sar file does not exist!')

        def plot_gps(ax, gpsfile, scale=1000):
            """Plot aftershocks if the file exists."""
            if os.path.exists(gpsfile):
                gps = np.genfromtxt(gpsfile)
                ax.quiver(gps[:,0], gps[:,1], np.zeros(len(gps)), gps[:,2], gps[:,3], np.zeros(len(gps)))
            else:
                print('Input gps file does not exist!')

        # Plot fault patches
        for i, row in felem.iterrows():
            if coordtype == 'llh':
                lon = row[['lt_lon', 'rt_lon', 'rb_lon', 'lb_lon']]
                lat = row[['lt_lat', 'rt_lat', 'rb_lat', 'lb_lat']]
                dep = row[['lt_dep', 'rt_dep', 'rb_dep', 'lb_dep']]
                verts = [list(zip(lon, lat, dep))]

            elif coordtype == 'enu':
                e = row[['lt_e', 'rt_e', 'rb_e', 'lb_e']]
                n = row[['lt_n', 'rt_n', 'rb_n', 'lb_n']]
                u = row[['lt_u', 'rt_u', 'rb_u', 'lb_u']]
                verts = [list(zip(e, n, u))]

            subfault = Poly3DCollection(verts, facecolors=plt.cm.jet(slip[i]), linewidths=1, alpha=0.3)
            subfault.set_linewidth(1)
            ax.add_collection3d(subfault)

            if show_idx:
                if coordtype == 'llh':
                    ax.text(np.mean(lon), np.mean(lat), np.mean(dep), str(i), color='red', fontsize=5)
                elif coordtype == 'enu':
                    ax.text(np.mean(n), np.mean(e), np.mean(u), str(i), color='red', fontsize=5)

        # Plot aftershocks
        plot_aftershocks(ax, afsfile)
        
        # Plot InSAR
        plot_sar(ax, sarfile)
        
        # Plot GPS
        plot_gps(ax, gpsfile)
        
        # Set labels and save plot
        if coordtype == 'llh':
            ax.set_xlabel('Longitude', fontsize=12, labelpad=12)
            ax.set_ylabel('Latitude', fontsize=12, labelpad=18)
        elif coordtype == 'enu':
            ax.set_xlabel('East (km)')
            ax.set_ylabel('North (km)')
        ax.set_zlabel('Depth (km)', fontsize=12)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        s  = ax.scatter(felem['lt_lon'], felem['lt_lat'], c=felem['total_slip']/1000.,cmap=plt.cm.jet, s=0.00001,lw=0)
        plt.colorbar(s, shrink=0.6, pad=0.1, label='Slip (m)')
        #   cb.set_label('Total slip (m)', fontsize=12)
        plt.show()

    def faultgeom2geojson(self, scale=1, filename="fault_patches.geojson"):
        """
        Export fault patches in DataFrame to a GeoJSON file.
        
        Parameters:
            felem : pd.DataFrame
                DataFrame containing at least the columns:
                ['lt_lon','lt_lat','rt_lon','rt_lat','rb_lon','rb_lat','lb_lon','lb_lat']
            filename : str
                Output GeoJSON filename
        """
        import json
        features = []
        felem = self.FaultGeom2AllVertex()
        for _, row in felem.iterrows():
            # Define polygon coordinates (GeoJSON expects [lon, lat])
            polygon = [
                [row['lt_lon'], row['lt_lat']],
                [row['rt_lon'], row['rt_lat']],
                [row['rb_lon'], row['rb_lat']],
                [row['lb_lon'], row['lb_lat']],
                [row['lt_lon'], row['lt_lat']],  # close the polygon
            ]

            feature = {
                "type": "Feature",
                "properties": {
                    "ss": row.get('strike_slip', 0.0),
                    "ds": row.get('dip_slip', 0.0),
                    "op": row.get('tensile_dip', 0.0),
                    "ts": row.get('total_slip', 0.0),
                    "rake": row.get('rake', 0.0),
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [polygon]
                }
            }
            features.append(feature)

        geojson_dict = {
            "type": "FeatureCollection",
            "features": features
        }

        # Write to file
        with open(filename, "w") as f:
            json.dump(geojson_dict, f, indent=2)

        print(f"GeoJSON saved to {filename}")

    def Change_Dip(self, top_dip, bottom_dip):
        """
        Change the dip of the fault
        """
        new_dip   = np.linspace(top_dip, bottom_dip, self.ndeps)
        dis_geom  = self.dis_geom_grid
        for i in range(self.ndeps):
            for j in range(self.nsegs):
                k    =  i*self.nsegs+j
                strk = self.dis_geom_grid[k, 4]
                dip  = new_dip[i]
                [x1, y1, z1] = self.dis_geom_grid[k, [5,6,2]]
                s     = np.array([np.sin(np.deg2rad(strk)), 
                                np.cos(np.deg2rad(strk)), 0.0])
                d     = np.array([-np.cos(np.deg2rad(strk))*np.cos(np.deg2rad(dip)), 
                                np.sin(np.deg2rad(strk))*np.cos(np.deg2rad(dip)), 
                                np.sin(np.deg2rad(dip))])

                x2    = x1 + self.dis_geom_grid[k,0]*s[0]
                y2    = y1 + self.dis_geom_grid[k,0]*s[1]
                z2    = z1

                x3    = x2 + self.dis_geom_grid[k,1]*d[0]
                y3    = y2 + self.dis_geom_grid[k,1]*d[1]
                z3    = z2 - self.dis_geom_grid[k,0]*d[2]

                x4    = x1 + self.dis_geom_grid[k,0]*d[0]
                y4    = y1 + self.dis_geom_grid[k,0]*d[1]
                z4    = z1 - self.dis_geom_grid[k,0]*d[2]

                self.dis_geom_grid[k,3] = dip

                if i+1<self.ndeps:
                    self.dis_geom_grid[k+self.nsegs, 5] = 0.5*(x3+x4)
                    self.dis_geom_grid[k+self.nsegs, 6] = 0.5*(y3+y4)
                    self.dis_geom_grid[k+self.nsegs, 2] = 0.5*(z3+z4)
        # Update self.geom_grid
        self.geom_grid = self.disgeom2geom(self.dis_geom_grid)        
        
    def Ouput_Faultgeom(self, slip=None, scale=1.0):
        """
        """
        fmt   = 7*'%11.3f'+3*'%13.3f'
        header  = " width_km   depth_km    dip_deg   lat0_deg   lon0_deg   lat1_deg   lon1_deg  slp_strk_mm  slp_ddip_mm  slp_tens_mm"
        if slip is None and self.nf>0:
            np.savetxt('output.faultgeom', self.geom_grid, fmt=fmt, header=header)
            logging.info('Output a new faultgeom file.')
        else:
            if self.nf == len(slip):
                out = np.column_stack((self.geom_grid[:,0:7], slip*scale))
                np.savetxt('output.faultgeom', out, fmt=fmt, header=header)
                logging.info('Output a new faultgeom file.')
        return

    def Merge_Faultgeoms(self, faultgeom_list):
        """
        """
        
        
if __name__ == '__main__':
    pass
