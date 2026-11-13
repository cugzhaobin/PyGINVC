#!/usr/bin/env python

import numpy as np
import pandas as pd
import os, logging
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
        self.nf            = 0
        self.doSubFault    = doSubFault
        self.origin        = np.asarray(origin)
        self.nsegs         = nsegs
        self.ndeps         = ndeps
        self.dis_geom_grid = np.array([])
        self.geom_grid     = np.array([])    

    def load_fault(self):
        # load the original fault geometry
        geom = self.LoadFaultGeom()

        # PatchFault subdivides one rectangular fault plane; broadcasting
        # nsegs/ndeps over several rows is not implemented.  Reject the
        # combination explicitly instead of failing later in PatchFault.
        if self.doSubFault and len(geom) > 1:
            raise ValueError(
                'doSubFault=True only supports a single-row fault geometry '
                'file, but %s contains %d rows. Split the rows into separate '
                'fault entries, or set doSubFault=False.'
                % (self.faultfile, len(geom)))

        # convert the fault into DISMODEL format
        if len(self.origin) != 2:
            self.origin = np.array([np.mean(geom[:,[3,5]]), np.mean(geom[:,[4,6]])])
        logging.info(f'origin is set to {self.origin[0]:10.4f} {self.origin[1]:10.4f}')
        dis_geom = self.SetOrigin(geom)

        # if we need subdivide the fault
        if not self.doSubFault:
            n_patch = int(self.nsegs) * int(self.ndeps)
            if n_patch != len(geom):
                raise ValueError(
                    'doSubFault=False requires nsegs*ndeps (%d*%d=%d) to '
                    'equal the number of rows in %s (%d).'
                    % (self.nsegs, self.ndeps, n_patch, self.faultfile,
                       len(geom)))
            self.dis_geom_grid = dis_geom
            self.geom_grid     = geom
            self.nf            = len(geom)
            logging.info(f'Fault geometry is loaded. There are {self.nf} patches in it.')
        else:
            pm = self.PatchFault(dis_geom)
            self.nf = len(pm)
            self.dis_geom_grid = np.column_stack((pm, zeros((self.nf,3))))
            self.geom_grid     = self.disgeom2geom(self.dis_geom_grid)
            logging.info(f'Fault geometry is loaded and subdivided into {self.nsegs} * {self.ndeps} patches.')
    
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
    
        # check the file exists
        if not os.path.isfile(self.faultfile):
            raise FileNotFoundError(
                'Fault geometry file does not exist: %s' % self.faultfile)

        # load the data
        geom = np.genfromtxt(self.faultfile, comments='#')
        if geom.size == 0:
            raise ValueError(
                'Fault geometry file is empty: %s' % self.faultfile)
        if geom.ndim == 1:
            geom = geom.reshape((1, len(geom)))
        if geom.shape[1] != 10:
            raise ValueError(
                'Fault geometry file must have 10 columns, got %d: %s'
                % (geom.shape[1], self.faultfile))
        if not np.isfinite(geom).all():
            raise ValueError(
                'Fault geometry contains non-numeric or non-finite values: %s'
                % self.faultfile)
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
        p         = np.outer(cos_dip*(jw*jjs - m[0,1]), ones(i))
        q         = np.outer(ones(j), (iw*iis) - 0.5 * (m[0,0] + iw))
        r         = np.outer(m[0,2] - jw * sin_dip * (j - jjs), ones(i))
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
            DataFrame, shape(m,33). The first nine columns are length, width,
            strike, dip, rake, slip components; then four corners in the order

                gll (lower-left), glr (lower-right),
                gur (upper-right), gul (upper-left)

            each as lon/lat/depth followed by east/north/up.  The FaultGeom
            input edge (SLat/SLon -> ELat/ELon) is the DOWN-DIP edge, so the
            corners derived from it are the lower ones and the up-dip corners
            (top) are one fault width away.
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
            ts    = sqrt(ss**2 + ds**2 + op**2)
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
        
        # felems[i] stores the corners in the order
        #   gll (lower-left), glr (lower-right), gur (upper-right),
        #   gul (upper-left) -- each as lon/lat/depth then east/north/up.
        columns = ["length", "width", "strike", "dip", "rake",
                   "strike_slip", "dip_slip", "tensile_slip", "total_slip",
                   "lb_lon", "lb_lat", "lb_dep",
                   "rb_lon", "rb_lat", "rb_dep",
                   "rt_lon", "rt_lat", "rt_dep",
                   "lt_lon", "lt_lat", "lt_dep",
                   "lb_e", "lb_n", "lb_u",
                   "rb_e", "rb_n", "rb_u",
                   "rt_e", "rt_n", "rt_u",
                   "lt_e", "lt_n", "lt_u"]
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

        # EDCMP source format:
        #   id slip north east dep length width strike dip rake
        # (north, east, dep) is EDCMP's "upper reference point for strike",
        # i.e. the upper-left (up-dip) corner, which is lt_* here.  The
        # FaultGeom dip convention is 180 degree larger than EDCMP's
        # (see docs), hence the modulo 180.  Distances are converted km -> m
        # and slip mm -> m.
        edcmp_flt = np.column_stack((
            np.arange(1, self.nf + 1),
            felem['total_slip'].to_numpy() / 1e3,
            felem['lt_n'].to_numpy() * 1e3,
            felem['lt_e'].to_numpy() * 1e3,
            -felem['lt_u'].to_numpy() * 1e3,
            felem['length'].to_numpy() * 1e3,
            felem['width'].to_numpy() * 1e3,
            felem['strike'].to_numpy() % 360,
            felem['dip'].to_numpy() % 180,
            felem['rake'].to_numpy()))

        # open output file
        with open('edcmp.fault', 'w') as fid:
            fid.write("# origin: {:10.4f} {:10.4f} \n".format(self.origin[0], self.origin[1]))
            fmt = "%4d %10.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3f %10.3f %10.3f"
            np.savetxt(fid, edcmp_flt, fmt=fmt)

        return edcmp_flt


    def FaultGeom2VISCO1D(self):
        '''
        Convert FaultGeom of DISMODEL program into the format of VISCO1D progam
        Written By Zhao Bin, Apr 8, 2016
        MOD by Zhao Bin, Nov 2, 2016. correct for Rp 
        MOD by Zhao Bin, July 8, 2017. call FaultGeom2AllVertex function
    
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
                  r['total_slip']/10, -r['lb_dep'], -r['lt_dep'], r['dip']))


    def FaultGeom2VISCO2PT5D(self):
        '''
        Convert FaultGeom of DISMODEL program into the format of VISCO2.5D program
        Written By Zhao Bin, June 10, 2016
        MOD by Zhao Bin, Nov 2, 2016. correct for Rp
    
        Output file:
            file in the format of VISCO2PT5D program, the following is a sample of the file:
    
            1                                   % number of fault plane
            15.00 0.0 30.0                      % lower and upper depth of fault in km, and dip angle in deg
            0.899361 0.233661 200. 0. 90. 100.  % latitude(deg) and longitude(deg) of the lowermost corner of the fault closest to the strike direction, length(km), strke angle(deg), rake angle(deg), slip(cm)
        '''
        felem  = self.FaultGeom2AllVertex()
        active = felem[felem['total_slip']>0] 
    
        print("{}".format(len(active)))
        for _, row in active.iterrows():
            if row['dip']<180:
                dip  = row['dip']
                rake = row['rake']
            else:
                dip  = row['dip']-180
                rake = row['rake']-90
            # same layout as VISCO1D: lower depth then upper depth, and the
            # lowermost corner closest to the strike direction (rb).
            lower_dep = -row['lb_dep']
            upper_dep = -row['lt_dep']
            print(f"{lower_dep:10.3f} {upper_dep:10.3f} {dip:5.2f}")
            print(f"{row['rb_lat']:10.3f} {row['rb_lon']:10.3f} {row['length']:10.3f} "
              f"{row['strike']:8.2f} {rake:8.2f} {row['total_slip']/10.:10.3f}")


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
                    # write corners: lower-left, lower-right, upper-right, upper-left
                    corners = [
                        ("lb_lon", "lb_lat", "lb_dep"),
                        ("rb_lon", "rb_lat", "rb_dep"),
                        ("rt_lon", "rt_lat", "rt_dep"),
                        ("lt_lon", "lt_lat", "lt_dep")
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

    def moment(self, slip, shear_modulus=3e10):
        '''
        '''
        area     = self.dis_geom_grid[:,0] * self.dis_geom_grid[:,1] * 1e6
        tol_slip = np.linalg.norm(slip, axis=1)
        Mo       = area * tol_slip * shear_modulus
        Mo_total = np.sum(Mo)/1000
        Mw_total = 2.0/3.0*np.log10(Mo_total) - 6.067

        return Mo_total, Mw_total

    def FaultGeom2VTK(self, filename='faultgeom.vtp', scale=1.0):
        '''
        Write the fault geometry and slip distribution to VTK file for Paraview.
        '''
        felem = self.FaultGeom2AllVertex()
        lines = []
        lines.extend(['<?xml version="1.0"?>',
                      '<VTKFile type="PolyData" version="0.1">',
                      '  <PolyData>'])

        for i, r in felem.iterrows():
            piece = f"""    <Piece NumberOfPoints="4" NumberOfPolys="1">
                   <Points>
                     <DataArray type="Float32" Name="Fault Patch" NumberOfComponents="3" format="ascii">
                       {r['lb_e']*scale} {r['lb_n']*scale} {r['lb_u']*scale}
                       {r['rb_e']*scale} {r['rb_n']*scale} {r['rb_u']*scale}
                       {r['rt_e']*scale} {r['rt_n']*scale} {r['rt_u']*scale}
                       {r['lt_e']*scale} {r['lt_n']*scale} {r['lt_u']*scale}
                     </DataArray>
                   </Points>
                   <Polys>
                     <DataArray type="Int32" Name="connectivity">0 1 2 3</DataArray>
                     <DataArray type="Int32" Name="offsets">4</DataArray>
                   </Polys>
                 <CellData Scalar="geometry">
                     <DataArray type="Float32" Name="strike">{r['strike']}</DataArray>
                     <DataArray type="Float32" Name="dip">{r['dip']}</DataArray>
                     <DataArray type="Float32" Name="strike_slip">{r['strike_slip']}</DataArray>
                     <DataArray type="Float32" Name="dip_slip">{r['dip_slip']}</DataArray>
                     <DataArray type="Float32" Name="total_slip">{r['total_slip']}</DataArray>
                 </CellData>
               </Piece>"""
            lines.append(piece)
        lines.extend(['  </PolyData>', '</VTKFile>'])
        
        with open(filename, 'w') as f:
            f.write('\n'.join(lines))

    def plot_faultgeom(self, value=None, afsfile="", coordtype="llh", show_idx=False,
                       azimuth=-40, elevation=40, sarfile="", gpsfile="",
                       show_map=True, scale=1000, show_slip=False, slip_scale=0.6,
                       slip_head=0.18, slip_head_angle=15.0, savefig=""):
        '''
        Plot the fault geometry with an optional cartopy map view.

        The map view is drawn on the left (coordtype='llh' only) and the 3D
        perspective view on the right. cartopy's GeoAxes is a 2-D projection
        and cannot be used as an mplot3d axis, so the two views are separate
        subplots. If cartopy is not installed, only the 3D view is drawn.

        show_slip : bool, optional
            Draw the slip vector of every patch: the tail of each arrow is at
            the patch centre and the whole arrow (shaft and head) lies in the
            fault plane. The magnitude is sqrt(strike_slip^2 + dip_slip^2) and
            the direction is built from the patch edges, so the arrows follow
            the fault surface exactly. Opening/tensile slip is ignored and
            zero-slip patches get no arrow. The default is False.
        slip_scale : float, optional
            Longest slip arrow as a fraction of the median patch size.
            The default is 0.6.
        slip_head : float, optional
            Arrow-head length as a fraction of the arrow length, applied to
            the 3D slip vectors. Smaller values give a more compact head.
            The default is 0.18.
        slip_head_angle : float, optional
            Half-angle (degrees) between each barb and the arrow shaft, applied
            to the 3D slip vectors. Smaller values make the barbs less splayed.
            The default is 15.0.
        savefig : str, optional
            If non-empty, save the figure to this path (at 300 dpi) before
            it is shown. The default is '', i.e. do not save.
        '''
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        felem = self.FaultGeom2AllVertex()
        ts    = np.asarray(felem['total_slip'], dtype=float)
        tmax0 = ts.max() if ts.size else 0.0
        slip  = ts / tmax0 if tmax0 > 0 else np.zeros_like(ts)

        if value is not None:
            v    = np.abs(np.asarray(value, dtype=float))
            vmax = v.max() if v.size else 0.0
            slip = v / vmax if vmax > 0 else np.zeros_like(v)
    
        if coordtype == 'llh':
            x = felem[['lb_lon', 'rb_lon', 'rt_lon', 'lt_lon']]
            y = felem[['lb_lat', 'rb_lat', 'rt_lat', 'lt_lat']]
            z = felem[['lb_dep', 'rb_dep', 'rt_dep', 'lt_dep']]
        elif coordtype == 'enu':
            x = felem[['lb_e', 'rb_e', 'rt_e', 'lt_e']]
            y = felem[['lb_n', 'rb_n', 'rt_n', 'lt_n']]
            z = felem[['lb_u', 'rb_u', 'rt_u', 'lt_u']]

        # cartopy is optional: without it (or for enu coordinates) only the
        # 3D view is drawn.
        ccrs = cfeature = None
        if show_map and coordtype == 'llh':
            try:
                import cartopy.crs as ccrs
                import cartopy.feature as cfeature
            except ImportError:
                print('cartopy is not installed; plotting the 3D view only.')

        def map_extent(lon, lat, margin=0.2):
            '''Return [lon_min, lon_max, lat_min, lat_max] with a margin.'''
            lon = np.asarray(lon, dtype=float).ravel()
            lat = np.asarray(lat, dtype=float).ravel()
            return [lon.min()-margin, lon.max()+margin,
                    lat.min()-margin, lat.max()+margin]

        def add_map_axes(fig, position, extent):
            '''Add a PlateCarree GeoAxes with coastlines, borders and gridlines.'''
            ax_map = fig.add_subplot(position, projection=ccrs.PlateCarree())
            ax_map.set_extent(extent, crs=ccrs.PlateCarree())
            # 10 m Natural Earth data: the default 110 m layers are far too
            # coarse for a fault map covering a fraction of a degree.
            ax_map.add_feature(cfeature.LAND.with_scale('10m'),
                               facecolor='#f7f7f2', zorder=0)
            ax_map.add_feature(cfeature.LAKES.with_scale('10m'),
                               facecolor='#dbe9f6', edgecolor='#8fa8bf',
                               linewidth=0.3, zorder=1)
            ax_map.add_feature(cfeature.COASTLINE.with_scale('10m'),
                               linewidth=0.5, zorder=2)
            ax_map.add_feature(cfeature.BORDERS.with_scale('10m'),
                               linewidth=0.4, linestyle=':', zorder=2)
            gl = ax_map.gridlines(draw_labels=True, linewidth=0.2, color='gray',
                                  alpha=0.4, linestyle='--')
            gl.top_labels   = False
            gl.right_labels = False
            return ax_map

        # Create the map view (left) and the 3D view (right)
        ax_map = None
        if ccrs is not None:
            fig    = plt.figure(figsize=(13, 6))
            ax_map = add_map_axes(fig, 121, map_extent(x.values, y.values))
            ax     = fig.add_subplot(122, projection='3d')
        else:
            fig = plt.figure()
            ax  = fig.add_subplot(projection='3d')
        ax.view_init(elev=elevation, azim=azimuth)
        ax.scatter(x, y, z, s=0.0)

        # True-scale 3D box.  The default box aspect is (4, 4, 3) and ignores
        # the data ranges, which maps ~0.4 deg of longitude and ~25 km of
        # depth onto the same box length and makes a dipping fault look like a
        # wide, flat ribbon.  Use the physical extents (km) instead, so one
        # kilometre has the same length on every axis.
        def _span(v):
            v = np.asarray(v, dtype=float).ravel()
            v = v[np.isfinite(v)]
            return float(v.max() - v.min()) if v.size else 0.0

        if coordtype == 'llh':
            lon_all = np.asarray(x, dtype=float).ravel()
            lat_all = np.asarray(y, dtype=float).ravel()
            coslat  = np.cos(np.deg2rad(np.nanmean(lat_all)))
            span_x  = _span(lon_all) * 111.32 * coslat
            span_y  = _span(lat_all) * 111.32
        else:
            span_x  = _span(x)
            span_y  = _span(y)
        span_z = _span(z)
        if min(span_x, span_y, span_z) > 0:
            ax.set_box_aspect((span_x, span_y, span_z))

        # Per-patch slip vectors, used by the optional vector overlay.
        #
        # Only the strike- and dip-slip components enter the vector; the
        # opening component is deliberately excluded so that every 3D arrow
        # lies exactly in the fault plane. The basis is built from the patch
        # edges:
        #   +strike = lt -> rt (right-hand strike),
        #   +dip    = lt -> lb (up-dip; DISMODEL stores the input edge as the
        #             deep edge, so lb is the shallow corner),
        # so felem['strike_slip'] (>0 left-lateral) and felem['dip_slip']
        # (>0 reverse) enter with positive signs.
        #
        # arrow_scale converts slip (mm) to a physical arrow length (km): the
        # largest-slip patch gets an arrow slip_scale * (median patch length).
        # The arrows start (tail) at the patch centre.
        slip_dir_enu = slip_mag = None
        arrow_scale  = 0.0
        if show_slip:
            def _unit(vec):
                norm = np.linalg.norm(vec)
                return vec / norm if norm > 0 else np.zeros_like(vec)

            ss = np.nan_to_num(np.asarray(felem['strike_slip'], dtype=float))
            ds = np.nan_to_num(np.asarray(felem['dip_slip'], dtype=float))
            slip_mag     = np.hypot(ss, ds)
            slip_dir_enu = np.zeros((len(felem), 3))
            sizes        = np.zeros(len(felem))
            for i, (_, row) in enumerate(felem.iterrows()):
                lt = np.array([row['lt_e'], row['lt_n'], row['lt_u']])
                rt = np.array([row['rt_e'], row['rt_n'], row['rt_u']])
                lb = np.array([row['lb_e'], row['lb_n'], row['lb_u']])
                t_vec = _unit(rt - lt)      # right-hand strike direction
                d_vec = _unit(lb - lt)      # up-dip direction
                slip_dir_enu[i] = _unit(ss[i]*t_vec + ds[i]*d_vec)
                sizes[i]        = np.linalg.norm(rt - lt)
            if slip_mag.max() > 0:
                arrow_scale = slip_scale * np.median(sizes) / slip_mag.max()

        def plot_aftershocks(ax, afsfile, transform=None):
            """Plot aftershocks if the file exists."""
            if os.path.exists(afsfile):
                afs = np.genfromtxt(afsfile)
                if transform is not None:
                    ax.scatter(afs[:, 0], afs[:, 1], s=4, marker='s',
                               transform=transform)
                else:
                    ax.scatter(afs[:, 0], afs[:, 1], afs[:, 2], s=4, marker='s')
            else:
                print('Input aftershock catalog does not exist!')

        def plot_sar(ax, sarfile, transform=None):
            """Plot InSAR LOS if the file exists."""
            if os.path.exists(sarfile):
                sar = np.genfromtxt(sarfile)
                if transform is not None:
                    s = ax.scatter(sar[:,0], sar[:,1], c=sar[:,2],
                                   cmap=plt.cm.rainbow, transform=transform)
                else:
                    s = ax.scatter(sar[:,0], sar[:,1], np.zeros(len(sar)),
                                   c=sar[:,2], cmap=plt.cm.rainbow)
                return s
            else:
                print('Input sar file does not exist!')

        def plot_gps(ax, gpsfile, scale=1000, transform=None):
            """Plot GPS vectors if the file exists."""
            if os.path.exists(gpsfile):
                gps = np.genfromtxt(gpsfile)
                if transform is not None:
                    ax.quiver(gps[:,0], gps[:,1], gps[:,2], gps[:,3],
                              transform=transform, scale_units='xy', scale=scale,
                              width=0.003)
                else:
                    ax.quiver(gps[:,0], gps[:,1], np.zeros(len(gps)),
                              gps[:,2], gps[:,3], np.zeros(len(gps)))
            else:
                print('Input gps file does not exist!')

        # Plot fault patches: 3D polygons on the right, map-view polygons on
        # the left. A PlateCarree axes uses lon/lat as data coordinates, so the
        # map polygons need no explicit transform.
        map_polys  = []
        map_colors = []
        arrows     = []   # (centre, slip vector, patch normal) in plot coords
        for i, row in felem.iterrows():
            if coordtype == 'llh':
                lon = row[['lb_lon', 'rb_lon', 'rt_lon', 'lt_lon']]
                lat = row[['lb_lat', 'rb_lat', 'rt_lat', 'lt_lat']]
                dep = row[['lb_dep', 'rb_dep', 'rt_dep', 'lt_dep']]
                verts = [list(zip(lon, lat, dep))]

            elif coordtype == 'enu':
                e = row[['lb_e', 'rb_e', 'rt_e', 'lt_e']]
                n = row[['lb_n', 'rb_n', 'rt_n', 'lt_n']]
                u = row[['lb_u', 'rb_u', 'rt_u', 'lt_u']]
                verts = [list(zip(e, n, u))]

            subfault = Poly3DCollection(verts, facecolors=plt.cm.jet(slip[i]), linewidths=1, alpha=0.3)
            subfault.set_linewidth(1)
            ax.add_collection3d(subfault)

            # Slip arrow with its tail at the patch centre. The displacement
            # is assembled directly from the patch edges in the coordinate
            # system being plotted:
            #   vec = (ss/Ls) * strike edge + (ds/Ld) * dip edge, scaled by
            #   arrow_scale (km/mm).
            # Both basis vectors are edges of the patch, so the arrow lies in
            # the fault plane exactly; no ENU -> lon/lat approximation is made
            # that could tilt it out of a steeply dipping plane. The result is
            # stored and drawn as explicit arrow segments after the loop.
            if (show_slip and slip_mag is not None and slip_mag[i] > 0
                    and arrow_scale > 0):
                if coordtype == 'llh':
                    p_lt = np.array([row['lt_lon'], row['lt_lat'],
                                     row['lt_dep']], dtype=float)
                    p_rt = np.array([row['rt_lon'], row['rt_lat'],
                                     row['rt_dep']], dtype=float)
                    p_lb = np.array([row['lb_lon'], row['lb_lat'],
                                     row['lb_dep']], dtype=float)
                    p_rb = np.array([row['rb_lon'], row['rb_lat'],
                                     row['rb_dep']], dtype=float)
                else:
                    p_lt = np.array([row['lt_e'], row['lt_n'], row['lt_u']],
                                    dtype=float)
                    p_rt = np.array([row['rt_e'], row['rt_n'], row['rt_u']],
                                    dtype=float)
                    p_lb = np.array([row['lb_e'], row['lb_n'], row['lb_u']],
                                    dtype=float)
                    p_rb = np.array([row['rb_e'], row['rb_n'], row['rb_u']],
                                    dtype=float)

                # Physical edge lengths (km) set the slip -> plot scale.
                e_s = np.array([row['rt_e'] - row['lt_e'],
                                row['rt_n'] - row['lt_n'],
                                row['rt_u'] - row['lt_u']], dtype=float)
                e_d = np.array([row['lb_e'] - row['lt_e'],
                                row['lb_n'] - row['lt_n'],
                                row['lb_u'] - row['lt_u']], dtype=float)
                L_s = np.linalg.norm(e_s)
                L_d = np.linalg.norm(e_d)

                p_s = p_rt - p_lt            # strike edge, plot coords
                p_d = p_lb - p_lt            # dip edge, plot coords
                centre = 0.25 * (p_lt + p_rt + p_lb + p_rb)
                if L_s > 0 and L_d > 0:
                    # Unit edge vectors in physical length.  In lon/lat/depth
                    # coordinates the horizontal components are degrees, so a
                    # geometric normal would mix units and make the barbs far
                    # too long.  The in-plane basis is therefore built from
                    # the two patch edges.
                    u_s   = p_s / L_s
                    u_d   = p_d / L_d
                    theta = np.arctan2(ds[i], ss[i])
                    v_hat = np.cos(theta)*u_s + np.sin(theta)*u_d
                    w_hat = -np.sin(theta)*u_s + np.cos(theta)*u_d
                    arrows.append((centre, slip_mag[i]*arrow_scale,
                                   v_hat, w_hat))

            if ax_map is not None:
                map_polys.append([(row['lb_lon'], row['lb_lat']),
                                  (row['rb_lon'], row['rb_lat']),
                                  (row['rt_lon'], row['rt_lat']),
                                  (row['lt_lon'], row['lt_lat'])])
                map_colors.append(plt.cm.jet(slip[i]))

            if show_idx:
                if coordtype == 'llh':
                    ax.text(np.mean(lon), np.mean(lat), np.mean(dep), str(i), color='red', fontsize=5)
                elif coordtype == 'enu':
                    ax.text(np.mean(e), np.mean(n), np.mean(u), str(i), color='red', fontsize=5)

        if ax_map is not None and map_polys:
            from matplotlib.collections import PolyCollection
            patches = PolyCollection(map_polys, facecolors=map_colors,
                                     edgecolors='k', linewidths=0.1, alpha=0.5)
            ax_map.add_collection(patches)

        # 3D slip arrows, drawn as explicit line segments.
        #
        # Axes3D.quiver rotates its arrowheads about a horizontal axis, so the
        # barbs leave the fault plane for a steeply dipping fault. Here the
        # shaft and both barbs are built from the two in-plane patch edges,
        # which keeps them inside the plane and at the correct length. The
        # shaft starts at the patch centre (tail pivot), not at the midpoint.
        if arrows:
            from mpl_toolkits.mplot3d.art3d import Line3DCollection

            segments = []
            # Head length as a fraction of the arrow and the half-angle
            # between each barb and the shaft (configurable, see docstring).
            head_frac = float(slip_head)
            head_ang  = np.deg2rad(float(slip_head_angle))
            for centre, vlen, v_hat, w_hat in arrows:
                if vlen <= 0:
                    continue

                tip     = centre + vlen * v_hat
                head    = head_frac * vlen
                back    = tip - head * np.cos(head_ang) * v_hat
                side    = head * np.sin(head_ang)

                segments.append([centre, tip])
                segments.append([tip, back + side * w_hat])
                segments.append([tip, back - side * w_hat])

            if segments:
                ax.add_collection3d(
                    Line3DCollection(segments, colors='k', linewidths=0.6))

        # Slip vectors on the map view: horizontal projection of the in-plane
        # slip direction, with lengths proportional to sqrt(ss^2 + ds^2).
        if ax_map is not None and show_slip and slip_mag is not None:
            lon_c  = felem[['lt_lon', 'rt_lon', 'rb_lon', 'lb_lon']].mean(axis=1).values
            lat_c  = felem[['lt_lat', 'rt_lat', 'rb_lat', 'lb_lat']].mean(axis=1).values
            coslat = np.maximum(np.cos(np.deg2rad(lat_c)), 1e-3)
            hx, hy = slip_dir_enu[:, 0], slip_dir_enu[:, 1]
            hnorm  = np.hypot(hx, hy)
            good   = (slip_mag > 0) & (hnorm > 0) & (arrow_scale > 0)
            if np.any(good):
                len_km = slip_mag * arrow_scale
                u = hx[good] / hnorm[good] * len_km[good] / (111.32*coslat[good])
                v = hy[good] / hnorm[good] * len_km[good] / 111.32
                ax_map.quiver(lon_c[good], lat_c[good], u, v, angles='xy',
                              scale_units='xy', scale=1, pivot='tail',
                              color='k', width=0.002,
                              transform=ccrs.PlateCarree(), zorder=5)

        # Plot aftershocks
        plot_aftershocks(ax, afsfile)
        if ax_map is not None:
            plot_aftershocks(ax_map, afsfile, transform=ccrs.PlateCarree())

        # Plot InSAR
        s_sar = plot_sar(ax, sarfile)
        s_sar_map = None
        if ax_map is not None:
            s_sar_map = plot_sar(ax_map, sarfile, transform=ccrs.PlateCarree())
        if s_sar_map is not None:
            fig.colorbar(s_sar_map, ax=ax_map, orientation='horizontal',
                         shrink=0.8, pad=0.12, label='LOS (mm)')
        elif s_sar is not None:
            fig.colorbar(s_sar, ax=ax, shrink=0.8, pad=0.1, label='LOS (mm)')

        # Plot GPS
        plot_gps(ax, gpsfile)
        if ax_map is not None:
            plot_gps(ax_map, gpsfile, scale=scale,
                     transform=ccrs.PlateCarree())

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

        # Shared slip colour bar for the map and 3D views
        tmax = felem['total_slip'].max() / 1000.
        sm   = plt.cm.ScalarMappable(norm=mcolors.Normalize(vmin=0, vmax=max(tmax, 1e-6)),
                                     cmap=plt.cm.jet)
        sm.set_array([])
        fig.colorbar(sm, ax=[ax_map, ax] if ax_map is not None else ax,
                     shrink=0.6, pad=0.12, label='Slip (m)')
        #   cb.set_label('Total slip (m)', fontsize=12)
        if savefig:
            fig.savefig(savefig, dpi=300, bbox_inches='tight')
        plt.show()

    def faultgeom2geojson(self, scale=1, filename="fault_patches.geojson"):
        """
        Export fault patches in DataFrame to a GeoJSON file.
        
        Parameters:
            felem : pd.DataFrame
                DataFrame containing at least the columns:
                ['lb_lon','lb_lat','rb_lon','rb_lat','rt_lon','rt_lat','lt_lon','lt_lat']
            filename : str
                Output GeoJSON filename
        """
        import json
        features = []
        felem = self.FaultGeom2AllVertex()
        for _, row in felem.iterrows():
            # Define polygon coordinates (GeoJSON expects [lon, lat])
            polygon = [
                [row['lb_lon'], row['lb_lat']],
                [row['rb_lon'], row['rb_lat']],
                [row['rt_lon'], row['rt_lat']],
                [row['lt_lon'], row['lt_lat']],
                [row['lb_lon'], row['lb_lat']],  # close the polygon
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

        
        
if __name__ == '__main__':
    pass
