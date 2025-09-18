#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Feb 14 2016
# Mod by Zhao Bin, Dec. 7, 2018. replace print() with print

import os, sys, logging
import numpy as np
from numpy import array, zeros
from pyginvc.libs import geotools as gt
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")


class Fault(object):
    '''
    Fault represents rectangular fault plane 
    '''
    
    nsegs         = 0
    ndeps         = 0
    dis_geom_grid = []
    geom_grid     = []
    origin        = []
    
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
        
        # load the original fault geometry
        geom = self.LoadFaultGeom(faultfile)
        if len(geom) == 0:
            logging.fatal('The input fault file {} is empty.'.format(faultfile))
            exit()
        
        # convert the fault into DISMODEL format
        if len(origin) != 2:
            origin = array([np.mean(geom[:,[3,5]]), np.mean(geom[:,[4,6]])])
        logging.info('origin is set to {:10.4f} {:10.4f}'.format(origin[0], origin[1]))
        [dis_geom, origin] = self.SetOrigin(geom, origin)
        
        # if we need subdivide the fault
        if doSubFault == True:
            pm = self.PatchFault(dis_geom, nsegs, ndeps)
            nf = nsegs*ndeps
            
            # init the output variables
            dis_geom_grid = np.column_stack((pm, zeros((nf,3))))            
            geom_grid     = self.disgeom2geom(dis_geom_grid, nf, origin)
    
            logging.info('Fault geometry is loaded and subdivided into %d * %d patches.' %(nsegs, ndeps))
        else:
            dis_geom_grid = dis_geom
            geom_grid     = geom
            logging.info('Fault geometry is loaded. There are %d * %d patches in it.' %(nsegs, ndeps))

        self.nf            = len(geom_grid)
        self.nsegs         = nsegs
        self.ndeps         = ndeps
        self.dis_geom_grid = dis_geom_grid
        self.geom_grid     = geom_grid
        self.origin        = origin
        
        flt_all = self.FaultGeom2AllVertex(self.geom_grid, self.dis_geom_grid, origin=self.origin)
        self.length        = flt_all[:,0]
        self.width         = flt_all[:,1]
        self.dip           = flt_all[:,2]
        self.strike        = flt_all[:,3]
        self.rake          = flt_all[:,8]
        self.ss            = flt_all[:,4]
        self.ds            = flt_all[:,5]
        self.op            = flt_all[:,6]
        self.ts            = flt_all[:,7]
        self.top_left_llh  = flt_all[:,9:12]
        self.top_right_llh = flt_all[:,12:15]
        self.bot_right_llh = flt_all[:,15:18]
        self.bot_left_llh  = flt_all[:,18:21]
        self.top_left_neu  = flt_all[:,21:24]
        self.top_right_neu = flt_all[:,24:27]
        self.bot_right_neu = flt_all[:,27:30]
        self.bot_left_neu  = flt_all[:,30:33]

        # Add by Zhao Bin. Aug 9, 2021
        self.center_llh      = np.zeros((self.nf,3))
        self.center_llh[:,0] = np.mean(flt_all[:,[9,12,15,18]], axis=1)
        self.center_llh[:,1] = np.mean(flt_all[:,[10,13,16,19]], axis=1)
        self.center_llh[:,2] = np.mean(flt_all[:,[11,14,17,20]], axis=1)
        self.center_neu      = np.zeros((self.nf,3))
        self.center_neu[:,0] = np.mean(flt_all[:,[21,24,27,30]], axis=1)
        self.center_neu[:,1] = np.mean(flt_all[:,[22,25,28,31]], axis=1)
        self.center_neu[:,2] = np.mean(flt_all[:,[23,26,29,32]], axis=1)
        self.felem           = flt_all
        return
    
    @staticmethod
    def LoadFaultGeom(filename):
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
        from numpy import array, genfromtxt
        # init parameter
        geom = array([])
    
        # check the file is exit
        isfile = os.path.isfile(filename)
        if ( isfile == False):
           logging.fatal('fault file does not exist! Please input another file!')
           sys.exit()
    
        # load the data
        geom = genfromtxt(filename, comments='#')
        if geom.ndim == 1: geom = geom.reshape((1, len(geom)))
        logging.info('Fault Geometry is loaded.')
        return geom

    @staticmethod
    def SetOrigin(geom, origin):
        '''
        Written by Zhao Bin when his is in UC Berkeley, Feb 20 2016
    
        Input:
            geom     = Fault Geometry 
                      Width  Depth  Dip  SLon SLat ELon Elat Strike-Slip Dip-slip Openning
            origin   = [lat, lon]
        Output:
            dis_geom = numpy array, containing Length Width Depth Dip Strike dE dN Strike-Slip Dip-slip Openning
            origin   = numpy array, containing latitude and longitude in degree
        '''
        from numpy import array, zeros, math, sqrt, rad2deg
        [nf, ncol] = geom.shape
            
        # initiate the ouput parameter dis_geom
        dis_geom   = zeros((nf,10))
       
        # for each subfaults
        for i in range(nf):
            
            # convert latitude and longitude to XY with respect to origin
            xy_1end  = gt.llh2utm([geom[i,3], geom[i,4]], origin)
            xy_2end  = gt.llh2utm([geom[i,5], geom[i,6]], origin)
            
            # Length of the fault
            leng     = sqrt(sum((xy_1end-xy_2end)**2))
            
            # strike angle of the fault
            strik = rad2deg(math.atan2((xy_2end[1]-xy_1end[1]), (xy_2end[0]-xy_1end[0])))
            strik = 90-strik
            
            # local coordinates with respect to origin
            delE = 0.5*(xy_1end[0] + xy_2end[0])
            delN = 0.5*(xy_1end[1] + xy_2end[1])
            
            dis_geom[i,:] = array([leng, geom[i,0], geom[i,1], geom[i,2], strik, delE, delN, geom[i,7], geom[i,8], geom[i,9]])    
        return dis_geom, origin
        
    @staticmethod
    def PatchFault(m, i, j):
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
        
        from numpy import sin, cos, arange, transpose, deg2rad, ones, dot, lexsort
        from numpy import zeros, mat, array
        dip     = deg2rad(m[0,3])
        strike  = deg2rad(-m[0,4])
        sin_dip = sin(dip)
        cos_dip = cos(dip)
        iw      = m[0,0]/i
        jw      = m[0,1]/j
        iis     = arange(1, i+1)
        jjs     = transpose(arange(1, j+1))
        n       = i*j
        
        # Calculate midpoints, depths of each patch    
        p         = mat(cos_dip*(jw*jjs - m[0,1])).T * ones((1, i))   
        q         = mat(ones((j,1)) * ((iw*iis) - 0.5 * (m[0,0] + iw))) 
        r         = mat((m[0,2] - jw * sin_dip * (j - jjs))).T * ones((1, i))
        p         = p.T.flatten().T
        q         = q.T.flatten().T
        r         = r.T.flatten().T
    
        mp        = array([p, q, r]).reshape(3, p.size).T
        
        # Adjust midpoints for strike
        R         = array([[cos(strike), -sin(strike), 0],
                           [sin(strike), cos(strike), 0],
                           [ 0, 0, 1]])
    
        mp        = dot(mp, R.T)
    
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
        idx = lexsort((pm[:,5], pm[:,2]))
        pm  = pm[idx]    
        return pm
        
    @staticmethod
    def disgeom2geom(dis_geom, nf, origin):
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
            nf          = number of faults
            origin      = numpy array, 1x10 vector, lat, lon in degree
        
        Output:
            geom        = numpy array, nf*7 vector defined as following
                            Width  Depth  Dip  SLon SLat ELon Elat
                            ( km    km    deg  deg  deg  deg  deg)
            
        '''
        from numpy import sin, cos, deg2rad, zeros
        geom = zeros((nf,7))
        
        if nf > 0:
            for i in range(0, nf):
                # convert middle of fault plane to upper corner (lat1, lon1)
                len_x     = dis_geom[i,0]*sin(deg2rad(dis_geom[i,4]))
                len_y     = dis_geom[i,0]*cos(deg2rad(dis_geom[i,4]))
                tmpx      = dis_geom[i,5] - 0.5*len_x
                tmpy      = dis_geom[i,6] - 0.5*len_y
                
                fltllh    = gt.xy2llh([tmpy, tmpx], origin)
                geom[i,3] = fltllh[0,0]
                geom[i,4] = fltllh[0,1]
                
                tmpx      = dis_geom[i,5] + 0.5*len_x
                tmpy      = dis_geom[i,6] + 0.5*len_y
                fltllh    = gt.xy2llh([tmpy, tmpx], origin)
                geom[i,5] = fltllh[0,0]
                geom[i,6] = fltllh[0,1]
                geom[i,0] = dis_geom[i,1]
                geom[i,1] = dis_geom[i,2]
                geom[i,2] = dis_geom[i,3]
        return geom
        
    @staticmethod
    def FaultGeom2AllVertex(geom, disgeom, origin=[]):
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
        from numpy import sin, cos, deg2rad, rad2deg, ones, pi, sqrt, math
        from numpy import array
    
        N            = 6378206.4
        e            = 1.0/298.257
    
        # get the shape
        [nf, ncol]   = geom.shape
        flt_elem_all = ones([nf, 33])

        if len(origin) != 2:
            origin = array([np.mean(geom[0,[3,5]]), np.mean(geom[0,[4,6]])])
    
        for i in range(nf):
        
            W     = sqrt(1-e**2*sin(deg2rad(geom[i,3])))
            Rp    = (N/W)*cos(deg2rad(geom[i,3]))
            Rm    = N*(1-e**2)/W**3
    
            dip   = disgeom[i,3]
            strk  = disgeom[i,4]
            dipr  = deg2rad(dip)
            strkr = deg2rad(strk)
            dipdir= pi - strkr
            ss    =-disgeom[i,7]
            ds    = disgeom[i,8]
            op    = disgeom[i,9]
            ts    = sqrt(ss**2 + ds**2)
            rake  = rad2deg(math.atan2(ds, ss))
            if ts == 0: rake = 0.0
    
            leng     = disgeom[i,0]
            wid      = disgeom[i,1]
            up       = -1*disgeom[i,2]
            bottom   = up
            top      = up + wid*sin(dipr)
            
            gll      = complex(geom[i,3], geom[i,4])
            glr      = complex(geom[i,5], geom[i,6])
    
            # factor of 1e3 converts km to meters
            horproj  = disgeom[i,1]*cos(dipr)*1.0e3
            offset   = rad2deg(horproj)*complex( sin(dipdir)/Rm, cos(dipdir)/Rp )
            
            # convert offset from radians to degrees
            gul      = gll + offset
            gur      = glr + offset
    
            # convert llh to neu
            en       = gt.llh2utm([gll.real, gll.imag], origin)
            gll_e    = en[0]
            gll_n    = en[1]
            gll_u    = bottom
            en       = gt.llh2utm([glr.real, glr.imag], origin)
            glr_e    = en[0]
            glr_n    = en[1]
            glr_u    = bottom
            en       = gt.llh2utm([gur.real, gur.imag], origin)
            gur_e    = en[0]
            gur_n    = en[1]
            gur_u    = top
            en       = gt.llh2utm([gul.real, gul.imag], origin)
            gul_e    = en[0]
            gul_n    = en[1]
            gul_u    = top
    
            
            # assign to flt_elem_all
            flt_elem_all[i,:] = array([leng, wid, dip, strk, ss, ds, op, ts, rake, gll.imag, gll.real, bottom, glr.imag, glr.real, bottom, gur.imag, gur.real, top, gul.imag, gul.real, top, gll_e, gll_n, gll_u, glr_e, glr_n, glr_u, gur_e, gur_n, gur_u, gul_e, gul_n, gul_u])
    
        return flt_elem_all

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
        from numpy import mat, array, zeros, sqrt, rad2deg, math
        # convert array to matrix
        geom = mat(geom)
        [nf, ncol] = geom.shape
        geom = array(geom)

        # initiate the ouput parameter x
        x = zeros((nf,7))
   
	    # for each subfaults
        for i in range(0,nf):
	        
	        # convert latitude and longitude to XY with respect to origin
	        xy_1end  = gt.llh2utm([geom[i,3], geom[i,4]], origin)
	        xy_2end  = gt.llh2utm([geom[i,5], geom[i,6]], origin)
	        
	        # Length of the fault
	        len  = sqrt( (xy_1end[0]-xy_2end[0])**2 +  (xy_1end[1]-xy_2end[1])**2);
	        
	        # strike angle of the fault
	        strik = rad2deg(math.atan2((xy_2end[1]-xy_1end[1]), (xy_2end[0]-xy_1end[0])))
	        strik = 90-strik
	        
	        # midpoint
	        midpointlat = 0.5*(geom[0,3] + geom[0,5])
	        midpointlon = 0.5*(geom[0,4] + geom[0,6])
	        
	        x[i,:] = array([len, geom[i,0], geom[i,1], geom[i,2], strik, midpointlat, midpointlon])    
	    
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

        from numpy import zeros, deg2rad, sin, cos, sqrt, rad2deg, array
        fmt = 7*"%10.4f"
        print(fmt %(x[0,0], x[0,1], x[0,2], x[0,3], x[0,4], x[0,5], x[0,6]))
        nf       = len(x)
        geom     = zeros((nf, 7))
        dis_geom = zeros((nf, 7))
    
        N        = 6378206.4
        e        = 1./298.257
    
        for i in range(nf):
    #        first    = 7*i
    #        tempgeom = x[first:first+7]
            tempgeom = x[i]
            theta    = deg2rad(90 - tempgeom[4])
            EOffset  = 0.5*tempgeom[0]*cos(theta)
            NOffset  = 0.5*tempgeom[0]*sin(theta)
    
            W        = sqrt(1-e**2*sin(deg2rad(tempgeom[5])))
            Rp       = (N/W)*cos(deg2rad(tempgeom[5]))
            Rm       = N*(1-e**2)/W**3
    
            NDeg     = rad2deg(1.0e3*NOffset/Rm)
            EDeg     = rad2deg(1.0e3*EOffset/Rp)
    
            geom[i,0]  = tempgeom[1]
            geom[i,1]  = tempgeom[2]
            geom[i,2]  = tempgeom[3]
            geom[i,3]  = tempgeom[5] - NDeg
            geom[i,4]  = tempgeom[6] - EDeg
            geom[i,5]  = tempgeom[5] + NDeg
            geom[i,6]  = tempgeom[6] + EDeg
    
            xy_1end = gt.llh2utm([geom[i,3], geom[i,4]], origin)
            xy_2end = gt.llh2utm([geom[i,5], geom[i,6]], origin)
    
            delE    = 0.5*(xy_1end[0] + xy_2end[0])
            delN    = 0.5*(xy_1end[1] + xy_2end[1])
    
            dis_geom[i,:] = array([tempgeom[0], geom[i,0], geom[i,1], geom[i,2], tempgeom[4], delE, delN])    
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
	    '''
	
	    # check the file is exit
	    isfile = os.path.isfile(filename)
	    if ( isfile == False):
	       logging.warning('{} is not exist! Please input another file!'.format(filename))
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
	    
	    return xlb, xub

    def FaultGeom2Coulomb(self):
        '''
        Convert FaultGeom file into the format Coulomb 
        By Zhao Bin, Mar 1st, 2016
        
        Notice that the in Coulomb x(east) and y(north)
        '''
        from numpy import deg2rad, sin
    
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
    
        edcmp_flt = np.zeros((self.nf, 10))
    
        # open output file
        with open('edcmp.fault', 'w') as fid:
            fid.write("# origin: {:10.4f} {:10.4f} ".format(self.origin[0], self.origin[1]))
            for i in range(self.nf):
                strik = self.strike[i] % 360
                dip   = self.dip[i] % 360
      
                # total slip, convert mm in FautlGeom to m in EDCMP
                slip  = self.ts/1e3
                
                leng  = self.length[i]*1e3
                width = self.width[i]*1e3
                north = self.top_left_neu[i,0]*1e3
                east  = self.top_left_neu[i,1]*1e3
                dep   =-self.top_left_neu[i,2]*1e3
            
    
                edcmp_flt[i,:] = np.array([i+1, slip[i], north, east, dep, leng, width, strik, dip, self.rake[i]])
                fid.write("%4d %10.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3f %10.3f %10.3f\n" 
                    %(i+1, slip[i], north, east, dep, leng, width, strik, dip, self.rake[i]))

        return edcmp_flt

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
    
        for i in range(self.nf):
            print("%10.3f %10.3f %10.3f %8.2f %8.2f %10.3f %10.3f %10.3f %8.2f" 
                %(self.bot_right_llh[i,1],
                  self.bot_right_llh[i,0],
                  self.width, self.strike, self.rake, self.ts/10.,
                  -self.top_left_llh[i,2],
                  -self.bot_left_llh[i,2],
                  self.dip))
                  
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
    
        idx = np.where(self.ts>0)[0]
        print("{}".format(len(idx)))
        for i in idx:
            if self.dip[i]<180:
                print("%10.3f %10.3f %8.2f"
                    %(-self.top_left_llh[i,2], 
                      -self.bot_left_llh[i,2], self.dip[i]))
                print("%10.3f %10.3f %10.3f %8.2f %8.2f %10.3f" 
                    %(self.bot_right_llh[i,1],
                      self.bot_right_llh[i,0],
                      self.length[i], self.strike[i], self.rake[i], self.ts[i]/10.))
            elif self.dip[i]>180:
                print("%10.3f %10.3f %8.2f"
                    %(-self.top_left_llh[i,2], 
                      -self.bot_left_llh[i,2], self.dip[i]-180))
                print("%10.3f %10.3f %10.3f %8.2f %8.2f %10.3f" 
                    %(self.bot_right_llh[i,1],
                      self.bot_right_llh[i,0],
                      self.length[i], self.strike[i], self.rake[i]-90, self.ts[i]/10.))
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
    
        # open output file
        with open('relax.fault', 'w') as fid:
            fid.write("# origin: {:10.4f} {:10.4f}\n".format(self.origin[0], self.origin[1]))
            # for each patches
            for i in range(self.nf):
                fid.write("%4d %10.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3f %10.3f %10.3f\n" 
                    %(i+1, self.ts[i]/1e3,
                      self.top_left_neu[i,0]*1e3,
                      self.top_left_neu[i,1]*1e3,
                     -self.top_left_neu[i,2]*1e3,
                      self.length[i]*1e3,
                      self.width[i]*1e3,
                      self.strike[i],
                      self.dip[i]-180,
                      self.rake[i]))
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
        
    
        n          = 1
        np_st      = 1
        np_di      = 1
        start_time = 0
        
        for i in range(self.nf):
            strik = self.strike[i] % 360
            dip   = self.dip[i] % 180
            print("%5d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %5d %5d %5d"
                %(n, self.top_left_llh[i,1],
                     self.top_left_llh[i,0],
                    -self.top_left_llh[i,2],
                     self.length[i],
                     self.width[i],
                     stirk,
                     dip,
                     np_st, np_di, start_time))
            print("     %10.3f %10.3f %10.3f %10.3f %10.3f" 
                %(self.length[i]/2.0, self.width[i]/2,
                  self.ss[i]/1e3, -self.ds[i]/1e3, self.op[i]/1e3))
            n = n+1
        return
            
    def FaultGeom2GMT(self):
        '''
        Convert fault slip model in standard FaultGeom to GMT plotting format
        Written by Zhao Bin, @ UC Berkeley April 22 2016
        '''        
    
        with open('ss_slip_llh.gmtlin', 'w') as fid:
            for i in range(self.nf):
                fid.write("> -Z%f\n" %(self.ss[i]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.top_left_llh[i,0], 
                                                     self.top_left_llh[i,1], 
                                                     self.top_left_llh[i,2]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.top_right_llh[i,0], 
                                                     self.top_right_llh[i,1], 
                                                     self.top_right_llh[i,2]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.bot_right_llh[i,0], 
                                                     self.bot_right_llh[i,1], 
                                                     self.bot_right_llh[i,2]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.bot_left_llh[i,0],
                                                     self.bot_left_llh[i,1],
                                                     self.bot_left_llh[i,2]))
        
        with open('ds_slip_llh.gmtlin', 'w') as fid:
            for i in range(self.nf):
                fid.write("> -Z%f\n" %(self.ds[i]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.top_left_llh[i,0], 
                                                     self.top_left_llh[i,1], 
                                                     self.top_left_llh[i,2]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.top_right_llh[i,0], 
                                                     self.top_right_llh[i,1], 
                                                     self.top_right_llh[i,2]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.bot_right_llh[i,0], 
                                                     self.bot_right_llh[i,1], 
                                                     self.bot_right_llh[i,2]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.bot_left_llh[i,0],
                                                     self.bot_left_llh[i,1],
                                                     self.bot_left_llh[i,2]))
                
        with open('ts_slip_llh.gmtlin', 'w') as fid:
            for i in range(self.nf):
                fid.write("> -Z%f\n" %(self.ts[i]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.top_left_llh[i,0], 
                                                     self.top_left_llh[i,1], 
                                                     self.top_left_llh[i,2]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.top_right_llh[i,0], 
                                                     self.top_right_llh[i,1], 
                                                     self.top_right_llh[i,2]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.bot_right_llh[i,0], 
                                                     self.bot_right_llh[i,1], 
                                                     self.bot_right_llh[i,2]))
                fid.write("%10.5f %10.5f %10.4f\n" %(self.bot_left_llh[i,0],
                                                     self.bot_left_llh[i,1],
                                                     self.bot_left_llh[i,2])) 
        return

    def Moment(self, shearmodulus=3e10):
        '''
        '''

        nf = len(self.dis_geom_grid)
        Mo = np.zeros(nf)
        
        for i in range(0, nf):
            Mo[i] = 1000*self.dis_geom_grid[i,0] * \
                    1000*self.dis_geom_grid[i,1] * \
                    self.ts[i] * shearmodulus
        
        Mo_total = np.sum(Mo)/1000
        Mw_total = 2.0/3.0*np.log10(Mo_total) - 6.067

        return Mo_total, Mw_total

    def FaultGeom2VTK(self, scale=1.0):
        '''
        Write the fault geometry and slip distribution to VTK file for Paraview.
        '''

        nf = len(self.dis_geom_grid)
        with open('faultgeom.vtp', 'w') as fid:
            fid.writelines(["<?xml version=\"1.0\"?>\n",
                       "<VTKFile type=\"PolyData\" version=\"0.1\">\n",
                       "  <PolyData>\n"])

            for i in range(0,nf):
                fid.write("    <Piece NumberOfPoints=\"4\" NumberOfPolys=\"1\">\n")
                fid.write("      <Points>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"Fault Patch\" NumberOfComponents=\"3\" format=\"ascii\">\n")
                fid.write("         {} {} {}\n".format(self.top_left_neu[i,0]*scale, self.top_left_neu[i,1]*scale, self.top_left_neu[i,2]*scale))
                fid.write("         {} {} {}\n".format(self.top_right_neu[i,0]*scale, self.top_right_neu[i,1]*scale, self.top_right_neu[i,2]*scale))
                fid.write("         {} {} {}\n".format(self.bot_right_neu[i,0]*scale, self.bot_right_neu[i,1]*scale, self.bot_right_neu[i,2]*scale))
                fid.write("         {} {} {}\n".format(self.bot_left_neu[i,0]*scale, self.bot_left_neu[i,1]*scale, self.bot_left_neu[i,2]*scale))
                fid.write("         </DataArray>\n")
                fid.write("      </Points>\n")
                fid.write("      <Polys>\n")
                fid.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"{}\">\n".format(nf-1))
                fid.write("0 1 2 3\n")
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"4\" RangeMax=\"4\">\n")
                fid.write("          4\n")
                fid.write("        </DataArray>\n")
                fid.write("      </Polys>\n")
                fid.write("      <CellData Scalar=\"geometry\">\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"strike\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(self.strike[i]))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"dip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(self.dip[i]))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"strike slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(self.ss[i]))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"dip slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(self.ds[i]))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(self.ts[i]))
                fid.write("        </DataArray>\n")
                fid.write("      </CellData>\n")
                fid.write("    </Piece>\n")
            fid.write("  </PolyData>\n")
            fid.write("</VTKFile>\n")
        
if __name__ == '__main__':
    print('Fault.py contain all functions related with FaultGeom')
