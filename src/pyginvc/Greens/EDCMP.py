#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Mar 7 2016
# Revised by Zhao Bin, Apr. 26, 2019. to Object Oriented Code.
import os
import shlex
import numpy as np
import subprocess, logging
from scipy import linalg
from pyginvc.libs import geotools as gt
from pyginvc.Greens.BaseGreen import BaseGreen
logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")


class EDCMP(BaseGreen):
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
        super(EDCMP, self).__init__(flt, data, dict_green)

    def build_greens(self):
        if not self.load_greens():
            self.generate_greens()
            self.save_greens()
        self.rotate_greens()
        logging.info("Building Green's functions finished.")


    def _grn_paths(self):
        '''返回 EDGRN 层状模型 Green 函数目录及三个分量文件，并检查其存在性。'''
        g      = self.dict_green
        grndir = g.get('grndir') or './'
        if grndir and not grndir.endswith(('/', '\\')):
            grndir += '/'
        files = [g.get('grn_ssfile', 'edgrnhs.ss'),
                 g.get('grn_dsfile', 'edgrnhs.ds'),
                 g.get('grn_clfile', 'edgrnhs.cl')]
        for f in files:
            path = os.path.join(grndir, f)
            if not os.path.exists(path):
                raise FileNotFoundError(
                    "EDGRN Green's function file %s does not exist. "
                    "Generate it with edgrn first, or set green.grndir / "
                    "grn_ssfile / grn_dsfile / grn_clfile in the YAML config."
                    % path)
        return grndir, files


    def _run_edcmp(self):
        '''运行 edcmp，输入文件名通过标准输入传入。'''
        exe = self.dict_green.get('edcmp_exe', 'edcmp')
        cmd = shlex.split(str(exe))
        if not cmd:
            raise ValueError('edcmp_exe is empty; set green.edcmp_exe.')
        # Direct exec (no shell): paths with spaces survive and no extra
        # /bin/sh process is spawned for every patch.
        with open('tmp.inp') as fin:
            subprocess.run(cmd, stdin=fin)


    @staticmethod
    def _read_disp(gdim):
        '''读取 edcmp 位移，并按反演要求的 [east, north, up] 顺序返回。

        EDCMP 输出的三列为 (Ux, Uy, Uz)，其中 x=north、y=east、
        z=downward（见 edcmain.f）。反演内部使用 [east, north, up]，
        因此交换前两列并对 Uz 取负。
        '''
        if gdim == 3:
            dat = np.genfromtxt('tmp.disp', usecols=[2, 3, 4], comments='#')
            dat = np.atleast_2d(dat)
            dat[:, 2] = -dat[:, 2]
            dat = dat[:, [1, 0, 2]]
        else:
            dat = np.genfromtxt('tmp.disp', usecols=[2, 3], comments='#')
            dat = np.atleast_2d(dat)
            dat = dat[:, [1, 0]]
        return dat.flatten()


    def generate_greens(self):
        '''
        Generate Green's Function.

        Output:
            self.G     = Green's function for GPS
            self.G_sar = Green's function for InSAR
        '''
        flt, data, dict_green = self.flt, self.data, self.dict_green
        llh_gps       = data.llh_gps
        llh_lev       = data.llh_lev
        llh_sar       = data.llh_sar
        ndim          = data.ndim
        unit          = data.unit
        origin        = flt.origin
        greentype     = dict_green['greentype']
        ss, ds, op    = greentype

        if op:
            raise NotImplementedError(
                "EDCMP (wang) does not support opening/tensile slip; "
                "set greentype[2] = 0.")
        if len(llh_lev) > 0:
            raise NotImplementedError(
                "Leveling data are not implemented for the EDCMP (wang) "
                "Green's function.")

        grndir, grnfiles = self._grn_paths()
        edgrn_ssfile, edgrn_dsfile, edgrn_clfile = grnfiles
        edcmp_flt     = flt.FaultGeom2EDCMP()
        if edcmp_flt is None or np.asarray(edcmp_flt).size == 0:
            raise ValueError(
                "Fault %s cannot be converted to EDCMP source geometry; "
                "FaultGeom2EDCMP() must return the 10-column source array."
                % (flt.faultfile if hasattr(flt, 'faultfile') else ''))
        edcmp_flt = np.atleast_2d(edcmp_flt)

        # convert the LLH to local coordinate
        xy_gps        = self._convert_to_local(llh_gps, origin)
        xy_sar        = self._convert_to_local(llh_sar, origin)

        n_gps = len(llh_gps) * ndim
        n_sar = len(llh_sar)
        nf3   = 3 * flt.nf
        self.G        = np.zeros((n_gps, nf3))
        self.G_sar    = np.zeros((n_sar, nf3))

        # if we have GPS data
        if xy_gps.size > 0 and n_gps > 0 and edcmp_flt is not None and edcmp_flt.size > 0:
            self.G = self.MakeGGPS(edcmp_flt, xy_gps, greentype, ndim,
                                   grndir, edgrn_ssfile, edgrn_dsfile, edgrn_clfile)
            logging.info('Green function for %d GPS stations are computed.' % (len(xy_gps)))

        # if we have SAR data
        if xy_sar.size > 0 and n_sar > 0 and edcmp_flt is not None and edcmp_flt.size > 0:
            self.G_sar = self.MakeGSAR(unit, edcmp_flt, xy_sar, greentype,
                                       grndir, edgrn_ssfile, edgrn_dsfile, edgrn_clfile)
            logging.info('Green function for %d InSAR points are computed.' % (len(xy_sar)))

        if dict_green.get('gps_ramp', False) and xy_gps.size > 0:
            self.G_gps_ramp = self.make_green_gps_ramp(xy_gps, ndim, method=1)
        else:
            self.G_gps_ramp = np.empty((self.G.shape[0], 0))
        if dict_green.get('sar_ramp', False) and xy_sar.size > 0:
            sizes = data.n_sar
            ramp  = []
            for i in range(len(sizes)):
                start = sum(sizes[:i])
                end   = sum(sizes[:i+1])
                ramp.append(self.make_green_sar_ramp(xy_sar[start:end]))
            self.G_sar_ramp = linalg.block_diag(*ramp)
        else:
            self.G_sar_ramp = np.empty((self.G_sar.shape[0], 0))

        logging.info("Green's function are computed using EDCMP (Wang's layered model).")


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
                gen_edcmp_inp(edcmp_flt[i], xy[:,1], xy[:,0], np.zeros(len(xy)),
                              grndir, grn_ssfile, grn_dsfile, grn_clfile)
                self._run_edcmp()
                G[:,3*i+0] = self._read_disp(gdim)

            # if we want to constrain dip slip
            if bool(ds) is True:
                edcmp_flt[i,1] = 1.0
                edcmp_flt[i,9] = 90
                gen_edcmp_inp(edcmp_flt[i], xy[:,1], xy[:,0], np.zeros(len(xy)),
                              grndir, grn_ssfile, grn_dsfile, grn_clfile)
                self._run_edcmp()
                G[:,3*i+1] = self._read_disp(gdim)

            # opening/tensile slip is not supported by EDCMP
            if bool(op) is True:
                raise NotImplementedError(
                    "EDCMP (wang) does not support opening/tensile slip; "
                    "set greentype[2] = 0.")

        return G
        
        
    def MakeGSAR(self, unit, edcmp_flt, xy, greentype, edgrndir, edgrn_ssfile, edgrn_dsfile, edgrn_clfile):
        '''
        Make Green's Function for InSAR data

        Input:
            unit         -  unit vector(s) of InSAR data, [east, north, up]
            edcmp_flt    -  numpy array, EDCMP source geometry
            xy           -  numpy array, xy coordinates of stations [east, north]
            greentype    -  [ss, ds, op] switch

        Output:
            G            -  numpy array, design matrix, shape (ndata, 3*nf)
        '''
        #
        ss = greentype[0]
        ds = greentype[1]
        op = greentype[2]
        if op:
            raise NotImplementedError(
                "EDCMP (wang) does not support opening/tensile slip; "
                "set greentype[2] = 0.")
        if edcmp_flt.ndim == 1:
            edcmp_flt = edcmp_flt.reshape(-1, 10)
        if xy.ndim == 1:
            xy = xy.reshape(-1, 2)

        nf     = edcmp_flt.shape[0]
        ndata  = xy.shape[0]

        # Descending track unit vector
        if len(unit) == 0:
            lookangle = 23
            track     = 13.9
            rad_lookangle = np.deg2rad(lookangle)
            rad_track     = np.deg2rad(track)
            unit  = np.array([-np.cos(rad_track)*np.sin(rad_lookangle),
                               np.sin(rad_track)*np.sin(rad_lookangle),
                              -np.cos(rad_lookangle)])
            unit  = -1.0*unit
        unit = np.asarray(unit)
        G = np.zeros((ndata, 3*nf))

        # for each subfaults
        for i in range(nf):
            for k, (flag, rake) in enumerate(((ss, 0.0), (ds, 90.0))):
                if not bool(flag):
                    continue
                edcmp_flt[i,1] = 1.0
                edcmp_flt[i,9] = rake
                gen_edcmp_inp(edcmp_flt[i], xy[:,1], xy[:,0], np.zeros(len(xy)),
                              edgrndir, edgrn_ssfile, edgrn_dsfile, edgrn_clfile)
                self._run_edcmp()
                u = self._read_disp(3).reshape(ndata, 3)
                G[:,3*i+k] = (u * unit).sum(axis=1)

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
    if len(n):
        np.savetxt(fid, np.column_stack((np.asarray(n, dtype=float) * 1000.0,
                                         np.asarray(e, dtype=float) * 1000.0)),
                   fmt="(%e, %e)")

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
    np.savetxt(fid, sfault,
               fmt="%d %f %e %e %e %e %e %f %f %f")

    # 4. earth model
    fid.write('# earth model\n')
    fid.write('1\n')
    line = "'{}' '{}' '{}' '{}'".format(edgrndir, ss, ds, cl)
    fid.write(line)

    # close the file
    fid.close()
