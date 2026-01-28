# Written by Zhao Bin, when he is at UC Berkeley. April 18 2016
# Mod by Zhao Bin, Dec. 7, 2018. We use print() now.
import numpy as np
import logging, h5py
import os, time, glob
from pandas import pandas as pd
from numpy import column_stack, hstack, vstack, transpose
#from pyginvc.Geometry.Fault import Fault
from pyginvc.Geometry.Patch import Fault
from pyginvc.Geometry.Triangle import Triangle
from pyginvc.Forward.TriForward import TriForward

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%m-%Y %H:%M:%S")

class Output(object):
    '''
    Output is a class representing output inversions, summary, figures
    '''

    def __init__(self, flt, data, green, sol=None, archdir=None):
        '''
        Constructor.

        Parameters
        ----------
        flt : instance
            an instance of Triangle or Fault.
        data : instance
            an instance of GeoData
        green : instance
            an instance of Okada, Nikkhoo, TriPoly3D
        sol : instance
            an instance of GeoInversion
        archdir : string, optional
            directory

        Returns
        -------
        None.

        '''

        self.flt       = flt
        self.data      = data
        self.green     = green
        self.sol       = sol
        self.archdir   = archdir
        
        self.flttype = 'rectangle'
        if isinstance(flt, Triangle):
            self.flttype = 'triangle'

        return

    def OutputSolution(self):
        '''
        Output solutions to hdf5 file and ascii files

        Returns
        -------
        None.

        '''

        # output hdf5 file
        with h5py.File('solutions.h5', 'w') as h5:
            # observation
            grp = h5.create_group('obs')
            grp.create_dataset('llh_gps', data = self.data.llh_gps, compression='gzip')
            grp.create_dataset('llh_lev', data = self.data.llh_lev, compression='gzip')
            grp.create_dataset('llh_sar', data = self.data.llh_sar, compression='gzip')
            grp.create_dataset('unit',  data = self.data.unit,  compression='gzip')
            grp.create_dataset('d_gps', data = self.data.d_gps, compression='gzip')
            grp.create_dataset('d_lev', data = self.data.d_lev, compression='gzip')
            grp.create_dataset('d_sar', data = self.data.d_sar, compression='gzip')
            grp['ndim']     = self.data.ndim

#           grp['gps_site'] = self.data.station_gps
            
            # fault geometry
            grp = h5.create_group('flt')
            grp.create_dataset('origin', data = self.flt.origin, compression='gzip')
            if self.flttype == 'rectangle':
                grp.create_dataset('flttype', data = 'rectangle')
                grp.create_dataset('dis_geom', data = self.flt.dis_geom_grid, compression='gzip')
                grp.create_dataset('geom', data = self.flt.geom_grid, compression='gzip')
            if self.flttype == 'triangle':
                grp.create_dataset('flttype', data = 'triangle')
                grp.create_dataset('vertex_enu', data = self.flt.vertex_enu, compression='gzip')
                grp.create_dataset('vertex_llh', data = self.flt.vertex_llh, compression='gzip')
                grp.create_dataset('element', data = self.flt.element, compression='gzip')

            # solution
            grp = h5.create_group('sol')
            grp.create_dataset('slip', data = self.sol.slip, compression='gzip')
#           grp.create_dataset('sig_slip', data = self.sol.sig_slip, compression='gzip')
            grp.create_dataset('r', data = self.sol.r, compression='gzip')
            grp.create_dataset('misfit', data = self.sol.misfit, compression='gzip')
            grp.create_dataset('dhat', data = self.sol.dhat, compression='gzip')

        # output ascii files
        if isinstance(self.sol, TriForward):
            self.WriteData()
            self.archive_outfile()
        else:
            self.WriteSummary()
            self.WriteSlipModel()
            self.WriteGMTSlip()
            self.WriteLCurve()
            self.WriteData()
            self.archive_outfile()
        return

    def WriteGMTSlip(self):
        '''
        Ouput slip model in GMT format

        Returns
        -------
        None.

        '''
        nf       = self.flt.nf
        slip     = self.sol.slip[-1,0:nf*3]
        sig_slip = self.sol.sig_slip[-1,0:nf*3]
        slip     = slip.reshape(nf, 3)
        sig_slip = sig_slip.reshape(nf, 3)
        slip_mag = np.sqrt(slip[:,0]**2 + slip[:,1]**2)
        
        if self.flttype == 'rectangle':
            self.write_rectangle_slip('ss_slip_llh.gmtlin', slip[:,0])
            self.write_rectangle_slip('ds_slip_llh.gmtlin', slip[:,1])
            self.write_rectangle_slip('ts_slip_llh.gmtlin', slip_mag)
            
        if self.flttype == 'triangle':
            self.write_triangle_slip('ss_slip_llh.gmtlin', slip[:,0])
            self.write_triangle_slip('ds_slip_llh.gmtlin', slip[:,1])
            self.write_triangle_slip('ts_slip_llh.gmtlin', slip_mag)
            
        # print the status
        logging.info('Fault slip distribution is ouput in geography system.')
    
    def write_rectangle_slip(self, filename, value):
        """
        Write slip distribution in GMT format
        """
        nf       = self.flt.nf
        all_data = []
        
        felem    = self.flt.FaultGeom2AllVertex()
        for i in range(nf):
            all_data.append(f'>-Z{value[i]:.3f}')
            #tl = self.flt.top_left_llh[i]
            #tr = self.flt.top_right_llh[i]
            #br = self.flt.bot_right_llh[i]
            #bl = self.flt.bot_left_llh[i]
            tl = [felem['lt_lon'][i], felem['lt_lat'][i], felem['lt_dep'][i]]
            tr = [felem['rt_lon'][i], felem['rt_lat'][i], felem['rt_dep'][i]]
            br = [felem['rb_lon'][i], felem['rb_lat'][i], felem['rb_dep'][i]]
            bl = [felem['lb_lon'][i], felem['lb_lat'][i], felem['lb_dep'][i]]
            
            all_data.append(f'{tl[0]:10.4f} {tl[1]:10.4f} {tl[2]:10.4f}')
            all_data.append(f'{tr[0]:10.4f} {tr[1]:10.4f} {tr[2]:10.4f}')
            all_data.append(f'{br[0]:10.4f} {br[1]:10.4f} {br[2]:10.4f}')
            all_data.append(f'{bl[0]:10.4f} {bl[1]:10.4f} {bl[2]:10.4f}')
        
        header = '# Rectangle fault slip model'
        np.savetxt(filename, all_data, fmt='%s', header=header)
    
    def write_triangle_slip(self, filename, value):
        """
        Write slip distribution in GMT format
        """
        nf = self.flt.nf
        v  = self.flt.vertex_llh[:,[1,0,2]]
        e  = self.flt.element-1
        all_data = []
        
        for i in range(nf):
            all_data.append(f'>-Z{value[i]:.3f}')
            all_data.append(f'{v[e[i,0],0]:10.3f} {v[e[i,0],1]:10.3f} {v[e[i,0],2]:10.3f}')
            all_data.append(f'{v[e[i,1],0]:10.3f} {v[e[i,1],1]:10.3f} {v[e[i,1],2]:10.3f}')
            all_data.append(f'{v[e[i,2],0]:10.3f} {v[e[i,2],1]:10.3f} {v[e[i,2],2]:10.3f}')
        
        header = '# Rectangle fault slip model'
        np.savetxt(filename, all_data, fmt='%s', header=header)
 
    def WriteSlipModel(self):
        '''
        Output slip model in FaultGeom or Triangular format

        Returns
        -------
        None.

        '''
        
        slip = self.sol.slip
        if slip.ndim == 1:
            slip = np.expand_dims(slip, axis=0)

        nscount, nslip = slip.shape
        
        if self.flttype == 'rectangle':
            for i in range(nscount):
                islip = slip[i].reshape((self.flt.nf, 3))
                patchFaultGeom = np.hstack((self.flt.geom_grid[:, 0:7], islip))
                fname = 'slipmodel_{:03d}.faultgeom'.format(i+1)
                fmt   = 7*'%11.4f'+3*'%13.3f'
                header  = " width_km   depth_km    dip_deg   lat0_deg   lon0_deg   lat1_deg   lon1_deg  slp_strk_mm  slp_ddip_mm  slp_tens_mm"
                np.savetxt(fname, patchFaultGeom, fmt=fmt, header=header)

                # Rotate slip back to normal fault coordinate system
                if self.green.rake_beta != 0.0:
                    beta = np.deg2rad(self.green.rake_beta)
                    R = np.array([[np.cos(beta), -np.sin(beta)],
                                  [np.sin(beta), np.cos(beta)]])
                    islip    = slip[i].reshape((self.flt.nf, 3))
                    new_slip = R.dot(islip[:,[0,1]].T).T
                    patchFaultGeom = np.column_stack((self.flt.geom_grid[:,0:7], new_slip, islip[:,2]))
                    fname = 'slipmodel_rot_{:03d}.faultgeom'.format(i+1)
                    np.savetxt(fname, patchFaultGeom, fmt=fmt, header=header)
        elif self.flttype == 'triangle':
            for i in range(nscount):
                islip = slip[i].reshape((self.flt.nf, 3))
                patchFaultGeom = np.hstack((self.flt.element, islip))
                fname = 'slipmodel_%03d.tri' % (i + 1)
                fmt   = '%5d %5d %5d %15.4f %15.4f %15.4f'
                np.savetxt(fname, patchFaultGeom, fmt=fmt)

                # Rotate slip back to normal fault coordinate system
                if self.green.rake_beta != 0.0:
                    beta = self.green.rake_beta
                    R = np.array([[np.cos(beta), -np.sin(beta)],
                                  [np.sin(beta), np.cos(beta)]])
                    islip = slip[i].reshape((self.flt.nf, 3))
                    new_slip = R.dot(islip[:,[0,1]].T).T
                    patchFaultGeom = np.column_stack((self.flt.element, new_slip, islip[:,2]))
                    fname = 'slipmodel_rot_%03d.tri' % (i + 1)
                    np.savetxt(fname, patchFaultGeom, fmt=fmt)
        return
 
    def WriteSummary(self):
        '''
        Ouput summary of Inversion.

        Returns
        -------
        None.

        '''
        with open('summary', 'a') as (fid):
            fid.write('#*#*#*#*#*#*#*#*#*#Geodetic Data Summary#*#*#*#*#*#*#*#*#*#*\n')
            fid.write('Dimension of GPS                 = %d\n' % self.data.ndim)
            fid.write('Number of GPS stations           = %d\n' % len(self.data.llh_gps))
            fid.write('Number of InSAR piexl            = %d\n' % len(self.data.llh_sar))
            fid.write('Number of Leveling data          = %d\n' % len(self.data.llh_lev))
            fid.write('\n\n')
            if self.flttype == 'rectangle':
                fid.write('#*#*#*#*#*#*#*#*#*#*##*Fault Summary#*#*#*#*#*#*#*#*#*#*#*#*\n')
                fid.write('Number of patches along strike   = %d\n' % self.flt.nsegs)
                fid.write('Number of patches along dip      = %d\n' % self.flt.ndeps)
                fid.write('\n\n')
            elif self.flttype == 'triangle':
                fid.write('#*#*#*#*#*#*#*#*#*#*##*Fault Summary#*#*#*#*#*#*#*#*#*#*#*#*\n')
                fid.write('Number of patches along strike   = %d\n' % len(self.flt.vertex_llh))
                fid.write('Number of patches along dip      = %d\n' % len(self.flt.element))
                fid.write('\n\n')
            fid.write('#*#*#*#*#*#*#*#*#*#*#Inversion Summary#*#*#*#*#*#*#*#*#*#*#*\n')
            for i in range(len(self.sol.smo_facts)):
                fid.write('Smoothing factor                 = %f\n' % self.sol.smo_facts[i])
                if len(self.sol.data.llh_gps) > 0:
                    fid.write('GPS WRSS/N                       = %f\n' 
                            %(self.sol.misfit_gps[i] / len(self.sol.data.d_gps)))
                if len(self.sol.data.llh_sar) > 0:
                    fid.write('InSAR WRSS/N                     = %f\n' 
                            %(self.sol.misfit_sar[i] / len(self.sol.data.d_sar)))
                fid.write('Mo                               = %e\n' % self.sol.moment[(i, 0)])
                fid.write('Mw                               = %4.2f\n' % self.sol.moment[(i, 1)])
        return
 
    def WriteLCurve(self):
        """
        Output roughness misfit and smooth factors, which can be used to plot tradeoff between roughness and fittness
        Written by Zhao Bin @ UC Berkeley April 18 2016
         
        Input:
            ruff      = array, shape(m,)   roughness
            misfit    = array, shape(m,)   misfit
            smo_facts = array, shape(m,)   smooth factor
        Output:
            ruff_misfit: file
        """
        
        if len(self.sol.smo_facts) < 3:
            logging.warning('Only %d smoothing factors.' % len(self.sol.smo_facts))
            return
        ruff_misfit = column_stack((self.sol.ruff, self.sol.misfit, self.sol.smo_facts, self.sol.moment))
        fname       = 'ruff_misfit'
        header      = " Ruffness\t Misfit\t Smoothing_factor\t Moment"
        np.savetxt(fname, ruff_misfit, header=header, fmt='%10.2f\t%10.2f\t%10.4f\t%10.2e\t%10.2f')
        logging.info('Roughness vs misfit is printed into ruff_misfit.')
        return
	
        
    def WriteSurfSlip(self):
        '''
        '''        
        
        geom_grid   = self.flt.geom_grid
        nsegs       = self.flt.nsegs
        patch       = 0.5*(geom_grid[0:nsegs,4] + geom_grid[0:nsegs,5])
        ss_slip     = self.sol.slip[0::3]
        ds_slip     = self.sol.slip[1::3]
        ss_sig_slip = self.sol.sig_slip[0::3]
        ds_sig_slip = self.sol.sig_slip[1::3]
	    
        patch = 0.5*(geom_grid[0:nsegs,4] + geom_grid[0:nsegs,5])
        ss_surface_slip = transpose(vstack((patch, ss_slip[0:nsegs], ss_sig_slip[0:nsegs])))
        ds_surface_slip = transpose(vstack((patch, ds_slip[0:nsegs], ds_sig_slip[0:nsegs])))
        
        fname = 'surface_slip_ss.txt'
        np.savetxt(fname, ss_surface_slip[0:3], fmt="%10.2f%10.2f%10.2f")
        fname = 'surface_slip_ds.txt'
        np.savetxt(fname, ds_surface_slip[0:3], fmt="%10.2f%10.2f%10.2f")
        return
        
        
    def WriteData(self, dhat=None, r=None):
        '''
        Output modeled and residual GPS/LEV/SAR data
        Written by Zhao Bin, Institute of Seismology @ UC Berkeley. April 2016
        Mod by Zhao Bin, Apr. 25, 2019. If the gfiletype is IOS3D
	
        Input:
            dhat       -   numpy array, shape(m,), modeled data
            r          -   numpy array, shape(m,), residual data
        Output:
	        gps_obs.gmtvec   -  observed GPS data
	        gps_mod.gmtvec   -  modeled GPS data
	        gps_res.gmtvec   -  residual GPS data
	        gps_obs_up.gmtvec -  vertical observed GPS data
	        gps_mod_up.gmtvec -  vertical modeled GPS data
	        gps_res_up.gmtvec -  vertical ressidual GPS data
	        lev_mod.gmtvec   -  modeled LEV data
	        lev_res.gmtvec   -  residual LEV data
	        sar_mod.gmtvec   -  modeled SAR data
	        sar_res.gmtvec   -  residual SAR data
	    '''
	
        
#         if dhat is None or r is None:
#             dhat = self.sol.dhat
#             r    = self.sol.r
#         len_gps = len(self.data.d_gps)
#         len_lev = len(self.data.d_lev)
#         len_sar = len(self.data.d_sar)
#         len_geod= len_gps+len_lev
#         len_all = len_geod+len_sar
#         if len(self.data.llh_gps) > 0:
#             nsta = len(self.data.llh_gps)
#             if self.data.ndim == 2:
#                 mod  = dhat[0:nsta*2].reshape((nsta,2))
#                 obs  = self.data.d_gps.reshape((nsta,2))
#                 sig  = np.sqrt(np.diag(self.data.cov_gps)).reshape((nsta,2))
#                 llh  = self.data.llh_gps[:,[1,0]]
#                 if len(self.data.station_gps) == nsta:
#                     # observed
#                     fmt  = "%10.4f\t%10.4f\t%10.2f\t%10.2f\t%5.2f\t%5.2f\t%5.2f\t%8s\n"
#                     fmt3 = "%10.4f\t%10.4f\t%10.2f\t%10.2f\t%10.2f\t%5.2f\t%5.2f\t%5.2f"
#                     with open('gps_obs.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt %(llh[i,0], llh[i,1], obs[i,0], obs[i,1], sig[i,0],
#                                     sig[i,1], 0.0, self.data.station_gps[i]))

#                     # modeled
#                     with open('gps_mod.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt %(llh[i,0], llh[i,1], mod[i,0], mod[i,1], 0.0,
#                                     0.0, 0.0, self.data.station_gps[i]))

#                     # residual
#                     r_gps = obs - mod
#                     with open('gps_res.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt %(llh[i,0], llh[i,1], r_gps[i,0], r_gps[i,1], 0.0, 0.0, 0.0,
#                                 self.data.station_gps[i]))
#                 else:
#                     # observed
#                     fmt  = "%10.4f\t%10.4f\t%10.2f\t%10.2f\t%5.2f\t%5.2f\t%5.2f"
#                     null = np.zeros((nsta,1))
#                     outmatrix = hstack((llh, obs, sig, null))
#                     np.savetxt("gps_obs.gmtvec", outmatrix, fmt=fmt)

#                     # modeled
#                     null = np.zeros((nsta,3))
#                     outmatrix = hstack((llh, mod, null))
#                     np.savetxt("gps_mod.gmtvec", outmatrix, fmt=fmt)

#                     # residual
#                     outmatrix = hstack((llh, obs-mod, null))
#                     np.savetxt("gps_res.gmtvec", outmatrix, fmt=fmt)
#             if self.data.ndim == 3:
#                 mod  = dhat[0:nsta*3].reshape((nsta,3))
#                 obs  = self.data.d_gps.reshape((nsta,3))
#                 sig  = np.sqrt(np.diag(self.data.cov_gps)).reshape((nsta,3))
#                 llh  = self.data.llh_gps[:,[1,0]]
#                 if len(self.data.station_gps) == nsta:
#                     # observed
#                     fmt  = "%10.4f\t%10.4f\t%10.2f\t%10.2f\t%5.2f\t%5.2f\t%5.2f\t%8s\n"
#                     fmt3 = "%10.4f\t%10.4f\t%10.2f\t%10.2f\t%10.2f\t%5.2f\t%5.2f\t%8s\t%10.2f\t%5.2f\n"
#                     with open('gps_obs.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt %(llh[i,0], llh[i,1], obs[i,0], obs[i,1], sig[i,0],
#                                     sig[i,1], 0.0, self.data.station_gps[i]))

#                     with open('gps_obs_up.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt %(llh[i,0], llh[i,1], 0.0, obs[i,2], 0.0, sig[i,2], 0.0,
#                                 self.data.station_gps[i]))

#                     # modeled
#                     with open('gps_mod.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt %(llh[i,0], llh[i,1], mod[i,0], mod[i,1], 0.0,
#                                     0.0, 0.0, self.data.station_gps[i]))

#                     with open('gps_mod_up.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt %(llh[i,0], llh[i,1], 0.0, mod[i,2], 0.0, 0.0, 0.0,
#                                 self.data.station_gps[i]))

#                     # 3D 
#                     with open('gps_mod_3d.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt3 %(llh[i,0], llh[i,1], mod[i,0], mod[i,1], 0.0, 0.0, 0.0,
#                                 self.data.station_gps[i], mod[i,2], 0.0))

#                     # residual
#                     r_gps = obs - mod
#                     with open('gps_res.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt %(llh[i,0], llh[i,1], r_gps[i,0], r_gps[i,1], 0.0, 0.0, 0.0,
#                                 self.data.station_gps[i]))
#                     with open('gps_res_up.gmtvec', 'w') as fid:
#                         for i in range(nsta):
#                             fid.write(fmt %(llh[i,0], llh[i,1], 0.0, r_gps[i,2], 0.0, 0.0, 0.0,
#                                 self.data.station_gps[i]))
#                 else:
#                     # observed
#                     fmt = "%10.4f\t%10.4f\t%10.2f\t%10.2f\t%5.2f\t%5.2f\t%5.2f"
#                     fmt3 = "%10.4f\t%10.4f\t%10.2f\t%10.2f\t%10.2f\t%5.2f\t%5.2f\t%5.2f"
#                     outmatrix = np.column_stack((llh, obs[:,0:2], sig[:,0:2], np.zeros(nsta)))
#                     np.savetxt("gps_obs.gmtvec", outmatrix, fmt=fmt)
#                     outmatrix = np.column_stack((llh, np.zeros(nsta), obs[:,2], sig[:,2], np.zeros((nsta,2))))
#                     np.savetxt("gps_obs_up.gmtvec", outmatrix, fmt=fmt)

#                     # modeled
#                     null = np.zeros((nsta,3))
#                     outmatrix = np.column_stack((llh, mod[:,0:2], null))
#                     np.savetxt("gps_mod.gmtvec", outmatrix, fmt=fmt)
#                     outmatrix = np.column_stack((llh, null[:,0], mod[:,2], null))
#                     np.savetxt("gps_mod_up.gmtvec", outmatrix, fmt=fmt)
#                     outmatrix = np.column_stack((llh, mod, null))
#                     np.savetxt("gps_mod_3d.gmtvec", outmatrix, fmt=fmt3)

# 	                # residual
#                     r_gps = obs - mod
#                     outmatrix = np.column_stack((llh, r_gps[:,0:2], null))
#                     np.savetxt("gps_res.gmtvec", outmatrix, fmt=fmt)
#                     outmatrix = np.column_stack((llh, null[:,0], r_gps[:,2], null))
#                     np.savetxt("gps_res_up.gmtvec", outmatrix, fmt=fmt)

# 	        # print the status

        if dhat is None:
            dhat = self.sol.dhat
            r    = self.sol.r

        self.Write_GPS_Data(dhat, r)
        self.Write_LEV_Data(dhat, r)
        self.Write_SAR_Data(dhat, r)

    
    def Write_GPS_Data(self, dhat, r):
        """Write observed, modeled and residual GPS data"""
        len_gps = len(self.data.d_gps)
        nsta    = len(self.data.llh_gps)
        if nsta == 0:
            return
        
        llh     = self.data.llh_gps
        station = self.data.station_gps
        ndim    = self.data.ndim
        
        obs     = self.data.d_gps.reshape(nsta, ndim)
        mod     = dhat[0:len_gps].reshape(nsta, ndim)
        sig     = np.sqrt(np.diag(self.data.cov_gps)).reshape(nsta,ndim)
        
        output_cfg = []
        if ndim == 2:
            output_cfg = [
                ('obs', obs[:,0], obs[:,1], sig[:,0], sig[:,1]),
                ('mod', mod[:,0], mod[:,1], 0, 0),
                ('res', obs[:,0]-mod[:,0], obs[:,1]-mod[:,1], 0, 0)
                ]
        if ndim == 3:
            output_cfg = [
                ('obs', obs[:,0], obs[:,1], sig[:,0], sig[:,1]),
                ('mod', mod[:,0], mod[:,1], 0, 0),
                ('res', obs[:,0]-mod[:,0], obs[:,1]-mod[:,1], 0, 0),
                ('obs_up', 0, obs[:,2], 0, sig[:,2]),
                ('mod_up', 0, mod[:,2], 0, 0),
                ('res_up', 0, obs[:,2]-mod[:,2], 0, 0)
                ]
        for name, east, north, sig_e, sig_n in output_cfg:
            df = pd.DataFrame({
                'lon': llh[:,1],
                'lat': llh[:,0],
                'east': east,
                'north': north,
                'sig_e': sig_e,
                'sig_n': sig_n,
                'corr': 0,
                'station': station
                })
            
            fname = f'gps_{name}.gmtvec'
            df.to_csv(fname, sep='\t', index=False, float_format='%10.4f', header=False)
        for name, _, _, _, _ in output_cfg:
            fname = f'gps_{name}.gmtvec'
            with open(fname) as f:
                content = f.read()

            if name[-2:] != 'up':
                with open(fname, 'w') as f:
                    f.write("#      Lon             Lat        East(mm)       North(mm)  SigE(mm)  SigN(mm)  Corr    Site\n" + content)
            else:
                with open(fname, 'w') as f:
                    f.write("#      Lon             Lat        NULL(mm)          Up(mm)  NULL(mm)  SigU(mm)  Corr    Site\n" + content)
                
        logging.info('Modeled and residual GPS points are output.')

    def Write_LEV_Data(self, dhat, r):
        """Write modeled and residual leveling data"""
        len_gps = len(self.data.d_gps)
        len_lev = len(self.data.d_lev)
        len_geod= len_gps+len_lev
        
        if len(self.data.llh_lev) > 0:
            llh_lev   = self.data.llh_lev[:,[1,0]]
            d_lev     = self.data.d_lev
            dhat_lev  = dhat[len_gps:len_geod]
            
            outmatrix = hstack((llh_lev[:,1], llh_lev[:,0], dhat_lev)).T
            np.savetxt("lev_mod.xyz", outmatrix, fmt="%10.4f\t%10.4f\t%10.2f")
            res_lev = d_lev - dhat_lev
            outmatrix = hstack((llh_lev[:,1], llh_lev[:,0], res_lev)).T
            np.savetxt("lev_res.xyz", outmatrix, fmt="%10.4f\t%10.4f\t%10.2f")
            
            # print the status
            logging.info('%d modeled and residual LEV points are output.' %(len_lev))
    
    def Write_SAR_Data(self, dhat, r):
        """Write modeled and residual SAR data"""
        len_gps = len(self.data.d_gps)
        len_lev = len(self.data.d_lev)
        len_sar = len(self.data.d_sar)
        len_geod= len_gps+len_lev
        len_all = len_geod+len_sar
        
        if len(self.data.d_sar) < 1:
            return
        
        wsar       = self.sol.WSAR
        llh_sar    = self.data.llh_sar[:,[1,0]]
        d_sar      = self.data.d_sar
        n_sar      = self.data.n_sar
        r_sar      = r[len_geod:len_all]/np.diag(wsar)
        
        # For the case there are multiple sar files
        for i in range(len(n_sar)):
            idx0 = sum(n_sar[:i])
            idx1 = sum(n_sar[:i+1])
            
            header = "Lon  Lat  LOS(mm)"
            # residual
            lon, lat  = llh_sar[idx0:idx1,0], llh_sar[idx0:idx1,1]
            rsar      = r_sar[idx0:idx1]
            outmatrix = np.column_stack((lon, lat, rsar))
            np.savetxt(f"sar_res_{i+1}.xyz", outmatrix, header=header, fmt="%10.4f\t%10.4f\t%10.2f")
	
            # model
            msar       = (d_sar[idx0:idx1] - rsar)/np.diag(wsar)[idx0:idx1]
            outmatrix  = np.column_stack((lon, lat, msar))
            np.savetxt(f"sar_mod_{i+1}.xyz", outmatrix, header=header, fmt="%10.4f\t%10.4f\t%10.2f")

	        # print the status
            logging.info('%d modeled and residual SAR points are output.' %(len_sar))
            
    def archive_outfile(self):
        '''
        '''
        dirname = self.archdir
        if dirname == "":
            logging.error('This input archive directory is null')
            self.dirname = 'test'
        if os.path.exists(dirname) == False:
            os.mkdir(dirname)

        # GPS
        if len(self.data.llh_gps) > 0:
            if len(glob.glob('gps*.gmtvec'))>0:
                os.popen('mv gps*.gmtvec '+dirname)
        # InSAR
        if len(self.data.llh_sar) > 0:
            if len(glob.glob('sar*.xyz'))>0:
                os.popen('mv sar*.xyz '+dirname)
        # Roughness_misfit curve
        if os.path.isfile('ruff_misfit'):
            os.popen('mv ruff_misfit '+dirname)
        # Solution
        if len(glob.glob('slipmodel*.tri'))>0:
            os.popen('mv slipmodel*.tri '+dirname)
        if len(glob.glob('slipmodel*.faultgeom'))>0:
            os.popen('mv slipmodel*.faultgeom '+dirname)
        if len(glob.glob('*h5'))>0:
            os.popen('mv *.h5 '+dirname)
        if len(glob.glob('*gmtlin'))>0:
            os.popen('mv *.gmtlin '+dirname)
        if len(glob.glob('summary'))>0:
            os.popen('mv summary '+dirname)
        time.sleep(10)
        return
