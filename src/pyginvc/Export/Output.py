# Written by Zhao Bin, when he is at UC Berkeley. April 18 2016
# Mod by Zhao Bin, Dec. 7, 2018. We use print() now.
#
# Output is a thin facade: it extracts arrays/paths from the fault, data,
# green and solution objects and delegates all file writing to the
# object-free helpers in pyginvc.Export.OutputUtils.
import numpy as np
import logging
from pyginvc.Geometry.Triangle import Triangle
from pyginvc.Export import OutputUtils as ou

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

    # ------------------------------------------------------------------ #
    # Main entry point
    # ------------------------------------------------------------------ #
    def OutputSolution(self):
        '''
        Output solutions to hdf5 file and ascii files

        Returns
        -------
        None.

        '''
        obs = {key: getattr(self.data, key) for key in
               ('llh_gps', 'llh_lev', 'llh_sar', 'unit',
                'd_gps', 'd_lev', 'd_sar')}
        obs['ndim'] = self.data.ndim

        fault = {'origin': self.flt.origin, 'flttype': self.flttype}
        if self.flttype == 'rectangle':
            fault['dis_geom'] = self.flt.dis_geom_grid
            fault['geom']     = self.flt.geom_grid
        else:
            fault['vertex_enu'] = self.flt.vertex_enu
            fault['vertex_llh'] = self.flt.vertex_llh
            fault['element']    = self.flt.element

        solution = {key: getattr(self.sol, key)
                    for key in ('slip', 'r', 'misfit', 'dhat')}

        ou.write_solution_h5('solutions.h5', obs, fault, solution)

        # output ascii files
        from pyginvc.Forward.Forward import Forward
        from pyginvc.Forward.TriForward import TriForward
        if isinstance(self.sol, (Forward, TriForward)):
            self.WriteData()
            self.archive_outfile()
        else:
            self.WriteSummary()
            self.WriteSlipModel()
            self.WriteGMTSlip()
            self.WriteLCurve()
            self.WriteData()
            self.archive_outfile()

    # ------------------------------------------------------------------ #
    # Slip in GMT format
    # ------------------------------------------------------------------ #
    def _rectangle_vertices(self):
        '''Return (nf, 4, 3) [lon, lat, depth] patch corners.'''
        felem = self.flt.FaultGeom2AllVertex()
        corners = [('lt_lon', 'lt_lat', 'lt_dep'),
                   ('rt_lon', 'rt_lat', 'rt_dep'),
                   ('rb_lon', 'rb_lat', 'rb_dep'),
                   ('lb_lon', 'lb_lat', 'lb_dep')]
        return np.stack([np.column_stack([felem[c[0]], felem[c[1]], felem[c[2]]])
                         for c in corners], axis=1)

    def _triangle_vertices(self):
        '''Return (nf, 3, 3) [lon, lat, depth] triangle corners.'''
        v = self.flt.vertex_llh[:, [1, 0, 2]]
        return v[self.flt.element - 1]

    def WriteGMTSlip(self):
        '''
        Ouput slip model in GMT format

        Returns
        -------
        None.

        '''
        nf       = self.flt.nf
        slip     = self.sol.slip[-1, 0:nf*3].reshape(nf, 3)
        slip_mag = np.sqrt(slip[:, 0]**2 + slip[:, 1]**2)

        if self.flttype == 'rectangle':
            vertices = self._rectangle_vertices()
            fmt      = '%10.4f %10.4f %10.4f'
            header   = '# Rectangle fault slip model'
        else:
            vertices = self._triangle_vertices()
            fmt      = '%10.3f %10.3f %10.3f'
            header   = '# Triangle fault slip model'

        ou.write_gmt_slip('ss_slip_llh.gmtlin', vertices, slip[:, 0],
                          header=header, fmt=fmt)
        ou.write_gmt_slip('ds_slip_llh.gmtlin', vertices, slip[:, 1],
                          header=header, fmt=fmt)
        ou.write_gmt_slip('ts_slip_llh.gmtlin', vertices, slip_mag,
                          header=header, fmt=fmt)

        logging.info('Fault slip distribution is ouput in geography system.')

    def write_rectangle_slip(self, filename, value):
        """
        Write slip distribution in GMT format
        """
        ou.write_gmt_slip(filename, self._rectangle_vertices(), value,
                          header='# Rectangle fault slip model',
                          fmt='%10.4f %10.4f %10.4f')
        logging.info('Rectangular fault slip distribution is ouput in GMT format.')

    def write_triangle_slip(self, filename, value):
        """
        Write slip distribution in GMT format
        """
        ou.write_gmt_slip(filename, self._triangle_vertices(), value,
                          header='# Triangle fault slip model',
                          fmt='%10.3f %10.3f %10.3f')
        logging.info('Triangular fault slip distribution is ouput in GMT format.')

    # ------------------------------------------------------------------ #
    # Slip models
    # ------------------------------------------------------------------ #
    def WriteSlipModel(self):
        '''
        Output slip model in FaultGeom or Triangular format

        Returns
        -------
        None.

        '''
        slip = np.asarray(self.sol.slip)
        rake_beta = getattr(self.green, 'rake_beta', 0.0)

        if self.flttype == 'rectangle':
            ou.write_slip_models('.', self.flt.geom_grid, slip,
                                 rake_beta=rake_beta, kind='rectangle')
        else:
            ou.write_slip_models('.', self.flt.element, slip,
                                 rake_beta=rake_beta, kind='triangle')

    # ------------------------------------------------------------------ #
    # Statistics
    # ------------------------------------------------------------------ #
    def WriteSummary(self):
        '''
        Ouput summary of Inversion.

        Returns
        -------
        None.

        '''
        if self.flttype == 'rectangle':
            nstrike, ndip = self.flt.nsegs, self.flt.ndeps
        else:
            nstrike, ndip = len(self.flt.vertex_llh), len(self.flt.element)

        ou.write_summary(
            'summary',
            ndim       = self.data.ndim,
            n_gps      = len(self.data.llh_gps),
            n_sar      = len(self.data.llh_sar),
            n_lev      = len(self.data.llh_lev),
            flttype    = self.flttype,
            nstrike    = nstrike,
            ndip       = ndip,
            smo_facts  = self.sol.smo_facts,
            misfit_gps = self.sol.misfit_gps,
            misfit_sar = self.sol.misfit_sar,
            n_d_gps    = len(self.data.d_gps),
            n_d_sar    = len(self.data.d_sar),
            moment     = self.sol.moment,
        )

    def WriteLCurve(self):
        """
        Output roughness misfit and smooth factors, which can be used to plot
        tradeoff between roughness and fittness
        Written by Zhao Bin @ UC Berkeley April 18 2016

        Input:
            ruff      = array, shape(m,)   roughness
            misfit    = array, shape(m,)   misfit
            smo_facts = array, shape(m,)   smooth factor
        Output:
            ruff_misfit: file
        """
        ok = ou.write_lcurve('ruff_misfit', self.sol.ruff, self.sol.misfit,
                             self.sol.smo_facts, self.sol.moment)
        if not ok:
            logging.warning('Only %d smoothing factors.'
                            % len(self.sol.smo_facts))
            return
        logging.info('Roughness vs misfit is printed into ruff_misfit.')

    # ------------------------------------------------------------------ #
    # Observed / modeled / residual data
    # ------------------------------------------------------------------ #
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

        ndim    = self.data.ndim
        obs     = self.data.d_gps.reshape(nsta, ndim)
        mod     = dhat[0:len_gps].reshape(nsta, ndim)
        sigma   = (self.data.get_sigma()
                   if hasattr(self.data, 'get_sigma') else None)
        station = getattr(self.data, 'station_gps', None)

        ou.write_gps_data(self.data.llh_gps, obs, mod,
                          sigma=sigma, station=station)
        logging.info('Modeled and residual GPS points are output.')

    def Write_LEV_Data(self, dhat, r):
        """Write modeled and residual leveling data"""
        len_gps = len(self.data.d_gps)
        len_lev = len(self.data.d_lev)

        if len(self.data.llh_lev) > 0:
            d_lev    = self.data.d_lev
            dhat_lev = dhat[len_gps:len_gps + len_lev]

            ou.write_lev_data(self.data.llh_lev, d_lev, dhat_lev)
            logging.info('%d modeled and residual LEV points are output.'
                         % len_lev)

    def Write_SAR_Data(self, dhat, r):
        """Write modeled and residual SAR data"""
        len_gps = len(self.data.d_gps)
        len_lev = len(self.data.d_lev)
        len_sar = len(self.data.d_sar)
        len_geod = len_gps + len_lev

        if len_sar < 1:
            return

        d_sar = self.data.d_sar
        dhat_sar = dhat[len_geod:len_geod + len_sar]

        ou.write_sar_data(self.data.llh_sar, d_sar, dhat_sar,
                          n_sar=self.data.n_sar)
        logging.info('%d modeled and residual SAR points are output.' % len_sar)

    # ------------------------------------------------------------------ #
    # Archiving
    # ------------------------------------------------------------------ #
    def archive_outfile(self):
        '''
        Move output files into the archive directory (self.archdir).
        '''
        ou.archive_files(self.archdir, ou.DEFAULT_ARCHIVE_PATTERNS)
