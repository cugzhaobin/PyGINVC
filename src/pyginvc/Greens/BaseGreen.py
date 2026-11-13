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
        self.flt        = flt
        self.data       = data
        self.dict_green = dict_green

        self.G          = np.empty((0,0))
        self.G_sar      = np.empty((0,0))
        self.G_gps_ramp = np.empty((0,0))
        self.G_sar_ramp = np.empty((0,0))

        self.modulus    = float(dict_green.get('modulus', 3e10))
        self.greenfile  = dict_green['greenfile']
        self.rake_beta  = dict_green.get('rake_beta', 0)


    def _cache_disabled(self):
        '''True when no Green's function cache should be read or written.'''
        path = self.greenfile
        return path is None or str(path).strip().upper() == 'NONE'

    def load_greens(self):
        '''
        Load Green's function.
        '''
        if self._cache_disabled():
            logging.info('Green function cache is disabled.')
            return False
        greenfile = self.greenfile
        if not os.path.exists(greenfile):
            logging.info('Green file {} does not exist.'.format(greenfile))
            return False
        
        ext = greenfile.split('.')[-1].lower()
        if ext == 'h5':
            with h5py.File(greenfile, 'r') as h5:
                G          = h5['G'][()]
                G_sar      = h5['G_sar'][()]
                G_gps_ramp = h5['G_gps_ramp'][()]
                G_sar_ramp = h5['G_sar_ramp'][()]
        elif ext in ('npy', 'npz'):
            # np.load returns a bare ndarray for .npy files, which cannot
            # store the four matrices of a cache; only .npz is a keyed
            # archive.  Anything incomplete is regenerated instead of
            # raising KeyError.
            dat = np.load(greenfile, allow_pickle=False)
            if isinstance(dat, np.ndarray):
                logging.warning(
                    'Green file {} is a bare .npy array without '
                    'G_sar/ramp terms. Regenerating.'.format(greenfile))
                return False
            try:
                missing = [k for k in ('G', 'G_sar', 'G_gps_ramp',
                                       'G_sar_ramp') if k not in dat.files]
                if missing:
                    logging.warning(
                        'Green file {} is missing {}. Regenerating.'.format(
                            greenfile, ', '.join(missing)))
                    return False
                G          = dat['G']
                G_sar      = dat['G_sar']
                G_gps_ramp = dat['G_gps_ramp']
                G_sar_ramp = dat['G_sar_ramp']
            finally:
                dat.close()
        else:
            logging.warning('Unsupported green file {}'.format(greenfile))
            return False

        # A cached green file must match the current data layout.  Otherwise
        # (e.g. the file was generated with InSAR points that are not used in
        # this run) reusing it gives matrices with the wrong number of rows.
        nf3    = 3 * self.flt.nf
        n_geod = len(self.data.llh_gps) * self.data.ndim + len(self.data.llh_lev)
        n_sar  = len(self.data.llh_sar)
        if (G.shape, G_sar.shape) != ((n_geod, nf3), (n_sar, nf3)):
            logging.warning(
                'Green file {} does not match the current data: '
                'G {} (expected {}), G_sar {} (expected {}). Regenerating.'.format(
                    greenfile, G.shape, (n_geod, nf3),
                    G_sar.shape, (n_sar, nf3)))
            return False

        self.G          = G
        self.G_sar      = G_sar
        self.G_gps_ramp = G_gps_ramp
        self.G_sar_ramp = G_sar_ramp
        logging.info('Load Greens function from {}'.format(greenfile))

        gps_ramp = self.dict_green.get('gps_ramp', False)
        sar_ramp = self.dict_green.get('sar_ramp', False)
        
        if not gps_ramp:
            len_gps         = self.G.shape[0]
            self.G_gps_ramp = np.empty((len_gps,0))
        if not sar_ramp:
            len_sar         = self.G_sar.shape[0]
            self.G_sar_ramp = np.empty((len_sar,0))
        return True

    def save_greens(self, greenfile=None):
        """Save Green's functions into hdf5 file"""
        if greenfile is None:
            greenfile = self.greenfile
        if greenfile is None or str(greenfile).strip().upper() == 'NONE':
            logging.info('Green function cache is disabled; not saving.')
            return
        if not greenfile or str(greenfile).strip().upper() == 'SAVE':
            greenfile = 'green.h5'
        if os.path.exists(greenfile):
            os.remove(greenfile)
        with h5py.File(greenfile, 'w') as h5:
            h5.create_dataset('G', data = self.G, compression='gzip')
            h5.create_dataset('G_sar', data = self.G_sar, compression='gzip')
            h5.create_dataset('G_gps_ramp', data = self.G_gps_ramp, compression='gzip')
            h5.create_dataset('G_sar_ramp', data = self.G_sar_ramp, compression='gzip')
            logging.info('Green function is saved to {}'.format(greenfile))

    def rotate_greens(self, greenfile='green_rot.h5'):
        '''
        Rotate the green function.
           ds'     ds      ss'
            \      |      /
              \    |    /
                \  |  /  deta
        __________\|/____________ss

        The strike/dip columns of every patch are rotated in place with one
        vectorised operation per matrix, instead of a Python loop over patches.
        '''
        beta = self.rake_beta
        if beta is None: return
        rbeta = np.deg2rad(beta)
        cosb, sinb = np.cos(rbeta), np.sin(rbeta)

        nf = self.flt.nf if self.flt is not None else None

        for name in ('G', 'G_sar'):
            M = getattr(self, name)
            if M.size == 0:
                continue
            ncol = 3 * nf if nf is not None else M.shape[1]
            u = M[:, 0:ncol:3]
            v = M[:, 1:ncol:3]
            M = np.array(M, dtype=float, copy=True)
            M[:, 0:ncol:3] = u * cosb + v * sinb
            M[:, 1:ncol:3] = -u * sinb + v * cosb
            setattr(self, name, M)

        # The rotated matrices are only needed in memory.  Writing them to the
        # current working directory used to litter ``./green_new.h5`` and
        # could fail on read-only directories or clash between concurrent runs.
        logging.info('Greens function rotated in memory (rake_beta=%.6g deg).'
                     % beta)

    @staticmethod
    def make_green_sar_ramp(xy):
        '''
        Per-track InSAR ramp: constant offset + east/north tilt.

        Input:
            xy: local coordinate of InSAR data, array (n,2) [east, north]
        Output:
            G: array (n,3), columns [1, east, north]
        '''
        xy = np.asarray(xy, dtype=float)
        if xy.size == 0:
            return np.zeros((0, 3))
        xy = xy.reshape(-1, 2)
        G      = np.zeros((xy.shape[0], 3))
        G[:,0] = 1.0          # per-track constant offset
        G[:,1] = xy[:,0]
        G[:,2] = xy[:,1]
        return G

    def make_green_gps_ramp(self, xy, dim, method=2):
        """Build the GPS ramp columns for local coordinates ``xy``.

        Input:
            xy     = (n,2) local [east, north] in km
            dim    = 2 or 3 GPS components per station
            method = 1 (translation + east/north term) or
                     2 (spherical rotation, using ``self.data.llh_gps``)
        """

        if len(xy) == 0:
            return np.zeros((0,3))

        xy = np.asarray(xy, dtype=float).reshape(-1, 2)
        try:
            method = int(method)
        except (TypeError, ValueError):
            raise ValueError('Unsupported GPS ramp method: %r' % (method,))

        n = len(xy)
        G = np.zeros((dim*n, 3))
        east_row = dim * np.arange(n)
        north_row = east_row + 1

        if method == 1:
            G[east_row, 0]  = 1.0
            G[north_row, 1] = 1.0
            G[east_row, 2]  = -xy[:, 0]
            G[north_row, 2] = xy[:, 1]
        elif method == 2:
            R    = 6370
            rllh = np.deg2rad(self.data.llh_gps)
            sin_lat = np.sin(rllh[:n, 0])
            cos_lat = np.cos(rllh[:n, 0])
            sin_lon = np.sin(rllh[:n, 1])
            cos_lon = np.cos(rllh[:n, 1])
            G[east_row, 0]  = R * (-sin_lat * cos_lon)
            G[east_row, 1]  = R * (-sin_lat * sin_lon)
            G[east_row, 2]  = R * cos_lat
            G[north_row, 0] = R * sin_lon
            G[north_row, 1] = R * (-cos_lon)
        else:
            raise ValueError('Unsupported GPS ramp method: %r' % (method,))
        return G

    @staticmethod
    def _convert_to_local(llh, origin):
        return gt.llh2localxy_batch(llh, origin)
