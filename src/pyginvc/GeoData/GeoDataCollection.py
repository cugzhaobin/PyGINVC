#!/usr/bin/env python
"""Aggregate view over a collection of per-file geodetic data objects."""

import numpy as np
from scipy import linalg

from pyginvc.GeoData.GPSData import GPSData
from pyginvc.GeoData.SARData import SARData
from pyginvc.GeoData.LEVData import LEVData


class GeoDataCollection(object):
    '''
    GeoDataCollection represents a collection of per-file geodetic data
    objects (GPSData, SARData and LEVData).

    Each input file is loaded into its own object by the caller and the
    objects are passed as a single flat list, e.g.

        GeoDataCollection([gps1, gps2, sar1, sar2, lev1],
                          [gps_cfg1, gps_cfg2, sar_cfg1, sar_cfg2, lev_cfg])

    The kind of every object is determined from the object itself. Aggregated
    attributes are concatenated on demand and mirror the GeoData interface,
    so the collection can be consumed anywhere a GeoData instance is expected.

    configs is an optional per-file configuration aligned with files. It is
    paired with the files by pairs(); the weight properties additionally read
    the optional 'weight' key from each config.
    '''

    def __init__(self, files, configs=None):
        '''
        Constructor.

        Input:
            files    = list of loaded GPSData / SARData / LEVData objects,
                       one object per input file
            configs  = optional list of per-file configurations, aligned
                       with files
        '''
        files = list(files)
        if configs is not None and len(files) != len(configs):
            raise ValueError('files and configs must have the same length')
        self.files   = files
        self.configs = configs
        for f in files:
            self.kind_of(f)

    # ------------------------------------------------------------------ #
    # Basic container protocol
    # ------------------------------------------------------------------ #
    def __iter__(self):
        return iter(self.files)

    def __len__(self):
        return len(self.files)

    @staticmethod
    def kind_of(obj):
        '''Return the kind ('gps', 'sar' or 'lev') of a file object.'''
        if isinstance(obj, GPSData):
            return 'gps'
        if isinstance(obj, SARData):
            return 'sar'
        if isinstance(obj, LEVData):
            return 'lev'
        raise TypeError('unsupported data object: %s' % type(obj).__name__)

    def pairs(self, kind=None):
        '''Return (file object, config) for each file, optionally by kind.'''
        cfgs = self.configs if self.configs is not None else [None] * len(self.files)
        for f, c in zip(self.files, cfgs):
            if kind is None or self.kind_of(f) == kind:
                yield f, c

    def gps_pairs(self):
        return self.pairs('gps')

    def sar_pairs(self):
        return self.pairs('sar')

    def lev_pairs(self):
        return self.pairs('lev')

    @property
    def gps_data(self):
        return [f for f in self.files if self.kind_of(f) == 'gps']

    @property
    def sar_data(self):
        return [f for f in self.files if self.kind_of(f) == 'sar']

    @property
    def lev_data(self):
        return [f for f in self.files if self.kind_of(f) == 'lev']

    # ------------------------------------------------------------------ #
    # Concatenation helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def _stack(objs, name, shape, axis=0):
        '''Concatenate the named array of every object, dropping empties.'''
        arrays = [np.asarray(getattr(o, name)) for o in objs]
        arrays = [a for a in arrays if a.size > 0]
        if not arrays:
            return np.empty(shape)
        if axis == 0:
            return np.vstack(arrays)
        return np.concatenate(arrays).ravel()

    @staticmethod
    def _merge_cov(covs):
        '''Merge 1-D variances / 2-D covariance matrices of several files.'''
        covs = [np.asarray(c) for c in covs]
        covs = [c for c in covs if c.size > 0]
        if not covs:
            return np.array([])
        if all(c.ndim == 1 for c in covs):
            return np.concatenate([c.ravel() for c in covs])
        return linalg.block_diag(*[np.diag(c) if c.ndim == 1 else c
                                   for c in covs])

    # ------------------------------------------------------------------ #
    # GPS aggregation
    # ------------------------------------------------------------------ #
    @property
    def llh_gps(self):
        return self._stack(self.gps_data, 'llh_gps', (0, 2))

    @property
    def d_gps(self):
        return self._stack(self.gps_data, 'd_gps', (0,), axis=1)

    @property
    def cov_gps(self):
        return self._merge_cov([g.cov_gps for g in self.gps_data])

    @property
    def n_gps(self):
        return np.array([len(g.llh_gps) for g in self.gps_data], dtype=int)

    @property
    def station_gps(self):
        gps = self.gps_data
        if not gps:
            return np.array([], dtype=object)
        stations = [np.asarray(g.station_gps).ravel() for g in gps]
        if any(len(s) != len(g.llh_gps) for s, g in zip(stations, gps)):
            return np.array([], dtype=object)
        return np.concatenate(stations).astype(object)

    @property
    def ndim(self):
        dims = {int(g.ndim) for g in self.gps_data}
        if len(dims) > 1:
            raise ValueError('GPS files have inconsistent ndim: %s'
                             % sorted(dims))
        return dims.pop() if dims else 2

    # ------------------------------------------------------------------ #
    # InSAR aggregation
    # ------------------------------------------------------------------ #
    @property
    def llh_sar(self):
        return self._stack(self.sar_data, 'llh_sar', (0, 2))

    @property
    def d_sar(self):
        return self._stack(self.sar_data, 'd_sar', (0,), axis=1)

    @property
    def unit(self):
        return self._stack(self.sar_data, 'unit', (0, 3))

    @property
    def cov_sar(self):
        return self._merge_cov([s.cov_sar for s in self.sar_data])

    @property
    def n_sar(self):
        return np.array([len(s.llh_sar) for s in self.sar_data], dtype=int)

    # ------------------------------------------------------------------ #
    # Leveling aggregation
    # ------------------------------------------------------------------ #
    @property
    def llh_lev(self):
        return self._stack(self.lev_data, 'llh_lev', (0, 2))

    @property
    def d_lev(self):
        return self._stack(self.lev_data, 'd_lev', (0,), axis=1)

    @property
    def n_lev(self):
        return np.array([len(l.llh_lev) for l in self.lev_data], dtype=int)

    # ------------------------------------------------------------------ #
    # Weights
    # ------------------------------------------------------------------ #
    @staticmethod
    def _var1d(cov):
        '''Extract per-observation variance from 1-D or 2-D covariance.'''
        if cov.size == 0:
            return np.array([])
        v = cov if cov.ndim == 1 else np.diag(cov)
        return np.where(v > 0, v, 1.0)

    @property
    def w_gps(self):
        w   = 1.0 / np.sqrt(self._var1d(self.cov_gps))
        row = 0
        for gps, cfg in self.gps_pairs():
            n = len(gps.llh_gps) * self.ndim
            w[row:row + n] *= (cfg or {}).get('weight', 1.0)
            row += n
        return w

    @property
    def w_lev(self):
        return np.ones(self.d_lev.size)

    @property
    def w_sar(self):
        w   = 1.0 / np.sqrt(self._var1d(self.cov_sar))
        row = 0
        for sar, cfg in self.sar_pairs():
            n = len(sar.llh_sar)
            w[row:row + n] *= (cfg or {}).get('weight', 1.0)
            row += n
        return w

    @property
    def w_geod(self):
        return np.concatenate([self.w_gps, self.w_lev])

    def get_weights(self):
        '''Return diagonal weights for GPS+leveling and InSAR data.'''
        return np.diag(self.w_geod), np.diag(self.w_sar)
