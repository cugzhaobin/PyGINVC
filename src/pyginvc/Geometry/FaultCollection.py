#!/usr/bin/env python
"""Aggregate view over a collection of rectangular fault patches."""

import numpy as np


class FaultCollection(object):
    '''
    FaultCollection represents a collection of Fault instances.

    It holds the individual faults and exposes aggregated geometry on demand,
    so no concatenation happens until an aggregated attribute is accessed.

    A sub-fault must provide the attributes nf, origin, dis_geom_grid,
    geom_grid, nsegs and ndeps.

    configs is an optional per-fault configuration. It is opaque to this
    class and is only paired with the faults by pairs().
    '''

    def __init__(self, flts, configs=None):
        '''
        Constructor.

        Input:
            flts     = list of Fault instances
            configs  = optional list of per-fault configurations
        '''
        if configs is not None and len(flts) != len(configs):
            raise ValueError('flts and configs must have the same length')
        self.flts    = flts
        self.configs = configs

    def __iter__(self):
        return iter(self.flts)

    def __len__(self):
        return len(self.flts)

    def pairs(self):
        '''Return (fault, config) for each sub-fault.'''
        if self.configs is None:
            return zip(self.flts, [None] * len(self.flts))
        return zip(self.flts, self.configs)

    def stack(self, name):
        '''Vertically concatenate the named 2-D array of every sub-fault.

        All sub-faults must have the same number of columns.  The usual cause
        of a mismatch is mixing doSubFault=True faults (10-column geom_grid)
        with doSubFault=False ones (7 columns) in the same inversion.
        '''
        arrays = [np.asarray(getattr(f, name)) for f in self.flts]
        ncols  = sorted({a.shape[1] for a in arrays if a.ndim == 2})
        if len(ncols) > 1:
            raise ValueError(
                "Cannot stack %r: sub-faults have different column counts "
                "%s; check that doSubFault is consistent for all faults."
                % (name, ncols))
        return np.vstack(arrays)

    @property
    def nf(self):
        return sum(f.nf for f in self.flts)

    @property
    def origin(self):
        return self.flts[0].origin

    @property
    def dis_geom_grid(self):
        return self.stack('dis_geom_grid')

    @property
    def geom_grid(self):
        return self.stack('geom_grid')

    @property
    def nsegs(self):
        return [f.nsegs for f in self.flts]

    @property
    def ndeps(self):
        return [f.ndeps for f in self.flts]

    def total_moment(self, slip, default_modulus=3e10):
        '''Sum the seismic moment of every fault and return it with Mw.'''
        slip = np.asarray(slip).reshape(-1)
        Mo   = []
        off  = 0
        for flt, cfg in self.pairs():
            islip = slip[off:off + 3 * flt.nf].reshape(-1, 3)
            green = (cfg or {}).get('green') or {}
            G     = float(green.get('modulus', default_modulus))
            Mo.append(flt.moment(islip, shear_modulus=G)[0])
            off += 3 * flt.nf
        Mo_total = float(np.sum(Mo))
        Mw = (2.0 / 3.0) * np.log10(Mo_total) - 6.067 if Mo_total > 0 else 0.0
        return Mo_total, Mw
