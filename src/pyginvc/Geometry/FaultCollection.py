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
        '''Vertically concatenate the named 2-D array of every sub-fault.'''
        return np.vstack([getattr(f, name) for f in self.flts])

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
