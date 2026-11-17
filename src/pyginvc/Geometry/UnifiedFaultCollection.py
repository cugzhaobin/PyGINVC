#!/usr/bin/env python
"""Aggregate rectangular and triangular faults behind one interface.

``UnifiedFaultCollection`` mirrors the subset of ``FaultCollection`` used by
the inversion and output layers, but it works with both ``Patch.Fault``
rectangular faults and ``Triangle`` meshes.  Every fault must be loaded with
the same reference origin; this is checked when the collection is built.
"""

import logging

import numpy as np

from pyginvc.Geometry.Patch import Fault
from pyginvc.Geometry.Triangle import Triangle

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%M-%Y %H:%M:%S")


def fault_kind(flt):
    """Return ``'rectangle'`` or ``'triangle'`` for a loaded fault object."""
    if isinstance(flt, Triangle):
        return 'triangle'
    if isinstance(flt, Fault):
        return 'rectangle'
    raise TypeError('unsupported fault object: %s' % type(flt).__name__)


def safe_name(name, default):
    """Return ``name`` reduced to characters that are safe in file names."""
    if not name:
        return default
    text = ''.join(ch if (ch.isalnum() or ch in '._-') else '_'
                   for ch in str(name))
    return text or default


class UnifiedFaultCollection(object):
    """Collection of rectangular and/or triangular faults.

    Parameters
    ----------
    flts : list
        Loaded ``Fault`` and/or ``Triangle`` instances.
    configs : list, optional
        Per-fault YAML configuration, aligned with ``flts``.
    """

    def __init__(self, flts, configs=None):
        flts = list(flts)
        if not flts:
            raise ValueError('At least one fault is required.')
        if configs is not None and len(flts) != len(configs):
            raise ValueError('flts and configs must have the same length')

        self.flts = flts
        self.configs = list(configs) if configs is not None else [None] * len(flts)

        # Mixed Green's functions are all computed relative to the same
        # origin.  A mismatch would silently shift one fault relative to the
        # others, so make it an explicit error.
        origins = [np.asarray(getattr(f, 'origin', []), dtype=float)
                   for f in self.flts]
        if any(o.size != 2 for o in origins):
            raise ValueError('Every fault must carry a two-element origin.')
        origin0 = origins[0]
        for flt, origin in zip(self.flts, origins):
            if not np.allclose(origin, origin0, rtol=0.0, atol=1e-9):
                raise ValueError(
                    'Faults use different origins (%s and %s). Load all '
                    'faults with one common origin.'
                    % (origin0.tolist(), origin.tolist()))

    def __iter__(self):
        return iter(self.flts)

    def __len__(self):
        return len(self.flts)

    def pairs(self):
        """Yield ``(fault, config)`` for every fault."""
        return zip(self.flts, self.configs)

    @property
    def nf(self):
        """Total number of fault patches."""
        return int(sum(int(f.nf) for f in self.flts))

    @property
    def origin(self):
        return np.asarray(self.flts[0].origin, dtype=float)

    def kinds(self):
        return [fault_kind(f) for f in self.flts]

    def names(self):
        names, seen = [], set()
        for i, (flt, cfg) in enumerate(self.pairs()):
            raw = (cfg or {}).get('name')
            name = safe_name(raw, 'fault%02d' % (i + 1))
            base, suffix = name, 1
            while name in seen:
                suffix += 1
                name = '%s_%d' % (base, suffix)
            seen.add(name)
            names.append(name)
        return names

    def groups(self):
        """Yield ``(index, fault, config, kind, name)`` in model order."""
        for i, ((flt, cfg), kind, name) in enumerate(
                zip(self.pairs(), self.kinds(), self.names())):
            yield i, flt, cfg, kind, name

    @property
    def offsets(self):
        """Start column of each fault in the flattened slip vector."""
        out, offset = [], 0
        for flt in self.flts:
            out.append(offset)
            offset += 3 * int(flt.nf)
        return out

    @property
    def nf3(self):
        return 3 * self.nf

    def rake_betas(self):
        """Per-fault ``rake_beta`` in degrees (zero when not configured)."""
        betas = []
        for flt, cfg in self.pairs():
            green = (cfg or {}).get('green') or {}
            betas.append(float(green.get('rake_beta', 0.0) or 0.0))
        return betas

    def total_moment(self, slip, default_modulus=3e10):
        """Sum the seismic moment of every fault and return ``(Mo, Mw)``."""
        slip = np.asarray(slip, dtype=float).reshape(-1)
        if slip.size != self.nf3:
            raise ValueError('slip has %d values but %d are expected.'
                             % (slip.size, self.nf3))

        moments = []
        for flt, cfg, offset in zip(self.flts, self.configs, self.offsets):
            nf = int(flt.nf)
            islip = slip[offset:offset + 3 * nf].reshape(nf, 3)
            green = (cfg or {}).get('green') or {}
            modulus = float(green.get('modulus', default_modulus))
            moments.append(float(flt.moment(islip, shear_modulus=modulus)[0]))

        mo_total = float(np.sum(moments))
        mw = (2.0 / 3.0) * np.log10(mo_total) - 6.067 if mo_total > 0 else 0.0
        return mo_total, mw

    def summary_rows(self, slip=None):
        """Per-fault metadata used by the summary writer and reports."""
        rows = []
        for i, flt, cfg, kind, name in self.groups():
            if kind == 'rectangle':
                counts = (int(flt.nsegs), int(flt.ndeps))
            else:
                counts = (len(flt.vertex_llh), len(flt.element))
            row = {
                'index': i,
                'name': name,
                'kind': kind,
                'nf': int(flt.nf),
                'nstrike': counts[0],
                'ndip': counts[1],
            }
            if slip is not None:
                offset = self.offsets[i]
                islip = np.asarray(slip).reshape(-1)[offset:offset + 3 * flt.nf]
                row['mean_slip'] = float(np.mean(
                    np.linalg.norm(islip.reshape(-1, 3), axis=1)))
            rows.append(row)
        return rows
