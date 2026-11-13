# Written by Zhao Bin, when he is at UC Berkeley. April 18 2016
# Mod by Zhao Bin, Dec. 7, 2018. We use print() now.
#
# Output is the shared output implementation.  It extracts arrays/paths from
# a fault, data, green and solution object, and delegates every file write to
# the object-free helpers in pyginvc.Export.OutputUtils.
#
# Multi-fault / multi-file output is implemented by
# pyginvc.Export.MultiOutput.MultiOutput, which subclasses Output and only
# overrides the small adapter hooks below.
import os
import logging
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np

from pyginvc.Geometry.Triangle import Triangle
from pyginvc.Export import OutputUtils as ou

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%m-%Y %H:%M:%S")


@dataclass
class ObsBlock:
    '''One output group of observations (one GPS file, level file or track).'''
    kind: str
    name: str
    llh: Any
    d_obs: Any
    ndim: int = 1
    unit: Optional[Any] = None
    sigma: Optional[Any] = None
    station: Optional[Any] = None
    n_sar: Optional[Any] = None


@dataclass
class FaultGroup:
    '''One fault to be written; ``kind`` is 'rectangle' or 'triangle'.'''
    flt: Any
    name: str = 'fault01'
    kind: str = 'rectangle'
    rake_beta: float = 0.0


@dataclass
class SolutionView:
    '''Common view of an inversion / forward solution.'''
    slip: Any
    dhat: Optional[Any] = None
    r: Optional[Any] = None
    misfit: Optional[Any] = None
    misfit_gps: Optional[Any] = None
    misfit_sar: Optional[Any] = None
    smo_facts: Optional[Any] = None
    moment: Optional[Any] = None
    ruff: Optional[Any] = None
    ramp_sar: Optional[Any] = None
    sar_ramp_flag: Optional[Any] = None
    forward: bool = False


class Output(object):
    '''
    Output is a class representing output inversions, summary, figures.

    The class works on normalized pieces (ObsBlock / FaultGroup /
    SolutionView).  Single-object users (Output itself) and collection users
    (MultiOutput) only differ in the hook implementations.
    '''

    def __init__(self, flt, data, green, sol=None, archdir=None, outdir='.'):
        '''
        Constructor.

        Parameters
        ----------
        flt : instance
            an instance of Triangle or Fault.
        data : instance
            an instance of GeoData.
        green : instance
            an instance of Okada, Nikkhoo, TriPoly3D.
        sol : instance
            an inversion / forward solution.
        archdir : string, optional
            archive directory.
        outdir : string, optional
            output directory (default '.', preserving the legacy behaviour).
        '''
        self.flt      = flt
        self.data     = data
        self.green    = green
        self.sol      = sol
        self.archdir  = archdir
        self.outdir   = outdir or '.'

        self.flttype = 'rectangle'
        if isinstance(flt, Triangle):
            self.flttype = 'triangle'

    # ------------------------------------------------------------------ #
    # Adapter hooks -- subclasses override these
    # ------------------------------------------------------------------ #
    def _obs_blocks(self):
        '''Return observations as a list of ObsBlock in GPS/LEV/SAR order.'''
        data = self.data
        blocks = []

        llh_gps = np.asarray(getattr(data, 'llh_gps', np.empty((0, 2))))
        if len(llh_gps):
            sigma = (data.get_sigma()
                     if hasattr(data, 'get_sigma') else None)
            blocks.append(ObsBlock(
                'gps', 'gps', llh_gps, np.asarray(data.d_gps),
                ndim=int(data.ndim), sigma=sigma,
                station=getattr(data, 'station_gps', None)))

        llh_lev = np.asarray(getattr(data, 'llh_lev', np.empty(0)))
        if len(llh_lev):
            blocks.append(ObsBlock(
                'lev', 'lev', llh_lev, np.asarray(data.d_lev), ndim=1))

        llh_sar = np.asarray(getattr(data, 'llh_sar', np.empty(0)))
        if len(llh_sar):
            blocks.append(ObsBlock(
                'sar', 'sar', llh_sar, np.asarray(data.d_sar), ndim=1,
                unit=getattr(data, 'unit', None), n_sar=self._n_sar()))

        return blocks

    def _fault_groups(self):
        '''Return faults as a list of FaultGroup (one for the single case).'''
        rake_beta = getattr(self.green, 'rake_beta', 0.0)
        return [FaultGroup(self.flt, 'fault01', self.flttype, rake_beta)]

    def _solution_view(self):
        '''Return a SolutionView, or None when no solution is available.'''
        sol = self.sol
        if sol is None:
            return None

        ramp_sar = self._sar_ramp_effect()
        sar_flags = getattr(sol, 'sar_ramp_flag', None)
        if sar_flags is None and ramp_sar is not None:
            sar_flags = np.asarray(self._n_sar()) > 0

        return SolutionView(
            slip=getattr(sol, 'slip', None),
            dhat=getattr(sol, 'dhat', None),
            r=getattr(sol, 'r', None),
            misfit=getattr(sol, 'misfit', None),
            misfit_gps=getattr(sol, 'misfit_gps', None),
            misfit_sar=getattr(sol, 'misfit_sar', None),
            smo_facts=getattr(sol, 'smo_facts', None),
            moment=getattr(sol, 'moment', None),
            ruff=getattr(sol, 'ruff', None),
            ramp_sar=ramp_sar,
            sar_ramp_flag=sar_flags,
            forward=self._is_forward())

    def _origin(self):
        return self.flt.origin

    def _flttype(self):
        return self.flttype

    def _summary_counts(self, groups):
        '''Return (nstrike, ndip) as scalars or per-fault lists.'''
        if not groups:
            return 0, 0
        if groups[0].kind == 'triangle':
            flt = groups[0].flt
            return len(flt.vertex_llh), len(flt.element)
        if len(groups) == 1:
            return groups[0].flt.nsegs, groups[0].flt.ndeps
        return ([g.flt.nsegs for g in groups],
                [g.flt.ndeps for g in groups])

    def _is_forward(self):
        # kept local to avoid an Output <-> Forward import cycle
        from pyginvc.Forward.Forward import Forward
        from pyginvc.Forward.TriForward import TriForward
        return isinstance(self.sol, (Forward, TriForward))

    # ------------------------------------------------------------------ #
    # Small helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def _as_array(value):
        '''None -> empty array, so empty datasets can still be written.'''
        return np.array([]) if value is None else np.asarray(value)

    @staticmethod
    def _slip_rows(slip):
        '''Return slip as a 2-D array of shape (nmodel, 3*nf).'''
        slip = np.asarray(slip)
        if slip.ndim == 1:
            slip = slip[None, :]
        return slip

    def _n_sar(self):
        '''Per-track InSAR point counts, inferred when not stored.'''
        llh_sar = np.asarray(getattr(self.data, 'llh_sar', np.empty(0)))
        if len(llh_sar) == 0:
            return np.array([], dtype=int)

        n_sar = getattr(self.data, 'n_sar', None)
        if n_sar is None or np.size(n_sar) == 0:
            return np.array([len(llh_sar)], dtype=int)

        n_sar = np.atleast_1d(n_sar).astype(int)
        if int(n_sar.sum()) != len(llh_sar):
            return np.array([len(llh_sar)], dtype=int)
        return n_sar

    def _path(self, filename, outdir=None):
        outdir = outdir or self.outdir
        if os.path.isabs(filename):
            return filename
        return os.path.join(outdir, filename)

    def _sar_ramp_effect(self):
        '''
        Contribution of the estimated ramp to the InSAR observations.

        Uses the last smoothing factor of an inversion solution. Returns
        None when no InSAR ramp was estimated.
        '''
        ramp        = getattr(self.sol, 'ramp', None)
        G_sar_ramp  = getattr(self.green, 'G_sar_ramp', None)
        if (ramp is None or np.size(ramp) == 0
                or G_sar_ramp is None or G_sar_ramp.shape[1] == 0
                or len(self.data.llh_sar) == 0):
            return None

        ramp_last  = np.asarray(ramp, dtype=float)
        ramp_last  = ramp_last[-1] if ramp_last.ndim > 1 else ramp_last
        G_gps_ramp = getattr(self.green, 'G_gps_ramp', np.empty((0, 0)))
        n_gps_ramp = G_gps_ramp.shape[1]
        return np.asarray(G_sar_ramp) @ ramp_last[n_gps_ramp:]

    # ------------------------------------------------------------------ #
    # HDF5 / observation dictionaries
    # ------------------------------------------------------------------ #
    def _obs_dict(self):
        data = self.data
        obs = {key: getattr(data, key) for key in
               ('llh_gps', 'llh_lev', 'llh_sar', 'unit',
                'd_gps', 'd_lev', 'd_sar')}
        obs['ndim'] = data.ndim
        obs['n_sar'] = self._n_sar()

        sv = self._solution_view()
        ramp  = None if sv is None else sv.ramp_sar
        flags = None if sv is None else sv.sar_ramp_flag
        if ramp is not None and np.size(ramp) == 0:
            ramp = None
        if flags is None and ramp is not None:
            flags = np.asarray(obs['n_sar']) > 0
        obs['d_sar_ramp'] = ramp
        obs['sar_ramp']   = flags
        return obs

    def _fault_dict(self):
        groups = self._fault_groups()
        if not groups:
            raise ValueError('No fault is available for HDF5 output.')

        kinds = {g.kind for g in groups}
        flttype = 'triangle' if kinds == {'triangle'} else 'rectangle'
        fault = {'origin': np.asarray(self._origin()), 'flttype': flttype}

        if flttype == 'rectangle':
            fault['dis_geom'] = np.vstack(
                [np.asarray(g.flt.dis_geom_grid) for g in groups])
            fault['geom'] = np.vstack(
                [np.asarray(g.flt.geom_grid) for g in groups])
        else:
            flt = groups[0].flt
            fault['vertex_enu'] = np.asarray(flt.vertex_enu)
            fault['vertex_llh'] = np.asarray(flt.vertex_llh)
            fault['element']    = np.asarray(flt.element)
        return fault

    def _solution_dict(self):
        sv = self._solution_view()
        if sv is None:
            raise ValueError('No solution is available for HDF5 output.')
        rows = self._slip_rows(self._as_array(sv.slip))
        return {
            'slip'  : self._as_array(sv.slip),
            'r'     : self._as_array(sv.r),
            'misfit': self._as_array(sv.misfit),
            'dhat'  : self._as_array(sv.dhat),
            # extra metadata (optional datasets in write_solution_h5)
            'smo_facts'     : self._as_array(sv.smo_facts),
            'ruff'          : self._as_array(sv.ruff),
            'ramp_sar'      : self._as_array(sv.ramp_sar),
            'sar_ramp_flag' : self._as_array(sv.sar_ramp_flag),
            'chosen_index'  : np.array([rows.shape[0] - 1], dtype=int),
        }

    def WriteHDF5(self, filename='solutions.h5'):
        '''Write observations, fault geometry and solution to HDF5.'''
        path = self._path(filename)
        os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
        ou.write_solution_h5(path, self._obs_dict(),
                             self._fault_dict(), self._solution_dict())
        logging.info('Solution is written to %s.' % path)
        return path

    # ------------------------------------------------------------------ #
    # Observed / modeled / residual data
    # ------------------------------------------------------------------ #
    @staticmethod
    def _block_ndata(block):
        n = len(block.llh)
        return n * int(block.ndim or 1) if block.kind == 'gps' else n

    def _write_blocks(self, dhat, outdir, kinds=None, parts=None):
        '''Write every observation block (or only the requested kinds).

        parts limits which of the obs/mod/res triplets are written.
        '''
        dhat = np.asarray(dhat).reshape(-1)
        written = []
        offset = 0
        for block in self._obs_blocks():
            n = len(block.llh)
            if n == 0:
                continue
            ndata = self._block_ndata(block)
            selected = kinds is None or block.kind in kinds

            if block.kind == 'gps':
                ndim = int(block.ndim or 1)
                obs = np.asarray(block.d_obs).reshape(n, ndim)
                mod = dhat[offset:offset + ndata].reshape(n, ndim)
                if selected:
                    written += ou.write_gps_data(
                        block.llh, obs, mod, sigma=block.sigma,
                        station=block.station, outdir=outdir,
                        prefix=block.name, parts=parts)
            elif block.kind == 'lev':
                obs = np.asarray(block.d_obs).reshape(-1)
                mod = dhat[offset:offset + ndata].reshape(-1)
                if selected:
                    written += ou.write_lev_data(
                        block.llh, obs, mod, outdir=outdir,
                        prefix=block.name, parts=parts)
            elif block.kind == 'sar':
                obs = np.asarray(block.d_obs).reshape(-1)
                mod = dhat[offset:offset + ndata].reshape(-1)
                if selected:
                    written += ou.write_sar_data(
                        block.llh, obs, mod, n_sar=block.n_sar,
                        outdir=outdir, prefix=block.name, parts=parts)

            offset += ndata
        return written

    def _write_sar_ramp(self, ramp, outdir, sv):
        if ramp is None or not np.size(ramp):
            return []
        sar = next((b for b in self._obs_blocks()
                    if b.kind == 'sar' and len(b.llh) > 0), None)
        if sar is None:
            return []
        flags = None if sv is None else sv.sar_ramp_flag
        return ou.write_sar_ramp_data(sar.llh, ramp, n_sar=sar.n_sar,
                                      outdir=outdir, flags=flags)

    def WriteData(self, dhat=None, r=None, outdir=None, ramp=None, parts=None):
        '''
        Output modeled and residual GPS/LEV/SAR data.

        Input:
            dhat  = modeled data, shape(m,); default self.sol.dhat
            r     = kept for backward compatibility (residuals are always
                    recomputed as obs - model, as in the historical code)
            outdir= output directory; default self.outdir
            ramp  = InSAR ramp contribution; default from the solution
            parts = subset of ('obs', 'mod', 'res'); default all three
        '''
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)

        sv = self._solution_view()
        if dhat is None:
            if sv is None or sv.dhat is None:
                raise ValueError('No dhat given and the solution has none.')
            dhat = sv.dhat

        written = self._write_blocks(dhat, outdir, parts=parts)

        # The ramp contribution is part of the modeled data, so it is only
        # written when mod/res are requested.
        if parts is None or any(p in ('mod', 'res') for p in parts):
            if ramp is None and sv is not None:
                ramp = sv.ramp_sar
            written += self._write_sar_ramp(ramp, outdir, sv)

        logging.info('%d data-fit file(s) written to %s.'
                     % (len(written), outdir))
        return written

    def Write_GPS_Data(self, dhat, r=None, outdir=None):
        '''Write observed/modeled/residual GPS data only.'''
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)
        return self._write_blocks(dhat, outdir, kinds={'gps'})

    def Write_LEV_Data(self, dhat, r=None, outdir=None):
        '''Write observed/modeled/residual leveling data only.'''
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)
        return self._write_blocks(dhat, outdir, kinds={'lev'})

    def Write_SAR_Data(self, dhat, r=None, outdir=None):
        '''Write observed/modeled/residual InSAR data and the ramp effect.'''
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)
        written = self._write_blocks(dhat, outdir, kinds={'sar'})
        sv = self._solution_view()
        ramp = None if sv is None else sv.ramp_sar
        written += self._write_sar_ramp(ramp, outdir, sv)
        return written

    # ------------------------------------------------------------------ #
    # Slip in GMT format
    # ------------------------------------------------------------------ #
    def _vertices(self, group):
        '''Return (nf, ncorner, 3) [lon, lat, depth] patch corners.'''
        if group.kind == 'triangle':
            v = np.asarray(group.flt.vertex_llh)[:, [1, 0, 2]]
            return v[np.asarray(group.flt.element) - 1]

        felem = group.flt.FaultGeom2AllVertex()
        # ring order: lower-left -> lower-right -> upper-right -> upper-left
        corners = [('lb_lon', 'lb_lat', 'lb_dep'),
                   ('rb_lon', 'rb_lat', 'rb_dep'),
                   ('rt_lon', 'rt_lat', 'rt_dep'),
                   ('lt_lon', 'lt_lat', 'lt_dep')]
        return np.stack(
            [np.column_stack([np.asarray(felem[c[0]]),
                              np.asarray(felem[c[1]]),
                              np.asarray(felem[c[2]])])
             for c in corners], axis=1)

    def _rectangle_vertices(self):
        return self._vertices(FaultGroup(self.flt, 'fault01', 'rectangle'))

    def _triangle_vertices(self):
        return self._vertices(FaultGroup(self.flt, 'fault01', 'triangle'))

    def WriteGMTSlip(self, slip=None, outdir=None):
        '''Write ss/ds/ts slip in GMT multi-segment format (last model).'''
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)

        if slip is None:
            sv = self._solution_view()
            if sv is None or sv.slip is None:
                raise ValueError('No slip given and the solution has none.')
            slip = sv.slip
        row = np.asarray(self._slip_rows(slip)[-1]).reshape(-1)

        groups = self._fault_groups()
        vertices, values, off = [], [], 0
        for group in groups:
            nf = group.flt.nf
            vertices.append(self._vertices(group))
            values.append(row[off:off + 3 * nf].reshape(nf, 3))
            off += 3 * nf

        vertices = np.vstack(vertices)
        slip3    = np.vstack(values)
        slip_mag = np.sqrt(slip3[:, 0] ** 2 + slip3[:, 1] ** 2)

        if any(g.kind == 'triangle' for g in groups):
            header = '# Triangle fault slip model'
            fmt    = '%10.3f %10.3f %10.3f'
        else:
            header = '# Rectangle fault slip model'
            fmt    = '%10.4f %10.4f %10.4f'

        written = []
        for fname, value in (('ss_slip_llh.gmtlin', slip3[:, 0]),
                             ('ds_slip_llh.gmtlin', slip3[:, 1]),
                             ('ts_slip_llh.gmtlin', slip_mag)):
            path = os.path.join(outdir, fname)
            ou.write_gmt_slip(path, vertices, value, header=header, fmt=fmt)
            written.append(path)

        logging.info('Fault slip distribution is output in geography system.')
        return written

    def write_rectangle_slip(self, filename, value):
        '''Write one rectangle slip field in GMT format.'''
        ou.write_gmt_slip(filename, self._rectangle_vertices(), value,
                          header='# Rectangle fault slip model',
                          fmt='%10.4f %10.4f %10.4f')
        logging.info('Rectangular fault slip distribution is ouput in GMT format.')

    def write_triangle_slip(self, filename, value):
        '''Write one triangle slip field in GMT format.'''
        ou.write_gmt_slip(filename, self._triangle_vertices(), value,
                          header='# Triangle fault slip model',
                          fmt='%10.3f %10.3f %10.3f')
        logging.info('Triangular fault slip distribution is ouput in GMT format.')

    # ------------------------------------------------------------------ #
    # Slip models
    # ------------------------------------------------------------------ #
    @staticmethod
    def _fault_table(group):
        if group.kind == 'triangle':
            return np.asarray(group.flt.element)
        return np.asarray(group.flt.geom_grid)

    @staticmethod
    def _write_model(path, table, slip2d, kind):
        if kind == 'triangle':
            ou.write_tri_model(path, table, slip2d)
        else:
            ou.write_faultgeom(path, table, slip2d)
        return path

    def WriteSlipModel(self, slip=None, k=None, outdir=None,
                       prefix='slipmodel'):
        '''
        Write slip model(s) in FaultGeom or triangular format.

        Input:
            slip   = (3*nf,) or (nsmooth, 3*nf); default self.sol.slip
            k      = write only the k-th smoothing model (0-based)
            outdir = output directory
            prefix = file name prefix
        '''
        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)

        if slip is None:
            sv = self._solution_view()
            if sv is None or sv.slip is None:
                raise ValueError('No slip given and the solution has none.')
            slip = sv.slip

        single = np.asarray(slip).ndim == 1
        rows = self._slip_rows(slip)
        if single:
            indices = [0 if k is None else k]
        elif k is None:
            indices = list(range(rows.shape[0]))
        else:
            rows = rows[k:k + 1]
            indices = [k]

        groups = self._fault_groups()
        if not groups:
            raise ValueError('No fault is available for slip output.')

        kinds = {g.kind for g in groups}
        if len(kinds) != 1:
            raise ValueError('Cannot mix rectangle and triangle faults.')
        kind = kinds.pop()

        nf_total = sum(int(g.flt.nf) for g in groups)
        nf3 = 3 * nf_total
        ext = '.tri' if kind == 'triangle' else '.faultgeom'
        if kind == 'rectangle':
            table = np.vstack([self._fault_table(g) for g in groups])
        else:
            table = self._fault_table(groups[0])

        written = []
        for row, idx in zip(rows, indices):
            row = np.asarray(row).reshape(-1)
            if row.size < nf3:
                raise ValueError('slip size %d does not match %d patches'
                                 % (row.size, nf_total))
            row  = row[:nf3]
            stem = '%s_%03d' % (prefix, idx + 1)

            # all faults combined
            fname = os.path.join(outdir, stem + ext)
            self._write_model(fname, table, row.reshape(nf_total, 3), kind)
            written.append(fname)

            if len(groups) == 1:
                beta = groups[0].rake_beta
                if beta:
                    rname = os.path.join(outdir, stem + '_rot' + ext)
                    self._write_model(
                        rname, table,
                        ou.rotate_slip(row.reshape(nf_total, 3), beta), kind)
                    written.append(rname)
                continue

            # per-fault files (multi-fault case)
            off = 0
            for group in groups:
                nf = int(group.flt.nf)
                islip = row[off:off + 3 * nf].reshape(nf, 3)
                fname = os.path.join(
                    outdir, '%s_%s%s' % (stem, group.name, ext))
                self._write_model(fname, self._fault_table(group), islip, kind)
                written.append(fname)

                if group.rake_beta:
                    rname = fname[:-len(ext)] + '_rot' + ext
                    self._write_model(rname, self._fault_table(group),
                                      ou.rotate_slip(islip, group.rake_beta),
                                      kind)
                    written.append(rname)
                off += 3 * nf

        return written

    # ------------------------------------------------------------------ #
    # Statistics
    # ------------------------------------------------------------------ #
    def WriteSummary(self, filename='summary', outdir=None):
        '''Write the data / fault / inversion summary.'''
        sv = self._solution_view()
        if sv is None:
            raise ValueError('No solution is available for summary output.')

        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)
        blocks = self._obs_blocks()
        gps_blocks = [b for b in blocks if b.kind == 'gps']
        ndim = int(gps_blocks[0].ndim) if gps_blocks else int(self.data.ndim)
        n_gps = sum(len(b.llh) for b in gps_blocks)
        n_lev = sum(len(b.llh) for b in blocks if b.kind == 'lev')
        n_sar = sum(len(b.llh) for b in blocks if b.kind == 'sar')

        groups = self._fault_groups()
        nstrike, ndip = self._summary_counts(groups)

        path = self._path(filename, outdir)
        ou.write_summary(
            path,
            ndim       = ndim,
            n_gps      = n_gps,
            n_sar      = n_sar,
            n_lev      = n_lev,
            flttype    = self._flttype(),
            nstrike    = nstrike,
            ndip       = ndip,
            smo_facts  = self._as_array(sv.smo_facts),
            misfit_gps = self._as_array(sv.misfit_gps),
            misfit_sar = self._as_array(sv.misfit_sar),
            n_d_gps    = np.asarray(self.data.d_gps).size,
            n_d_sar    = np.asarray(self.data.d_sar).size,
            moment     = self._as_array(sv.moment),
        )
        logging.info('Summary is written to %s.' % path)
        return path

    def WriteLCurve(self, filename='ruff_misfit', outdir=None):
        '''Write the roughness / misfit / smoothing-factor trade-off table.'''
        sv = self._solution_view()
        if sv is None:
            raise ValueError('No solution is available for L-curve output.')

        outdir = outdir or self.outdir
        os.makedirs(outdir, exist_ok=True)
        path = self._path(filename, outdir)
        smo_facts = self._as_array(sv.smo_facts)
        ok = ou.write_lcurve(path,
                             self._as_array(sv.ruff),
                             self._as_array(sv.misfit),
                             smo_facts,
                             self._as_array(sv.moment))
        if not ok:
            logging.warning('Only %d smoothing factors.' % len(smo_facts))
            return None
        logging.info('Roughness vs misfit is printed into %s.' % path)
        return path

    # ------------------------------------------------------------------ #
    # One smoothing factor
    # ------------------------------------------------------------------ #
    def WriteSmoothing(self, k, result, outdir=None):
        '''
        Write the k-th smoothing result.

        Input:
            k      = smoothing index (0-based)
            result = dict from MultiFaultInversion.solve(), containing at
                     least 'slip' and 'dhat'
        '''
        outdir = outdir or self.outdir
        self.WriteSlipModel(result['slip'], k=k, outdir=outdir)

        # 观测数据不随平滑因子变化：只在 k == 0 时写一份到 outdir，
        # 各 sm%03d 子目录只保留模型与残差，避免重复写同一份观测。
        if k == 0:
            self.WriteData(result['dhat'], outdir=outdir, parts=('obs',))

        subdir = os.path.join(outdir, 'sm%03d' % k)
        os.makedirs(subdir, exist_ok=True)
        self.WriteData(result['dhat'], outdir=subdir,
                       ramp=result.get('ramp_sar'), parts=('mod', 'res'))
        return subdir

    # ------------------------------------------------------------------ #
    # Main entry point / archiving
    # ------------------------------------------------------------------ #
    def OutputSolution(self):
        '''Write HDF5, summary, models, GMT slip, data fit and archive.'''
        os.makedirs(self.outdir, exist_ok=True)
        self.WriteHDF5()

        if self._is_forward():
            self.WriteData()
            self.archive_outfile()
            return

        self.WriteSummary()
        self.WriteSlipModel()
        self.WriteGMTSlip()
        self.WriteLCurve()
        self.WriteData()
        self.archive_outfile()

    def archive_outfile(self):
        '''Move output files into the archive directory (self.archdir).'''
        if not self.archdir:
            logging.error('The input archive directory is null')
            return []
        return ou.archive_files(self.archdir,
                                ou.DEFAULT_ARCHIVE_PATTERNS,
                                src_dir=self.outdir)
