#!/usr/bin/env python
"""
OutputUtils -- reusable, object-free writers for inversion / forward results.

Every function in this module takes plain numpy arrays, paths, scalars or
dictionaries.  Nothing here knows about Fault/GeoData/solution classes; the
higher-level ``Export.Output`` facade extracts those attributes and delegates
to these helpers.

File formats and headers are kept identical to the historical ``Output``
implementation so downstream GMT/plotting scripts keep working.
"""

import os
import glob
import shutil
import logging
import numpy as np
import h5py
import pandas as pd


# Default files produced by the inversion output path that should be archived.
DEFAULT_ARCHIVE_PATTERNS = (
    'gps*.gmtvec',
    'sar*.xyz',
    'lev*.xyz',
    'ruff_misfit',
    'slipmodel*.tri',
    'slipmodel*.faultgeom',
    # Only the solution file is archived; Green-function caches and other
    # intermediate HDF5 files must stay next to the run that produced them.
    'solutions.h5',
    '*.gmtlin',
    'summary',
)

_FAULTGEOM_HEADER = (
    " width_km   depth_km    dip_deg   lat0_deg   lon0_deg   lat1_deg   "
    "lon1_deg  slp_strk_mm  slp_ddip_mm  slp_tens_mm"
)


# --------------------------------------------------------------------------- #
# HDF5 solution file
# --------------------------------------------------------------------------- #
def write_solution_h5(filename, obs, fault, solution):
    '''
    Write observations, fault geometry and solution to an HDF5 file.

    Input:
        filename = output file name
        obs      = dict with llh_gps, llh_lev, llh_sar, unit, d_gps, d_lev,
                   d_sar and ndim; optional n_sar holds the point count of
                   each InSAR file (track), used to split plots per file
        fault    = dict with origin, flttype, plus rectangle keys
                   (dis_geom, geom) or triangle keys
                   (vertex_enu, vertex_llh, element)
        solution = dict with slip, r, misfit and dhat
    '''
    with h5py.File(filename, 'w') as h5:
        # observation
        grp = h5.create_group('obs')
        for key in ('llh_gps', 'llh_lev', 'llh_sar', 'unit',
                    'd_gps', 'd_lev', 'd_sar'):
            grp.create_dataset(key, data=obs[key], compression='gzip')
        grp['ndim'] = obs['ndim']
        # 逐 InSAR 文件的点数，绘图时据此拆分 track
        if obs.get('n_sar') is not None:
            grp.create_dataset('n_sar',
                               data=np.asarray(obs['n_sar'], dtype=int))
        # 估计的 ramp 对 InSAR 观测的贡献及逐 track 标记
        if obs.get('d_sar_ramp') is not None:
            grp.create_dataset('d_sar_ramp',
                               data=np.asarray(obs['d_sar_ramp']),
                               compression='gzip')
        if obs.get('sar_ramp') is not None:
            grp.create_dataset('sar_ramp',
                               data=np.asarray(obs['sar_ramp'], dtype=bool))

        # fault geometry
        grp = h5.create_group('flt')
        grp.create_dataset('origin', data=fault['origin'], compression='gzip')
        flttype = fault['flttype']
        grp.create_dataset('flttype', data=flttype)
        if flttype == 'rectangle':
            for key in ('dis_geom', 'geom'):
                grp.create_dataset(key, data=fault[key], compression='gzip')
        elif flttype == 'triangle':
            for key in ('vertex_enu', 'vertex_llh', 'element'):
                grp.create_dataset(key, data=fault[key], compression='gzip')
        else:
            raise ValueError('unsupported flttype: %s' % flttype)

        # solution
        grp = h5.create_group('sol')
        for key in ('slip', 'r', 'misfit', 'dhat'):
            grp.create_dataset(key, data=solution[key], compression='gzip')
        # optional metadata: smoothing factors, final model index, ramp, etc.
        for key in ('smo_facts', 'ruff', 'chosen_index',
                    'ramp_sar', 'sar_ramp_flag'):
            value = solution.get(key)
            if value is None:
                continue
            value = np.asarray(value)
            if value.size == 0:
                continue
            grp.create_dataset(key, data=value, compression='gzip')


# --------------------------------------------------------------------------- #
# Observed / modeled / residual data
# --------------------------------------------------------------------------- #
def write_gps_data(llh, d_obs, d_mod, sigma=None, station=None,
                   outdir='.', prefix='gps', parts=None):
    '''
    Write observed, modeled and residual GPS data as GMT vectors.

    Input:
        llh      = (n, 2) array of [lat, lon]
        d_obs    = observed components, (n, ndim)
        d_mod    = modeled components, (n, ndim)
        sigma    = optional formal uncertainties, (n, ndim)
        station  = optional station names, (n,)
        outdir   = output directory
        prefix   = file name prefix
        parts    = subset of ('obs', 'mod', 'res'), plus the *_up variants
                   when ndim == 3; default all

    Output files (n names as above):
        {prefix}_obs.gmtvec, {prefix}_mod.gmtvec, {prefix}_res.gmtvec
        and *_up variants when ndim == 3
    '''
    llh = np.asarray(llh)
    if llh.size == 0:
        return []

    nsta = len(llh)
    d_obs = np.asarray(d_obs).reshape(nsta, -1)
    d_mod = np.asarray(d_mod).reshape(nsta, -1)
    ndim = d_obs.shape[1]

    if sigma is None:
        sigma = np.zeros_like(d_obs)
    else:
        sigma = np.asarray(sigma).reshape(nsta, ndim)

    if station is not None:
        station = np.asarray(station).ravel()
        if station.size != nsta:
            station = None

    if ndim == 2:
        output_cfg = [
            ('obs', d_obs[:, 0], d_obs[:, 1], sigma[:, 0], sigma[:, 1]),
            ('mod', d_mod[:, 0], d_mod[:, 1], 0, 0),
            ('res', d_obs[:, 0] - d_mod[:, 0],
                    d_obs[:, 1] - d_mod[:, 1], 0, 0),
        ]
    else:
        output_cfg = [
            ('obs', d_obs[:, 0], d_obs[:, 1], sigma[:, 0], sigma[:, 1]),
            ('mod', d_mod[:, 0], d_mod[:, 1], 0, 0),
            ('res', d_obs[:, 0] - d_mod[:, 0],
                    d_obs[:, 1] - d_mod[:, 1], 0, 0),
            ('obs_up', 0, d_obs[:, 2], 0, sigma[:, 2]),
            ('mod_up', 0, d_mod[:, 2], 0, 0),
            ('res_up', 0, d_obs[:, 2] - d_mod[:, 2], 0, 0),
        ]

    site = "    Site\n" if station is not None else "\n"
    header_2d = ("#      Lon             Lat        East(mm)       North(mm)  "
                 "SigE(mm)  SigN(mm)  Corr" + site)
    header_up = ("#      Lon             Lat        NULL(mm)          Up(mm)  "
                 "NULL(mm)  SigU(mm)  Corr" + site)

    parts = None if parts is None else {str(p) for p in parts}

    files = []
    for name, east, north, sig_e, sig_n in output_cfg:
        # 'obs_up'/'mod_up'/'res_up' belong to the same base part as 'obs'
        # etc., so selecting ('obs',) also writes the vertical component.
        base = name[:-3] if name.endswith('_up') else name
        if parts is not None and base not in parts:
            continue
        columns = {
            'lon': llh[:, 1],
            'lat': llh[:, 0],
            'east': east,
            'north': north,
            'sig_e': sig_e,
            'sig_n': sig_n,
            'corr': 0,
        }
        if station is not None:
            columns['station'] = station

        fname = os.path.join(outdir, '%s_%s.gmtvec' % (prefix, name))
        header = header_up if name.endswith('up') else header_2d
        # Write the header first and append the body: one pass instead of
        # writing the body, reading it back and rewriting header + body.
        with open(fname, 'w') as fid:
            fid.write(header)
        pd.DataFrame(columns).to_csv(fname, sep='\t', index=False,
                                     float_format='%10.4f', header=False,
                                     mode='a')
        files.append(fname)

    return files


def write_lev_data(llh, d_obs, d_mod, outdir='.', prefix='lev', parts=None):
    '''
    Write modeled and residual leveling data.

    Input:
        llh    = (n, 2) array of [lat, lon]
        d_obs  = observed values, (n,)
        d_mod  = modeled values, (n,)

    Output files:
        {prefix}_obs.xyz, {prefix}_mod.xyz, {prefix}_res.xyz with columns lon lat value
    '''
    llh = np.asarray(llh)
    if llh.size == 0:
        return []

    d_obs = np.asarray(d_obs).reshape(-1)
    d_mod = np.asarray(d_mod).reshape(-1)
    lon, lat = llh[:, 1], llh[:, 0]

    parts = None if parts is None else {str(p) for p in parts}

    files = []
    for name, value in (('obs', d_obs), ('mod', d_mod), ('res', d_obs - d_mod)):
        if parts is not None and name not in parts:
            continue
        fname = os.path.join(outdir, '%s_%s.xyz' % (prefix, name))
        np.savetxt(fname, np.column_stack((lon, lat, value)),
                   fmt='%10.4f\t%10.4f\t%10.2f')
        files.append(fname)
    return files


def write_sar_data(llh, d_obs, d_mod, n_sar=None, outdir='.', prefix='sar',
                   parts=None):
    '''
    Write modeled and residual InSAR data, one file per track.

    Input:
        llh    = (n, 2) array of [lat, lon]
        d_obs  = observed LOS values, (n,)
        d_mod  = modeled LOS values, (n,)  (raw model, e.g. sol.dhat)
        n_sar  = optional number of points in each SAR track

    Output files (one triplet per non-empty track):
        {prefix}_obs_{track}.xyz, {prefix}_mod_{track}.xyz,
        {prefix}_res_{track}.xyz
    '''
    llh = np.asarray(llh)
    if llh.size == 0:
        return []

    d_obs = np.asarray(d_obs).reshape(-1)
    d_mod = np.asarray(d_mod).reshape(-1)
    n = len(llh)

    if n_sar is None:
        counts = [n]
    elif np.isscalar(n_sar):
        counts = [int(n_sar)]
    else:
        counts = [int(v) for v in np.asarray(n_sar).reshape(-1)]

    lon, lat = llh[:, 1], llh[:, 0]
    header = 'Lon  Lat  LOS(mm)'
    parts = None if parts is None else {str(p) for p in parts}

    files = []
    start = 0
    for i, count in enumerate(counts):
        end = start + count
        # 空 track 不写文件，保证观测/模拟文件数一致
        if count > 0:
            for name, value in (('obs', d_obs[start:end]),
                                ('mod', d_mod[start:end]),
                                ('res', d_obs[start:end] - d_mod[start:end])):
                if parts is not None and name not in parts:
                    continue
                fname = os.path.join(outdir, '%s_%s_%d.xyz' % (prefix, name, i + 1))
                np.savetxt(fname,
                           np.column_stack((lon[start:end], lat[start:end], value)),
                           header=header, fmt='%10.4f\t%10.4f\t%10.2f')
                files.append(fname)
        start = end

    return files


def write_sar_ramp_data(llh, d_ramp, n_sar=None, outdir='.',
                        prefix='sar_ramp', flags=None):
    '''
    Write the estimated ramp contribution, one file per InSAR track.

    Input:
        llh    = (n, 2) array of [lat, lon]
        d_ramp = estimated ramp value per point, (n,)
        n_sar  = optional number of points in each SAR track
        flags  = optional per-track booleans; only True tracks are written

    Output files:
        {prefix}_{track}.xyz with columns lon lat ramp
    '''
    llh = np.asarray(llh)
    if llh.size == 0:
        return []

    d_ramp = np.asarray(d_ramp).reshape(-1)
    if n_sar is None:
        counts = [len(llh)]
    elif np.isscalar(n_sar):
        counts = [int(n_sar)]
    else:
        counts = [int(v) for v in np.asarray(n_sar).reshape(-1)]

    if flags is None:
        flags = [True] * len(counts)
    else:
        flags = [bool(v) for v in np.asarray(flags).reshape(-1)]

    lon, lat = llh[:, 1], llh[:, 0]
    header = 'Lon  Lat  Ramp(mm)'

    files = []
    start = 0
    for i, count in enumerate(counts):
        end = start + count
        # 只输出真正估计了 ramp 的非空 track
        if count > 0 and i < len(flags) and flags[i]:
            fname = os.path.join(outdir, '%s_%d.xyz' % (prefix, i + 1))
            np.savetxt(fname,
                       np.column_stack((lon[start:end], lat[start:end],
                                        d_ramp[start:end])),
                       header=header, fmt='%10.4f\t%10.4f\t%10.2f')
            files.append(fname)
        start = end

    return files


# --------------------------------------------------------------------------- #
# Slip models
# --------------------------------------------------------------------------- #
def write_faultgeom(filename, geom, slip, fmt=None, header=None):
    '''
    Write a rectangular slip model in FaultGeom format.

    Input:
        geom   = (nf, >=7) array whose first seven columns are written
        slip   = (nf, 3) array (or the flattened equivalent)
    '''
    geom = np.asarray(geom)
    slip = np.asarray(slip).reshape(len(geom), 3)
    out = np.hstack((geom[:, 0:7], slip))
    if fmt is None:
        fmt = 7 * '%11.4f' + 3 * '%13.3f'
    if header is None:
        header = _FAULTGEOM_HEADER
    np.savetxt(filename, out, fmt=fmt, header=header)
    return filename


def write_tri_model(filename, element, slip, fmt=None):
    '''
    Write a triangular slip model: node indices (1-based) plus slip.
    '''
    element = np.asarray(element)
    slip = np.asarray(slip).reshape(len(element), 3)
    out = np.hstack((element, slip))
    if fmt is None:
        fmt = '%5d %5d %5d %15.4f %15.4f %15.4f'
    np.savetxt(filename, out, fmt=fmt)
    return filename


def rotate_slip(slip, rake_beta):
    '''
    Rotate strike/dip slip back from the constrained rake frame.

    Input:
        slip      = (..., 3) array
        rake_beta = rotation angle in degrees (0 means no rotation)
    '''
    slip = np.asarray(slip)
    if not rake_beta:
        return slip.copy()

    beta = np.deg2rad(rake_beta)
    cos_b, sin_b = np.cos(beta), np.sin(beta)
    out = slip.copy()
    ss = slip[..., 0].copy()
    ds = slip[..., 1].copy()
    out[..., 0] = cos_b * ss - sin_b * ds
    out[..., 1] = sin_b * ss + cos_b * ds
    return out


def write_slip_models(outdir, table, slip, rake_beta=0.0,
                      prefix='slipmodel', kind='rectangle'):
    '''
    Write one slip-model file per smoothing factor.

    Input:
        outdir    = output directory
        table     = (nf, 7+) geom_grid for 'rectangle',
                    or (nf, 3) element table for 'triangle'
        slip      = (nf*3,) or (nsmooth, nf*3) slip array
        rake_beta = rake rotation in degrees; non-zero also writes *_rot_* files
        kind      = 'rectangle' or 'triangle'
    '''
    table = np.asarray(table)
    slip = np.asarray(slip)
    if slip.ndim == 1:
        slip = slip[None, :]

    nf = len(table)
    if slip.shape[1] != 3 * nf:
        raise ValueError('slip size %d does not match %d patches'
                         % (slip.shape[1], nf))

    ext = '.faultgeom' if kind == 'rectangle' else '.tri'
    writer = write_faultgeom if kind == 'rectangle' else write_tri_model

    files = []
    for i in range(slip.shape[0]):
        islip = slip[i].reshape(nf, 3)
        fname = os.path.join(outdir, '%s_%03d%s' % (prefix, i + 1, ext))
        writer(fname, table, islip)
        files.append(fname)

        if rake_beta:
            rname = os.path.join(outdir, '%s_rot_%03d%s' % (prefix, i + 1, ext))
            writer(rname, table, rotate_slip(islip, rake_beta))
            files.append(rname)

    return files


def write_gmt_slip(filename, vertices, values, header='# Fault slip model',
                   fmt='%10.4f %10.4f %10.4f'):
    '''
    Write slip polygons in GMT multi-segment format.

    Input:
        vertices = (nf, ncorner, 3) array of [lon, lat, depth]
        values   = (nf,) slip value used for the -Z segment header
    '''
    vertices = np.asarray(vertices)
    values = np.asarray(values).reshape(-1)

    lines = []
    for i in range(len(vertices)):
        lines.append('>-Z%.3f' % values[i])
        for vertex in vertices[i]:
            lines.append(fmt % (vertex[0], vertex[1], vertex[2]))

    np.savetxt(filename, lines, fmt='%s', header=header)
    return filename


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #
def write_summary(filename, ndim, n_gps, n_sar, n_lev, flttype,
                  nstrike, ndip, smo_facts, misfit_gps, misfit_sar,
                  n_d_gps, n_d_sar, moment, mode='w'):
    '''
    Write a human-readable inversion summary.

    ``mode='w'`` (default) overwrites any previous summary; pass ``'a'`` to
    append explicitly.
    '''
    with open(filename, mode) as fid:
        fid.write('#*#*#*#*#*#*#*#*#*#Geodetic Data Summary#*#*#*#*#*#*#*#*#*#*\n')
        fid.write('Dimension of GPS                 = %d\n' % ndim)
        fid.write('Number of GPS stations           = %d\n' % n_gps)
        fid.write('Number of InSAR piexl            = %d\n' % n_sar)
        fid.write('Number of Leveling data          = %d\n' % n_lev)
        fid.write('\n\n')

        if flttype in ('rectangle', 'triangle'):
            fid.write('#*#*#*#*#*#*#*#*#*#*##*Fault Summary#*#*#*#*#*#*#*#*#*#*#*#*\n')
            if np.ndim(nstrike) == 0:
                fid.write('Number of patches along strike   = %d\n' % nstrike)
                fid.write('Number of patches along dip      = %d\n' % ndip)
            else:
                for i, (ns, nd) in enumerate(zip(np.ravel(nstrike),
                                                 np.ravel(ndip))):
                    fid.write('Fault %2d: patches along strike  = %d\n'
                              % (i + 1, ns))
                    fid.write('Fault %2d: patches along dip     = %d\n'
                              % (i + 1, nd))
            fid.write('\n\n')

        fid.write('#*#*#*#*#*#*#*#*#*#*#Inversion Summary#*#*#*#*#*#*#*#*#*#*#*\n')
        for i in range(len(smo_facts)):
            fid.write('Smoothing factor                 = %f\n' % smo_facts[i])
            if n_gps > 0:
                fid.write('GPS WRSS/N                       = %f\n'
                          % (misfit_gps[i] / n_d_gps))
            if n_sar > 0:
                fid.write('InSAR WRSS/N                     = %f\n'
                          % (misfit_sar[i] / n_d_sar))
            fid.write('Mo                               = %e\n' % moment[(i, 0)])
            fid.write('Mw                               = %4.2f\n' % moment[(i, 1)])
    return filename


def write_lcurve(filename, ruff, misfit, smo_facts, moment):
    '''
    Write roughness / misfit / smoothing-factor / moment trade-off table.

    Returns False when there are fewer than two smoothing factors.
    '''
    smo_facts = np.asarray(smo_facts)
    if len(smo_facts) < 2:
        return False

    ruff_misfit = np.column_stack((ruff, misfit, smo_facts, moment))
    header = " Ruffness\t Misfit\t Smoothing_factor\t Mo\t Mw"
    np.savetxt(filename, ruff_misfit, header=header,
               fmt='%10.2f\t%10.2f\t%10.4f\t%10.2e\t%10.2f')
    return True


# --------------------------------------------------------------------------- #
# Archiving
# --------------------------------------------------------------------------- #
def archive_files(dest_dir, patterns=DEFAULT_ARCHIVE_PATTERNS, src_dir='.',
                  overwrite=True):
    '''
    Move output files matching ``patterns`` into ``dest_dir``.

    Returns the list of moved source paths.  When ``dest_dir`` is empty or is
    the current directory, nothing is moved.  Existing files are overwritten
    (with a warning) when ``overwrite`` is true; otherwise they are kept and
    the new source file stays in ``src_dir``.
    '''
    if not dest_dir:
        logging.error('The input archive directory is null')
        return []

    os.makedirs(dest_dir, exist_ok=True)

    if os.path.abspath(dest_dir) == os.path.abspath(src_dir):
        logging.info('Archive directory is current directory, skipping file move.')
        return []

    moved = []
    for pattern in patterns:
        for path in glob.glob(os.path.join(src_dir, pattern)):
            target = os.path.join(dest_dir, os.path.basename(path))
            if os.path.exists(target):
                if not overwrite:
                    logging.warning(
                        'Archive target %s already exists; leaving %s in place.',
                        target, path)
                    continue
                logging.warning(
                    'Archive target %s already exists; overwriting it with %s.',
                    target, path)
            shutil.move(path, target)
            moved.append(path)

    logging.info('Archived %d file(s) into %s.' % (len(moved), dest_dir))
    return moved
