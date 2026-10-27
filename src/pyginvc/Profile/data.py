#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
Raw GPS and InSAR readers.

GPSData reads .vel/.gmtvec into PyGINVC's GPSData fields; InSARData
reads point tables, GMT grids and GeoTIFFs into SARData fields.  Neither
class selects a profile, weights or offsets; profile.py does that.

'''
import glob
import os

import numpy as np

from pyginvc.GeoData.GPSData import GPSData as _GPSData
from pyginvc.GeoData.SARData import SARData as _SARData


__all__ = ['GPSData', 'InSARData', 'FORMATS', 'GLOBK_COLS',
           'GMT_COLS', 'read_gpsfile']

GLOBK_COLS = dict(lon=0, lat=1, ve=2, vn=3, se=6, sn=7, rho=8, name=12)
GMT_COLS   = dict(lon=0, lat=1, ve=2, vn=3, se=4, sn=5, rho=6, name=7)

FORMATS = {
    '.vel':    ('globk', GLOBK_COLS),
    '.gmtvec': ('gmt',   GMT_COLS),
}

_TWO_PI = 2.0*np.pi


def _rows(path):
    '''
    Yield (line number, whitespace-split fields) for every non-blank,
    non-comment line.
    '''
    with open(path, encoding='utf-8') as f:
        for lineno, line in enumerate(f, 1):
            s = line.strip()
            if s and not s.startswith('#'):
                yield lineno, s.split()


def _num(tok):
    try:
        return float(tok)
    except ValueError:
        return np.nan


def read_gpsfile(path):
    '''
    Read a GPS velocity field.  The format is chosen by the file suffix.

    Input:
        path = velocity file
    Output:
        dict of arrays: lon, lat, ve, vn, se, sn, rho, name, plus
        'format', 'path' and 'n_dropped' (sites dropped for a non-numeric
        field, e.g. a missing value)
    '''
    suffix = os.path.splitext(path)[1].lower()
    if suffix not in FORMATS:
        raise ValueError('%s: unknown velocity-file suffix %r, supported: %s'
                         % (path, suffix, ' '.join(sorted(FORMATS))))
    fmt, cols = FORMATS[suffix]

    need = max(cols[k] for k in ('lon', 'lat', 've', 'vn', 'se', 'sn', 'name'))
    out  = {k: [] for k in ('lon', 'lat', 've', 'vn', 'se', 'sn', 'rho')}
    names, ndrop = [], 0

    for lineno, f in _rows(path):
        if len(f) <= need:
            raise ValueError('%s:%d: only %d columns, %s needs at least %d'
                             % (path, lineno, len(f), fmt.upper(), need + 1))
        v = {k: _num(f[i]) for k, i in cols.items() if k != 'name'}
        if not all(np.isfinite(v[k]) for k in ('lon', 'lat', 've', 'vn',
                                               'se', 'sn')):
            ndrop += 1
            continue
        if 'rho' not in v or not np.isfinite(v['rho']):
            v['rho'] = 0.0
        if v['se'] <= 0. or v['sn'] <= 0.:
            raise ValueError('%s:%d: sigma_e=%.4g sigma_n=%.4g, must be '
                             'positive' % (path, lineno, v['se'], v['sn']))
        if abs(v['rho']) > 1.:
            raise ValueError('%s:%d: correlation rho=%.4g is outside [-1, 1]'
                             % (path, lineno, v['rho']))
        for k in out:
            out[k].append(v[k])
        i = cols['name']
        names.append(f[i] if len(f) > i else 'site%d' % lineno)

    if not out['lon']:
        raise ValueError('%s: no usable sites' % path)

    res = {k: np.asarray(v, dtype=float) for k, v in out.items()}
    res['name']      = np.asarray(names, dtype=object)
    res['format']    = fmt
    res['path']      = path
    res['n_dropped'] = ndrop
    return res


class GPSData(_GPSData):
    '''
    Raw GPS velocities (.vel GLOBK or .gmtvec GMT table).

    Reading only.  Inherits PyGINVC's GPSData: `llh_gps` ([lat, lon]),
    `d_gps`, `cov_gps` (variances), `station_gps`, `nsta`, `ndim`,
    `W_gps` and `plotGPSVector`.  Profile selection, weights and
    offsets belong to profile.py.
    '''
    kind = 'gps'

    def __init__(self, path):
        _GPSData.__init__(self, path, 'GMT2D')
        self.load()             # parent load() -> our load_gps_data()

    def load_gps_data(self):
        dat = read_gpsfile(self.gpsfile)
        self.nsta        = len(dat['lon'])
        self.ndim        = 2
        self.llh_gps     = np.column_stack((dat['lat'], dat['lon']))
        self.d_gps       = np.column_stack((dat['ve'], dat['vn'])).ravel()
        self.cov_gps     = np.column_stack((dat['se']**2,
                                            dat['sn']**2)).ravel()
        self.station_gps = dat['name']
        self.rho         = dat['rho']
        self.format      = dat['format']
        self.n_dropped   = dat['n_dropped']

    @property
    def lon(self):
        return self.llh_gps[:, 1]

    @property
    def lat(self):
        return self.llh_gps[:, 0]

    @property
    def ve(self):
        return self.d_gps[0::2]

    @property
    def vn(self):
        return self.d_gps[1::2]

    @property
    def se(self):
        return np.sqrt(self.cov_gps[0::2])

    @property
    def sn(self):
        return np.sqrt(self.cov_gps[1::2])

    @property
    def name(self):
        return self.station_gps


def _expand_paths(paths):
    '''
    A string (possibly a glob), a list of strings or a list of globs.
    '''
    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]
    out = []
    for p in paths:
        p = os.fspath(p)
        out.extend(sorted(glob.glob(p)) if any(c in p for c in '*?[')
                   else [p])
    if not out:
        raise ValueError('no InSAR files matched %r' % (paths,))
    for p in out:
        if not os.path.isfile(p):
            raise ValueError('InSAR file does not exist: %s' % p)
    return out


def _insar_format(paths, fmt):
    suffixes = {os.path.splitext(p)[1].lower() for p in paths}
    if fmt is None:
        if suffixes <= {'.tif', '.tiff'}:
            return 'tif'
        if suffixes <= {'.grd', '.nc'}:
            return 'grd'
        return 'point'
    fmt = str(fmt).lower()
    if fmt not in ('point', 'grd', 'tif'):
        raise ValueError('format must be point, grd or tif, got %r' % fmt)
    return fmt


def _unit_from_angles(inc, head, shape):
    '''
    LOS unit vector (ground -> satellite), right-looking radar:

        u_E =  sin(i) cos(h)
        u_N = -sin(i) sin(h)
        u_U =  cos(i)

    `head` is the satellite heading, clockwise from north.
    '''
    i = np.broadcast_to(np.radians(inc), shape)
    h = np.broadcast_to(np.radians(head), shape)
    return np.stack((np.sin(i)*np.cos(h),
                     -np.sin(i)*np.sin(h),
                     np.cos(i)), axis=-1)


def _unit_ok(unit):
    n = np.linalg.norm(unit, axis=1)
    if not np.all(np.isfinite(n)) or np.any(np.abs(n - 1.) > 0.02):
        raise ValueError('LOS unit vectors are not unit length (%.3f .. %.3f)'
                         % (np.nanmin(n), np.nanmax(n)))
    return unit/n[:, None]


def _read_grd(path):
    '''
    Read a GMT/NetCDF grid: x, y and z (z shaped (len(y), len(x))).
    '''
    def get(nc, *names):
        for n in names:
            if n in nc.variables:
                return np.array(nc.variables[n][:], dtype=float)
        raise KeyError('none of %s' % (', '.join(names),))

    try:
        from scipy.io import netcdf_file
        with netcdf_file(path, 'r', mmap=False) as nc:
            x = get(nc, 'x', 'lon', 'longitude')
            y = get(nc, 'y', 'lat', 'latitude')
            z = get(nc, 'z', 'los', 'band1')
    except Exception as e:
        try:
            from netCDF4 import Dataset
        except ImportError:
            raise ImportError(
                '%s: cannot read as NetCDF-3 (%s).  If it is NetCDF-4, '
                'install netCDF4: pip install pysliprate[grid]'
                % (path, e)) from e
        with Dataset(path) as nc:
            x = get(nc, 'x', 'lon', 'longitude')
            y = get(nc, 'y', 'lat', 'latitude')
            z = get(nc, 'z', 'los', 'band1')

    z = np.squeeze(z)
    if z.ndim != 2:
        raise ValueError('%s: z must be a 2-D grid, got shape %r'
                         % (path, z.shape))
    if z.shape != (len(y), len(x)) and z.shape == (len(x), len(y)):
        z = z.T
    if z.shape != (len(y), len(x)):
        raise ValueError('%s: z shape %r does not match x(%d) y(%d)'
                         % (path, z.shape, len(x), len(y)))
    if x[0] > x[-1]:
        x, z = x[::-1], z[:, ::-1]
    if y[0] > y[-1]:
        y, z = y[::-1], z[::-1, :]
    return x, y, z


def _crop(axis, lo, hi):
    i0 = int(np.searchsorted(axis, lo, side='left'))
    i1 = int(np.searchsorted(axis, hi, side='right'))
    return i0, i1


class InSARData(_SARData):
    '''
    Raw InSAR LOS observations: point tables, GMT grids or GeoTIFFs.

    Reading only.  Inherits PyGINVC's SARData (`llh_sar`, `d_sar`,
    `unit`, `cov_sar`).  Selection, binning, weights and offsets belong
    to profile.py.
    '''

    kind = 'insar'

    def __init__(self, paths, format=None, incidence=None, heading=None,
                 unit_bands=None, sigma=3.0, scale=1.0, overlap='first'):
        '''
        Input:
            paths       = file, list of files or glob, one track
            format      = point | grd | tif, default from the suffix
            incidence   = look angle (deg), scalar, per-file list or grid path
            heading     = satellite heading (deg), same forms
            unit_bands  = [e, n, u] band numbers if the tif carries them
            sigma       = LOS uncertainty (scalar or per point); ignored for
                          point files that carry their own sigma column
            scale       = multiply the LOS values (e.g. m -> mm)
            overlap     = first | mean, for overlapping tif tiles
        '''
        paths = _expand_paths(paths)

        # PyGINVC's __init__ does not call super(); sigma/lamda stay None
        # because build_spatial_cov only needs ones.
        _SARData.__init__(self, paths[0], scale=float(scale),
                          sigma=None, lamda=None)

        self.paths = paths
        self.format = _insar_format(self.paths, format)
        self.incidence = incidence
        self.heading = heading
        self.unit_bands = unit_bands
        self.sigma = sigma
        self.scale = float(scale)
        self.overlap = str(overlap).lower()
        if self.overlap not in ('first', 'mean'):
            raise ValueError('overlap must be first or mean, got %r'
                             % (overlap,))
        if self.format != 'point':
            s = np.asarray(sigma, dtype=float)
            if s.ndim != 0 or not np.isfinite(s) or s <= 0.:
                raise ValueError('%s: gridded/tif InSAR needs a positive '
                                 'scalar sigma' % self._name)
            self.sigma = float(s)
        if self.format == 'point':
            self.load()         # parent load() -> our load_sar_data()

    @property
    def _name(self):
        return os.path.basename(self.paths[0])

    def load_sar_data(self, scale=None):
        '''
        Read point LOS files into the PyGINVC fields.  `scale` defaults
        to the constructor's value (PyGINVC silently dropped it).

        The columns follow PyGINVC's GetSARdata (`lon lat LOS u_e u_n u_u`)
        with an optional seventh sigma column; parsing here also fixes
        GetSARdata's failure on a one-row file.
        '''
        scale = self.scale if scale is None else float(scale)
        llh, val, unit, sig = [], [], [], []
        have = []
        for p in self.paths:
            dat = np.genfromtxt(p, comments='#')
            if dat.ndim == 1:
                dat = dat.reshape(1, -1)
            if dat.shape[1] < 6:
                raise ValueError('%s: point InSAR needs lon lat los ue un uu '
                                 '[sigma], got %d columns'
                                 % (p, dat.shape[1]))
            llh.append(dat[:, [1, 0]])          # lat, lon
            val.append(dat[:, 2])
            unit.append(dat[:, 3:6])
            has = dat.shape[1] >= 7
            have.append(has)
            sig.append(dat[:, 6] if has else None)

        if any(have) and not all(have):
            raise ValueError('%s: either every point file has a sigma column '
                             'or none does' % self._name)
        sigma = np.concatenate(sig) if all(have) else self.sigma
        if sigma is None:
            raise ValueError('%s: no sigma column and no sigma set'
                             % self._name)

        self.llh_sar = np.vstack(llh)
        self.d_sar   = np.concatenate(val)*scale
        self.unit    = _unit_ok(np.vstack(unit))
        self.sigma   = np.broadcast_to(np.asarray(sigma, dtype=float),
                                       len(self.d_sar)).copy()
        if np.any(~np.isfinite(self.sigma)) or np.any(self.sigma <= 0.):
            raise ValueError('%s: sigma must be positive and finite'
                             % self._name)
        return self

    def build_spatial_cov(self):
        self.cov_sar = np.ones(len(self.llh_sar))

    def _angle(self, name, i, nfiles):
        a = getattr(self, name)
        if a is None:
            return None
        if isinstance(a, (str, os.PathLike)):
            return os.fspath(a)
        if np.isscalar(a):
            return float(a)
        a = list(a)
        if len(a) == nfiles and all(np.isscalar(v) for v in a):
            return float(a[i])
        if len(a) == 1 and np.isscalar(a[0]):
            return float(a[0])
        return np.asarray(a, dtype=float)

    def _read_grd_box(self, bounds):
        lon0, lat0, lon1, lat1 = bounds
        lons, lats, vals, units, sigs = [], [], [], [], []

        for i, p in enumerate(self.paths):
            x, y, z = _read_grd(p)
            ix0, ix1 = _crop(x, lon0, lon1)
            iy0, iy1 = _crop(y, lat0, lat1)
            xs = x[ix0:ix1]
            ys = y[iy0:iy1]
            zz = z[iy0:iy1, ix0:ix1]

            inc = self._angle('incidence', i, len(self.paths))
            hd  = self._angle('heading', i, len(self.paths))
            if inc is None or hd is None:
                raise ValueError('%s: gridded InSAR needs incidence and '
                                 'heading' % self._name)
            inc = self._angle_grid(inc, x, y, ix0, ix1, iy0, iy1)
            hd  = self._angle_grid(hd, x, y, ix0, ix1, iy0, iy1)
            unit = _unit_from_angles(inc, hd, zz.shape)

            lon2, lat2 = np.meshgrid(xs, ys)
            ok = np.isfinite(zz)
            lons.append(lon2[ok])
            lats.append(lat2[ok])
            vals.append(zz[ok]*self.scale)
            units.append(unit[ok])
            sigs.append(np.full(ok.sum(), self.sigma, dtype=float))
        return (np.concatenate(lons), np.concatenate(lats),
                np.concatenate(vals), _unit_ok(np.vstack(units)),
                np.concatenate(sigs))

    def _angle_grid(self, a, x, y, ix0, ix1, iy0, iy1):
        if np.isscalar(a):
            return float(a)
        if isinstance(a, (str, os.PathLike)):
            ax, ay, az = _read_grd(os.fspath(a))
            if len(ax) != len(x) or len(ay) != len(y):
                raise ValueError('%s: angle grid shape differs from the LOS '
                                 'grid' % os.fspath(a))
            az = np.squeeze(az)
            if az.shape != (len(y), len(x)):
                az = az.T
            return az[iy0:iy1, ix0:ix1]
        a = np.asarray(a, dtype=float)
        if a.shape == (len(y), len(x)):
            return a[iy0:iy1, ix0:ix1]
        return a

    def _read_tif_box(self, bounds):
        try:
            import rasterio
            from rasterio.warp import transform_bounds, transform as rio_transform
        except ImportError:
            raise ImportError('GeoTIFF support needs rasterio: '
                              'pip install pysliprate[grid]')
        from rasterio.windows import from_bounds, Window

        lon0, lat0, lon1, lat1 = bounds

        lons, lats, vals, units, sigs = [], [], [], [], []
        for i, p in enumerate(self.paths):
            inc = self._angle('incidence', i, len(self.paths))
            hd  = self._angle('heading', i, len(self.paths))
            with rasterio.open(p) as src:
                if src.crs is None:
                    raise ValueError('%s: GeoTIFF has no CRS' % p)
                b = transform_bounds('EPSG:4326', src.crs, lon0, lat0,
                                     lon1, lat1, densify_pts=21)
                win = from_bounds(*b, transform=src.transform)
                win = win.intersection(Window(0, 0, src.width, src.height))
                if win.width < 1 or win.height < 1:
                    continue
                win = Window(win.col_off, win.row_off,
                             int(win.width), int(win.height))
                data = src.read(1, window=win, masked=True).astype(float)
                shape = data.shape
                if self.unit_bands is not None:
                    eb, nb, ub = self.unit_bands
                    unit = np.stack([src.read(k, window=win).astype(float)
                                     for k in (eb, nb, ub)], axis=-1)
                else:
                    if inc is None or hd is None:
                        raise ValueError('%s: GeoTIFF needs incidence and '
                                         'heading or unit_bands' % self._name)
                    inc = self._tif_angle(inc, src, win, shape)
                    hd  = self._tif_angle(hd, src, win, shape)
                    unit = _unit_from_angles(inc, hd, shape)

                rows, cols = np.meshgrid(
                    np.arange(shape[0]) + win.row_off,
                    np.arange(shape[1]) + win.col_off, indexing='ij')
                xs, ys = rasterio.transform.xy(src.transform, rows, cols,
                                               offset='center')
                xs = np.asarray(xs, dtype=float)
                ys = np.asarray(ys, dtype=float)
                if src.crs.to_string() != 'EPSG:4326':
                    lo, la = rio_transform(src.crs, 'EPSG:4326',
                                           xs.ravel(), ys.ravel())
                    lons_ = np.asarray(lo).reshape(shape)
                    lats_ = np.asarray(la).reshape(shape)
                else:
                    lons_, lats_ = xs, ys

                zz = np.ma.getdata(data)
                ok = (~np.ma.getmaskarray(data)) & np.isfinite(zz)
                # Many InSAR GeoTIFFs store the nodata fill (e.g. -9999) as
                # ordinary pixels without setting the nodata tag, so the
                # masked read above leaves them in.  Drop anything far below
                # any physical LOS velocity.
                ok &= zz > -999.0
                lons.append(lons_[ok])
                lats.append(lats_[ok])
                vals.append(zz[ok]*self.scale)
                units.append(unit[ok])
                sigs.append(np.full(ok.sum(), self.sigma, dtype=float))
        if not lons:
            raise ValueError('%s: no GeoTIFF tile intersects the box'
                             % self._name)
        return (np.concatenate(lons), np.concatenate(lats),
                np.concatenate(vals), _unit_ok(np.vstack(units)),
                np.concatenate(sigs))

    def _tif_angle(self, a, src, win, shape):
        if np.isscalar(a):
            return float(a)
        if isinstance(a, (str, os.PathLike)):
            import rasterio
            with rasterio.open(os.fspath(a)) as ang:
                arr = ang.read(1, window=win).astype(float)
            if arr.shape != shape:
                raise ValueError('%s: angle raster shape differs from LOS'
                                 % os.fspath(a))
            return arr
        a = np.asarray(a, dtype=float)
        if a.shape != shape:
            raise ValueError('%s: angle array shape differs from LOS'
                             % self._name)
        return a

    def read_box(self, bounds):
        '''
        Read a gridded/tif dataset inside (lon0, lat0, lon1, lat1) and
        fill the PyGINVC fields.  Point files are already loaded; this
        is a no-op for them.
        '''
        if self.format == 'point':
            return self
        if self.format == 'grd':
            lon, lat, value, unit, sigma = self._read_grd_box(bounds)
        else:
            lon, lat, value, unit, sigma = self._read_tif_box(bounds)

        if self.format == 'tif' and self.overlap == 'first' \
                and len(self.paths) > 1 and len(lon):
            key = (np.rint(lon*1e6).astype(np.int64)*10**9
                   + np.rint(lat*1e6).astype(np.int64))
            _, keep = np.unique(key, return_index=True)
            sel = np.zeros(len(lon), dtype=bool)
            sel[keep] = True
            lon, lat = lon[sel], lat[sel]
            value, unit, sigma = value[sel], unit[sel], sigma[sel]

        self.llh_sar = np.column_stack((lat, lon))
        self.d_sar   = value
        self.unit    = unit
        self.sigma   = sigma
        self.cov_sar = np.ones(len(value))
        return self
