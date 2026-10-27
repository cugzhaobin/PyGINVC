#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
The profile layer: geometry, projection and the flat observation table.

    raw data (data.py)
        │  Profile.cut_gps / cut_insar
        ▼
    ProfileData
        value = cp*v_par + cn*v_norm + cu*v_up + bias
        sigma, along, normal, label, component, ...

GPS is rotated to fault-parallel / fault-normal; InSAR keeps its LOS and
stores the three projection coefficients.  A list of datasets is stacked
by concat(), so GPS + ascending + descending share one parameter vector.

Coordinates are the PyGINVC local projection by default
(`frame='pyginvc'`); `frame='localenu'` reproduces the old geo.py frame.
'''
import numpy as np
from dataclasses import dataclass

from pyginvc.libs import geotools as _gt

__all__ = ['Profile', 'ProfileData', 'concat', 'to_profile', 'from_profile',
           'lonlat_to_enu', 'enu_to_lonlat', 'vel_to_profile', 'rot_sigma',
           'write_profile', 'read_profile']

FRAMES = ('pyginvc', 'localenu')

# WGS84
_A  = 6378.137
_F  = 1./298.257223563
_E2 = 2*_F - _F*_F


def _radii(lat):
    w = 1. - _E2*np.sin(np.radians(lat))**2
    return _A*(1. - _E2)/w**1.5, _A/np.sqrt(w)


def _localenu_to_enu(lon, lat, origin):
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    lon0, lat0 = (float(v) for v in origin)
    _, n0 = _radii(lat0)
    m0, _ = _radii(lat0)
    E = n0*np.cos(np.radians(0.5*(lat0 + lat)))*np.radians(lon - lon0)
    N = m0*np.radians(lat - lat0)
    return E, N


def _localenu_to_lonlat(E, N, origin):
    E = np.asarray(E, dtype=float)
    N = np.asarray(N, dtype=float)
    lon0, lat0 = (float(v) for v in origin)
    _, n0 = _radii(lat0)
    m0, _ = _radii(lat0)
    lat = lat0 + np.degrees(N/m0)
    lon = lon0 + np.degrees(E/(n0*np.cos(np.radians(0.5*(lat0 + lat)))))
    return lon, lat


def lonlat_to_enu(lon, lat, origin, frame='pyginvc'):
    '''
    Local east/north of sites relative to origin.

    Input:
        lon, lat = degree, array or scalar
        origin   = (lon0, lat0)
        frame    = 'pyginvc' (default) or 'localenu'
    Output:
        E, N     = km
    '''
    if frame == 'localenu':
        return _localenu_to_enu(lon, lat, origin)
    if frame != 'pyginvc':
        raise ValueError('frame must be one of %s, got %r'
                         % (', '.join(FRAMES), frame))
    lon = np.atleast_1d(np.asarray(lon, dtype=float))
    lat = np.atleast_1d(np.asarray(lat, dtype=float))
    lon, lat = np.broadcast_arrays(lon, lat)
    shape = lon.shape
    lon0, lat0 = (float(v) for v in origin)
    llh = np.column_stack((lat.ravel(), lon.ravel(), np.zeros(lon.size)))
    out = _gt.llh2local(llh, [lat0, lon0, 0.])
    return out[:, 0].reshape(shape), out[:, 1].reshape(shape)


def enu_to_lonlat(E, N, origin, frame='pyginvc'):
    '''
    Inverse of lonlat_to_enu.
    '''
    if frame == 'localenu':
        return _localenu_to_lonlat(E, N, origin)
    if frame != 'pyginvc':
        raise ValueError('frame must be one of %s, got %r'
                         % (', '.join(FRAMES), frame))
    E = np.atleast_1d(np.asarray(E, dtype=float))
    N = np.atleast_1d(np.asarray(N, dtype=float))
    E, N = np.broadcast_arrays(E, N)
    shape = E.shape
    lon0, lat0 = (float(v) for v in origin)
    xyz = np.column_stack((E.ravel(), N.ravel(), np.zeros(E.size)))
    out = _gt.local2llh(xyz, [lat0, lon0, 0.])
    return out[:, 0].reshape(shape), out[:, 1].reshape(shape)


def to_profile(E, N, strike):
    '''
    East/north -> (along, normal).  Vectorized equivalent of
    geotools.rotate2d([E, N], -strike).
    '''
    phi = np.radians(strike)
    return (np.sin(phi)*E + np.cos(phi)*N,
            np.cos(phi)*E - np.sin(phi)*N)


def from_profile(along, normal, strike):
    '''
    (along, normal) -> east/north.
    '''
    phi = np.radians(strike)
    return (np.sin(phi)*along + np.cos(phi)*normal,
            np.cos(phi)*along - np.sin(phi)*normal)


def vel_to_profile(ve, vn, strike):
    return to_profile(ve, vn, strike)


def rot_sigma(se, sn, rho, strike):
    '''
    Rotate E/N uncertainties to fault-parallel / fault-normal.

        sig_par^2  = s^2 se^2 + 2 s c rho se sn + c^2 sn^2
        sig_norm^2 = c^2 se^2 - 2 c s rho se sn + s^2 sn^2

    The par/norm cross-covariance is dropped; the two components are
    treated as independent (the old profile method did the same).
    '''
    phi = np.radians(strike)
    s, c = np.sin(phi), np.cos(phi)
    se  = np.asarray(se, dtype=float)
    sn  = np.asarray(sn, dtype=float)
    rho = np.asarray(rho, dtype=float)
    sp2 = s*s*se*se + 2*s*c*rho*se*sn + c*c*sn*sn
    sn2 = c*c*se*se - 2*c*s*rho*se*sn + s*s*sn*sn
    return (np.sqrt(np.maximum(sp2, 0.)),
            np.sqrt(np.maximum(sn2, 0.)))


@dataclass
class ProfileData(object):
    '''
    Flat observation table for one profile (one or several datasets).

        value = cp*v_par + cn*v_norm + cu*v_up + bias @ par
    '''
    label:       np.ndarray
    component:   np.ndarray
    site:        np.ndarray
    lon:         np.ndarray
    lat:         np.ndarray
    E:           np.ndarray
    N:           np.ndarray
    along:       np.ndarray
    normal:      np.ndarray
    value:       np.ndarray
    sigma:       np.ndarray
    cp:          np.ndarray
    cn:          np.ndarray
    cu:          np.ndarray
    count:       np.ndarray
    bias:        np.ndarray
    bias_names:  tuple
    bias_bounds: dict
    strike:      float
    weight:      np.ndarray

    def __len__(self):
        return len(self.value)


def concat(data_list):
    '''
    Stack ProfileData objects; bias columns become the union of names.
    '''
    data_list = list(data_list)
    if not data_list:
        raise ValueError('no profile data to concatenate')

    names = []
    bounds = {}
    for d in data_list:
        for k, n in enumerate(d.bias_names):
            if n in bounds and bounds[n] != d.bias_bounds[n]:
                raise ValueError('%s: inconsistent bias bounds' % n)
            bounds[n] = d.bias_bounds[n]
            if n not in names:
                names.append(n)

    strike = data_list[0].strike
    for d in data_list:
        if d.strike != strike:
            raise ValueError('cannot concatenate profiles with different '
                             'strike')

    def cat(name):
        return np.concatenate([getattr(d, name) for d in data_list])

    n = sum(len(d) for d in data_list)
    bias = np.zeros((n, len(names)))
    off = 0
    for d in data_list:
        for k, nm in enumerate(d.bias_names):
            bias[off:off + len(d), names.index(nm)] = d.bias[:, k]
        off += len(d)

    return ProfileData(
        label=cat('label'), component=cat('component'), site=cat('site'),
        lon=cat('lon'), lat=cat('lat'), E=cat('E'), N=cat('N'),
        along=cat('along'), normal=cat('normal'),
        value=cat('value'), sigma=cat('sigma'),
        cp=cat('cp'), cn=cat('cn'), cu=cat('cu'),
        count=cat('count'), bias=bias, bias_names=tuple(names),
        bias_bounds=bounds, strike=strike, weight=cat('weight'))


class Profile(object):
    '''
    Geometry of one cross-fault profile and the cut/projection methods.
    '''

    def __init__(self, origin, strike, along_bounds, normal_bounds,
                 name='profile', frame='pyginvc'):
        origin = tuple(float(v) for v in origin)
        if len(origin) != 2:
            raise ValueError('origin must be (lon, lat), got %r' % (origin,))
        lon0, lat0 = origin
        if not -180. <= lon0 <= 180. or not -90. <= lat0 <= 90.:
            raise ValueError('origin lon=%.4f lat=%.4f out of range '
                             '(is the order (lon, lat)?)' % (lon0, lat0))

        along  = tuple(float(v) for v in along_bounds)
        normal = tuple(float(v) for v in normal_bounds)
        for what, b in (('along', along), ('normal', normal)):
            if len(b) != 2:
                raise ValueError('%s_bounds must be (min, max), got %r'
                                 % (what, b))
            if not np.isfinite(b).all():
                raise ValueError('%s_bounds contains a non-finite value %r'
                                 % (what, b))
            if b[0] >= b[1]:
                raise ValueError('%s_bounds min=%g must be smaller than max=%g'
                                 % (what, b[0], b[1]))
        if frame not in FRAMES:
            raise ValueError('frame must be one of %s, got %r'
                             % (', '.join(FRAMES), frame))

        self.origin        = origin
        self.strike        = float(strike)
        self.along_bounds  = along
        self.normal_bounds = normal
        self.name          = str(name)
        self.frame_name    = str(frame)

    def frame(self, lon, lat):
        '''
        Project sites: returns E, N, along, normal (km or degree).
        '''
        E, N = lonlat_to_enu(lon, lat, self.origin, self.frame_name)
        along, normal = to_profile(E, N, self.strike)
        return E, N, along, normal

    def contains(self, lon, lat):
        '''
        Boolean mask of the sites inside the profile box.
        '''
        _, _, along, normal = self.frame(lon, lat)
        a0, a1 = self.along_bounds
        n0, n1 = self.normal_bounds
        return (along >= a0) & (along <= a1) & \
               (normal >= n0) & (normal <= n1)

    def bounds_lonlat(self):
        '''
        (lon0, lat0, lon1, lat1) of the profile box, for grid readers.
        '''
        lon, lat = self.box()
        return (float(np.min(lon)), float(np.min(lat)),
                float(np.max(lon)), float(np.max(lat)))

    def cut_gps(self, gps, label='GPS', weight=1.0,
                components=('par', 'norm'), offsets=None):
        '''
        Rotate GPS velocities to fault-parallel / fault-normal and return
        a ProfileData (one row per site and component).
        '''
        comps = tuple(components)
        if not comps or set(comps) - {'par', 'norm'}:
            raise ValueError('components must be a subset of (par, norm), '
                             'got %r' % (components,))
        w = float(weight)
        if not np.isfinite(w) or w <= 0.:
            raise ValueError('%s: weight must be positive, got %r'
                             % (label, weight))
        offsets = dict(offsets or {})
        for k, b in offsets.items():
            if k not in comps:
                raise ValueError('%s: offset %r is not in components %r'
                                 % (label, k, comps))
            lo, hi = (float(v) for v in b)
            if lo > hi:
                raise ValueError('%s: offset %r min=%g > max=%g'
                                 % (label, k, lo, hi))
            offsets[k] = (lo, hi)

        E, N, along, normal = self.frame(gps.lon, gps.lat)
        a0, a1 = self.along_bounds
        n0, n1 = self.normal_bounds
        m = ((along >= a0) & (along <= a1) &
             (normal >= n0) & (normal <= n1))
        idx = np.nonzero(m)[0]
        if idx.size == 0:
            raise ValueError('%s: no sites inside the profile box' % label)

        phi = np.radians(self.strike)
        s, c = np.sin(phi), np.cos(phi)
        vpar  = s*gps.ve[idx] + c*gps.vn[idx]
        vnorm = c*gps.ve[idx] - s*gps.vn[idx]
        spar, snorm = rot_sigma(gps.se[idx], gps.sn[idx],
                                gps.rho[idx], self.strike)
        vals = {'par': vpar, 'norm': vnorm}
        sigs = {'par': spar, 'norm': snorm}

        nc = len(comps)
        ns = len(idx)
        n  = nc*ns
        value = np.empty(n)
        sigma = np.empty(n)
        cp = np.zeros(n)
        cn = np.zeros(n)
        cu = np.zeros(n)
        comp_arr = np.empty(n, dtype=object)
        for j, comp in enumerate(comps):
            value[j::nc] = vals[comp]
            sigma[j::nc] = sigs[comp]/w
            if comp == 'par':
                cp[j::nc] = 1.
            else:
                cn[j::nc] = 1.
            comp_arr[j::nc] = comp

        bias_keys = [k for k in comps if k in offsets]
        bias_names = tuple('%s.%s' % (label, k) for k in bias_keys)
        bias = np.zeros((n, len(bias_names)))
        for j, comp in enumerate(comps):
            if comp in bias_keys:
                bias[j::nc, bias_keys.index(comp)] = 1.

        return ProfileData(
            label=np.full(n, label, dtype=object),
            component=comp_arr,
            site=np.repeat(np.asarray(gps.name)[idx], nc),
            lon=np.repeat(np.asarray(gps.lon)[idx], nc),
            lat=np.repeat(np.asarray(gps.lat)[idx], nc),
            E=np.repeat(E[idx], nc), N=np.repeat(N[idx], nc),
            along=np.repeat(along[idx], nc),
            normal=np.repeat(normal[idx], nc),
            value=value, sigma=sigma, cp=cp, cn=cn, cu=cu,
            count=np.ones(n, dtype=int),
            bias=bias, bias_names=bias_names,
            bias_bounds={'%s.%s' % (label, k): offsets[k]
                         for k in bias_keys},
            strike=self.strike, weight=np.full(n, w))

    def cut_insar(self, insar, label='insar', weight=1.0,
                  along_bin=0.0, normal_bin=0.0, offset=None, ramp=None):
        '''
        Cut and (optionally) bin an InSAR dataset, keeping the LOS and the
        projection coefficients cp/cn/cu.

        `along_bin` / `normal_bin` are the along-strike and fault-normal bin
        widths (km); 0 = no binning on that axis (pixels averaged across it).
        '''
        w = float(weight)
        if not np.isfinite(w) or w <= 0.:
            raise ValueError('%s: weight must be positive, got %r'
                             % (label, weight))
        if insar.format != 'point':
            insar.read_box(self.bounds_lonlat())

        llh = np.asarray(insar.llh_sar, dtype=float)
        if len(llh) == 0:
            raise ValueError('%s: no InSAR points' % label)
        lon, lat = llh[:, 1], llh[:, 0]
        E, N, along, normal = self.frame(lon, lat)

        a0, a1 = self.along_bounds
        n0, n1 = self.normal_bounds
        along_ok = (along >= a0) & (along <= a1)
        sig_all = np.asarray(insar.sigma, dtype=float)
        if sig_all.ndim == 0:
            sig_all = np.full(len(lon), float(sig_all))
        value_all = np.asarray(insar.d_sar, dtype=float)
        m = (along_ok & (normal >= n0) & (normal <= n1) &
             np.isfinite(value_all) & (sig_all > 0.))
        if m.sum() == 0:
            raise ValueError('%s: no InSAR points inside the band and box'
                             % label)

        lon, lat = lon[m], lat[m]
        E, N = E[m], N[m]
        along, normal = along[m], normal[m]
        value = value_all[m]
        unit = np.asarray(insar.unit, dtype=float)[m]
        sigma = sig_all[m]
        count = np.ones(len(value), dtype=int)

        if along_bin > 0. or normal_bin > 0.:
            (lon, lat, along, normal, value, unit, sigma,
             count) = _bin(lon, lat, along, normal, value, unit, sigma,
                           a0, float(along_bin), n0, float(normal_bin))
            E, N, _, _ = self.frame(lon, lat)

        phi = np.radians(self.strike)
        s, c = np.sin(phi), np.cos(phi)
        cp = unit[:, 0]*s + unit[:, 1]*c
        cn = unit[:, 0]*c - unit[:, 1]*s
        cu = unit[:, 2]

        names, cols, bounds = [], [], {}
        if offset is not None:
            o = tuple(float(v) for v in offset)
            if o[0] > o[1]:
                raise ValueError('%s: offset min > max' % label)
            names.append('%s.offset' % label)
            cols.append(np.ones(len(value)))
            bounds[names[-1]] = o
        for k, axis in (('e', E), ('n', N)):
            if ramp and k in ramp:
                b = tuple(float(v) for v in ramp[k])
                if b[0] > b[1]:
                    raise ValueError('%s: ramp.%s min > max' % (label, k))
                names.append('%s.ramp_%s' % (label, k))
                cols.append(np.asarray(axis, dtype=float))
                bounds[names[-1]] = b
        bias = (np.column_stack(cols) if cols
                else np.zeros((len(value), 0)))

        n = len(value)
        return ProfileData(
            label=np.full(n, label, dtype=object),
            component=np.full(n, 'los', dtype=object),
            site=np.array(['%s_%05d' % (label, i) for i in range(n)],
                          dtype=object),
            lon=lon, lat=lat, E=E, N=N, along=along, normal=normal,
            value=value, sigma=sigma/w, cp=cp, cn=cn, cu=cu,
            count=count, bias=bias, bias_names=tuple(names),
            bias_bounds=bounds, strike=self.strike, weight=np.full(n, w))

    def cut(self, specs):
        '''
        Cut a list of dataset specs and concatenate them.

            {'type': 'gps',   'data': GPSData,  'label': ..., ...}
            {'type': 'insar', 'data': InSARData, 'label': ..., ...}

        Every configured GPS component is kept; selecting the rows a
        model can predict is the estimator's job, not the cut's.  InSAR
        LOS is kept because it mixes both components.
        '''
        out = []
        for spec in specs:
            typ = spec.get('type')
            data = spec['data']
            if typ == 'gps':
                comps = tuple(spec.get('components', ('par', 'norm')))
                offsets = {k: v for k, v in (spec.get('offsets') or {}).items()
                           if k in comps}
                out.append(self.cut_gps(
                    data,
                    label=spec.get('label', 'GPS'),
                    weight=spec.get('weight', 1.0),
                    components=comps,
                    offsets=offsets))
            elif typ == 'insar':
                out.append(self.cut_insar(
                    data,
                    label=spec.get('label', 'insar'),
                    weight=spec.get('weight', 1.0),
                    along_bin=spec.get('along_bin', 0.0),
                    normal_bin=spec.get('normal_bin', 0.0),
                    offset=spec.get('offset'),
                    ramp=spec.get('ramp')))
            else:
                raise ValueError('unknown dataset type %r' % (typ,))
        return concat(out)

    def sample(self, normal, along=0.0):
        '''
        Coordinates on the profile trace at given normal distances.
        Returns a dict of arrays.
        '''
        normal = np.asarray(normal, dtype=float)
        along  = np.broadcast_to(np.asarray(along, dtype=float), normal.shape)
        E, N = from_profile(along, normal, self.strike)
        lon, lat = enu_to_lonlat(E, N, self.origin, self.frame_name)
        return {'lon': lon, 'lat': lat, 'E': E, 'N': N,
                'along': along, 'normal': normal, 'strike': self.strike}

    def box(self):
        a0, a1 = self.along_bounds
        n0, n1 = self.normal_bounds
        a = np.array([a0, a1, a1, a0, a0], dtype=float)
        n = np.array([n0, n0, n1, n1, n0], dtype=float)
        E, N = from_profile(a, n, self.strike)
        return enu_to_lonlat(E, N, self.origin, self.frame_name)

    def trace(self):
        a = np.array(self.along_bounds, dtype=float)
        E, N = from_profile(a, np.zeros(2), self.strike)
        return enu_to_lonlat(E, N, self.origin, self.frame_name)

    def __repr__(self):
        return ('Profile(%r, strike=%g, along=(%g, %g), normal=(%g, %g), '
                'frame=%s)'
                % (self.origin, self.strike, self.along_bounds[0],
                   self.along_bounds[1], self.normal_bounds[0],
                   self.normal_bounds[1], self.frame_name))


def _bin(lon, lat, along, normal, value, unit, sigma,
         a0, along_width, n0, normal_width):
    '''
    Inverse-variance weighted bins on the along/normal grid.  A width > 0
    bins along its axis; a width of 0 keeps that axis unbinned (pixels
    averaged across it).
    '''
    i = (np.floor((normal - n0)/normal_width + 1e-3).astype(np.int64)
         if normal_width > 0. else np.zeros(len(normal), dtype=np.int64))
    j = (np.floor((along - a0)/along_width + 1e-3).astype(np.int64)
         if along_width > 0. else np.zeros(len(along), dtype=np.int64))
    if along_width > 0. and normal_width > 0.:
        idx = i * (j.max() + 1) + j
    elif along_width > 0.:
        idx = j
    else:
        idx = i
    _, inv = np.unique(idx, return_inverse=True)
    wgt = 1.0/sigma**2
    sw  = np.bincount(inv, weights=wgt)

    def mean(a):
        return np.bincount(inv, weights=wgt*a)/sw

    lon_b, lat_b = mean(lon), mean(lat)
    along_b, normal_b = mean(along), mean(normal)
    val_b = mean(value)
    u = np.column_stack([mean(unit[:, k]) for k in range(3)])
    u /= np.linalg.norm(u, axis=1)[:, None]
    sig_b = 1.0/np.sqrt(sw)
    count = np.bincount(inv).astype(int)
    return lon_b, lat_b, along_b, normal_b, val_b, u, sig_b, count


def _fmt(v):
    return '%.15g' % v


def write_profile(path, profile, data):
    '''
    Write one flat profile table (all datasets) as self-describing text.
    '''
    h = ['pysliprate-profile 2',
         'name %s' % profile.name,
         'frame %s' % profile.frame_name,
         'origin %s %s' % (_fmt(profile.origin[0]), _fmt(profile.origin[1])),
         'strike %s' % _fmt(profile.strike),
         'along %s %s' % (_fmt(profile.along_bounds[0]),
                          _fmt(profile.along_bounds[1])),
         'normal %s %s' % (_fmt(profile.normal_bounds[0]),
                           _fmt(profile.normal_bounds[1]))]
    for n in data.bias_names:
        lo, hi = data.bias_bounds[n]
        h.append('bias %s %s %s' % (n, _fmt(lo), _fmt(hi)))

    cols = ['label', 'component', 'site', 'lon', 'lat', 'along', 'normal',
            'value', 'sigma', 'cp', 'cn', 'cu', 'n', 'weight']
    cols += ['b:' + n for n in data.bias_names]
    h.append('columns ' + ' '.join(cols))

    with open(path, 'w', encoding='utf-8') as f:
        for line in h:
            f.write('# ' + line + '\n')
        for i in range(len(data)):
            row = [str(data.label[i]), str(data.component[i]),
                   str(data.site[i]),
                   _fmt(data.lon[i]), _fmt(data.lat[i]),
                   _fmt(data.along[i]), _fmt(data.normal[i]),
                   _fmt(data.value[i]), _fmt(data.sigma[i]),
                   _fmt(data.cp[i]), _fmt(data.cn[i]), _fmt(data.cu[i]),
                   '%d' % int(data.count[i]), _fmt(data.weight[i])]
            row += [_fmt(data.bias[i, k])
                    for k in range(len(data.bias_names))]
            f.write(' '.join(row) + '\n')


def read_profile(path):
    '''
    Read a table written by write_profile; returns (Profile, ProfileData).
    E/N are recomputed from lon/lat and the profile frame.
    '''
    hdr, bias, cols, rows = {}, [], None, []
    with open(path, encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith('#'):
                p = s[1:].strip().split()
                if not p:
                    continue
                if p[0] == 'bias':
                    bias.append((p[1], float(p[2]), float(p[3])))
                elif p[0] == 'columns':
                    cols = p[1:]
                elif p[0] == 'name':
                    hdr['name'] = ' '.join(p[1:])
                else:
                    hdr[p[0]] = p[1:]
            else:
                rows.append(s.split())
    if cols is None:
        raise ValueError('%s: no columns header' % path)
    if not rows:
        raise ValueError('%s: no data rows' % path)

    profile = Profile(
        (float(hdr['origin'][0]), float(hdr['origin'][1])),
        float(hdr['strike'][0]),
        (float(hdr['along'][0]), float(hdr['along'][1])),
        (float(hdr['normal'][0]), float(hdr['normal'][1])),
        name=hdr.get('name', 'profile'),
        frame=hdr.get('frame', ['pyginvc'])[0])

    ix = {c: i for i, c in enumerate(cols)}
    arr = {c: np.array([r[ix[c]] for r in rows]) for c in cols
           if c not in ('label', 'component', 'site')}
    arr = {c: v.astype(float) for c, v in arr.items()}

    lon, lat = arr['lon'], arr['lat']
    E, N, _, _ = profile.frame(lon, lat)
    bias_names = tuple(b[0] for b in bias)
    bias_bounds = {b[0]: (b[1], b[2]) for b in bias}
    bias_arr = (np.column_stack([arr['b:' + n] for n in bias_names])
                if bias_names else np.zeros((len(lon), 0)))

    data = ProfileData(
        label=np.array([r[ix['label']] for r in rows], dtype=object),
        component=np.array([r[ix['component']] for r in rows], dtype=object),
        site=np.array([r[ix['site']] for r in rows], dtype=object),
        lon=lon, lat=lat, E=E, N=N,
        along=arr['along'], normal=arr['normal'],
        value=arr['value'], sigma=arr['sigma'],
        cp=arr['cp'], cn=arr['cn'], cu=arr['cu'],
        count=arr['n'].astype(int), bias=bias_arr,
        bias_names=bias_names, bias_bounds=bias_bounds,
        strike=profile.strike, weight=arr['weight'])
    return profile, data
