#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
Forward models for one profile.

Every model has a `params` dict of default values and a
`predict(par, data)` method returning

    (v_par, v_norm, v_up)

at the rows of a ProfileData.  The 2D models use `data.normal`; Okada3D
uses `data.E`, `data.N` and `data.strike`.  There is no model base class
and no parameter registry: the estimator only needs the defaults and
the method.

`components` lists the fault components a model fits.  StrikeSlip2D
fits only `par`, DipSlip2D only `norm`, and Oblique2D/Okada3D fit both
unless `sense_ss` or `sense_ds` is None, which switches that component
off and removes its slip parameter.

Sign conventions:

    strike-slip  sense_ss='left'  -> s > 0 (left-lateral)
    dip-slip     sense_ds='normal' -> legacy kernel s > 0
                 sense_ds='thrust' -> Okada ds > 0

The legacy dipslip2d returns a vertical component positive downward;
the models negate it so v_up is up-positive.
'''
import numpy as np

from .DipSlip2D import dipslip2d
from .StrikeSlip2D import strikeslip2d

__all__ = ['StrikeSlip2D', 'DipSlip2D', 'Oblique2D', 'Okada3D',
           'model_from_config']

_SS = {'left': 1.0, 'right': -1.0}
_DS2D = {'normal': 1.0, 'thrust': -1.0}      # legacy kernel: s > 0 normal
_DS3D = {'normal': -1.0, 'thrust': 1.0}      # Okada: ds > 0 thrust


class StrikeSlip2D(object):
    params = {'slip_ss': 5.0, 'lock': 15.0}
    labels = {'slip_ss': 'slip rate (mm/yr)', 'lock': 'locking depth (km)'}
    name = 'strike2d'
    components = ('par',)

    def __init__(self, sense_ss='left'):
        if sense_ss not in _SS:
            raise ValueError('sense_ss must be left or right, got %r'
                             % (sense_ss,))
        self.sense_ss = sense_ss
        self.params = dict(type(self).params)

    def predict(self, par, data):
        x = np.atleast_1d(np.asarray(data.normal, dtype=float))
        s = _SS[self.sense_ss]*par['slip_ss']
        vpar = strikeslip2d(s, np.inf, par['lock'], 0, 0.).forward(x)
        zero = np.zeros_like(x)
        return vpar, zero, zero


class DipSlip2D(object):
    params = {'slip_ds': 5.0, 'dip': 45.0, 'lock': 15.0}
    labels = {'slip_ds': 'dip-slip rate (mm/yr)', 'dip': 'dip (degree)',
              'lock': 'locking depth (km)'}
    name = 'dip2d'
    components = ('norm',)

    def __init__(self, dip=45., sense_ds='normal'):
        if sense_ds not in _DS2D:
            raise ValueError('sense_ds must be normal or thrust, got %r'
                             % (sense_ds,))
        if not 0. < float(dip) <= 90.:
            raise ValueError('dip must be in (0, 90], got %r' % (dip,))
        self.sense_ds = sense_ds
        self.params = dict(type(self).params)
        self.params['dip'] = float(dip)

    def predict(self, par, data):
        x = np.atleast_1d(np.asarray(data.normal, dtype=float))
        s = _DS2D[self.sense_ds]*par['slip_ds']
        u, w = dipslip2d(s, par['dip'], par['lock'], 0.).forward(x)
        return np.zeros_like(x), u, -w


class Oblique2D(object):
    params = {'slip_ss': 5.0, 'slip_ds': 5.0, 'dip': 45.0, 'lock': 15.0}
    labels = {'slip_ss': 'strike-slip rate (mm/yr)',
              'slip_ds': 'dip-slip rate (mm/yr)',
              'dip': 'dip (degree)', 'lock': 'locking depth (km)'}
    name = 'oblique2d'

    def __init__(self, dip=45., sense_ss='left', sense_ds='normal'):
        if sense_ss is not None and sense_ss not in _SS:
            raise ValueError('sense_ss must be left or right, got %r'
                             % (sense_ss,))
        if sense_ds is not None and sense_ds not in _DS2D:
            raise ValueError('sense_ds must be normal or thrust, got %r'
                             % (sense_ds,))
        if sense_ss is None and sense_ds is None:
            raise ValueError('at least one of sense_ss and sense_ds must be '
                             'given')
        if not 0. < float(dip) <= 90.:
            raise ValueError('dip must be in (0, 90], got %r' % (dip,))
        self.sense_ss = sense_ss
        self.sense_ds = sense_ds
        self.components = tuple(
            c for c, s in (('par', sense_ss), ('norm', sense_ds))
            if s is not None)
        self.params = dict(type(self).params)
        self.params['dip'] = float(dip)
        if sense_ss is None:
            del self.params['slip_ss']
        if sense_ds is None:
            # The 2-D strike-slip kernel does not use the dip.
            del self.params['slip_ds']
            del self.params['dip']

    def predict(self, par, data):
        x = np.atleast_1d(np.asarray(data.normal, dtype=float))
        zero = np.zeros_like(x)
        if self.sense_ss is None:
            vpar = zero
        else:
            sss = _SS[self.sense_ss]*par['slip_ss']
            vpar = strikeslip2d(sss, np.inf, par['lock'], 0, 0.).forward(x)
        if self.sense_ds is None:
            vnorm = vup = zero
        else:
            sds = _DS2D[self.sense_ds]*par['slip_ds']
            u, w = dipslip2d(sds, par['dip'], par['lock'], 0.).forward(x)
            vnorm, vup = u, -w
        return vpar, vnorm, vup


class Okada3D(object):
    '''
    One uniform-slip rectangle below the locking depth, dipping toward
    +normal.  Uses PyGINVC's Okada.disloc through OkadaKernel (only the
    numpy-2 np.mat line is changed).
    '''
    params = {'slip_ss': 5.0, 'slip_ds': 1.0, 'lock': 15.0,
              'dip': 90.0, 'width': 200.0, 'length': 1000.0}
    labels = {'slip_ss': 'strike-slip rate (mm/yr)',
              'slip_ds': 'dip-slip rate (mm/yr)',
              'lock': 'top depth (km)', 'dip': 'dip (degree)',
              'width': 'down-dip width (km)',
              'length': 'along-strike length (km)'}
    name = 'okada3d'

    def __init__(self, nu=0.25, dip=90., width=200., length=1000.,
                 sense_ss='left', sense_ds='normal'):
        if sense_ss is not None and sense_ss not in _SS:
            raise ValueError('sense_ss must be left or right, got %r'
                             % (sense_ss,))
        if sense_ds is not None and sense_ds not in _DS3D:
            raise ValueError('sense_ds must be normal or thrust, got %r'
                             % (sense_ds,))
        if sense_ss is None and sense_ds is None:
            raise ValueError('at least one of sense_ss and sense_ds must be '
                             'given')
        try:
            from .okada_kernel import OkadaKernel
        except ImportError as e:
            raise ImportError(
                'Okada3D needs PyGINVC with a compiled okada85.  Install '
                'PyGINVC and build okada85 for this platform (see the '
                'README), then try again.') from e
        self._okada = OkadaKernel
        self.nu = float(nu)
        self.sense_ss = sense_ss
        self.sense_ds = sense_ds
        self.components = tuple(
            c for c, s in (('par', sense_ss), ('norm', sense_ds))
            if s is not None)
        self.params = dict(type(self).params)
        self.params.update(dip=float(dip), width=float(width),
                           length=float(length))
        if sense_ss is None:
            del self.params['slip_ss']
        if sense_ds is None:
            del self.params['slip_ds']

    def predict(self, par, data):
        ss = (0. if self.sense_ss is None
              else _SS[self.sense_ss]*par['slip_ss'])
        ds = (0. if self.sense_ds is None
              else _DS3D[self.sense_ds]*par['slip_ds'])
        dip, wid = float(par['dip']), float(par['width'])
        leng, lock = float(par['length']), float(par['lock'])

        dipr = np.radians(dip)
        dep  = lock + wid*np.sin(dipr)          # lower edge depth
        d    = wid*np.cos(dipr)
        phi  = np.radians(data.strike)
        delE = d*np.cos(phi)
        delN = -d*np.sin(phi)

        disgeom = np.array([leng, wid, dep, dip, data.strike,
                            delE, delN, ss, ds, 0.])
        x = np.column_stack((np.atleast_1d(np.asarray(data.E, float)),
                             np.atleast_1d(np.asarray(data.N, float))))
        u = self._okada.disloc(self.nu, disgeom, x)
        vpar  = np.sin(phi)*u[:, 0] + np.cos(phi)*u[:, 1]
        vnorm = np.cos(phi)*u[:, 0] - np.sin(phi)*u[:, 1]
        return vpar, vnorm, u[:, 2]


def model_from_config(spec):
    '''
    Build a model from the `model:` section of a YAML config.
    '''
    if not isinstance(spec, dict):
        raise ValueError('model must be a mapping')
    typ = str(spec.get('type', '')).lower()
    if typ in ('strike2d', 'strike'):
        return StrikeSlip2D(sense_ss=spec.get('sense_ss') or 'left')
    if typ in ('dip2d', 'dip'):
        return DipSlip2D(dip=spec.get('dip', 45.),
                         sense_ds=spec.get('sense_ds') or 'normal')
    if typ in ('oblique2d', 'oblique'):
        # One sense alone fits that component only; with neither given,
        # keep the historical default of both components.
        ss, ds = spec.get('sense_ss'), spec.get('sense_ds')
        if ss is None and ds is None:
            ss, ds = 'left', 'normal'
        return Oblique2D(dip=spec.get('dip', 45.), sense_ss=ss, sense_ds=ds)
    if typ in ('okada3d', 'okada'):
        ss, ds = spec.get('sense_ss'), spec.get('sense_ds')
        if ss is None and ds is None:
            ss, ds = 'left', 'normal'
        return Okada3D(nu=spec.get('nu', 0.25),
                       dip=spec.get('dip', 90.),
                       width=spec.get('width', 200.),
                       length=spec.get('length', 1000.),
                       sense_ss=ss, sense_ds=ds)
    raise ValueError('model.type must be strike2d, dip2d, oblique2d or '
                     'okada3d, got %r' % (spec.get('type'),))
