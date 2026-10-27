#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
Read and validate the run configuration.

    config.yaml
        profile:   geometry
        model:     mechanism and fixed geometry
        datasets:  GPS / InSAR files plus cut options
        parameters: free parameters and uniform bounds; nuisance terms
                    are named <dataset label>.<component>, e.g.
                    GPS.par, GPS.norm, asc.offset, asc.ramp_e
        estimator: leastsq or mcmc options
        output:    directory and overwrite

In batch mode `catalog:` replaces profile/model/parameters and points at
a fault table, one row per profile.  All relative paths are resolved
against the directory holding the YAML.
'''
import csv
import os
import sys

import numpy as np
import yaml

from .data import GPSData, InSARData
from .models import model_from_config
from .profile import Profile, FRAMES


class ConfigError(Exception):
    '''Anything wrong in the YAML or the fault table.'''


_MODEL_KEYS = {
    'strike':  {'type', 'sense_ss'},
    'dip':     {'type', 'dip', 'sense_ds'},
    'oblique': {'type', 'dip', 'sense_ss', 'sense_ds'},
    'okada3d': {'type', 'nu', 'dip', 'width', 'length',
                'sense_ss', 'sense_ds'},
}

_GPS_BIAS   = ('par', 'norm')
_INSAR_BIAS = ('offset', 'ramp_e', 'ramp_n')

_EST_KEYS = {
    'leastsq': {'type', 'loss', 'f_scale', 'max_nfev', 'ndraws', 'seed'},
    'mcmc':    {'type', 'nwalkers', 'nsteps', 'discard', 'thin', 'seed',
                'jitter', 'jitter_bounds', 'progress'},
}

_EST_DEFAULTS = {
    'leastsq': {'loss': 'linear', 'f_scale': 1.0, 'max_nfev': 4000,
                'ndraws': 4000, 'seed': 0},
    'mcmc':    {'nwalkers': 32, 'nsteps': 5000, 'discard': 'auto',
                'thin': 'auto', 'seed': 0, 'jitter': True,
                'jitter_bounds': (-5., 2.), 'progress': False},
}


def _err(path, msg):
    raise ConfigError('%s: %s' % (path, msg))


def _keys(path, what, spec, allowed):
    extra = sorted(set(spec) - set(allowed))
    if extra:
        _err(path, '%s: unknown key(s) %s; allowed: %s'
             % (what, ', '.join(extra), ', '.join(sorted(allowed))))


def _pair(path, what, v):
    if not isinstance(v, (list, tuple)) or len(v) != 2:
        _err(path, '%s must be [min, max], got %r' % (what, v))
    lo, hi = float(v[0]), float(v[1])
    if not np.isfinite([lo, hi]).all():
        _err(path, '%s contains a non-finite value %r' % (what, v))
    if lo > hi:
        _err(path, '%s: min=%g > max=%g' % (what, lo, hi))
    return lo, hi


def _abspath(base, p):
    p = os.fspath(p)
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(base, p))


def model_key(typ):
    t = str(typ).lower()
    if t in ('strike', 'strike2d'):
        return 'strike'
    if t in ('dip', 'dip2d'):
        return 'dip'
    if t in ('oblique', 'oblique2d'):
        return 'oblique'
    if t in ('okada', 'okada3d'):
        return 'okada3d'
    raise ConfigError('unknown model type %r' % (typ,))


def _model(path, spec):
    if not isinstance(spec, dict):
        _err(path, 'model must be a mapping')
    try:
        mk = model_key(spec.get('type'))
    except ConfigError as e:
        _err(path, str(e))
    _keys(path, 'model', spec, _MODEL_KEYS[mk])
    try:
        return model_from_config(spec)
    except ValueError as e:
        _err(path, 'model: %s' % e)


def _angles(base, v):
    if isinstance(v, str):
        return _abspath(base, v)
    if isinstance(v, (list, tuple)):
        return [_abspath(base, x) if isinstance(x, str) else x for x in v]
    return v


def _dataset(path, base, spec):
    if not isinstance(spec, dict) or 'type' not in spec:
        _err(path, 'each dataset needs a type')
    typ = str(spec['type']).lower()

    if typ == 'gps':
        _keys(path, 'gps dataset', spec,
              {'type', 'path', 'label', 'weight', 'components'})
        if 'path' not in spec:
            _err(path, 'gps dataset needs a path')
        fp = _abspath(base, spec['path'])
        if not os.path.isfile(fp):
            _err(path, 'gps file does not exist: %s' % fp)

        comps = tuple(spec.get('components', ('par', 'norm')))
        if not comps or set(comps) - {'par', 'norm'}:
            _err(path, 'gps components must be a subset of (par, norm), '
                       'got %r' % (spec.get('components'),))
        label = str(spec.get('label') or
                    os.path.splitext(os.path.basename(fp))[0])
        return {'type': 'gps', 'data': GPSData(fp), 'label': label,
                'weight': float(spec.get('weight', 1.0)),
                'components': comps, 'offsets': {}}

    if typ == 'insar':
        _keys(path, 'insar dataset', spec,
              {'type', 'path', 'paths', 'format', 'label', 'weight',
               'incidence', 'heading', 'unit_bands', 'sigma', 'scale',
               'along_bin', 'normal_bin', 'overlap'})
        paths = spec.get('paths', spec.get('path'))
        if paths is None:
            _err(path, 'insar dataset needs path or paths')
        if isinstance(paths, (list, tuple)):
            paths = [_abspath(base, p) for p in paths]
            first = paths[0]
        else:
            paths = _abspath(base, paths)
            first = paths
        label = str(spec.get('label') or
                    os.path.splitext(os.path.basename(first))[0])
        data = InSARData(
            paths,
            format=spec.get('format'),
            incidence=_angles(base, spec.get('incidence')),
            heading=_angles(base, spec.get('heading')),
            unit_bands=spec.get('unit_bands'),
            sigma=spec.get('sigma', 3.0),
            scale=spec.get('scale', 1.0),
            overlap=spec.get('overlap', 'first'))
        return {'type': 'insar', 'data': data, 'label': label,
                'weight': float(spec.get('weight', 1.0)),
                'along_bin': float(spec.get('along_bin', 0.0)),
                'normal_bin': float(spec.get('normal_bin', 0.0)),
                'offset': None, 'ramp': None}

    _err(path, 'dataset type must be gps or insar, got %r' % (spec['type'],))


def _bias_names(spec):
    '''
    Every nuisance parameter a dataset can take, named by its label.
    Bounds live in the top-level `parameters:` block.
    '''
    label = spec['label']
    suffix = _GPS_BIAS if spec['type'] == 'gps' else _INSAR_BIAS
    return ['%s.%s' % (label, k) for k in suffix]


def _all_bias_names(datasets):
    names = set()
    for d in datasets:
        names.update(_bias_names(d))
    return names


def _bias_targets(datasets):
    '''
    {'GPS.par': (dataset spec, 'par'), ...}
    '''
    out = {}
    for d in datasets:
        for nm in _bias_names(d):
            out[nm] = (d, nm.rsplit('.', 1)[1])
    return out


def _apply_bias_parameters(path, datasets, params):
    '''
    Turn `parameters:` entries such as `GPS.par: [-50, 50]` or
    `asc.offset: [-50, 50]` into offsets on the matching dataset.
    '''
    targets = _bias_targets(datasets)
    for name, b in params.items():
        if name not in targets:
            if '.' in name:
                _err(path, 'parameters: unknown nuisance parameter %r; '
                          'known: %s' % (name, ', '.join(sorted(targets))))
            continue
        d, suffix = targets[name]
        if d['type'] == 'gps':
            if suffix not in d['components']:
                _err(path, 'parameters.%s: component %r is not in the '
                          'dataset components %r'
                     % (name, suffix, d['components']))
            d['offsets'][suffix] = b
        elif suffix == 'offset':
            d['offset'] = b
        else:
            if d['ramp'] is None:
                d['ramp'] = {}
            d['ramp'][suffix[-1]] = b


def _active_slips(mk, senses, options=None):
    '''
    Slip parameters a catalog row fits.  A model with only one sense
    fits only that component; with neither given both are kept.
    '''
    if mk == 'strike':
        return ('slip_ss',)
    if mk == 'dip':
        return ('slip_ds',)
    options = options or {}
    ss = senses.get('sense_ss', options.get('sense_ss'))
    ds = senses.get('sense_ds', options.get('sense_ds'))
    if ss is None and ds is None:
        ss, ds = 'left', 'normal'
    return tuple(n for n, s in (('slip_ss', ss), ('slip_ds', ds))
                 if s is not None)


def _model_params(mk, slips):
    '''
    Free parameters a model type can take, restricted to the active
    slip components.
    '''
    if mk == 'strike':
        return {'slip_ss', 'lock'}
    if mk == 'dip':
        return {'slip_ds', 'dip', 'lock'}
    if mk == 'oblique':
        out = {'lock'} | set(slips)
        if 'slip_ds' in slips:
            out.add('dip')
        return out
    return {'dip', 'lock', 'width', 'length'} | set(slips)


def _all_model_params(mk):
    '''
    Every parameter the model type could take, before senses narrow it.
    Bounds for the inactive ones are accepted and then ignored, which
    lets one `parameters:` block serve a batch of mixed rows.
    '''
    return _model_params(mk, ('slip_ss', 'slip_ds'))


def _estimator(path, spec):
    if spec is None:
        spec = {}
    if not isinstance(spec, dict):
        _err(path, 'estimator must be a mapping')
    typ = str(spec.get('type', 'mcmc')).lower()
    if typ in ('leastsq', 'least_squares', 'lsq'):
        typ = 'leastsq'
    elif typ == 'mcmc':
        typ = 'mcmc'
    else:
        _err(path, 'estimator.type must be leastsq or mcmc, got %r'
                   % (spec.get('type'),))
    _keys(path, 'estimator', spec, _EST_KEYS[typ])

    out = dict(_EST_DEFAULTS[typ])
    out.update({k: v for k, v in spec.items() if k != 'type'})
    out['type'] = typ

    if typ == 'mcmc':
        for k in ('nwalkers', 'nsteps'):
            if not isinstance(out[k], (int, np.integer)) or out[k] <= 0:
                _err(path, 'estimator.%s must be a positive integer, got %r'
                           % (k, out[k]))
        if not isinstance(out['seed'], (int, np.integer)) or out['seed'] < 0:
            _err(path, 'estimator.seed must be a non-negative integer')
        for k in ('discard', 'thin'):
            if out[k] != 'auto' and not (isinstance(out[k], int)
                                         and out[k] > 0):
                _err(path, "estimator.%s must be 'auto' or a positive "
                           'integer, got %r' % (k, out[k]))
        out['jitter_bounds'] = _pair(path, 'estimator.jitter_bounds',
                                     out['jitter_bounds'])
    return out


def _parameters(path, spec, known):
    if spec is None:
        return {}
    if not isinstance(spec, dict):
        _err(path, 'parameters must be a mapping')
    out = {}
    for k, v in spec.items():
        if known is not None and k not in known:
            _err(path, 'parameters: unknown parameter %r; the model and '
                       'datasets have %s'
                       % (k, ', '.join(sorted(known))))
        out[k] = _pair(path, 'parameters.%s' % k, v)
    return out


def load_config(path):
    '''
    Read a single-profile or batch configuration.
    '''
    try:
        with open(path, encoding='utf-8') as f:
            raw = yaml.safe_load(f)
    except OSError as e:
        raise ConfigError('cannot read %s: %s' % (path, e))
    except yaml.YAMLError as e:
        raise ConfigError('failed to parse %s as YAML: %s' % (path, e))
    if not isinstance(raw, dict):
        _err(path, 'the top level must be a mapping')

    base = os.path.dirname(os.path.abspath(path))
    _keys(path, 'top level', raw,
          {'profile', 'model', 'datasets', 'parameters', 'estimator',
           'output', 'catalog', 'frame'})

    out = raw.get('output') or {}
    if not isinstance(out, dict):
        _err(path, 'output must be a mapping')
    _keys(path, 'output', out, {'dir', 'overwrite'})
    outdir = _abspath(base, out.get('dir', 'out'))

    specs = raw.get('datasets')
    if not isinstance(specs, list) or not specs:
        _err(path, 'datasets must be a non-empty list')
    datasets = [_dataset(path, base, d) for d in specs]

    labels = [d['label'] for d in datasets]
    dup = sorted({x for x in labels if labels.count(x) > 1})
    if dup:
        _err(path, 'duplicate dataset label(s): %s' % ', '.join(dup))
    bias_names = _all_bias_names(datasets)

    estimator = _estimator(path, raw.get('estimator'))
    common = {'base': base, 'datasets': datasets, 'estimator': estimator,
              'output_dir': outdir,
              'overwrite': bool(out.get('overwrite', True))}

    if 'catalog' in raw:
        cat = _abspath(base, raw['catalog'])
        if not os.path.isfile(cat):
            _err(path, 'no such catalog file: %s' % cat)
        opts = raw.get('model') or {}
        if not isinstance(opts, dict):
            _err(path, 'model must be a mapping')
        params = _parameters(path, raw.get('parameters'), None)
        _apply_bias_parameters(path, datasets, params)
        frame = str(raw.get('frame', 'pyginvc'))
        rows = load_catalog(cat, params, frame, bias_names=bias_names,
                            model_options=opts)
        return dict(common, batch=True, catalog=cat, parameters=params,
                    model_options=opts, rows=rows)

    for k in ('profile', 'model'):
        if k not in raw:
            _err(path, 'missing %r (or use catalog: for batch mode)' % k)
    ps = raw['profile']
    if not isinstance(ps, dict):
        _err(path, 'profile must be a mapping')
    _keys(path, 'profile', ps,
          {'name', 'origin', 'strike', 'along', 'normal', 'frame'})
    for k in ('origin', 'strike', 'along', 'normal'):
        if k not in ps:
            _err(path, 'profile.%s is required' % k)
    frame = str(ps.get('frame', raw.get('frame', 'pyginvc')))
    if frame not in FRAMES:
        _err(path, 'profile.frame must be one of %s, got %r'
             % (', '.join(FRAMES), frame))
    try:
        profile = Profile(ps['origin'], ps['strike'], ps['along'],
                          ps['normal'], name=ps.get('name', 'profile'),
                          frame=frame)
    except (TypeError, ValueError) as e:
        _err(path, 'profile: %s' % e)

    model = _model(path, raw['model'])
    known = (set(model.params) | bias_names |
             _all_model_params(model_key(model.name)))
    params = _parameters(path, raw.get('parameters'), known)
    _apply_bias_parameters(path, datasets, params)
    return dict(common, batch=False, profile=profile, model=model,
                model_spec=raw['model'], parameters=params)


GEOM_COLS = ('name', 'lon', 'lat', 'strike',
             'along_min', 'along_max', 'normal_min', 'normal_max', 'model')


def load_catalog(path, defaults, frame='pyginvc', bias_names=(),
                 model_options=None):
    '''
    Read a fault table, one row per profile.

    `slip_min/max` is mapped to the slip parameters the row's senses
    activate; optional sense_ss/sense_ds columns choose the direction.
    One sense alone fits only that component.
    '''
    header, hln, rows = None, None, []
    try:
        with open(path, newline='', encoding='utf-8') as f:
            for lineno, line in enumerate(f, 1):
                s = line.strip()
                if not s or s.startswith('#'):
                    continue
                fields = next(csv.reader([line]))
                if header is None:
                    header, hln = [x.strip() for x in fields], lineno
                else:
                    rows.append((lineno, fields))
    except OSError as e:
        raise ConfigError('cannot read %s: %s' % (path, e))
    if header is None:
        _err(path, 'no header line')

    idx = {c: i for i, c in enumerate(header)}
    for c in GEOM_COLS:
        if c not in idx:
            _err(path, '%s:%d: missing required column %r'
                 % (path, hln, c))

    out, seen = [], set()
    for lineno, fields in rows:
        def cell(col):
            i = idx[col]
            return fields[i].strip() if i < len(fields) else ''

        name = cell('name')
        if not name:
            _err(path, '%s:%d: empty name' % (path, lineno))
        if name in seen:
            _err(path, '%s:%d: duplicate name %r' % (path, lineno, name))
        seen.add(name)

        model_type = cell('model') or 'strike2d'
        # An empty cell explicitly disables that sense, overriding the
        # shared model block; an absent column falls back to it.
        senses = {k: (cell(k) or None) for k in ('sense_ss', 'sense_ds')
                  if k in idx}
        try:
            mk = model_key(model_type)
            slips = _active_slips(mk, senses, model_options)
            allowed = _model_params(mk, slips)
        except ConfigError as e:
            _err(path, '%s:%d: %s' % (path, lineno, e))

        try:
            profile = Profile((float(cell('lon')), float(cell('lat'))),
                              float(cell('strike')),
                              (float(cell('along_min')),
                               float(cell('along_max'))),
                              (float(cell('normal_min')),
                               float(cell('normal_max'))),
                              name=name, frame=frame)
        except (TypeError, ValueError) as e:
            _err(path, '%s:%d: %s' % (path, lineno, e))

        params = {k: v for k, v in defaults.items()
                  if k in allowed or k in bias_names}
        for col in header:
            if not col.endswith('_min'):
                continue
            p = col[:-4]
            if p in ('along', 'normal'):
                continue
            lo_s, hi_s = cell(p + '_min'), cell(p + '_max')
            if lo_s == '' and hi_s == '':
                continue
            if p == 'slip':
                targets = list(slips)
            elif p == 'constant':
                print('warning: %s:%d: constant is now per-dataset; ignored'
                      % (path, lineno), file=sys.stderr)
                continue
            elif p in allowed:
                targets = [p]
            else:
                print('warning: %s:%d: model=%s does not use %s; ignored'
                      % (path, lineno, model_type, p), file=sys.stderr)
                continue
            if lo_s == '' or hi_s == '':
                _err(path, '%s:%d: %s_min and %s_max must both be given or '
                           'both empty' % (path, lineno, p, p))
            try:
                b = _pair(path, '%s:%d %s' % (path, lineno, p),
                          (float(lo_s), float(hi_s)))
            except ValueError:
                _err(path, '%s:%d: %s bounds are not numbers: %r %r'
                     % (path, lineno, p, lo_s, hi_s))
            for t in targets:
                params[t] = b

        out.append({'name': name, 'profile': profile,
                    'model_type': model_type, 'parameters': params,
                    'senses': senses, 'line': lineno})
    if not out:
        _err(path, 'the fault table has no data rows')
    return out
