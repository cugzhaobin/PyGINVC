#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
The profile workflow, called by the scripts in scripts/.

    prepare(cfgfile)   read data, cut the profile, save data/profile.csv
    invert(cfgfile)    prepare + fit + fit/result/post files
    batch(cfgfile)     one catalog row per profile

Figures are produced by report.py (next step); this module only writes
the tables, the posterior samples and the summary.
'''
import csv
import json
import os
import shutil
import sys

import numpy as np

from . import report
from .config import ConfigError, load_config, model_key
from .models import model_from_config
from .profile import write_profile, read_profile
from .inversion import Result, fit, fit_components, fit_data, predict

SUMMARY_COLS = ('profile', 'datasets', 'estimator', 'param', 'mean', 'std',
                'p16', 'p50', 'p84', 'acceptance', 'tau', 'ndata',
                'status', 'message')


def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, np.ndarray):
        return [_jsonable(x) for x in v.tolist()]
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, (np.integer, int)):
        return int(v)
    return v


def _mkdirs(pdir):
    for d in ('data', 'fit', 'post', 'figures'):
        os.makedirs(os.path.join(pdir, d), exist_ok=True)


def _data_figures(figdir, profile, datasets, data):
    report.map_figure(os.path.join(figdir, 'map.png'), profile, datasets,
                      data)
    report.units_figure(os.path.join(figdir, 'units.png'), data)
    report.profile_figure(os.path.join(figdir, 'profile.png'), data)
    report.bins_figure(os.path.join(figdir, 'bins.png'), data)


def _fit_figures(figdir, data, model, result, pred):
    report.fit_figure(os.path.join(figdir, 'fit.png'), data, pred, model)
    report.joint_figure(os.path.join(figdir, 'joint.png'), data, pred, model)
    report.posterior_figure(os.path.join(figdir, 'posterior.png'), data,
                            model, result)
    report.corner_figure(os.path.join(figdir, 'corner.png'), result)


def _check_outdir(outdir, overwrite):
    if os.path.isdir(outdir) and os.listdir(outdir) and not overwrite:
        raise ConfigError('%s is not empty and overwrite is false' % outdir)
    os.makedirs(outdir, exist_ok=True)


def _fmt(v):
    if isinstance(v, float):
        return '' if not np.isfinite(v) else '%.6g' % v
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return v


def _model_for(cfg, row=None):
    if row is None:
        return cfg['model'], cfg['model_spec']
    spec = dict(cfg['model_options'])
    senses = row.get('senses', {})
    if model_key(row['model_type']) in ('strike', 'dip'):
        # These models have a single sense that cannot be switched off,
        # so an empty catalog cell falls back to the shared model block
        # instead of forcing the model's default sign.
        senses = {k: v for k, v in senses.items() if v is not None}
    spec.update(senses)
    spec['type'] = row['model_type']
    return model_from_config(spec), spec


def _write_json(path, profile, model, model_spec, estimator, parameters,
                datasets, result=None):
    d = {'name': profile.name,
         'origin': list(profile.origin),
         'strike': profile.strike,
         'along_bounds': list(profile.along_bounds),
         'normal_bounds': list(profile.normal_bounds),
         'frame': profile.frame_name,
         'model_spec': model_spec,
         'model_name': model.name,
         'datasets': [{'label': s['label'], 'type': s['type'],
                       'weight': s['weight']} for s in datasets],
         'parameters': {k: list(v) for k, v in parameters.items()},
         'estimator': estimator,
         'result': result.meta if result is not None else None}
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(_jsonable(d), f, indent=2, ensure_ascii=False)


def _write_result(path, profile, model, model_spec, result, data):
    d = {'profile': profile.name, 'frame': profile.frame_name,
         'model': model.name, 'model_spec': model_spec,
         'datasets': sorted({str(x) for x in data.label}),
         'labels': result.labels, 'best': result.best,
         'fixed': result.fixed, 'stats': result.stats(),
         'meta': result.meta}
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(_jsonable(d), f, indent=2, ensure_ascii=False)


def _write_fit(path, data, pred, model=None):
    fitted = (fit_components(data, model) if model is not None
              else np.ones(len(data), dtype=bool))
    with open(path, 'w', encoding='utf-8') as f:
        f.write('label component site lon lat along normal '
                'observed predicted residual sigma fitted\n')
        for i in range(len(data)):
            if fitted[i]:
                p, r = pred[i], data.value[i] - pred[i]
            else:
                p = r = float('nan')
            f.write('%s %s %s %.8f %.8f %.6f %.6f '
                    '%.6f %.6f %.6f %.6f %d\n'
                    % (data.label[i], data.component[i], data.site[i],
                       data.lon[i], data.lat[i], data.along[i],
                       data.normal[i], data.value[i], p, r, data.sigma[i],
                       int(fitted[i])))


def _datasets_label(data):
    return '+'.join(sorted({str(x) for x in data.label}))


def _summary_rows(profile, data, result):
    st = result.stats()
    common = {'profile': profile.name,
              'datasets': _datasets_label(data),
              'estimator': result.meta.get('estimator', ''),
              'acceptance': ('' if result.acceptance is None
                             else result.acceptance),
              'tau': ('' if result.tau is None
                      else float(np.max(result.tau))),
              'ndata': result.meta.get('ndata', ''), 'status': 'ok',
              'message': ''}
    rows = []
    for lab in result.labels:
        d = st[lab]
        rows.append(dict(common, param=lab, mean=d['mean'], std=d['std'],
                         p16=d['p16'], p50=d['p50'], p84=d['p84']))
    return rows


def _prepared_row(profile, data):
    return {'profile': profile.name, 'datasets': _datasets_label(data),
            'estimator': '', 'param': '', 'mean': '', 'std': '', 'p16': '',
            'p50': '', 'p84': '', 'acceptance': '', 'tau': '', 'ndata': '',
            'status': 'prepared', 'message': ''}


def _failed_row(name, datasets, msg):
    labels = '+'.join(d['label'] for d in datasets) if datasets else ''
    return {'profile': name, 'datasets': labels, 'estimator': '',
            'param': '', 'mean': '', 'std': '', 'p16': '', 'p50': '',
            'p84': '', 'acceptance': '', 'tau': '', 'ndata': '',
            'status': 'failed', 'message': str(msg)}


def _write_summary(path, rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(SUMMARY_COLS)
        for r in rows:
            w.writerow([_fmt(r.get(c, '')) for c in SUMMARY_COLS])


def _used_parameters(profile, parameters, model, data):
    '''
    Keep only parameters this model and these data actually use.  Batch
    rows share one `parameters` block, so a component or slip parameter
    switched off by the row's model must be dropped here.  The data
    still hold every component; `fit_data` decides which rows and bias
    columns enter the fit.
    '''
    fdata = fit_data(data, model)
    known = set(model.params) | set(fdata.bias_bounds)
    used, dropped = {}, []
    for k, v in parameters.items():
        if k in known:
            used[k] = v
        else:
            dropped.append(k)
    if dropped:
        print('warning: %s: ignoring parameter(s) not used by this model '
              'and data: %s' % (profile.name, ', '.join(sorted(dropped))),
              file=sys.stderr)
    return used


def _fit_and_save(cfg, profile, model, model_spec, data, pdir, parameters):
    parameters = _used_parameters(profile, parameters, model, data)
    _mkdirs(pdir)
    write_profile(os.path.join(pdir, 'data', 'profile.csv'), profile, data)
    _write_json(os.path.join(pdir, 'profile.json'), profile, model,
                model_spec, cfg['estimator'], parameters, cfg['datasets'])

    figdir = os.path.join(pdir, 'figures')
    _data_figures(figdir, profile, cfg['datasets'], data)

    est = dict(cfg['estimator'])
    method = est.pop('type')
    result = fit(data, model, parameters, method=method, **est)
    pred = predict(data, model, result.pars())

    _fit_figures(figdir, data, model, result, pred)
    _write_fit(os.path.join(pdir, 'fit', 'profile_fit.csv'), data, pred,
               model)
    if result.draws is not None:
        np.save(os.path.join(pdir, 'post', 'flat.npy'), result.draws)
    if result.chain is not None:
        np.save(os.path.join(pdir, 'post', 'chain.npy'), result.chain)
    _write_result(os.path.join(pdir, 'result.json'), profile, model,
                  model_spec, result, data)
    _write_json(os.path.join(pdir, 'profile.json'), profile, model,
                model_spec, cfg['estimator'], parameters, cfg['datasets'],
                result)
    return result


def prepare(cfgfile):
    '''
    Read data, cut the profile and save data/profile.csv (no fit).
    '''
    cfg = load_config(cfgfile)
    if cfg['batch']:
        raise ConfigError('prepare needs a single-profile config, not a '
                          'catalog; use invert/batch for catalog runs')
    outdir = cfg['output_dir']
    _check_outdir(outdir, cfg['overwrite'])
    shutil.copyfile(cfgfile, os.path.join(outdir, 'config.resolved.yaml'))

    profile = cfg['profile']
    data = profile.cut(cfg['datasets'])
    parameters = _used_parameters(profile, cfg['parameters'], cfg['model'],
                                  data)
    pdir = os.path.join(outdir, profile.name)
    _mkdirs(pdir)
    write_profile(os.path.join(pdir, 'data', 'profile.csv'), profile, data)
    _write_json(os.path.join(pdir, 'profile.json'), profile, cfg['model'],
                cfg['model_spec'], cfg['estimator'], parameters,
                cfg['datasets'])
    _data_figures(os.path.join(pdir, 'figures'), profile, cfg['datasets'],
                  data)
    _write_summary(os.path.join(outdir, 'summary.csv'),
                   [_prepared_row(profile, data)])
    print('prepared %s: %d rows (%s) -> %s'
          % (profile.name, len(data), _datasets_label(data), pdir))
    return 0


def invert(cfgfile):
    '''
    Prepare and fit one profile; save fit/result/post and summary.csv.
    '''
    cfg = load_config(cfgfile)
    if cfg['batch']:
        return batch(cfgfile)
    outdir = cfg['output_dir']
    _check_outdir(outdir, cfg['overwrite'])
    shutil.copyfile(cfgfile, os.path.join(outdir, 'config.resolved.yaml'))

    profile = cfg['profile']
    data = profile.cut(cfg['datasets'])
    pdir = os.path.join(outdir, profile.name)
    result = _fit_and_save(cfg, profile, cfg['model'], cfg['model_spec'],
                           data, pdir, cfg['parameters'])
    rows = _summary_rows(profile, data, result)
    _write_summary(os.path.join(outdir, 'summary.csv'), rows)
    print('done %s: %d rows, %s -> %s'
          % (profile.name, len(data), _datasets_label(data), pdir))
    return 0


def batch(cfgfile):
    '''
    One catalog row per profile.
    '''
    cfg = load_config(cfgfile)
    if not cfg['batch']:
        raise ConfigError('batch needs a config with a catalog: key')
    outdir = cfg['output_dir']
    _check_outdir(outdir, cfg['overwrite'])
    shutil.copyfile(cfgfile, os.path.join(outdir, 'config.resolved.yaml'))

    rows = []
    for i, row in enumerate(cfg['rows'], 1):
        profile = row['profile']
        try:
            model, spec = _model_for(cfg, row)
            data = profile.cut(cfg['datasets'])
            pdir = os.path.join(outdir, profile.name)
            result = _fit_and_save(cfg, profile, model, spec, data, pdir,
                                   row['parameters'])
            rows.extend(_summary_rows(profile, data, result))
            print('[%d/%d] %s (%s): %d rows'
                  % (i, len(cfg['rows']), profile.name, row['model_type'],
                     len(data)))
        except (ConfigError, ValueError, RuntimeError) as e:
            print('  failed: %s' % e, file=sys.stderr)
            rows.append(_failed_row(profile.name, cfg['datasets'], e))

    _write_summary(os.path.join(outdir, 'summary.csv'), rows)
    print('done: %d profiles -> %s' % (len(cfg['rows']), outdir))
    return 0


def plot(outdir):
    '''
    Redraw the profile/fit/posterior figures from data/profile.csv,
    profile.json and the saved samples; the raw velocity files are not
    read.
    '''
    with open(os.path.join(outdir, 'profile.json'), encoding='utf-8') as f:
        meta = json.load(f)
    profile, data = read_profile(
        os.path.join(outdir, 'data', 'profile.csv'))
    model = model_from_config(meta['model_spec'])

    with open(os.path.join(outdir, 'result.json'), encoding='utf-8') as f:
        rj = json.load(f)
    flat_p = os.path.join(outdir, 'post', 'flat.npy')
    chain_p = os.path.join(outdir, 'post', 'chain.npy')
    flat = np.load(flat_p) if os.path.isfile(flat_p) else None
    chain = np.load(chain_p) if os.path.isfile(chain_p) else None
    result = Result(labels=rj['labels'], best=rj['best'],
                    fixed=rj['fixed'], draws=flat, chain=chain,
                    meta=rj.get('meta', {}))
    pred = predict(data, model, result.pars())

    figdir = os.path.join(outdir, 'figures')
    os.makedirs(figdir, exist_ok=True)
    report.profile_figure(os.path.join(figdir, 'profile.png'), data)
    report.bins_figure(os.path.join(figdir, 'bins.png'), data)
    report.fit_figure(os.path.join(figdir, 'fit.png'), data, pred, model)
    report.joint_figure(os.path.join(figdir, 'joint.png'), data, pred, model)
    report.posterior_figure(os.path.join(figdir, 'posterior.png'), data,
                            model, result)
    report.corner_figure(os.path.join(figdir, 'corner.png'), result)
    print('redrew figures in %s' % figdir)
    return 0
