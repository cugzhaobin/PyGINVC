#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
Figures for the flat ProfileData.

    map        raw sites + profile box + fault trace
    units      InSAR incidence / heading / up component
    profile    value against fault-normal distance
    bins       binned InSAR values and pixel counts
    fit        observations, model and residuals
    joint      all datasets on the same normal axis
    posterior  200 posterior predictions, traces, marginals
    corner     parameter pair grid

Every function takes an output path, returns the Figure (path=None keeps
it open for a notebook) and writes nothing when there is no applicable
data.
'''
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from .profile import lonlat_to_enu
from .inversion import predict

C = plt.get_cmap('tab10')


def _colors(n):
    return [C(i % C.N) for i in range(n)]


def _save(fig, path, show=False):
    if path:
        fig.savefig(path, dpi=150, bbox_inches='tight')
    if show:
        plt.show()
    elif path:
        plt.close(fig)
    return fig


def _groups(data):
    '''
    [(label, component, boolean mask)] in order of first appearance.
    '''
    keys = []
    for lab, comp in zip(data.label, data.component):
        k = (str(lab), str(comp))
        if k not in keys:
            keys.append(k)
    return [(lab, comp, (data.label == lab) & (data.component == comp))
            for lab, comp in keys]


def _name(lab, comp):
    return str(lab) if comp == 'los' else '%s %s' % (lab, comp)


def _is_fitted(model, comp):
    '''
    True when the model predicts this component.  ProfileData keeps
    every component, so unfitted groups are still plotted as data.
    '''
    if model is None:
        return True
    return str(comp) == 'los' or str(comp) in getattr(model, 'components',
                                                     ())


def physical_labels(result):
    return [l for l in result.labels if not l.startswith('ln_jitter')]


def map_figure(path, profile, datasets, data, title=None, show=False):
    '''
    Raw sites, the selected sites, the profile box and the fault trace.
    '''
    fig, ax = plt.subplots(figsize=(7, 7))
    colors = _colors(max(len(datasets), 1))

    for i, spec in enumerate(datasets):
        d = spec['data']
        lab = spec['label']
        c = colors[i]
        lon = lat = None
        if spec['type'] == 'gps':
            lon, lat = d.lon, d.lat
        elif getattr(d, 'format', '') == 'point':
            lon, lat = d.llh_sar[:, 1], d.llh_sar[:, 0]

        if lon is not None and len(lon):
            E, N = lonlat_to_enu(lon, lat, profile.origin, profile.frame_name)
            inside = profile.contains(lon, lat)
            ax.plot(E[~inside], N[~inside], 'x', color='0.8', ms=3, mew=.6)
            ax.plot(E[inside], N[inside], 'o', color=c, ms=4,
                    label='%s (%d/%d)' % (lab, inside.sum(), len(inside)))
        else:
            m = data.label == lab
            if m.any():
                E, N = lonlat_to_enu(data.lon[m], data.lat[m],
                                     profile.origin, profile.frame_name)
                ax.plot(E, N, '.', color=c, ms=4,
                        label='%s (%d)' % (lab, m.sum()))

    lon, lat = profile.box()
    E, N = lonlat_to_enu(lon, lat, profile.origin, profile.frame_name)
    ax.plot(E, N, '-', color='k', lw=1.2, zorder=5)
    lon, lat = profile.trace()
    E, N = lonlat_to_enu(lon, lat, profile.origin, profile.frame_name)
    ax.plot(E, N, '-', color='r', lw=2, zorder=4, label='fault trace')
    ax.plot(0., 0., 'k+', ms=10, mew=1.5, zorder=6)

    h = max(abs(profile.normal_bounds[0]), abs(profile.normal_bounds[1]),
            max(abs(profile.along_bounds[0]),
                abs(profile.along_bounds[1])), 10.)
    ax.set_xlim(-1.1*h, 1.1*h)
    ax.set_ylim(-1.1*h, 1.1*h)
    ax.set_aspect('equal')
    ax.set_xlabel('East (km)')
    ax.set_ylabel('North (km)')
    ax.set_title(title or profile.name)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=8, loc='best')
    fig.tight_layout()
    return _save(fig, path, show)


def units_figure(path, data, show=False):
    '''
    InSAR incidence, heading and up component derived from cp/cn/cu.
    '''
    groups = [g for g in _groups(data) if g[1] == 'los']
    if not groups:
        return None
    fig, axes = plt.subplots(3, 1, figsize=(7, 6), sharex=True)
    colors = _colors(len(groups))
    phi = np.radians(data.strike)
    s, c = np.sin(phi), np.cos(phi)
    for (lab, comp, m), col in zip(groups, colors):
        cp, cn, cu = data.cp[m], data.cn[m], data.cu[m]
        ue = cp*s + cn*c
        un = cp*c - cn*s
        inc = np.degrees(np.arccos(np.clip(cu, -1., 1.)))
        head = np.degrees(np.arctan2(-un, ue)) % 360.
        for ax, y in ((axes[0], inc), (axes[1], head), (axes[2], cu)):
            ax.plot(data.normal[m], y, '.', ms=3, color=col, label=lab)
    for ax, ylab in zip(axes, ('incidence (deg)', 'heading (deg)',
                               'up component')):
        ax.set_ylabel(ylab)
        ax.grid(alpha=.25)
    axes[-1].set_xlabel('Distance from the fault (km)')
    axes[0].legend(fontsize=8)
    fig.suptitle('LOS unit vector', fontsize=10)
    fig.tight_layout()
    return _save(fig, path, show)


def profile_figure(path, data, show=False):
    groups = _groups(data)
    if not groups:
        return None
    fig = plt.figure(figsize=(9.5, 2.4 + 2.0*len(groups)),
                     constrained_layout=True)
    gs = fig.add_gridspec(len(groups), 2, width_ratios=[2.2, 1])
    for i, (lab, comp, m) in enumerate(groups):
        ax = fig.add_subplot(gs[i, 0])
        ax.errorbar(data.normal[m], data.value[m], data.sigma[m],
                    fmt='o', ms=3.5, lw=.8, capsize=2, color='C0')
        ax.axvline(0., color='0.6', lw=.8, ls=':')
        ax.axhline(0., color='0.85', lw=.8)
        ax.grid(alpha=.25)
        ax.set_ylabel('%s (mm/yr)' % _name(lab, comp))

        ax = fig.add_subplot(gs[i, 1])
        sc = ax.scatter(data.along[m], data.normal[m], c=data.value[m],
                        s=8, cmap='RdBu_r')
        ax.axhline(0., color='0.6', lw=.8, ls=':')
        ax.set_ylabel('normal (km)')
        ax.set_xlabel('along (km)')
        fig.colorbar(sc, ax=ax, pad=.02)
    fig.suptitle('profile data', fontsize=10)
    return _save(fig, path, show)


def bins_figure(path, data, show=False):
    groups = [g for g in _groups(data) if np.any(data.count[g[2]] > 1)]
    if not groups:
        return None
    fig, ax = plt.subplots(figsize=(8, 4.5))
    colors = _colors(len(groups))
    for (lab, comp, m), col in zip(groups, colors):
        ax.errorbar(data.normal[m], data.value[m], data.sigma[m],
                    fmt='o', ms=4, lw=.8, capsize=2, color=col,
                    label=_name(lab, comp))
    ax.set_xlabel('Distance from the fault (km)')
    ax.set_ylabel('LOS (mm/yr)')
    ax2 = ax.twinx()
    for (lab, comp, m), col in zip(groups, colors):
        n = np.sort(data.normal[m])
        d = np.diff(n)
        w = float(np.median(d)) if len(d) else 1.
        ax2.bar(data.normal[m], data.count[m], width=.8*w, alpha=.15,
                color=col)
    ax2.set_ylabel('pixels per bin')
    ax2.grid(False)
    ax.set_title('binned InSAR profile', fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    return _save(fig, path, show)


def fit_figure(path, data, pred, model=None, show=False):
    groups = _groups(data)
    if not groups:
        return None
    fig = plt.figure(figsize=(9.5, 2.8*len(groups)), constrained_layout=True)
    gs = fig.add_gridspec(len(groups), 2, width_ratios=[2.2, 1])
    for i, (lab, comp, m) in enumerate(groups):
        x = data.normal[m]
        fitted = _is_fitted(model, comp)
        ax = fig.add_subplot(gs[i, 0])
        ax.errorbar(x, data.value[m], data.sigma[m], fmt='o', ms=3.5,
                    lw=.8, capsize=2, color='C0', label='observed')
        if fitted:
            order = np.argsort(x)
            ax.plot(x[order], pred[m][order], '-', color='k', lw=1.5,
                    label='model')
        else:
            ax.set_title('not fitted', fontsize=8, color='0.4', loc='left')
        ax.axvline(0., color='0.6', lw=.8, ls=':')
        ax.grid(alpha=.25)
        ax.set_ylabel('%s (mm/yr)' % _name(lab, comp))
        ax.legend(fontsize=8)

        ax = fig.add_subplot(gs[i, 1])
        if fitted:
            r = (data.value[m] - pred[m])/data.sigma[m]
            ax.hist(r[np.isfinite(r)], 30, color='C1')
        else:
            ax.text(.5, .5, 'not fitted', transform=ax.transAxes,
                    ha='center', va='center', color='0.4')
        ax.set_xlabel('(obs - model) / sigma')
        ax.set_ylabel('count')
    fig.suptitle('fit', fontsize=10)
    return _save(fig, path, show)


def joint_figure(path, data, pred, model=None, show=False):
    groups = _groups(data)
    if not groups:
        return None
    fig, axes = plt.subplots(len(groups), 1, figsize=(9.5, 2.3*len(groups)),
                             sharex=True, squeeze=False)
    colors = _colors(len(groups))
    for ax, (lab, comp, m), col in zip(axes[:, 0], groups, colors):
        x = data.normal[m]
        fitted = _is_fitted(model, comp)
        ax.errorbar(x, data.value[m], data.sigma[m], fmt='o', ms=3, lw=.7,
                    capsize=2, color=col, alpha=.75, label='observed')
        if fitted:
            order = np.argsort(x)
            ax.plot(x[order], pred[m][order], '-', color=col, lw=1.2)
        else:
            ax.set_title('not fitted', fontsize=7, color='0.4', loc='left')
        ax.axvline(0., color='0.6', lw=.8, ls=':')
        ax.grid(alpha=.25)
        ax.set_ylabel('%s (mm/yr)' % _name(lab, comp))
        ax.legend(fontsize=7)
    axes[-1, 0].set_xlabel('Distance from the fault (km)')
    fig.suptitle('joint fit', fontsize=11)
    fig.tight_layout()
    return _save(fig, path, show)


def posterior_figure(path, data, model, result, title=None, ncurves=200,
                     show=False):
    if result.draws is None:
        return None
    labs = physical_labels(result)
    ndim = max(len(labs), 1)
    fig = plt.figure(figsize=(max(6., 3.2*ndim), 9.5),
                     constrained_layout=True)
    gs = fig.add_gridspec(3, ndim, height_ratios=[1.5, 1, 1])

    ax = fig.add_subplot(gs[0, :])
    groups = _groups(data)
    colors = _colors(len(groups))
    draws = result.draws
    rng = np.random.default_rng(0)
    idx = rng.choice(len(draws), min(ncurves, len(draws)), replace=False)
    for (lab, comp, m), col in zip(groups, colors):
        fitted = _is_fitted(model, comp)
        label = _name(lab, comp) + ('' if fitted else ' (not fitted)')
        ax.errorbar(data.normal[m], data.value[m], data.sigma[m],
                    fmt='o', ms=3, lw=.7, capsize=2, color=col,
                    label=label)
        if fitted:
            order = np.argsort(data.normal[m])
            for k in idx:
                par = dict(result.fixed,
                           **dict(zip(result.labels, draws[k])))
                p = predict(data, model, par)
                ax.plot(data.normal[m][order], p[m][order], '-', color=col,
                        alpha=.05, lw=.6)
    ax.axvline(0., color='0.6', lw=.8, ls=':')
    ax.set_xlabel('Distance from the fault (km)')
    ax.set_ylabel('mm/yr')
    ax.set_title(title or result.meta.get('estimator', ''), fontsize=10)
    ax.legend(fontsize=8)

    for i, lab in enumerate(labs):
        j = result.labels.index(lab)
        if result.chain is not None:
            a = fig.add_subplot(gs[1, i])
            a.plot(result.chain[:, :, j], color='k', lw=.2, alpha=.35)
            a.set_ylabel(lab)
            a.set_xlabel('step')
            a.tick_params(labelsize=7)
        a = fig.add_subplot(gs[2, i])
        a.hist(draws[:, j], 35)
        st = result.stats()[lab]
        for p, ls, lw in ((st['p16'], '--', .9), (st['p50'], '-', 1.2),
                          (st['p84'], '--', .9)):
            if np.isfinite(p):
                a.axvline(p, color='C1', ls=ls, lw=lw)
        a.set_xlabel(lab)
        a.tick_params(labelsize=7)
    return _save(fig, path, show)


def corner_figure(path, result, show=False):
    labs = physical_labels(result)
    if not labs or result.draws is None:
        return None
    draws = result.draws
    if len(draws) > 3000:
        draws = draws[np.random.default_rng(0).choice(
            len(draws), 3000, replace=False)]
    nd = len(labs)
    fig, axes = plt.subplots(nd, nd, figsize=(2.2*nd, 2.2*nd),
                             squeeze=False, constrained_layout=True)
    for i, a in enumerate(labs):
        for j, b in enumerate(labs):
            ax = axes[i, j]
            xi, xj = result.labels.index(a), result.labels.index(b)
            if i == j:
                ax.hist(result.draws[:, xi], 35)
                st = result.stats()[a]
                if np.isfinite(st['p50']):
                    ax.axvline(st['p50'], color='C1', lw=1.2)
            else:
                ax.scatter(draws[:, xj], draws[:, xi], s=1, alpha=.2,
                           color='C0')
            if i == nd - 1:
                ax.set_xlabel(b, fontsize=8)
            if j == 0:
                ax.set_ylabel(a, fontsize=8)
            ax.tick_params(labelsize=7)
    return _save(fig, path, show)
