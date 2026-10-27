#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
Least-squares and MCMC estimation for one flat ProfileData.

    pred = cp*v_par + cn*v_norm + cu*v_up + bias @ par
    r    = (value - pred) / sigma

    fit(data, model, free, fixed, method='leastsq'|'mcmc')

ProfileData keeps every component.  `fit_data` masks the rows a model
cannot predict (a strike-slip model does not fit GPS norm rows) before
estimation, so the cut stays model-independent.  No estimator classes
and no orchestration layer: `parameter_space` merges model defaults and
dataset bias bounds, the two `fit_*` functions do the work, and both
return the same Result.
'''
import sys
from dataclasses import dataclass, field

import numpy as np

from .profile import ProfileData

_EMCEE_HINT = ('MCMC needs emcee: pip install pysliprate[mcmc] '
               'or pip install emcee')


def require_emcee():
    try:
        import emcee
    except ImportError:
        raise ImportError(_EMCEE_HINT)
    return emcee


@dataclass
class Result(object):
    '''
        labels      free parameter names (jitter included for MCMC)
        best        best-fit value of every label + fixed values
        fixed       parameters held fixed
        cov         linearized covariance (LSQ), else None
        draws       (nsample, ndim) posterior samples
        chain       (nsteps, nwalkers, ndim) MCMC chain, physical units
        stats()     mean/std/percentiles for every label
        pars()      parameter dict for predict()
    '''
    labels:     list
    best:       dict
    fixed:      dict = field(default_factory=dict)
    cov:        np.ndarray = None
    draws:      np.ndarray = None
    chain:      np.ndarray = None
    lnprob:     np.ndarray = None
    tau:        np.ndarray = None
    acceptance: float = None
    meta:       dict = field(default_factory=dict)

    def stats(self):
        out = {}
        if self.draws is not None:
            for i, lab in enumerate(self.labels):
                c = self.draws[:, i]
                out[lab] = {'mean': float(c.mean()),
                            'std':  float(c.std()),
                            'p16':  float(np.percentile(c, 16)),
                            'p50':  float(np.percentile(c, 50)),
                            'p84':  float(np.percentile(c, 84))}
            return out

        sd = (np.sqrt(np.diag(self.cov))
              if self.cov is not None else np.full(len(self.labels), np.nan))
        for i, lab in enumerate(self.labels):
            out[lab] = {'mean': float(self.best[lab]),
                        'std':  float(sd[i]),
                        'p16':  float('nan'),
                        'p50':  float(self.best[lab]),
                        'p84':  float('nan')}
        return out

    def pars(self):
        s = self.stats()
        p = {k: float(v) for k, v in self.fixed.items()}
        for lab in self.labels:
            v = s[lab]['p50']
            p[lab] = float(self.best[lab]) if not np.isfinite(v) else v
        return p


def predict(data, model, par):
    '''
    Model prediction for every row of data, bias included.

    A bias column that only touches components the model does not fit
    may be absent from `par`; it contributes zero.  This lets the full
    data (all components) be passed to plotting after a component-wise
    fit, while an absent bias for a fitted row is still an error.
    '''
    vpar, vnorm, vup = model.predict(par, data)
    pred = data.cp*vpar + data.cn*vnorm + data.cu*vup
    names = data.bias_names
    if names:
        if all(n in par for n in names):
            b = np.array([par[n] for n in names], dtype=float)
        else:
            fitted = fit_components(data, model)
            b = np.zeros(len(names))
            for k, n in enumerate(names):
                if n in par:
                    b[k] = par[n]
                elif np.any(data.bias[fitted, k] != 0.):
                    raise KeyError('missing bias parameter %r' % n)
        pred = pred + data.bias @ b
    return pred


def residual(data, model, par):
    '''
    Whitened residual, one value per data row.
    '''
    return (data.value - predict(data, model, par))/data.sigma


def fit_components(data, model):
    '''
    Rows the model can predict: its active GPS components plus every
    LOS row, which mixes those components through cp/cn/cu.
    '''
    comps = set(getattr(model, 'components', ('par', 'norm')))
    return np.asarray([str(c) == 'los' or str(c) in comps
                       for c in data.component], dtype=bool)


def fit_data(data, model):
    '''
    View of `data` restricted to the rows and bias columns the model
    uses.  ProfileData itself keeps every component; only the estimator
    and the fit figures see this narrower view.
    '''
    mask = fit_components(data, model)
    if not mask.any():
        raise ValueError(
            'the model fits %s but the data contain only %s'
            % (', '.join(getattr(model, 'components', ())) or 'nothing',
               ', '.join(sorted({str(c) for c in data.component}))))
    if mask.all():
        return data
    keep = [k for k, n in enumerate(data.bias_names)
            if np.any(data.bias[mask, k] != 0.)]
    names = tuple(data.bias_names[k] for k in keep)
    bias = data.bias[mask][:, keep]
    return ProfileData(
        label=data.label[mask], component=data.component[mask],
        site=data.site[mask], lon=data.lon[mask], lat=data.lat[mask],
        E=data.E[mask], N=data.N[mask], along=data.along[mask],
        normal=data.normal[mask], value=data.value[mask],
        sigma=data.sigma[mask], cp=data.cp[mask], cn=data.cn[mask],
        cu=data.cu[mask], count=data.count[mask], bias=bias,
        bias_names=names,
        bias_bounds={n: data.bias_bounds[n] for n in names},
        strike=data.strike, weight=data.weight[mask])


def _par(labels, values, fixed):
    p = {k: float(v) for k, v in fixed.items()}
    p.update({n: float(v) for n, v in zip(labels, values)})
    return p


def parameter_space(data, model, free, fixed):
    '''
    Merge user bounds, model defaults and dataset bias bounds.

    Parameters not in `free` are held at the model default or at the
    dataset's bias value; a bias with lo == hi is fixed.
    '''
    free = {k: tuple(float(x) for x in v) for k, v in (free or {}).items()}
    fixed = {k: float(v) for k, v in (fixed or {}).items()}

    both = set(free) & set(fixed)
    if both:
        raise ValueError('%s are both free and fixed'
                         % ', '.join(sorted(both)))

    known = set(model.params)
    for n, p in model.params.items():
        if n not in free and n not in fixed:
            fixed[n] = float(p)

    bias = data.bias_bounds
    for n, b in bias.items():
        b = tuple(float(x) for x in b)
        if n in free or n in fixed:
            continue
        if b[0] == b[1]:
            fixed[n] = b[0]
        else:
            free[n] = b

    unknown = sorted((set(free) | set(fixed)) - known - set(bias))
    if unknown:
        raise ValueError('unknown parameter(s): %s; the model has %s, the '
                         'datasets have %s'
                         % (', '.join(unknown), ', '.join(model.params),
                            ', '.join(bias) or 'none'))
    for n, b in free.items():
        if len(b) != 2 or not np.isfinite(b).all():
            raise ValueError('%s: bounds must be two finite numbers, got %r'
                             % (n, b))
        if b[0] > b[1]:
            raise ValueError('%s: min=%g > max=%g' % (n, b[0], b[1]))
        if n in ('slip_ss', 'slip_ds') and b[0] < 0.:
            raise ValueError('%s is a magnitude and must have non-negative '
                             'bounds; use the model sense to choose the '
                             'direction' % n)
    for n in ('slip_ss', 'slip_ds'):
        if n in fixed and fixed[n] < 0.:
            raise ValueError('%s is a magnitude and cannot be negative' % n)
    return free, fixed


def _start_point(free, model):
    p = []
    for n, (lo, hi) in free.items():
        v = model.params.get(n, 0.5*(lo + hi))
        eps = 1e-6*(hi - lo)
        p.append(np.clip(v, lo + eps, hi - eps))
    return np.asarray(p, dtype=float)


def _flat_residual(data, model, labels, fixed):
    def fun(p):
        r = residual(data, model, _par(labels, p, fixed))
        return np.where(np.isfinite(r), r, 1e6)
    return fun


def fit_leastsq(data, model, free, fixed=None, loss='linear', f_scale=1.0,
                max_nfev=4000, ndraws=4000, seed=0):
    '''
    Bounded nonlinear least squares with a linearized covariance.
    '''
    from scipy.optimize import least_squares

    data = fit_data(data, model)
    free, fixed = parameter_space(data, model, free, fixed)
    if not free:
        raise ValueError('no free parameters to invert')
    labels = list(free)
    lo = np.array([free[n][0] for n in labels])
    hi = np.array([free[n][1] for n in labels])
    p0 = _start_point(free, model)
    fun = _flat_residual(data, model, labels, fixed)

    sol = least_squares(fun, p0, bounds=(lo, hi), loss=loss,
                        f_scale=f_scale, max_nfev=max_nfev)

    ndata = len(fun(sol.x))
    dof = ndata - len(labels)
    cov = None
    if dof > 0:
        J = sol.jac
        try:
            cov = np.linalg.inv(J.T @ J)*(2.0*sol.cost/dof)
        except np.linalg.LinAlgError:
            cov = None

    draws = None
    if cov is not None and ndraws:
        rng = np.random.default_rng(seed)
        vals, vecs = np.linalg.eigh(cov)
        L = vecs*np.sqrt(np.maximum(vals, 0.))
        draws = sol.x + rng.standard_normal((ndraws, len(labels))) @ L.T
        draws = np.clip(draws, lo, hi)

    return Result(
        labels=labels, best=_par(labels, sol.x, fixed), fixed=fixed,
        cov=cov, draws=draws,
        meta={'estimator': 'least_squares', 'loss': loss,
              'nfev': int(sol.nfev), 'cost': float(sol.cost),
              'dof': int(dof), 'ndata': int(ndata),
              'bounds': {n: list(free[n]) for n in labels}})


def _map_init(data, model, free, fixed):
    from scipy.optimize import least_squares

    labels = list(free)
    lo = np.array([free[n][0] for n in labels])
    hi = np.array([free[n][1] for n in labels])
    p0 = _start_point(free, model)
    fun = _flat_residual(data, model, labels, fixed)
    try:
        sol = least_squares(fun, p0, bounds=(lo, hi), max_nfev=2000)
        p = np.clip(sol.x, lo, hi)
    except Exception:
        p = p0

    span = hi - lo
    at_edge = (p - lo < 1e-3*span) | (hi - p < 1e-3*span)
    names = [n for n, e in zip(labels, at_edge) if e]
    if names:
        print('warning: the least-squares starting point is at a bound '
              '(%s); the MCMC will still explore the box'
              % ', '.join(names), file=sys.stderr)
    return p


def fit_mcmc(data, model, free, fixed=None, nwalkers=32, nsteps=5000,
             discard='auto', thin='auto', seed=0, jitter=True,
             jitter_bounds=(-5., 2.), progress=False):
    '''
    emcee with a uniform box prior; one ln_jitter per dataset label.
    '''
    emcee = require_emcee()

    data = fit_data(data, model)
    free, fixed = parameter_space(data, model, free, fixed)
    if not free:
        raise ValueError('no free parameters to invert')
    labels = list(free)
    lo = np.array([free[n][0] for n in labels], dtype=float)
    hi = np.array([free[n][1] for n in labels], dtype=float)
    nfree = len(labels)

    groups = np.unique(data.label)
    jlabels = (['ln_jitter.%s' % g for g in groups] if jitter else [])
    all_labels = labels + jlabels
    lo_all = np.append(lo, np.full(len(jlabels), jitter_bounds[0]))
    hi_all = np.append(hi, np.full(len(jlabels), jitter_bounds[1]))
    ndim = len(all_labels)

    if nwalkers < 2*ndim:
        raise ValueError('nwalkers=%d is too few for %d parameters, need '
                         'at least %d' % (nwalkers, ndim, 2*ndim))
    if nwalkers < 4*ndim:
        print('warning: nwalkers=%d is less than 4*ndim=%d; the chain may '
              'be poorly mixed' % (nwalkers, 4*ndim), file=sys.stderr)

    p_model = _map_init(data, model, free, fixed)
    theta_map = (p_model - lo)/np.maximum(hi - lo, np.finfo(float).tiny)
    if jitter:
        theta_map = np.append(theta_map, np.full(len(jlabels), 0.5))

    step = np.maximum(1e-3*np.abs(lo_all + (hi_all - lo_all)*theta_map),
                      1e-6)
    rng = np.random.default_rng(seed)
    theta0 = theta_map + (step/(hi_all - lo_all))*rng.standard_normal(
        (nwalkers, ndim))
    theta0 = np.clip(theta0, 1e-4, 1. - 1e-4)

    def loglike(theta):
        if np.any(theta < 0.) or np.any(theta > 1.):
            return -np.inf
        p = lo_all + (hi_all - lo_all)*theta
        par = dict(fixed, **dict(zip(labels, p[:nfree])))
        r = data.value - predict(data, model, par)
        ll = 0.0
        for k, g in enumerate(groups):
            m = data.label == g
            err2 = data.sigma[m]**2
            if jitter:
                err2 = err2 + np.exp(2.0*p[nfree + k])
            ll += -0.5*np.sum(r[m]**2/err2 + np.log(2.0*np.pi*err2))
        return ll if np.isfinite(ll) else -np.inf

    sampler = emcee.EnsembleSampler(nwalkers, ndim, loglike)
    sampler.random_state = np.random.RandomState(seed).get_state()
    sampler.run_mcmc(theta0, nsteps, progress=progress)

    try:
        tau = sampler.get_autocorr_time(tol=0)
    except emcee.autocorr.AutocorrError:
        tau = None
    taumax = float(np.max(tau)) if tau is not None else None

    if discard == 'auto':
        discard = min(int(5*taumax), nsteps//2) if taumax else nsteps//2
    else:
        discard = int(discard)
    if thin == 'auto':
        thin = max(int(0.5*taumax), 1) if taumax else 1
    else:
        thin = int(thin)
    if discard >= nsteps or (nsteps - discard)//thin < 2:
        raise ValueError('nsteps=%d is too few: no samples left after '
                         'discard=%d thin=%d, increase nsteps'
                         % (nsteps, discard, thin))

    theta_chain = sampler.get_chain()
    chain = lo_all + (hi_all - lo_all)*theta_chain
    flat  = chain[discard::thin].reshape(-1, ndim)
    lnp   = sampler.get_log_prob()
    best_i = np.unravel_index(np.argmax(lnp[discard::thin]),
                              lnp[discard::thin].shape)
    best = _par(all_labels, chain[discard::thin][best_i], fixed)

    return Result(
        labels=all_labels, best=best, fixed=fixed,
        cov=np.cov(flat.T) if len(flat) > ndim else None,
        draws=flat, chain=chain, lnprob=lnp, tau=tau,
        acceptance=float(np.mean(sampler.acceptance_fraction)),
        meta={'estimator': 'mcmc', 'nwalkers': nwalkers, 'nsteps': nsteps,
              'discard': discard, 'thin': thin, 'seed': seed,
              'ndata': int(len(data)), 'jitter': jitter,
              'jitter_labels': list(jlabels),
              'bounds': {n: list(free[n]) for n in labels}})


def fit(data, model, free, fixed=None, method='mcmc', **kw):
    '''
    Dispatch to least squares or MCMC.

    Extra keyword arguments are filtered to the chosen method, so one
    estimator block in the config can be passed unchanged.
    '''
    m = str(method).lower()
    if m in ('leastsq', 'least_squares', 'lsq'):
        allowed = ('loss', 'f_scale', 'max_nfev', 'ndraws', 'seed')
        return fit_leastsq(data, model, free, fixed,
                           **{k: v for k, v in kw.items() if k in allowed})
    if m == 'mcmc':
        allowed = ('nwalkers', 'nsteps', 'discard', 'thin', 'seed',
                   'jitter', 'jitter_bounds', 'progress')
        return fit_mcmc(data, model, free, fixed,
                        **{k: v for k, v in kw.items() if k in allowed})
    raise ValueError('method must be leastsq or mcmc, got %r' % (method,))
