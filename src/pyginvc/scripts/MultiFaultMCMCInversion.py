#!/usr/bin/env python
"""
    python MultiFaultMCMCInversion.py --cfgfile config/multifault_mcmc_cfg.yaml
"""

import os
import sys
import yaml
import logging
import argparse
import numpy as np
from pyginvc.Inversion.MultiFaultInversion import MultiFaultInversion
from pyginvc.Inversion.MCMCSlipInversion import MCMCSlipInversion
from pyginvc.Export import View
from pyginvc.scripts.MultiFaultSlipInversion import (
    build_origin,
    build_faults,
    build_data,
    resolve_gfiletype,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')

__all__ = ['run', 'main']


def run(cfgfile, nsteps=None, burnin=None):
    with open(cfgfile, 'r') as fid:
        cfg = yaml.safe_load(fid)

    G      = cfg.get('global', {}) or {}
    FAULTS = cfg.get('faults', []) or []
    GPSC   = cfg.get('gps', []) or []
    SARC   = cfg.get('insar', []) or []
    LEVC   = cfg.get('leveling', []) or []
    EXPC   = cfg.get('export', {}) or {}
    MCMC   = cfg.get('mcmc', {}) or {}

    if not FAULTS:
        logging.error('No faults defined in config.')
        sys.exit(1)

    # 命令行覆盖只用于快速调试；正式结果请写进 YAML
    if nsteps is not None or burnin is not None:
        sampler = dict(MCMC.get('sampler') or {})
        if nsteps is not None:
            sampler['nsteps'] = int(nsteps)
        if burnin is not None:
            sampler['burnin'] = int(burnin)
        MCMC['sampler'] = sampler

    gfiletype = resolve_gfiletype(G, GPSC)
    origin    = G.get('origin') or build_origin(FAULTS)
    if len(origin) != 2:
        origin = build_origin(FAULTS)

    faults = build_faults(FAULTS, origin)
    data   = build_data(GPSC, SARC, LEVC, gfiletype)

    inversion = MultiFaultInversion(data, faults, origin=origin,
                                    config=G, export=EXPC)
    mc        = MCMCSlipInversion(inversion, MCMC, export=EXPC)
    result    = mc.run()
    mc.write_outputs(result)

    print('\n=== MCMC posterior (lambda, Mw, Mo, acceptance) ===')
    print('lambda 16/50/84 : %.6g  %.6g  %.6g'
          % (result['lambda_p16'], result['lambda_median'],
             result['lambda_p84']))
    print('Mw 16/50/84     : %.3f  %.3f  %.3f'
          % tuple(np.percentile(result['Mw_samples'], [16, 50, 84])))
    print('Mo mean         : %.6e' % result['Mo'])
    print('acceptance      : %.3f' % result['acceptance'])
    print('MCMC outputs    : %s/' % mc.outdir)

    if EXPC.get('view', False) or MCMC.get('view', False):
        h5file = os.path.join(mc.outdir, 'mcmc_solutions.h5')
        View.plot_sar_ramp(h5file)
        View.plot_obs_mod(h5file, scale=EXPC.get('vecscale'))
        View.plot_slip_3d(
            h5file,
            depth_exaggeration=EXPC.get('depth_exaggeration', 2.0))

    return mc, result


def main():
    parser = argparse.ArgumentParser(
        description='Multi-fault / multi-file Bayesian (MCMC) slip inversion.')
    parser.add_argument('--cfgfile', type=str, required=True,
                        help='path to YAML config')
    parser.add_argument('--nsteps', type=int, default=None,
                        help='override mcmc.sampler.nsteps (debugging only)')
    parser.add_argument('--burnin', type=int, default=None,
                        help='override mcmc.sampler.burnin (debugging only)')
    args = parser.parse_args()
    run(args.cfgfile, nsteps=args.nsteps, burnin=args.burnin)


if __name__ == '__main__':
    main()
