#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA.
'''
One catalog row per profile.

    ProfileBatch.py --cfgfile example/config_batch.yaml
'''
import argparse
import sys

from pyginvc.Profile import run
from pyginvc.Profile.config import ConfigError


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='batch cross-fault profiles from a catalog')
    ap.add_argument('--cfgfile', required=True, help='YAML configuration')
    args = ap.parse_args(argv)
    try:
        return run.batch(args.cfgfile)
    except (ConfigError, ValueError, OSError, ImportError) as e:
        print('error: %s' % e, file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
