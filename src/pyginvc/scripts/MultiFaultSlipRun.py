#!/usr/bin/env python
"""Deprecated entry point for the legacy multi-fault inversion script.

This script used the pre-2020 module-level API (``GreenFunction``,
``Inverse``, ``LoadData``, ``Output``), which no longer exists.  The
supported interface is the YAML-driven command line tool::

    python MultiFaultSlipInversion.py --cfgfile config/multifault_cfg.yaml

The full legacy source is available in the git history.
"""

import sys

__all__ = ['main']


def main():
    sys.stderr.write(
        'MultiFaultSlipRun.py is deprecated and no longer functional.\n'
        'Use MultiFaultSlipInversion.py with a YAML config instead, e.g.\n'
        '  python MultiFaultSlipInversion.py --cfgfile <config.yaml>\n')
    return 2


if __name__ == '__main__':
    sys.exit(main())
