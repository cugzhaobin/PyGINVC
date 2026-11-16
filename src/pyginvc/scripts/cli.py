#!/usr/bin/env python
"""Unified ``pyginvc`` command-line dispatcher.

This first version only exposes the ``faultgeom2*`` conversion commands. Each
sub-command is imported lazily, so a missing optional dependency in one tool
never affects the others.

Usage::

    pyginvc list
    pyginvc faultgeom2coulomb --help
    pyginvc faultgeom2edcmp faultgeom

The arguments after the sub-command are passed through unchanged to the
target script's ``main()``.
"""
import sys
import difflib
import argparse

# group key -> (display title, one-line description)
GROUPS = {
    'convert': ('convert', 'Format conversion of fault geometry / slip models'),
    'invert':  ('invert',  'Slip / geometry inversion'),
    'forward': ('forward', 'Forward modelling'),
    'plot':    ('plot',    'Plotting'),
    'profile': ('profile', 'Profile analysis'),
    'model':   ('model',   'Stress / other modelling'),
}

# name -> (module, callable, one-line help, group)
COMMANDS = {
    'faultgeom2coulomb': (
        'pyginvc.scripts.faultgeom2coulomb', 'main',
        'Convert FaultGeom to Coulomb input', 'convert'),
    'faultgeom2curvedfault': (
        'pyginvc.scripts.faultgeom2curvedfault', 'main',
        'Convert FaultGeom to curved-fault format', 'convert'),
    'faultgeom2curvedfault_profile': (
        'pyginvc.scripts.faultgeom2curvedfault_profile', 'main',
        'Convert FaultGeom to curved-fault profile', 'convert'),
    'faultgeom2edcmp': (
        'pyginvc.scripts.faultgeom2edcmp', 'main',
        'Convert FaultGeom to EDCMP input', 'convert'),
    'faultgeom2geojson': (
        'pyginvc.scripts.faultgeom2geojson', 'main',
        'Convert FaultGeom patches to GeoJSON', 'convert'),
    'faultgeom2gmtlin': (
        'pyginvc.scripts.faultgeom2gmtlin', 'main',
        'Convert FaultGeom to GMT line format', 'convert'),
    'faultgeom2info': (
        'pyginvc.scripts.faultgeom2info', 'main',
        'Print FaultGeom summary information', 'convert'),
    'faultgeom2moment': (
        'pyginvc.scripts.faultgeom2moment', 'main',
        'Compute seismic moment from FaultGeom', 'convert'),
    'faultgeom2planefault': (
        'pyginvc.scripts.faultgeom2planefault', 'main',
        'Convert FaultGeom to a planar fault', 'convert'),
    'faultgeom2poly3din': (
        'pyginvc.scripts.faultgeom2poly3din', 'main',
        'Convert FaultGeom to Poly3D input', 'convert'),
    'faultgeom2pscmp': (
        'pyginvc.scripts.faultgeom2pscmp', 'main',
        'Convert FaultGeom to PSCMP input', 'convert'),
    'faultgeom2relax': (
        'pyginvc.scripts.faultgeom2relax', 'main',
        'Convert FaultGeom to RELAX input', 'convert'),
    'faultgeom2visco1d': (
        'pyginvc.scripts.faultgeom2visco1d', 'main',
        'Convert FaultGeom to VISCO1D input', 'convert'),
    'faultgeom2visco2pt5d': (
        'pyginvc.scripts.faultgeom2visco2pt5d', 'main',
        'Convert FaultGeom to VISCO2.5D input', 'convert'),
    'make-faultgeom': (
        'pyginvc.scripts.make_faultgeom', 'main',
        'Build a FaultGeom from a focal mechanism', 'convert'),
    'faultgeom2vtk': (
        'pyginvc.scripts.faultgeom2vtk', 'main',
        'Convert FaultGeom to VTK for ParaView', 'convert'),

    'slip': (
        'pyginvc.scripts.SlipInversion', 'main',
        'Rectangular-patch slip inversion from YAML', 'invert'),
    'tri-slip': (
        'pyginvc.scripts.TriSlipInversion', 'main',
        'Triangular-mesh slip inversion from YAML', 'invert'),
    'multi-slip': (
        'pyginvc.scripts.MultiFaultSlipInversion', 'main',
        'Multi-fault / multi-file slip inversion from YAML', 'invert'),
    'multi-slip-legacy': (
        'pyginvc.scripts.MultiFaultSlipRun', 'main',
        'Legacy multi-fault slip inversion (hard-coded script)', 'invert'),
    'mcmc': (
        'pyginvc.scripts.MultiFaultMCMCInversion', 'main',
        'Bayesian (MCMC) slip inversion from YAML', 'invert'),
    'checkerboard': (
        'pyginvc.scripts.CheckerBoardRun', 'main',
        'Checkerboard resolution test', 'invert'),
    'resolution': (
        'pyginvc.scripts.faultgeomResolution', 'main',
        'Fault-patch resolution test', 'invert'),
    'block': (
        'pyginvc.scripts.block_inversion', 'main',
        'Block-motion inversion from GPS (with optional InSAR)', 'invert'),
    'block-backslip': (
        'pyginvc.scripts.block_backslip_inversion', 'main',
        'Block-model backslip interseismic inversion from YAML', 'invert'),
    'geom-opt': (
        'pyginvc.scripts.GeomOptRun', 'main',
        'Optimal fault geometry estimation from YAML', 'invert'),
    'block': (
        'pyginvc.scripts.block_inversion', 'main',
        'Block-model motion inversion from GPS', 'invert'),
    'block-backslip': (
        'pyginvc.scripts.block_backslip_inversion', 'main',
        'Block-model backslip (coupling) inversion from YAML', 'invert'),

    'plot-faultgeom': (
        'pyginvc.scripts.plot_faultgeom', 'main',
        'Plot fault geometry', 'plot'),
    'plot-trimesh': (
        'pyginvc.scripts.plot_trimesh', 'main',
        'Plot a triangular mesh', 'plot'),
    'plot-surface': (
        'pyginvc.scripts.plot_recdislocation_surface', 'main',
        'Plot surface displacements of a slip model', 'plot'),
    'plot-sar': (
        'pyginvc.scripts.plot_sar', 'cli',
        'Plot InSAR LOS data (multi-file, cartopy)', 'plot'),
    'view': (
        'pyginvc.scripts.InversionView', 'main',
        'View inversion results (obs vs prediction, slip)', 'plot'),

    'forward': (
        'pyginvc.scripts.ForwardRun', 'main',
        'Rectangular-patch forward modelling from YAML', 'forward'),
    'tri-forward': (
        'pyginvc.scripts.TriForwardRun', 'main',
        'Triangular-mesh forward modelling from YAML', 'forward'),
    'vol-forward': (
        'pyginvc.scripts.VolForwardRun', 'main',
        'Volcano (Mogi) forward modelling from YAML', 'forward'),

    'profile-prepare': (
        'pyginvc.scripts.profile_prepare', 'main',
        'Prepare profile data from a config file', 'profile'),
    'profile-invert': (
        'pyginvc.scripts.profile_inversion', 'main',
        'Profile (2-D) inversion', 'profile'),
    'profile-plot': (
        'pyginvc.scripts.profile_plot', 'main',
        'Plot profile observations and model', 'profile'),
    'profile-batch': (
        'pyginvc.scripts.profile_batch', 'main',
        'Batch profile processing', 'profile'),

    'stress': (
        'pyginvc.scripts.calc_stress', 'main',
        'Coulomb stress change on receiver faults', 'model'),
}


def _names_in_group(group):
    return sorted(k for k, v in COMMANDS.items() if v[3] == group)


def _print_group(group):
    title = GROUPS[group][0]
    width = max(len(k) for k in COMMANDS)
    print('%s -- %s' % (title, GROUPS[group][1]))
    for name in _names_in_group(group):
        print('  %-*s  %s' % (width, name, COMMANDS[name][2]))


def _print_list():
    width = max(len(k) for k in COMMANDS)
    for group, (title, desc) in GROUPS.items():
        names = _names_in_group(group)
        if not names:
            continue
        print('%s -- %s' % (title, desc))
        for name in names:
            print('  %-*s  %s' % (width, name, COMMANDS[name][2]))
        print()


def _usage():
    print('usage: pyginvc [<group>] <command> [args...]\n')
    _print_list()
    print('Run "pyginvc <command> --help" for the arguments of one command.')


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)

    if not argv or argv[0] in ('-h', '--help', 'help'):
        _usage()
        return 0
    if argv[0] in ('list', 'ls'):
        if len(argv) > 1 and argv[1] in GROUPS:
            _print_group(argv[1])
        else:
            _print_list()
        return 0

    # group-qualified form: "pyginvc <group> <command> ..."
    if argv[0] in GROUPS:
        if len(argv) == 1:
            _print_group(argv[0])
            return 0
        name, rest = argv[1], argv[2:]
        if name not in COMMANDS or COMMANDS[name][3] != argv[0]:
            sys.stderr.write("pyginvc: '%s' is not a command in group '%s'\n"
                             % (name, argv[0]))
            _print_group(argv[0])
            return 2
    else:
        name, rest = argv[0], argv[1:]

    if name not in COMMANDS:
        close = difflib.get_close_matches(name, COMMANDS, n=3)
        sys.stderr.write("pyginvc: unknown command '%s'\n" % name)
        if close:
            sys.stderr.write('Did you mean: %s\n' % ', '.join(close))
        sys.stderr.write('Run "pyginvc list" to see all commands.\n')
        return 2

    module_name, func_name = COMMANDS[name][0], COMMANDS[name][1]
    import importlib
    module = importlib.import_module(module_name)
    func = getattr(module, func_name)

    # Let the target script's argparse see "<command> <args...>".
    sys.argv = [name] + rest
    func()
    return 0


if __name__ == '__main__':
    sys.exit(main())
