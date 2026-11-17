# PyGINVC
PyGINVC is a python program used for fault slip inversion from geodetic observations (GNSS, InSAR, leveling) and forward simulation. It supports two kinds of fault geometeies, triangular and rectangular patches. It also can invert for the optimal fault geometry using a uniform rectangular patch.

## Download and install
Run the following commands on the Terimal, set the environments for PYTHONPATH and PATH.

    git clone git@github.com:cugzhaobin/PyGINVC.git

    cd PyGINVC

    pip install -e .

## Compile the fortran code.

    cd PyGINVC/pyginvc/Greens/okada

    f2py -c okadaio.f -m okadaio

    ln -s *.so ..

## Run
Edit yaml configure files for slip inversion or forward simulation

    SlipInversion.py --cfgfile config.yaml

### yaml template file is in scricpts.

## Mixed rectangular/triangular inversion

The unified runner inverts any combination of rectangular and triangular
faults in one solve.  Every fault entry has an `enable` switch, so a fault can
stay in the configuration without taking part in the inversion.

    pyginvc unified-slip --cfgfile config/unified_cfg.yaml
    pyginvc unified-slip --cfgfile config/unified_cfg.yaml --outpath ./run01

Results are written to `solutions.h5` with an explicit per-fault layout.
Use the unified viewer for both new and old solution files:

    pyginvc unified-view --h5file run01/solutions.h5 --plot all
