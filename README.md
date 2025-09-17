# PyGINVC
PyGINVC is a python program used for fault slip inversion from geodetic observations (GNSS, InSAR, leveling) and forward simulation. It supports two kinds of fault geometeies, triangular and rectangular patches. It also can invert for the optimal fault geometry using a uniform rectangular patch.

## To get started:
Run the following commands on the Terimal, set the environments for PYTHONPATH and PATH.
export PYTHONPATH="/the/path/of/the/code/PyGINVC/:$PYTHONPATH"
export PATH="/the/path/of/the/code/PyGINVC/pyginvc/bin:$PATH"

Compile the fortran code.
cd /the/path/of/the/code/PyGINVC/pyginvc/Greens/okada
f2py -c okadaio.f -m okadaio
ln -s *.so ..

Edit yaml configure files for slip inversion or forward simulation


