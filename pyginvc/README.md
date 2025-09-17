# Considering the DISMODEL in matlab is hard to read and maintain, I decide to rewrite the program
# using Python language. We don't want to redefine new data structure, but we combine several signle
# functions into a model file. The new name for the program is ginvpy (geodetic inversion python)
# Main Models:
# 1. Fault.py    - handle with fault related functions, such as load fault data, subdivide the fault
#                - plot the fault, and so on.
# 2. LoadData.py - Load different geodetic data from files, including GPS, InSAR, level data

# 3. GreenFunction.py - Compute green functions using Okada model and multilayered program
#
# 4. geotools.py - some useful functions to convert coordinates from one to other


# April 21 2016
# Now I use okadaio.f from DISMODEL to compute green function, other than Okada's DC3D.

# left slip is positive
# right slip is negative

# Dec 27, 2017
We compile PSGRN/PSCMP to obtain pstress.so, which can calculate coulomb stress changes.

# https://github.com/treverhines/BVLS

# From April 2019. I decide to encapsulate the program to oject-oriented program.
