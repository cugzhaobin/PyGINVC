#/usr/bin/env pyhton
# estimate variance-covarnance matrix from InSAR interferogram
# Written by Zhao Bin, Institute of Seismology, Jul. 30, 2019.
import h5py
import numpy as np
# Please install GIAnT first
from tsinsar import tscov

XJorg = h5py.File('XJorg.mat')
lon   = XJorg['XJorg']['lonv'].value
lat   = XJorg['XJorg']['latv'].value

XJ    = h5py.File('XJ.mat')
infA  = XJ['XJ']['infA'].value
infD  = XJ['XJ']['infD'].value

idx1  = np.where(np.logical_and(lon>74.65, lon<75.0))[0]
idx2  = np.where(np.logical_and(lat>39.00, lat<39.5))[0]
infA  = infA[idx1][:,idx2]

par,qual,x,ecv = tscov.phi2covar(infA, 0.5, 0.1, plot=True)
print par
