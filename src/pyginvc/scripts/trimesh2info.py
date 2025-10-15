#!/usr/bin/env python
# This script depends on PyUnicycle package
from pyginvc.Geometry.Triangle import Triangle
import sys
import numpy as np

if len(sys.argv) != 4:
    print('trimesh2info.py <vertexfile> <elementfile> <coordtype>')
    sys.exit()

vertexfile  = sys.argv[1]
elementfile = sys.argv[2]
coordtype   = sys.argv[3]

tri = Triangle(vertexfile, elementfile)

if coordtype == 'enu':
    x  = tri.vertex_enu
    v  = tri.element-1
    xc = np.zeros((len(tri.element),3))
    xc[:,0] = (x[v[:,0],0]+x[v[:,1],0]+x[v[:,2],0])/3
    xc[:,1] = (x[v[:,0],1]+x[v[:,1],1]+x[v[:,2],1])/3
    xc[:,2] = (x[v[:,0],2]+x[v[:,1],2]+x[v[:,2],2])/3
if coordtype == 'llh':
    x  = tri.vertex_llh
    v  = tri.element-1
    xc = np.zeros((len(tri.element),3))
    xc[:,0] = (x[v[:,0],0]+x[v[:,1],0]+x[v[:,2],0])/3
    xc[:,1] = (x[v[:,0],1]+x[v[:,1],1]+x[v[:,2],1])/3
    xc[:,2] = (x[v[:,0],2]+x[v[:,1],2]+x[v[:,2],2])/3

dat = np.column_stack((xc, tri.slip/1e3, tri.rake))
dat = dat[:,[1,0,2,3,4,5,6]]

if coordtype == 'enu':
    header = '  east_km    north_km    dep_km strike_slip_m dip_slip_m op_slip_m rake_deg'
if coordtype == 'llh':
    header = '  lon_deg    lat_deg     dep_km strike_slip_m dip_slip_m op_slip_m rake_deg'
np.savetxt('tri.slip', dat, fmt='%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f', header=header)
