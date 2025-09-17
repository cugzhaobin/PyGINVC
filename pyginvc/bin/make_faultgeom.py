#!/usr/bin/env python
# Make FaultGeom file based on focal mechanism
# Written by Zhao Bin, Institute of Seismology, CEA. Oct. 26, 2017
# MOD by Zhao Bin, Institute of Seismology, CEA. Nov. 20, 2017
# MOD by Zhao Bin, Institute of Seismology, CEA. Feb. 24, 2018. Adding parameter eqdep for the depth of epicenter

import sys
import numpy as np
import geotools as gt

def gen_fautgeom(origin, eqdep, leng, width, topdep, strike, dip, ss, ds, op):
    '''
    Generate faultgeom element using input parameter
    Written by Zhao Bin, Institute of Seismology, CEA. Nov 20, 2017
    MOD     by Zhao Bin, Institute of Seismology, CEA. Feb 24, 2018. Adding parameter eqdep for the depth of epicenter
    IN:
        origin    : [lat, lon] epicenter of earthquake
        eqdep     : depth of earthquake, downward is positive
        leng      : length of fault in km
        width     : width of fauly in km
        topdep    : depth of top of fault plane in km
        strike    : strike angle of fault
        dip       : dip angle of fault
        ss        : strike slip in mm
        ds        : dip slip in mm
        op        : openning slip in mm
    OUT:
        geom      : array of faultgeom format
    '''
    rstrike = np.deg2rad(strike)
    deltaN  = eqdep/np.tan(np.deg2rad(dip))
    ddn     = deltaN * np.cos(rstrike+np.pi/2)
    dde     = deltaN * np.sin(rstrike+np.pi/2)

    dn      = leng * np.cos(rstrike)/2.0 - ddn
    de      = leng * np.sin(rstrike)/2.0 - dde
    elatlon = gt.xy2llh([dn, de], origin)

    rstrike = np.deg2rad(strike+180)
    dn      = leng * np.cos(rstrike)/2.0 - ddn
    de      = leng * np.sin(rstrike)/2.0 - dde
    slatlon = gt.xy2llh([dn, de], origin)
    geom    = np.array([width, topdep, dip+180, slatlon[0,0], slatlon[0,1], elatlon[0,0], elatlon[0,1], ss, ds, op])
    
    return geom

if len(sys.argv) != 10:
    print(sys.argv[0] + ' <eqlon> <eqlat> <eqdep> <top-dep> <strike> <dip> <leng> <wid> <Mw>')
    print('   deg   deg  km  deg deg  km km Mw')
    sys.exit()
    

eqlon   = float(sys.argv[1])
eqlat   = float(sys.argv[2])
eqdep   = float(sys.argv[3])
topdep  = float(sys.argv[4])
strike  = float(sys.argv[5])
dip     = float(sys.argv[6])
leng    = float(sys.argv[7])
wid     = float(sys.argv[8])
mw      = float(sys.argv[9])
Mo      = gt.mw2mo(mw)
slip    = Mo/leng/wid/1e3/3e10
origin  = [eqlat, eqlon]
geom    = gen_fautgeom(origin, eqdep, leng, wid, topdep, strike, dip, slip, slip, 0.0)

fmt     = 7*'%10.4f'+3*'%12.4f'
print("#    width   top-dep       dip     s-lat     s-lon     e-lat     e-lon    s-slip        d-slip     opening")
print(fmt %(geom[0], geom[1], geom[2], geom[3], geom[4], geom[5], geom[6], geom[7], geom[8], geom[9]))
