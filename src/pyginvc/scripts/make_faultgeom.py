#!/usr/bin/env python
# Make FaultGeom file based on focal mechanism
# Written by Zhao Bin, Institute of Seismology, CEA. Oct. 26, 2017
# MOD by Zhao Bin, Institute of Seismology, CEA. Nov. 20, 2017
# MOD by Zhao Bin, Institute of Seismology, CEA. Feb. 24, 2018. Adding parameter eqdep for the depth of epicenter

import sys
import argparse
import numpy as np
from pyginvc.libs import geotools as gt

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

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Make a single-patch FaultGeom from a focal mechanism.')
    parser.add_argument('eqlon',  type=float, help='epicenter longitude (deg)')
    parser.add_argument('eqlat',  type=float, help='epicenter latitude (deg)')
    parser.add_argument('eqdep',  type=float, help='earthquake depth (km)')
    parser.add_argument('topdep', type=float, help='top of fault plane (km)')
    parser.add_argument('strike', type=float, help='strike (deg)')
    parser.add_argument('dip',    type=float, help='dip (deg)')
    parser.add_argument('leng',   type=float, help='fault length (km)')
    parser.add_argument('wid',    type=float, help='fault width (km)')
    parser.add_argument('mw',     type=float, help='moment magnitude Mw')
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    Mo     = gt.mw2mo(args.mw)
    slip   = Mo/args.leng/args.wid/1e3/3e10
    origin = [args.eqlat, args.eqlon]
    geom   = gen_fautgeom(origin, args.eqdep, args.leng, args.wid,
                          args.topdep, args.strike, args.dip,
                          slip, slip, 0.0)

    fmt = 7*'%10.4f'+3*'%12.4f'
    print("#    width   top-dep       dip     s-lat     s-lon     e-lat     e-lon"
          "    s-slip        d-slip     opening")
    print(fmt % (geom[0], geom[1], geom[2], geom[3], geom[4],
                 geom[5], geom[6], geom[7], geom[8], geom[9]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
