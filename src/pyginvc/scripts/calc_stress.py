#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Calculate stress tensor and coulomb stress changes at points using EDCMP program
Written by Zhao Bin, Institute of Seismology, CEA.
Created on Tue Jan  9 08:07:48 2018
Mod by Zhao Bin, Feb. 6, 2018. debug for input calfile
Mod by Zhao Bin, June 24, 2018. remove temple files

"""


from pyginvc.pscmp import pstress
from pyginvc.libs import geotools as gt
from pyginvc.libs import edcmp 
from pyginvc.libs import stress as strss
import subprocess
import numpy as np
import os, sys



def faultgeom2origin(faultfile):
    '''
    get origin of FaultGeom
    INPUT:
        faultfile  : FaultGeom file in format of GINVPY/DISMODEL
    OUTPUT:
        origin     : array [lat, lon]
    '''
    geom    = np.genfromtxt(faultfile)
    origin  = np.array([np.mean(geom[:,[3,5]]), np.mean(geom[:,[4,6]])])
    return origin

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(' calc_stress.py <faultgeom> <calfile>')
        print(' calfile : lon lat dep(km) strike dip rake friction')
        print('            dep:   positive in down-ward direction ')
        sys.exit()
    
    cosfault = sys.argv[1]
    calfile  = sys.argv[2]
    
    # source fault origin
    origin = faultgeom2origin(cosfault)
    sfault = edcmp.gen_edcmp_fault(cosfault, origin)
    
    # receiver fault/points
    llh    = np.genfromtxt(calfile, usecols=[0,1,2])
    strike = np.genfromtxt(calfile, usecols=[3])
    dip    = np.genfromtxt(calfile, usecols=[4])
    rake   = np.genfromtxt(calfile, usecols=[5])
    fric   = np.genfromtxt(calfile, usecols=[6])
    
    comb_stress = np.zeros((len(llh), 9))
    coulomb     = np.zeros(len(llh))
    shear       = np.zeros(len(llh))
    normal      = np.zeros(len(llh))
    for i in range(len(llh)):
        print("Calculating %d point ..." %(i+1))
        en = gt.llh2localxy([llh[i,1], llh[i,0]] , origin)
        edcmp.gen_edcmp_inp(sfault, en[1], en[0], llh[i,2])
#       os.popen('echo tmp.inp | edcmp')
        subprocess.call("echo tmp.inp | edcmp", shell=True)
        stress = np.genfromtxt('tmp.strss', comments='#') 
        os.popen('cat tmp.strss >> comb.strss')
        comb_stress[i] = np.array([stress[0], stress[1], stress[2], stress[3], stress[4], stress[5], stress[6], stress[7], llh[i,2]])
        snn  = stress[2]
        sne  = stress[5]
        snz  = stress[7]
        see  = stress[3]
        sez  = stress[6]
        szz  = stress[4]
        cou,sh,no = pstress.cmbfix(snn, see, szz, sne, sez, snz, 0.0, fric[i], strike[i], dip[i], rake[i])
        shear[i]   = sh
        normal[i]  = -no
        coulomb[i] = cou
        
    strss.output_coulomb(llh[:,0], llh[:,1], llh[:,2], shear, normal, coulomb, coordtype='llh')
    
    np.savetxt('stress.dat', comb_stress, fmt = 9*"%15.5e")
#   if os.path.isfile('tmp.inp'):
#       os.remove('tmp.inp')
    if os.path.isfile('tmp.strss'):
        os.remove('tmp.strss')
