#!/usr/bin/env python
# Some functions used to calcualte stress and coulomb stress changes.
# By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley Oct 5, 2016 
# Mod by Zhao Bin, Nov 10, 2016. We abandon functions from Meade. It does not work.
# Mod by Zhao Bin, Apr 12, 2018. We change the output file name to coulomb
# We rewrite the calc_coulomb function according to website: 
#           https://searchcode.com/codesearch/view/13644410/
# Mod Zhao Bin, Institute of Seismology, CEA. Dec 26, 2017.
# import lib pstress, which contain functions of calculating coulomb stress.
# Meade's program: X-EW  Y-NS  Z-DWN
# Wang's  program: X-NS  Y-EW  Z-DWN
# Barbot's program: X-NS Y-EW  Z-DWN

# WE VERFIED THE RESULTS OF CALA_COULOMB2 PERFORM EXACTLY SAME AS WANG'S RESULT
# ON DEC 28, 2017

# import libs
import numpy as np


def tprod(n1, tensor, n2):
    t = 0.0
    for i in range(3):
        for j in range(3):
            t = t + n1[i]*tensor[i,j]*n2[j]
    return t

def axb(a,b):
    c = np.zeros(3)
    c[0] = a[1]*b[2] - b[1]*a[2]
    c[1] = b[0]*a[2] - a[0]*b[2]
    c[2] = a[0]*b[1] - b[0]*a[1]

    return c


def calc_coulomb(strike, dip, rake, fric, see, snn, szz, sne, sez, snz):
    '''
    Calculate coulomb stress on a fault plane
    strike, dip and rake are in Aki and Richards convention
    +X is east and +Y is north
    We rewrite the calc_coulomb function according to website: 
             https://searchcode.com/codesearch/view/13644410/
             Author is: Syed Tabrez Ali, University of Wisconsin.
    
    INPUT:
        strike   : strike angle of fault plane
        dip      : dip angle of fault plane
        rake     : rake angle of slip
        fric     : frictional coefficient
        s??      : stress tensor in Pa
    
    OUTPUT:
        sh       : shear stress
        no       : normal stress
        cff      : coulomb stress
        
    '''
    phi    = np.deg2rad(strike)
    sdip   = np.sin(np.deg2rad(dip))
    cdip   = np.cos(np.deg2rad(dip))
    srake  = np.sin(np.deg2rad(rake))
    crake  = np.cos(np.deg2rad(rake))
    nhat   = np.array([sdip*np.cos(phi), -sdip*np.sin(phi), cdip])
    vhat   = np.array([cdip*np.cos(phi), -cdip*np.sin(phi), -sdip])
    that   = axb(nhat, vhat)
    sij    = np.array([[see, sne, sez], [sne, snn, snz], [sez,snz,szz] ])
    tnn    = tprod(nhat, sij, nhat)
    tnv    = tprod(nhat, sij, vhat)
    tvv    = tprod(vhat, sij, vhat)
    ttt    = tprod(that, sij, that)
    tnt    = tprod(nhat, sij, that)
    sh     = tnt*crake - tnv*srake
    no     = tnn
    cff    = sh+fric*no
    return -cff, -sh, -no


def calc_coulomb2(strike, dip, rake, fric, sxx, syy, szz, sxy, sxz, syz):
    '''
    Calculate Coulomb stress changes. Please refer Manual of RELAX program for detail.
    Written by Zhao Bin, Instiute of Seismology, CEA. Dec 27, 2017.
    Mod by Zhao Bin, Jan 3, 2018. correct bugs in "b"
    
    INPUT:
        strike   : strike angle of fault plane
        dip      : dip angle of fault plane
        rake     : rake angle of slip
        fric     : frictional coefficient
        s??      : stress tensor in Pa
                   x, y, z - north, east, up 
    OUTPUT:
        tc       : Coulomb stress in Pa 
        tn       : normal stress in Pa
        tss      : shear stress in Pa

    '''
    from numpy import sin, cos
    rstr   = np.deg2rad(strike)
    rdip   = np.deg2rad(dip)
    rrak   = np.deg2rad(rake)
    n      = np.array([-sin(rstr)*sin(rdip), cos(rstr)*sin(rdip), -cos(rdip)])
    s      = np.array([cos(rstr), sin(rstr), 0])
    d      = np.array([sin(rstr)*cos(rdip), -cos(rstr)*cos(rdip), -sin(rdip)])
    b      = np.array([cos(rrak)*s[0]+sin(rrak)*d[0], cos(rrak)*s[1]+sin(rrak)*d[1], cos(rrak)*s[2]+sin(rrak)*d[2]])
    t1     = sxx*n[0] + sxy*n[1] + sxz*n[2]
    t2     = sxy*n[0] + syy*n[1] + syz*n[2]
    t3     = sxz*n[0] + syz*n[1] + szz*n[2]

    tn     = t1*n[0] + t2*n[1] + t3*n[2]
    ts     = t1*s[0] + t2*s[1]
    td     = t1*d[0] + t2*d[1] + t3*d[2]
    tss    = t1*b[0] + t2*b[1] + t3*b[2]
    tc     = tss + fric*tn
#    tc     = tss - fric*tn

    return tc, tss, tn


def strain2stress(strain, lamb, mu):
    '''
    convert strain to stress
    Written by Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley, Nov 10, 2016
    Mod by Zhao Bin, Dec 26, 2017. data type of strain and stress are set as array
    
    Please refer subsroutine isotropicstressstrain of elastic3d.f90 in RELAX program
    Please refer route of edcoutput.f in EDCMP program
    It has been tested on Jan 3, 2018 by Zhao Bin
    

    Input:
        strain    : array [exx, eyy, ezz, exy, eyz, ezx]
                    array [eee, enn, ezz, ene, enz, eez]
                    x is east; y is north
        lambda    : Lame constant in pascal
        mu        : Lame constant in pascal

    Output:
        stress    : array [sxx, syy, szz, sxy, syz, szx]
                          [see, snn, szz, sne, snz, sez]
    '''
    stress     = np.copy(strain)
    dilatation = strain[:,0]+strain[:,1]+strain[:,2]
    stress[:,0]  = lamb*dilatation + 2.0*mu*strain[:,0]
    stress[:,1]  = lamb*dilatation + 2.0*mu*strain[:,1]
    stress[:,2]  = lamb*dilatation + 2.0*mu*strain[:,2]
    stress[:,3]  = 2.0*mu*strain[:,3]
    stress[:,4]  = 2.0*mu*strain[:,4]
    stress[:,5]  = 2.0*mu*strain[:,5]

    return stress



def output_coulomb(lon, lat, up, tau, nml, cfs, coordtype='llh'):
# NOTICE: THE SIGN OF NORMAL  STRESS IN MEADE'S PROGRAM IS DIFFERENT FROM COULOMB PROGRAM
# NOTICE: THE SIGN OF COULOMB STRESS IN MEADE'S PROGRAM IS DIFFERENT FROM COULOMB PROGRAM
    '''
    Output shear/normal/coulomb stress in Coulomb program's format (dcff.cou)
    Written by Zhao Bin, Insitutute of Seismology, CEA. Ocb 10 2016
    Mod by Zhao Bin, Dec 26, 2017. adding input parameter coordtype
    Mod by Zhao Bin, Apr 12, 2018. change output file format using scientific format

    Input:
        lon      : longitude or east coordinate in meter
        lat      : latitude  or north coordinate in meter
        up       : downward in positive in meter
        tau      : shear stress in Pa
        nml      : normal stress in Pa
        cfs      : coulomb stress in Pa
    Output:
        coulomb.out
        unit of stress is bar = 1e5 Pa
    '''
    # open file
    fid = open('coulomb.out', 'w')
    
    if coordtype == 'neu':
        lon = lon/1e3
        lat = lat/1e3
        # write header
        fid.write('#Coulomb stress changes\n')
        fid.write('# EW      NS      UP      Coulomb   Shear    Normal\n')
        fid.write('# (km)    (km)    (km)    (bar)     (bar)    (bar)\n')       
    elif coordtype == 'llh':
        # write header
        fid.write('#Coulomb stress changes\n')
        fid.write('# Lon     Lat      UP      Coulomb   Shear    Normal\n')
        fid.write('# (deg)   (deg)    (km)    (bar)     (bar)    (bar)\n')

    # convert to bar and change ths sign
    # need tested
    cfs = cfs/1.0e5
    tau = tau/1.0e5
    nml = -nml/1.0e5

    fmt = 6*"%12.5e  "+"\n"
    # write the file
    for i in range(0,len(lon)):
        line = fmt %(lon[i], lat[i], up[i], cfs[i], tau[i], nml[i]) 
        fid.write(line)
    print(' coulomb.out file are output.')


def plot_coulomb(coufile, afsfile, coordtype='llh'):
    '''
    plot shear/normal/coulomb stress
    By Zhao Bin, Insitutute of Seismology, CEA. Ocb 10 2016
    MOD by Zhao Bin. Apr 12, 2018. change the input parameter to coufile
    MOD by Zhao Bin. July 3, 2018. Adding afshockfile

    Input:
        coufile    : coulomb stress file from Coulomb, calc_stress.py
        afshockfile: aftershock file, lon, lat, dep
        ctype      : coordinate type
    Output:

    '''
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    import os

    dat = np.genfromtxt(coufile, comments='#')
    lon = dat[:,0]
    lat = dat[:,1]
    cfs = dat[:,3]
    tau = dat[:,4]
    nml = dat[:,5]



    if (coordtype == "neu"):
        x = lon/1.0e3
        y = lat/1.0e3
    else:
        x = lon
        y = lat
    xi = np.linspace(np.min(x), np.max(x), 200)
    yi = np.linspace(np.min(y), np.max(y), 200)
    X,Y= np.meshgrid(xi, yi)

    # shear stress
    plt.figure()
    zi = griddata((x, y), tau, (X, Y), method='linear')
    plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
    plt.contourf(xi, yi, zi, 15, cmap=plt.cm.jet,
                     norm=plt.Normalize(vmax=abs(zi).max(), vmin=-abs(zi).max()))
    cb = plt.colorbar() 
    cb.set_label('Shear Stress(bar)', rotation=90)
    if os.path.isfile(afsfile) is True:
        afs = np.genfromtxt(afsfile)
        plt.scatter(afs[:,0], afs[:,1], s=5)
    plt.savefig('shear_stress.pdf', format='pdf')
    plt.savefig('shear_stress.jpg', format='jpg')

    # normal stress
    plt.figure()
    zi = griddata((x, y), nml, (X, Y), method='linear')
    plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
    plt.contourf(xi, yi, zi, 15, cmap=plt.cm.seismic,
                     norm=plt.Normalize(vmax=abs(zi).max(), vmin=-abs(zi).max()))
    cb = plt.colorbar() 
    cb.set_label('Normal Stress(bar)', rotation=90)
    if os.path.isfile(afsfile) is True:
        afs = np.genfromtxt(afsfile)
        plt.scatter(afs[:,0], afs[:,1], s=5)
    plt.savefig('normal_stress.pdf', format='pdf')

    # coulomb stress
    plt.figure()
    zi = griddata((x, y), cfs, (X, Y), method='linear')
    plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
    plt.contourf(xi, yi, zi, 15, cmap=plt.cm.seismic,
                     norm=plt.Normalize(vmax=abs(zi).max(), vmin=-abs(zi).max()))
    cb = plt.colorbar()
    cb.set_label('Coulomb Stress(bar)', rotation=90)
    if os.path.isfile(afsfile) is True:
        afs = np.genfromtxt(afsfile)
        plt.scatter(afs[:,0], afs[:,1], s=5)
    plt.savefig('coulomb_stress.pdf', format='pdf')


def plot_coulomb_2d(coufile, stresstype, maxstress=5):
    '''
    Plot coulomb/shear/normal stress in 2D.
    Written by Zhao Bin, Institute of Seismology, CEA. Jan 11 2018
    MOD     by Zhao Bin, April 12, 2018. fix debug of read file

    INPUT:
        coufile   : coulomb file output from convert2faultplane
                    length  width  coulomb shear normal
        stresstype: coulomb/shear/normal
    OUTPUT:
        figure
    '''
    from matplotlib import pyplot as plt

    dat = np.genfromtxt(coufile, comments='#')
    
    if stresstype == 'coulomb' or stresstype == 'Coulomb':
        sdat = dat[:,3]
    if stresstype == 'shear' or stresstype == 'Shear':
        sdat = dat[:,4]
    if stresstype == 'normal' or stresstype == 'Normal':
        sdat = dat[:,5]
    
    indx = np.where(sdat > maxstress)[0]
    sdat[indx] = maxstress
    indx = np.where(sdat < -maxstress)[0]
    sdat[indx] = -maxstress
 

    plt.scatter(dat[:,0], dat[:,1], c=sdat, cmap='jet', s=50)
    plt.xlabel('Distance along strike(km)')
    plt.ylabel('Distance along dip(km)')
    cb = plt.colorbar()
    cb.set_label(stresstype+' stress (bar)')

    plt.show()



def plot_coulomb_3d(coufile, stresstype, maxstress=150):
    '''
    Plot coulomb/shear/normal stress in 3D.
    Written by Zhao Bin, Institute of Seismology, CEA. Jan 10 2018

    INPUT:
        coufile   : coulomb file in the format of Coulomb software
        stresstype: coulomb/shear/normal
    OUTPUT:
        figure
    '''
    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    dat = np.genfromtxt(coufile, comments='#')

    if stresstype == 'coulomb' or stresstype == 'Coulomb':
        sdat = dat[:,3]*100
    if stresstype == 'shear' or stresstype == 'Shear':
        sdat = dat[:,4]*100
    if stresstype == 'normal' or stresstype == 'Normal':
        sdat = dat[:,5]*100

    indx = np.where(sdat > maxstress)[0]
    sdat[indx] = maxstress
    indx = np.where(sdat < -maxstress)[0]
    sdat[indx] = -maxstress

    p = ax.scatter(dat[:,0], dat[:,1], -dat[:,2], c=sdat, cmap='jet', s=50)
    ax.set_xlabel('east(km)')
    ax.set_ylabel('north(km)')
    ax.set_zlabel('vertical(km)')
    cb=fig.colorbar(p)
    if stresstype == 'coulomb' or stresstype == 'Coulomb':
        cb.set_label('Coulomb stress (kPa)')
    if stresstype == 'shear' or stresstype == 'Shear':
        cb.set_label('Shear stress (kPa)')
    if stresstype == 'normal' or stresstype == 'Normal':
        cb.set_label('Normal stress (kPa)')

    plt.show()

def convert2local(coufile, lon0, lat0):
    '''
    Convert LLH coordinate system of coulombs file to local ENU coordinate system
    Written by Zhao Bin, Institute of Seismology, CEA. June 30, 2018

    Input:
        coufile    : coulomb file output from calc_stress.py
        lon0,lat0  : origin
    '''
    import geotools as gt
    cou    = np.genfromtxt(coufile, comments='#')
    for i in range(len(cou)):
        cou[i,[0,1]] = gt.llh2utm(cou[i,[1,0]], [lat0, lon0])
    header = ' East     North      UP      Coulomb   Shear    Normal'
    np.savetxt('coulomb_utm.out', cou, fmt=6*'%12.3e', header=header)

    
def convert2faultplane(coufile, strike, dip):
    '''
    Convert ENU coordinate system of coulomb file to fault plane coordinate system
    Written by Zhao Bin, Institute of Seismology, CEA. Jan 11, 2018
    
    INPUT:
        coufile    : coulomb file output from calc_coulomb_edcmp.py/Coulomb/calc_stress.py
                     E N U coulomb shear normal
        strike     : strike angle of fault plane
        dip        : dip angle of fault plane
    OUTPUT:
        coufile    : in fault plane coordinate system
    '''
    
    dat    = np.genfromtxt(coufile)
    east   = dat[:,0]
    north  = dat[:,1]
    dep    = dat[:,2]
    rstrk  = np.deg2rad(strike)
    rdip   = np.deg2rad(dip)
    en     = np.vstack((east, north)).T
    RS     = np.array([[np.cos(rstrk), -np.sin(rstrk)],[np.sin(rstrk), np.cos(rstrk)]])
    lw     = en.dot(RS.T)  # leng-width

    w      = dep/np.sin(rdip)
    xy     = np.vstack((lw[:,1], w)).T
    stress = dat[:,[2,3,4,5]]
    newdat = np.hstack((xy, stress))
    header = ' Leng  Depth  UP   Coulomb Shear  Normal\n km    km     km   bar     bar    bar'
    np.savetxt('coulomb_profile.out', newdat, fmt=6*"%12.3e", header=header)

if __name__ == '__main__':
    plot_coulomb_3d('/home/zhao/earthquake/2016Aktao/coulomb/stress.out', 'shear')
