from pyginvc.Geometry.Fault import Fault
import pyginvc.libs.geotools as gt
import numpy as np
import math
def gen_edcmp_fault(filename, origin):
    '''
    Convert FaultGeom of DISMODEL program into the format of EDCMP progam
    Written By Zhao Bin, Mar 1, 2016
    Mod by Zhao Bin, Mar 3, 2016 : output to file

    Notice:
         In EDCMP x = north  and y = east
         The unit is meter

    Input file:
        file name of FaultGeom
    Output file:
        file in the format of EDCMP program
    '''


    # we first load the FaultGeom file
    flt = Fault(filename, 1, 1, False, origin=origin)
    edcmp_flt = flt.FaultGeom2EDCMP()

    return edcmp_flt


def gen_edcmp_inp(edcmp_sfault, n, e, u):
    '''
    make input file for EDCMP program when we want to calculate the stress changes at a single point
    By Zhao Bin, Institute of Seismology, CEA. jan  9 2018

    Input:
        sfault:   source fault rupture in EDCMP format, output from faultgeom2edcmp.py
        n     :   north component in km   
        e     :   east component in km   
        u     :   up component in km. downward is positive
    Output:
    '''
    # input file name for EDCMP
    infile = 'tmp.inp'
    fid = open(infile, 'w')

    # write the file
    # 1. receiver pms
    fid.write("#receiver points\n")
    fid.write('0\n')
    fid.write("1\n")
    line = "(%e, %e)\n" %(n*1000.0, e*1000.0)
    fid.write(line)

    # 2. output file
    fid.write('# output\n')
    fid.write("'./'\n")
    fid.write('0               0              1                0\n')
    line =  "'tmp.disp'    'tmp.strn'  'tmp.strss'  'tmp.tilt'\n"
    fid.write(line)

    # 3. source fault
    fid.write('# source fault\n')
    fid.write(str(len(edcmp_sfault))+"\n")
    sfault = edcmp_sfault
    for j in range(len(sfault)):
        line =  "%d %f %e %e %e %e %e %f %f %f\n" %(sfault[j,0], sfault[j,1], sfault[j,2], sfault[j,3], sfault[j,4], sfault[j,5], sfault[j,6], sfault[j,7], sfault[j,8], sfault[j,9])
        fid.write(line)

    # 4. earth model
    fid.write('# earth model\n')
    fid.write('0\n')
    line = "%e %e %e\n" %(u*1000, 3.2074E+10, 3.2704E+10)
    fid.write(line)

    # close the file
    fid.close()

def write_edcmp_inp(edcmp_sfault, pm):
    '''
    make input file for EDCMP program when we want to calculate the stress changes on a fault plane rather than
    stress changes at a certain depth
    By Zhao Bin, Institute of Seismology, CEA. Oct 10 2016

    Input:
    sfault:   source fault rupture in EDCMP format, output from faultgeom2edcmp.py
    pm:       disgeom receiver fault geometries

    Output:
    INP files for EDCMP program
    '''

    # read in rupture source in EDCMP source format
    sfault = edcmp_sfault

    # get unique depth of receiver fault
    depths = pm[2] 


    # input file name for EDCMP
    inpfile = "tmp.inp"
    # open the file
    fid = open(inpfile, 'w')
    # write the file
    # 1. receiver pms
    fid.write("#receiver points\n")
    fid.write('0\n')
    fid.write("1\n")
    line = "(%e, %e)\n" %(pm[6]*1000, pm[5]*1000)
    fid.write(line)

    # 2. output file
    fid.write('# output\n')
    fid.write("'./'\n")
    fid.write('0               0              1                0\n')
    line =  "'tmp.disp'    'tmp.strn'  'tmp.strss'  'tmp.tilt'\n"
    fid.write(line)

    # 3. source fault
    fid.write('# source fault\n')
    fid.write(str(len(sfault))+"\n")

    for j in range(len(sfault)):
        line =  "%d %f %e %e %e %e %e %f %f %f\n" %(sfault[j,0], sfault[j,1], sfault[j,2], sfault[j,3], sfault[j,4], sfault[j,5], sfault[j,6], sfault[j,7], sfault[j,8], sfault[j,9])
        fid.write(line)

    # 4. earth model
    fid.write('# earth model\n')
    fid.write('0\n')
    line = "%e %e %e\n" %(depths*1000, 3.2074E+10, 3.2704E+10)
    fid.write(line)

    # close the file
    fid.close()


def read_stress(stressfile):
    '''
    Read stress file and store variables in dict
    Written by Zhao Bin, Institute of Seismology, CEA. Dec 27, 2017

    INPUT:
        stressfile    : stressfile output by EDCMP
    OUTPUT:
        tensor        : stress tensor
        east          : east coordinate 
        north         : north coordinate

    '''
    stress    = np.genfromtxt(stressfile)
    north     = stress[:,0]
    east      = stress[:,1]
    tensor    = {}
    tensor['snn'] = stress[:,2]
    tensor['sne'] = stress[:,5]
    tensor['snz'] = stress[:,7]
    tensor['see'] = stress[:,3]
    tensor['sez'] = stress[:,6]
    tensor['szz'] = stress[:,4]

    return north, east, tensor

def read_strain(strnfile):
    '''
    Read strain file
    Written by Zhao Bin, Institute of Seismology, CEA. Jan 3, 2018
    
    INPUT:
        strnfile     : strain file output by EDCMP
    OUTPUT:
        north        : north coordinate
        east         : east coordinate
        strain       : strain tensor
    '''
    straindata  = np.genfromtxt(strnfile, comments='#')
    north       = straindata[:,0]
    east        = straindata[:,1]
    strain      = straindata[:,[3,2,4,5,7,6]]
    strain      = straindata[:,[2,3,4,5,6,7]]
    
    return north, east, strain
