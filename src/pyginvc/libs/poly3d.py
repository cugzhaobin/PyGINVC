#!/usr/bin/env python
# Written by Zhao Bin, Institute of Seismology, CEA. Nov 8, 2017
# Modified by Zhao Bin, Oct. 26, 2020
# This contain several functions to manipulate POLY3D input and ouput file

# import libs
import numpy as np
import sys, math
import pyginvc.libs.geotools as gt

def gen_poly3d_input(vertex, element, bc, bctype=[], vcs='FaultCS', ecs='elocal'):
    '''
    Make input file for POLY3D using vertex and displacement element
    Written by Zhao Bin, Institute of Seismolgy, CEA. Nov 13, 2017

    IN:
        vertex       : list or array
        element      : lisr or array
        bc           : [ds, ss, op]
        bctype       : ttt/bbb
        vcs          : vertex coordicate system
        ecs          : element coordinate system
    OUT:
        poly3d.in

    '''

    # open file for writing
    with open('poly3d.in', 'w') as fid:

        # write header of POLY3D input file
        write_header(fid, 2000)
        # write vertex of POLY3D input file
        write_vertex(fid, vertex, vcs=vcs)
        # write element of POLY3D input file
        write_element(fid, element, bc, bctype=bctype, ecs=ecs)


def write_header(fid, maxslip):
    '''
    Write header for input file of POLY3D
    Written by Zhao Bin, Institute of Seismology, CEA. Nov 13, 2017

    IN:
        fid    : file ID
        maxslip: maximum slip in mm
    OUT:
        file
    '''

    fid.write('*******************************************************************************\n')
    fid.write('*                            Section 1: CONSTANTS                             *\n')
    fid.write('*******************************************************************************\n')
    fid.write('* Titles\n')
    fid.write('*--------\n')
    fid.write('title1 = "Nepal EQ Stress Driven Afterslip"\n')
    line    = 'title2 = "dip slip more than %5.2f meters"\n' %(maxslip/1000.0)
    fid.write(line)
    fid.write('\n')
    fid.write('* Elastic Constants\n')
    fid.write('* Specify any two. Leave the rest blank.\n')
    fid.write('*----------------------------------------\n')
    fid.write('shear_mod         = 30000.0\n')
    fid.write('psn_ratio         = 0.25\n')
    fid.write('youngs_mod        =\n')
    fid.write('bulk_mod          =\n')
    fid.write('lame_lambda       =\n')
    fid.write('\n')
    fid.write('* Remote Stresses/Strains\n')
    fid.write('*-------------------------\n')
    fid.write('rem_bc_type       = stress           *(stress/strain)\n')
    fid.write('s11r              = 0.0\n')
    fid.write('s22r              = 0.0\n')
    fid.write('s33r              = 0.0\n')
    fid.write('s12r              = 0.0\n')
    fid.write('s13r              = 0.0\n')
    fid.write('s23r              = 0.0\n')
    fid.write('\n')
    fid.write('* Options\n')
    fid.write('*---------\n')
    fid.write('half_space        = yes               *(yes/no)\n')
    fid.write('check_cond_num    = yes               *(yes/no)\n')
    fid.write('print_elt_geom    = no                *(yes/no)\n')
    fid.write('elt_geom_csys     = global\n')
    fid.write('null_value        = -999.0\n')
    fid.write('\n')
    fid.write('end *(CONSTANTS)\n')
	# write section two
    fid.write('*******************************************************************************\n')
    fid.write('*                        Section 2: USER COORDINATE SYSTEMS                   *\n')
    fid.write('*******************************************************************************\n')
    fid.write('* (1)        (2)     (3)    (4)    (5)     (6)     (7)     (8)       (9)\n')
    fid.write('*name      parent    x1o    x2o    x3o    rot1    rot2    rot3    rot order\n')
    fid.write('*--------- --------- ------ ------ ------ ------- ------- ------- ---------\n')
    fid.write('FaultCS   global    0.0    0.0    0.0    0.0001     00.0    0.0     321\n')
    fid.write('\n')
    fid.write('end *(USER COORDINATE SYSTEMS)\n')
    fid.write('\n')
    # write section three
    fid.write('*******************************************************************************\n')
    fid.write('*                         Section 3: OBSERVATION GRIDS                        *\n')
    fid.write('*******************************************************************************\n')
    fid.write('* (1)    (2) (3)      (4)        (5)        (6)     (7)    (8)    (9)    (10)   (11)   (12) (13) (14) (15)\n')
    fid.write('*name    dim outp  endpt csys obspt csys outp csys x1beg  x2beg  x3beg  x1end  x2end  x3end  nx1 nx2 nx3\n')
    fid.write('*------- --- ----- ---------- ---------- --------- ------ ------ ------ ------ ------ ------ --- --- ---\n')
    fid.write('ObsGrid  0    sd    global     global     global    -50.5  -20.5  0.0\n')
    fid.write('ObsPlane 2     d    global     global     global    -300     -300      0 300     300    0 30   30  1\n')
    fid.write('\n')
    fid.write('end *(OBSERVATION GRIDS)\n')
    # write section four
    fid.write('*******************************************************************************\n')
    fid.write('*                      Section 4: OBJECTS/ELEMENTS/VERTICES                   *\n')
    fid.write('*                      (o = object, e = element, v = vertex)                  *\n')
    fid.write('*******************************************************************************\n')
    fid.write('\n')
    fid.write('*(1) (2)            (3)        (4)\n')
    fid.write('*o name             outp    eltc csys\n')
    fid.write('*- ---------------- ------- ----------\n')
    fid.write('o "MFT_NEPAL"          bt      FaultCS\n')
    fid.write('\n')
    fid.write('*(1) (2)     (3)       (4)     (5)     (6)\n')
    fid.write('*v name     csys       x1      x2      x3\n')
    fid.write('*- -------- ---------  ------  ------  ------\n')


def write_vertex(fid, vertex, vcs='FaultCS', sid=0):
    '''
    Write vertex section into POLY3D input file
    Written by Zhao Bin, Institute of Seismology, CEA. Nov. 13, 2017
    Modified by Zhao Bin, Dec. 22, 2020

    IN
        fid      : file ID for poly3d.in
        vertex   : vertex of elements
        vcs      : coordinate system of vertex
        sid      : start index of vertex
    '''

    # write vertex
    for i in range(len(vertex)):
        line = "v %4d %s %10.3f %10.3f %10.3f\n" %(i+1+sid, vcs, vertex[i,0], vertex[i,1], vertex[i,2] )
        fid.write(line)


def write_element(fid, element, bc, bctype=[], ecs='elocal'):
    '''
    Write element section into POLY3D input file
    Written by Zhao Bin, Institute of Seismology, CEA. Nov. 13, 2017

    Mod by Zhao Bin, Mar. 24, 2020. adding another parameter bc

    IN
        fid         : file ID for poly3d.in
        element     : indx of vertex
        bc          : boundary condition of displacement, [ds, ss, op]
        bctype      : ttt/bbb
        ecs         : elocal/FaultCS
    '''

    fid.write('\n')
    fid.write('*(1) (2)   (3)      (4)    (5)       (6)       (7)      (8) (9) (10) (...)\n')
    fid.write('*e #vert BC csys   BC type BC1       BC2       BC3       v1 v2 v3 ...\n')
    fid.write(' - ----- --------- ------- --------- --------- --------- -- -- -- --\n')
    newline = ""

    
    if len(bctype) == 0:
        for i in range(len(bc)):
            bctype.append('bbb')

    if len(element) != len(bc):
        print('write_element: make sure element and disp_element are with same length')
        sys.exit()

    for i in range(len(bc)):
        if element.shape[1] == 4:
            newline =  'e  4  %s  %s  %10.2f  %10.2f  %10.2f  %d %d %d %d\n' \
                    %(ecs, bctype[i], bc[i,0], bc[i,1], bc[i,2], \
                    element[i,0], element[i,1], element[i,2], element[i,3])
        elif element.shape[1] == 3:
            newline =  'e  3  %s  %s  %10.2f  %10.2f  %10.2f  %d %d %d\n' \
                    %(ecs, bctype[i], bc[i,0], bc[i,1], bc[i,2], \
                    element[i,0], element[i,1], element[i,2])
        fid.write(newline)
    fid.write("end *(OBJECTS/ELEMENTS/VERTICES)")


def read_element_vertex(infile):
    '''
    read element and vertex info from POLY3D input file
    Modified by Zhao Bin, Institute of Seismology. Nov. 13, 2017
    IN
        infile      : POLY3D input file
    OUT
        vertex      : vertex array,  [x1, x2, x3]
        element     : element array, [element index]
        bctype      : element array, [element index]
        bc          : element array, [element index]
    '''

    # read all lines into content
    with open(infile, 'r') as fid:
        content = fid.readlines()

    vertex  = []
    element = []
    bc      = []
    bctype  = []
    belem   = False
    bvert   = False
    for i in range(1, len(content)):
        line =  content[i]
        if line[0:3] == 'end':
            belem = False
        if line[0:2] == '*e':
            belem = True
            bvert = False
        if line[0:2] == '*v':
            bvert = True
        if line[0:1] == "v" and bvert:
            x1    = float(line.split()[3])
            x2    = float(line.split()[4])
            x3    = float(line.split()[5])
            vertex.append([x1, x2, x3])
        if line[0:1] == 'e' and belem:
            nelem = int(line.split()[1])
            bct   = line.split()[3]
            bc1   = float(line.split()[4])
            bc2   = float(line.split()[5])
            bc3   = float(line.split()[6])
            bctype.append(bct)
            bc.append([bc1, bc2, bc3])
            elem  = []
            for i in range(nelem):
                elem.append(int(line.split()[7+i]))
            element.append(elem)

    vertex  = np.array(vertex)
    element = np.array(element)
    bc      = np.array(bc)   
    return vertex, element, bctype, bc

def read_disp_element(outfile):
    '''
    read displacement for each element
    right now this function can not process the case of multi-object

    Input:
    outfile   : output file of POLY3D
    Output:
    element_disp: list of displacement
    '''
    with open(outfile, 'r') as fid:
        content = fid.readlines()

    element_disp = []
    bobj   = False
    bdisp  = False
    bdash  = False
    for i in range(len(content)):
        line = content[i]
        if line[14:20] == "OBJECT":
            bobj  = True
        elif line[0:8] == "STRESS":
            bdisp = False
        elif line[0:4] == "DISP":
            bdisp = True
        elif line[0:1] == "-" and bdisp:
            bdash = True
            continue
        elif bobj and bdisp and bdash:
            if len(line.split())<1 and len(element_disp)>0:
                bobj, bdisp, bdash = False, False, False
                continue
            dipslip    =  float(line.split()[4])
            strikeslip =  float(line.split()[7])
            openslip   =  float(line.split()[10])
            element_disp.append([dipslip, strikeslip, openslip])

    element_disp = np.array(element_disp)
    return element_disp



def read_stress_element(outfile):
    '''
    read stress for each element
    right now this function can not process the case of multi-object

    Input:
    outfile   : output file of POLY3D
    Output:
    element_strss: list of stress
    '''
    with open(outfile, 'r') as fid:
        content = fid.readlines()

    element_strss = []
    bobj   = False
    bdash  = False
    bstress = False
    for i in range(len(content)):
        line = content[i]
        if line[14:20] == "OBJECT":
            bobj  = True
        elif line[0:8] == "STRESSES":
            bstress = True
        elif line[0:1] == "-" and bstress:
            bdash = True
            continue
        elif bobj and bstress and bdash:
            if len(line.split())<1 and len(element_strss)>0:
                bobj, bstress, bdash = False, False, False
                continue
            dipstrss    =  float(line.split()[4])
            strikestrss =  float(line.split()[5])
            openstrss   =  float(line.split()[6])
            element_strss.append([dipstrss, strikestrss, openstrss])
    
    # convert to array
    element_strss = np.array(element_strss)
    return element_strss


def element_depth(vertex, element):
    deplist = []
    deplist.append(vertex[int(element[7])-1][2])
    deplist.append(vertex[int(element[8])-1][2])
    deplist.append(vertex[int(element[9])-1][2])
    deplist.append(vertex[int(element[10])-1][2])

    return np.min(deplist)


#######################################################################################
# copy from faultgeom2node_element.py                                                 #
#######################################################################################
def unique2D(a):
    '''
    unique the row records of a 2D array
    Input:
    a      - 2D array
    '''
    order  = np.lexsort(a.T)
    a      = a[order]
    diff   = np.diff(a, axis=0)
    ui     = np.ones(len(a), 'bool')
    ui[1:] = (diff != 0).any(axis=1) 
    return a[ui]

def distance(x1, y1, z1, x2, y2, z2):
    '''
    calculate distance between two points
    '''
    dx = x1-x2
    dy = y1-y2
    dz = z1-z2
    dist = np.sqrt(dx**2+dy**2+dz**2)
    return dist


def mergelist(depth, tolerance):
    '''
    get uniqe list 'depth' when inter-distance between two records is less than the tolerence
    depth     - one dimension array
    tolerance - value
    '''
    new = []
    for i in range(len(depth)):
        for j in range(i+1,len(depth)):
            if np.abs(depth[i]-depth[j]) < tolerance:
                depth[j] = depth[i]
    new = np.unique(depth)
    return new

def merge_array(vertex, tolerance):
    '''
    merge vertex togeather if the distance between two points is less than the tolerance
    INPUT:
        vertax    - 2D array, contain lon, lat and depth
        tolerance - distance of two points, unit of km
    OUTPUT:
        newvertex - 2D array, contain lon, lat and depth
    '''

    for i in range(len(vertex)):
        for j in range(i+1, len(vertex)):
            dist = gt.distance3d([vertex[i,0], vertex[i,1], vertex[i,2]], [vertex[j,0], vertex[j,1], vertex[j,2]])
            if dist < tolerance:
                vertex[j,0] = vertex[i,0]
                vertex[j,1] = vertex[i,1]
                vertex[j,2] = vertex[i,2]

    newvertex = unique2D(vertex)

    return newvertex


def uniq_vertex(vertex):
    '''
    get unique vertexes of finite fault elements
    Input:
    vertex   - 2D array, contain lon, lat and depth
    '''
    # get unique depth
    depth = np.unique(vertex[:,2])
    depth = mergelist(depth, 0.1)

    # for each depth
    uniqv = []
    for i in range(len(depth)):
        indx = np.where((np.abs(vertex[:,2]-depth[i])<1))
        hvertex = vertex[indx]
        hvertex = merge_array(hvertex, 5)
        uniqv.append(hvertex)
#       print '***************************'
#       print len(hvertex), depth[i]
        
    # merge the vertexes at certain depths together
    tmp = uniqv[0]
    for i in range(1, len(uniqv)):
        tmp = np.vstack((tmp, uniqv[i]))

    newvertex =  tmp
    return newvertex

def getindx(node, vertex):
    dist = np.zeros(len(vertex))
    for i in range(len(vertex)):
        dist[i] = distance(node[0], node[1], node[2], vertex[i,0], vertex[i,1], vertex[i,2])
    indx = np.where(dist==dist.min())[0][0]
    return indx+1


def vertex_element_from_faultgeom(faultfile, origin):
    '''
    Read vertex and element from FaultGeom file
    Written by Zhao Bin, Institute of Seismology, CEA. Nov. 14, 2017
    IN:
        faultfile : GINPY file
    OUT:
        vertex    :
        elem      :
    '''
    from pyginvc.Geometry.Fault import Fault
    flt               = Fault(faultfile, 1, 1, False, origin=origin)
    felem             = flt.felem
    vertex            = np.zeros((4*len(felem), 3))
    elements          = []
    patch             = []
    disp_element      = felem[:,[5,4,6]] # ds, ss, opa # unit is mm

    for i in range(len(felem)):
        vertex[4*i,:]   = felem[i,21:24]
        vertex[4*i+1,:] = felem[i,24:27]
        vertex[4*i+2,:] = felem[i,27:30]
        vertex[4*i+3,:] = felem[i,30:33]

        patch.append(felem[i,21:24])
        patch.append(felem[i,24:27])
        patch.append(felem[i,27:30])
        patch.append(felem[i,30:33])
        elements.append(patch)
        patch = []

    # get unique vertexs
    newvertex = uniq_vertex(vertex)

    elem = np.zeros((len(felem),4))
    for i in range(len(elements)):
        patch = elements[i] 
        node1 = patch[0]
        indx1 = getindx(node1, newvertex)
        node2 = patch[1]
        indx2 = getindx(node2, newvertex)
        node3 = patch[2]
        indx3 = getindx(node3, newvertex)
        node4 = patch[3]
        indx4 = getindx(node4, newvertex)
        elem[i,:] = np.array([indx1, indx2, indx3, indx4])

    return newvertex, elem, disp_element
    
    
def write_gmt_disp(vertex, element, disp_element, disp_type="ds"):
    '''
    Write POLY3D output to GMT format
    Written by Zhao Bin, Institute of Seismology, CEA. Nov. 14, 2017
    IN:
        vertex      : narray, [m, 3], east, north, depth
        element     : narray, [m, 4], index of vertex for rectangular
        disp_element: narray, [m, 3], ds, ss, op
        disp_type   : ds/ss/op
    OUT:
        disp.gmtlin
    '''
    
    outfile = "disp_"+disp_type+".gmtlin"
    fmt     = 3*"%10.3f "+"\n"
    with open(outfile, 'w') as fid:
        for i in range(len(element)):
            if disp_type == 'ds':
                line = '> -Z%f\n' %(disp_element[i,0])
            elif disp_type == 'ss':
                line = '> -Z%f\n' %(disp_element[i,1])
            elif disp_type == 'op':
                line = '> -Z%f\n' %(np.sqrt(disp_element[i,0]**2 + disp_element[i,1]**2))
            fid.write(line)

            elem = element[i]
            for j in range(len(elem)):
                index = elem[j]-1
                line  = fmt %(vertex[index, 0], vertex[index, 1], vertex[index, 2])
                fid.write(line)


def write_gmt_stress(vertex, element, stress_element, stress_type="ds", mu=0.4, rake=0, scale=1.0):
    '''
    Write POLY3D output to GMT format
    Written by Zhao Bin, Institute of Seismology, CEA. Nov. 14, 2017
    Modified by Zhao Bin, Oct. 26, 2020, output Coulomb stress
    IN:
        vertex        : narray, [m, 3], east, north, depth
        element       : narray, [m, 4], index of vertex for rectangular
        stress_element: narray, [m, 3], ds, ss, op
        stress_type   : ds/ss/op/ts/cfs
        mu            : frictional coefficient
        rake          : rake angle
    OUT:
        disp.gmtlin
    '''
    
    # convert stress to MPa
    stress_element = stress_element*scale
    outfile = "stress_"+stress_type+".gmtlin"
    fmt     = 3*"%10.3f "+"\n"
    rake    = np.deg2rad(rake)
    with open(outfile, 'w') as fid:
        for i in range(len(element)):
            if stress_type == 'ds':
                line = '> -Z%f\n' %(stress_element[i,0])
            elif stress_type == 'ss':
                line = '> -Z%f\n' %(stress_element[i,1])
            elif stress_type == 'op':
                line = '> -Z%f\n' %(stress_element[i,2])
            elif stress_type == 'ts':
                line = '> -Z%f\n' %(np.linalg.norm(stress_element[i,0:2]))
            elif stress_type == 'cfs':
                cfs  = stress_element[i,0]*np.cos(rake) +\
                       stress_element[i,1]*np.sin(rake) -\
                       stress_element[i,2]*mu
                line = '> -Z%f\n' %(cfs)
            fid.write(line)

            elem = element[i]
            for j in range(len(elem)):
                index = elem[j]-1
                line  = fmt %(vertex[index, 0], vertex[index, 1], vertex[index, 2])
                fid.write(line)


def read_disp_surface(outfile):
    '''
    read displacement for each element
    right now this function can not process the case of multi-object

    Input:
    outfile   : output file of POLY3D
    Output:
    element_disp: list of displacement
    '''

    with open(outfile, 'r') as fid:
        content = fid.readlines()
    bdisp = False
    disp  = []
    for i in range(len(content)):
        if content[i].find('X1       X2       X3         U1         U2         U3')>0:
            bdisp = True
        if content[i].find('-') >0 and bdisp:
            if len(content[i].split())==6:
                dat  = [float(j) for j in content[i].split()]
                disp.append(dat)
            else:
                bdisp = False
    disp = np.array(disp)
    return disp


def gen_faultgeom(infile, outfile, origin):
    '''
    Output fault slip model in FaultGeom format
    
    Input:
        infile   = poly3d.in
        outfile  = poly3d.out
        origin   = [lat, lon]
    '''

    # read vertex, element, bctype and bc from input file
    vertex, elements, bctype, bc = read_element_vertex(infile)

    # read displacement on each element from output file from POLY3D
    elem_slip = read_disp_element(outfile)
    for i in range(len(elements)):
        patch = elements[i]-1
        v1    = vertex[patch[0]]
        v2    = vertex[patch[1]]
        v3    = vertex[patch[2]]
        v4    = vertex[patch[3]]
        width1  = np.sqrt((v1[0]-v4[0])**2 + (v1[1]-v4[1])**2 + (v1[2]-v4[2])**2)
        width2  = np.sqrt((v2[0]-v3[0])**2 + (v2[1]-v3[1])**2 + (v2[2]-v3[2])**2)
        width   = (width1+width2)/2

        # dip angle of the fault
        dh   = v3[2]-v2[2]
        dip  = np.rad2deg(math.asin(dh/width))
        dip  = np.abs(dip)+180

        llh1 = gt.xy2llh([v1[1], v1[0]], origin)
        llh2 = gt.xy2llh([v2[1], v2[0]], origin)

        fmt  = 3*"%10.3f"+4*"%10.4f"+3*"%15.3f"
        print(fmt %(width, np.abs(v1[2]), dip, llh1[0,0], llh1[0,1], llh2[0,0], llh2[0,1], elem_slip[i,0], elem_slip[i,1], elem_slip[i,2]))

def write_trielement(element, data, dtype='disp'):
    '''
    Write displacement/stress to .mesh file
    ss ds op
    '''
    if element.shape[0] != data.shape[0]:
        print(' check the dimension')
        sys.exit()
    else:
        ofile= dtype+'.mesh'
        with open(ofile, 'w') as fid:
            for i in range(len(element)):
                if len(element[i]) == 3:
                    fid.write("{:5d} {:5d} {:5d} {:15.3f} {:15.3f} {:15.3f}\n".format(
                        element[i][0], element[i][1], element[i][2], data[i,1], data[i,0], data[i,2]))
                elif len(element[i]) == 4:
                    fid.write("{:5d} {:5d} {:5d} {:5d} {:15.3f} {:15.3f} {:15.3f}\n".format(
                        element[i][0], element[i][1], element[i][2], element[i][3], data[i,1], data[i,0], data[i,2]))


def vertex_enu2llh(vertex, origin):
    '''
    Convert vertex coordinate from enu format to llh
    Written by Zhao Bin, Dec. 22, 2020

    Input:
        vertex = [e, n, u]
        origin = [lat, lon]
    Output:
        enu
    '''

    llh = np.copy(vertex)
    for i in range(len(llh)):
        llh[i,0:2] = gt.utm2llh(vertex[i,[1,0]], origin)

    return llh



def ExportVTK(vertex, element, slip, traction, disp):
    '''
    Write the fault geometry and slip distribution to VTK file for Paraview.
    '''

    nf  = len(element)
    npt = element.shape[1]
    with open('poly3d.vtp', 'w') as fid:
        fid.writelines(["<?xml version=\"1.0\"?>\n",
                   "<VTKFile type=\"PolyData\" version=\"0.1\">\n",
                   "  <PolyData>\n"])

        for i in range(0,nf):
            enu = vertex[element[i]-1]
            fid.write("    <Piece NumberOfPoints=\"{}\" NumberOfPolys=\"1\">\n".format(npt))
            fid.write("      <Points>\n")
            fid.write("        <DataArray type=\"Float32\" Name=\"Fault Patch\" NumberOfComponents=\"3\" format=\"ascii\">\n")
            for j in range(element.shape[1]):
                fid.write("         {} {} {}\n".format(enu[j,0], enu[j,1], enu[j,2]))
            fid.write("         </DataArray>\n")
            fid.write("      </Points>\n")
            fid.write("      <Polys>\n")
            fid.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"{}\">\n".format(nf-1))
            if npt == 3:
                fid.write("0 1 2\n")
            elif npt == 4:
                fid.write("0 1 2 3\n")
            else:
                pass
            fid.write("        </DataArray>\n")
            fid.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"{}\" RangeMax=\"{}\">\n".format(npt, npt))
            fid.write("          {}\n".format(npt))
            fid.write("        </DataArray>\n")
            fid.write("      </Polys>\n")
            fid.write("      <CellData Scalar=\"geometry\">\n")
            fid.write("        <DataArray type=\"Float32\" Name=\"dip slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
            fid.write("{}\n".format(slip[i,0]))
            fid.write("        </DataArray>\n")
            fid.write("        <DataArray type=\"Float32\" Name=\"strike slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
            fid.write("{}\n".format(slip[i,1]))
            fid.write("        </DataArray>\n")
            fid.write("        <DataArray type=\"Float32\" Name=\"slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
            fid.write("{}\n".format(np.sqrt(slip[i,0]**2+slip[i,1]**2)))
            fid.write("        </DataArray>\n")
            fid.write("        <DataArray type=\"Float32\" Name=\"dip traction\" NumberOfComponents=\"1\" format=\"ascii\">\n")
            fid.write("{}\n".format(traction[i,1]))
            fid.write("        </DataArray>\n")
            fid.write("        <DataArray type=\"Float32\" Name=\"strike traction\" NumberOfComponents=\"1\" format=\"ascii\">\n")
            fid.write("{}\n".format(traction[i,0]))
            fid.write("        </DataArray>\n")
            fid.write("      </CellData>\n")
            fid.write("    </Piece>\n")
        fid.write("  </PolyData>\n")
        fid.write("</VTKFile>\n")
