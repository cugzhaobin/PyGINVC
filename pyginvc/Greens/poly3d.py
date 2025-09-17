#!/usr/bin/env python
# By Zhao Bin, Institute of Seismology, CEA. Nov 8, 2017
# This contain several functions to manipulate POLY3D input and ouput file

# import libs
import numpy as np
import Fault as flt
import sys, math
import pyginvc.libs.geotools as gt

def gen_poly3d_input(vertex, element, bc, bctype=[], vcs='FaultCS', ecs='elocal', obspt=[],
                     shear_mod=3e4):
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
    fid  = open('poly3d.in', 'w')

    if shear_mod > 1e8:
        shear_mod = shear_mod/1e6
    # write header of POLY3D input file
    write_header(fid, 2000, obspt=obspt, shear_mod=shear_mod)
    # write vertex of POLY3D input file
    write_vertex(fid, vertex, vcs=vcs)
    # write element of POLY3D input file
    write_element(fid, element, bc, bctype=bctype, ecs=ecs)

    fid.close()



def write_header(fid, maxslip, obspt=[], shear_mod=3e4):
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
    fid.write('shear_mod         = {}\n'.format(shear_mod))
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
    for i in range(len(obspt)):
        fid.write('ObsGrid   0   sd   global  global  global {} {} {}\n'.format(obspt[i,0], obspt[i,1], 0.0))
#   fid.write('ObsGrid  0    sd    global     global     global    -50.5  -20.5  0.0\n')
#   fid.write('ObsPlane 2     d    global     global     global    -100  -200   0 300  200 0 30   30  1\n')
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


def write_vertex(fid, vertex, vcs='FaultCS'):
    '''
    Write vertex section into POLY3D input file
    Written by Zhao Bin, Institute of Seismology, CEA. Nov. 13, 2017

    IN
        fid      : file ID for poly3d.in
        vertex   : vertex of elements
    '''

    # write vertex
    for i in range(len(vertex)):
        line = "v %4d %s %10.3f %10.3f %10.3f\n" %(i+1, vcs, vertex[i,0], vertex[i,1], vertex[i,2] )
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
            bctype.append('ttt')

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
                    %(ecs, bctype, bc[i,0], bc[i,1], bc[i,2], \
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
    fid = open(infile, 'r')

    # read all lines into content
    content = fid.readlines()

    # close the file
    fid.close()

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
    fid = open(outfile, 'r')
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
                break 
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
    fid = open(outfile, 'r')
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
                break
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
        indx = np.where((np.abs(vertex[:,2]-depth[i])<0.01))
        hvertex = vertex[indx]
        hvertex = merge_array(hvertex, 1.5)
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


def vertex_element_from_faultgeom(faultfile):
    '''
    Read vertex and element from FaultGeom file
    Written by Zhao Bin, Institute of Seismology, CEA. Nov. 14, 2017
    IN:
        faultfile : GINPY file
    OUT:
        vertex    :
        elem      :
    '''
    geom              = flt.LoadFaultGeom(faultfile)
    [disgeom, origin] = flt.SetOrigin(geom)
    felem             = flt.FaultGeom2AllVertex(geom, disgeom)
    vertex            = np.zeros((4*len(felem), 3))
    elements          = []
    patch             = []
    disp_element      = felem[:,[5,4,6]] # ds, ss, op

    for i in range(len(felem)):
        vertex[4*i,:]   = felem[i,9:12]
        vertex[4*i+1,:] = felem[i,12:15]
        vertex[4*i+2,:] = felem[i,15:18]
        vertex[4*i+3,:] = felem[i,18:21]

        patch.append(felem[i,9:12])
        patch.append(felem[i,12:15])
        patch.append(felem[i,15:18])
        patch.append(felem[i,18:21])
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

    return newvertex, elem, disp_element, origin
    
    
def write_gmt_disp(vertex, element, disp_element, disp_type="ds"):
    '''
    Write POLY3D output to GMT format
    Written by Zhao Bin, Institute of Seismology, CEA. Nov. 14, 2017
    IN:
        vertex      : narray, [m, 3], east, north, depth
        element     : narray, [m, 4], index of vertex for rectangular
        disp_element: narray, [m, 3], ds, ss, op
        disp_type   : ds/ss/ts
    OUT:
        disp.gmtlin
    '''
    
    outfile = "disp_"+disp_type+".gmtlin"
    fid     = open(outfile, 'w')
    fmt     = 3*"%10.3f "+"\n"
    for i in range(len(element)):
        if disp_type == 'ds':
            line = '> -Z%f\n' %(disp_element[i,0])
        elif disp_type == 'ss':
            line = '> -Z%f\n' %(disp_element[i,1])
        elif disp_type == 'ts':
            line = '> -Z%f\n' %(np.sqrt(disp_element[i,0]**2 + disp_element[i,1]**2))
        fid.write(line)

        elem = element[i]
        for j in range(len(elem)):
            index = elem[j]-1
            line  = fmt %(vertex[index, 0], vertex[index, 1], vertex[index, 2])
            fid.write(line)


def write_gmt_stress(vertex, element, stress_element, stress_type="ds"):
    '''
    Write POLY3D output to GMT format
    Written by Zhao Bin, Institute of Seismology, CEA. Nov. 14, 2017
    IN:
        vertex        : narray, [m, 3], east, north, depth
        element       : narray, [m, 4], index of vertex for rectangular
        stress_element: narray, [m, 3], ds, ss, op
        stress_type   : ds/ss/ts
    OUT:
        disp.gmtlin
    '''
    
    # convert stress to MPa
    stress_element = stress_element
    outfile = "stress_"+stress_type+".gmtlin"
    fid     = open(outfile, 'w')
    fmt     = 3*"%10.3f "+"\n"
    for i in range(len(element)):
        if stress_type == 'ds':
            line = '> -Z%f\n' %(stress_element[i,0])
        elif stress_type == 'ss':
            line = '> -Z%f\n' %(stress_element[i,1])
        elif stress_type == 'ts':
            line = '> -Z%f\n' %(np.sqrt(stress_element[i,0]**2 + stress_element[i,1]**2))
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
    fid = open(outfile, 'r')
    content = fid.readlines()
    bdisp = False
    disp  = []
    for i in range(len(content)):
        if content[i].find('DISPLACEMENTS') != -1:
            bdisp = True
        if len(content[i].split())==6:
            if content[i].find('e') >0 and bdisp:
                dat  = [float(j) for j in content[i].split()]
                disp.append(dat)
    disp = np.array(disp)
    disp = disp[:,[3,4,5]]
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
        length1 = np.sqrt((v1[0] - v2[0])**2 + (v1[1] - v2[1])**2)
        length2 = np.sqrt((v3[0] - v4[0])**2 + (v3[1] - v4[1])**2)
        length  = (length1+length2)/2
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
