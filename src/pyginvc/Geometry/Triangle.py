#!/usr/bin/env python
# Written by Zhao Bin, when he is at UC Berkeley. Feb 14 2016
# Mod by Zhao Bin, Dec. 7, 2018. replace print() with print

import h5py
import numpy as np
import sys, os, logging
from pyginvc.libs import geotools as gt

logging.basicConfig(
                    level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt="%d-%M-%Y %H:%M:%S")

class Triangle(object):
    '''
    Triangle is a class representing triangular fault plane
    '''

    def __init__(self, vertexfile, elementfile, origin=[]):
    
        '''
        Load vertex and elements for triangualr dislocation method
        Written by Zhao Bin, Institute of Seismology, CEA. Dec 19, 2016
        
        Input:
        vertexfile  : file name of vertex file
                    : index lat lon dep(up positive)
        elementfile : file name of elements
                    : indx1 indx2 indx3 strike-slip dip-slip open-slip
        
        Output:
        origin      : [lat0, lon0]
        vertex_enu  : local coordinates for vertex [east north] in km
        vertex_llh  : geophysical coordinates for vertex [latitutde longitude]
        element     : index of vertex for each triangular element
        '''
        # check the file is exit
        isfile = os.path.isfile(vertexfile)
        if (isfile == False):
           logging.warning ('{} is not exist! Please input another file!'.format(vertexfile))
           sys.exit()
     
        # check the file is exit
        isfile = os.path.isfile(elementfile)
        if (isfile == False):
           logging.warning('{} is not exist! Please input another file!'.format(elementfile))
           sys.exit()
        
        # read vertex, element and slip
        vertex_llh = np.genfromtxt(vertexfile, comments='#', usecols=[1,2,3])
        element    = np.genfromtxt(elementfile, comments='#', usecols = [0,1,2], dtype = 'int')
        slip       = np.genfromtxt(elementfile, comments='#', usecols = [3,4,5])
        
        # calculate the origin
        if len(origin) != 2:
            lat0   = np.mean(vertex_llh[:,0])
            lon0   = np.mean(vertex_llh[:,1])
            origin = np.array([lat0, lon0])
        
        # convert geophysical coordinate into local coordinate
        vertex_enu = np.zeros((len(vertex_llh),3))
        for i in range(len(vertex_llh)):
#           en = gt.llh2localxy(vertex_llh[i,0:2], origin)
            en = gt.llh2utm(vertex_llh[i,0:2], origin)
            vertex_enu[i,:] = np.array([en[0], en[1], vertex_llh[i,2]])

        # IMPORTANT: ENFORCE CLOCKWISE CIRCULATION.
        # ADD by Zhao Bin, Dec. 10, 2021
        tri    = vertex_enu[element-1]
        crossd = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0])
        idx    = np.where(crossd[:,2]<0)[0]
        logging.info('{} elements will be enforced to clockwise'.format(len(idx)))
        element[idx] = element[idx][:,[0,2,1]]
        
        # return results
        self.vertex_enu = vertex_enu
        self.vertex_llh = vertex_llh
        self.element    = element
        self.origin     = origin
        self.nf         = len(element)
        self.slip       = slip
        self.rake       = np.rad2deg(np.arctan2(self.slip[:,1], self.slip[:,0]))
        return

    def Moment(self, shearmodulus=3e10):
        '''
        Calaulate moment
        '''
        vertex  = self.vertex_enu
        element = self.element
        idx1    = element[:,0]-1
        idx2    = element[:,1]-1
        idx3    = element[:,2]-1
        v1      = vertex[idx1]
        v2      = vertex[idx2]
        v3      = vertex[idx3]

        vec1    = v2 - v1
        vec2    = v3 - v1
        cross   = np.cross(vec1, vec2)
        area    = 0.5 * np.linalg.norm(cross, axis=1)
        ss      = self.slip[:,0]
        ds      = self.slip[:,1]
        tslip   = np.sqrt(ss**2 + ds**2)
        Mo      = 1e6 * area * tslip * shearmodulus
        Mo_total  = np.sum(Mo)/1000.0
        Mw_total  = 2.0/3.0*np.log10(Mo_total) - 6.067
        return Mo_total, Mw_total


    def Triangle2VTK(self):
        '''
        Write the fault geometry and slip distribution to VTK file for Paraview.
        '''

        nf = len(self.element)
        with open('faultgeom.vtp', 'w') as fid:
            fid.writelines(["<?xml version=\"1.0\"?>\n",
                       "<VTKFile type=\"PolyData\" version=\"0.1\">\n",
                       "  <PolyData>\n"])

            for i in range(0,nf):
                enu = self.vertex_enu[self.element[i]-1]
                fid.write("    <Piece NumberOfPoints=\"3\" NumberOfPolys=\"1\">\n")
                fid.write("      <Points>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"Fault Patch\" NumberOfComponents=\"3\" format=\"ascii\">\n")
                fid.write("         {} {} {}\n".format(enu[0,0], enu[0,1], -enu[0,2]))
                fid.write("         {} {} {}\n".format(enu[1,0], enu[1,1], -enu[1,2]))
                fid.write("         {} {} {}\n".format(enu[2,0], enu[2,1], -enu[2,2]))
                fid.write("         </DataArray>\n")
                fid.write("      </Points>\n")
                fid.write("      <Polys>\n")
                fid.write("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"{}\">\n".format(nf-1))
                fid.write("0 1 2\n")
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"3\" RangeMax=\"3\">\n")
                fid.write("          3\n")
                fid.write("        </DataArray>\n")
                fid.write("      </Polys>\n")
                fid.write("      <CellData Scalar=\"geometry\">\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"strike slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(self.slip[i,0]))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"dip slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(self.slip[i,1]))
                fid.write("        </DataArray>\n")
                fid.write("        <DataArray type=\"Float32\" Name=\"slip\" NumberOfComponents=\"1\" format=\"ascii\">\n")
                fid.write("{}\n".format(np.sqrt(self.slip[i,0]**2+self.slip[i,1]**2)))
                fid.write("        </DataArray>\n")
                fid.write("      </CellData>\n")
                fid.write("    </Piece>\n")
            fid.write("  </PolyData>\n")
            fid.write("</VTKFile>\n")


    def DumpFault(self):

        with h5py.File('fault.h5', 'w') as h5:
            h5['vertex_enu'] = self.vertex_enu
            h5['vertex_llh'] = self.vertex_llh
            h5['element']    = self.element
            h5['slip']       = self.slip
            
