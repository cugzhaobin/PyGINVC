#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 00:26:56 2020

@author: zhao
"""

import h5py
import numpy as np
import poly3d
import random

h5 = h5py.File('luhuo.e', 'r')
e  = h5['coordx'][()]*1e3
n  = h5['coordy'][()]*1e3
#u  = h5['coordz'][()]
u  = e*0.0
for i in range(len(e)):
    print(random.random())
#   e[i] = e[i]+random.random()/20
#   n[i] = n[i]+random.random()/20
#   u[i] = u[i]+random.random()/20
vertex = np.vstack((e,u,n)).T
mapelem = h5['connect1'][()]
disp    = np.zeros((len(mapelem),3))
disp[:,1] = disp[:,1]-30000

poly3d.gen_poly3d_input(vertex, mapelem, disp, bctype=[], vcs='FaultCS', ecs='elocal')
