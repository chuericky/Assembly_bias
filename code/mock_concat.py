### Concatenate subhalo and host halo mass and concentrations to the galaxy file.
### Updated: May 26, 2016

import numpy as np
import os, re, sys

mock_list = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/mock/"

os.chdir (mock_list)

infile = np.loadtxt("sm_1.5_0.5_3.0_catalog_pre.dat", unpack = True, dtype = {'names': ('id','x','y','z','vx','vy','vz','m_sub','vpeak','sm','gr','mhost','host_id'), 'formats': ('i8','float64','float64','float64','float64','float64','float64','float64','float64','float64','float64','float64','i8')})

halofile = np.loadtxt("sm_halo_catalog_withconc.dat", dtype = {'names': ('id','rs','m200c','conc'), 'formats': ('i8','float64','float64','float64')})

bc = np.sort(halofile, order = 'id', kind = 'mergesort')

halo_name = [0] * halofile.size

for j in range(halofile.size):
    halo_name[j] = bc[j][0]

un = np.searchsorted(halo_name, infile[0])

c = bc[un]

d = np.zeros((infile[0].size, 3))

for i in range(infile[0].size):
    for j in range(1,4):
        d[i][j - 1] = c[i][j]


mat = np.vstack((infile, d.T)).T

np.savetxt("sm_1.5_0.5_3.0_catalog.dat", mat, delimiter = "\t", fmt = "%lu\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%lu\t%g\t%g\t%g", newline = "\n")


#for j in range(infile[0].size):
#    print j, c[j][0], c[j][1], c[j][2], c[j][3]

#print c




