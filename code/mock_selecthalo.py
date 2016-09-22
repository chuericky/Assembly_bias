### Select the subsamples of halos to cut down computational time.
### Updated: May 26, 2016

import numpy as np
import os, sys, re

mock_list = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/mock/"
halo_list = "/Users/rickyccy/Documents/Research_assemblybias/mock/halo_catalog/"

os.chdir (mock_list)

infile = np.loadtxt("sm_1.5_0.5_3.0_catalog.dat", unpack = True, usecols = (0,12), dtype = {'names': ('id','host_id'), 'formats': ('i8','i8')})

host_id = infile[1][np.where((infile[1] != -1))[0]]

ids = np.hstack((infile[0], host_id))

### Unique ids
ids = np.unique(ids)

print ids.size
print ids

del infile

os.chdir (halo_list)

halofile = np.loadtxt("hlist_1.00030.list", usecols = (1,10,12,37), dtype = {'names': ('id','mvir','rs','m200c'), 'formats': ('i8','float64','float64','float64')})

bc = np.sort(halofile, order = 'id', kind = 'mergesort')

halo_name = [0] * halofile.size

for j in range(halofile.size):
    halo_name[j] = bc[j][0]

un = np.searchsorted(halo_name, ids)


#print un

#print bc[un]

np.savetxt(mock_list + "sm_halo_catalog.dat", bc[un], delimiter = "\t", fmt = "%lu\t%g\t%g\t%g", newline = "\n")


#print np.where((un != 0))[0]

#print un[np.where((un != 0))[0]]

del halofile