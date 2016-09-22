### This program selects the galaxies with halo tag = -1, and picks out the clusters.
### Updated: Nov 1, 2015

import numpy as np
import os, glob, re, itertools

flink = "/Users/rickyccy/Documents/Research_assemblybias/mock/"
os.chdir(flink)

grid = np.arange(0, 251, 25)

infile = np.loadtxt("sm_gr_fiducial_mock.dat", usecols =(0,1,2,3,4,5,6,7,8,9,10,11,12), unpack = True, dtype = {'names': ('halo_id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'm_vir', 'gg', 'M_r', 'g_r', 'm_host', 'host_id'), 'formats':('S10', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'S10')})

halo_id = infile[0].astype(np.int)
x = infile[1]
y = infile[2]
z = infile[3]
vx = infile[4]
vy = infile[5]
vz = infile[6]
m_vir = infile[7]
gg = infile[8]
M_r = infile[9]
g_r = infile[10]
m_host = infile[11]
host_id = infile[12].astype(np.int)

b = np.where(host_id == -1)[0]
a = host_id[b]

print np.max(a), np.min(a)

outfile = open("sm_nosat.dat", "w")

count = 0

for j in range(a.size):
    if (j % 50000 == 0): print str(j) + " out of " + str(halo_id.size)
    outfile.write(str(halo_id[b[j]]) + "\t" + str(x[b[j]]) + "\t" + str(y[b[j]]) + "\t" + str(z[b[j]]) + "\t" + str(vx[b[j]]) + "\t" + str(vy[b[j]]) + "\t" + str(vz[b[j]]) + "\t" + str(m_vir[b[j]]) + "\t" + str(gg[b[j]]) + "\t" + str(M_r[b[j]]) + "\t" + str(g_r[b[j]]) + "\t" + str(m_host[b[j]]) + "\t" + str(host_id[b[j]]) + "\n")
    count += 1

outfile.close()