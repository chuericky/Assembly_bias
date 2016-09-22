### Check the fraction of satellites being removed by different exclusion schemes.
### Updated: Nov 11, 2015

import numpy as np
import os, glob, re, itertools
import matplotlib.pyplot as plt

mock_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/"
os.chdir(mock_link)

infile = np.loadtxt("sm_gr_fiducial_mock.dat", usecols = (12,), unpack = True, dtype = {'names': ('host_id',), 'formats': ('S10',)})


host_id = np.array(infile).astype(np.int)

sat_num_all = np.where(host_id != -1)[0].size

infile1 = np.loadtxt("sm_probab_1.5_0.3_3.0.dat", usecols = (12,), unpack = True, dtype = {'names': ('host_id',), 'formats': ('S10',)})
infile2 = np.loadtxt("sm_probab_1.5_0.4_3.0.dat", usecols = (12,), unpack = True, dtype = {'names': ('host_id',), 'formats': ('S10',)})
infile3 = np.loadtxt("sm_probab_1.5_0.5_3.0.dat", usecols = (12,), unpack = True, dtype = {'names': ('host_id',), 'formats': ('S10',)})
infile4 = np.loadtxt("sm_probab_1.5_0.5_5.0.dat", usecols = (12,), unpack = True, dtype = {'names': ('host_id',), 'formats': ('S10',)})
infile5 = np.loadtxt("sm_probab_1.2_0.5_3.0.dat", usecols = (12,), unpack = True, dtype = {'names': ('host_id',), 'formats': ('S10',)})
infile6 = np.loadtxt("sm_probab_1.0_0.5_3.0.dat", usecols = (12,), unpack = True, dtype = {'names': ('host_id',), 'formats': ('S10',)})

host_id1 = np.array(infile1).astype(np.int)
host_id2 = np.array(infile2).astype(np.int)
host_id3 = np.array(infile3).astype(np.int)
host_id4 = np.array(infile4).astype(np.int)
host_id5 = np.array(infile5).astype(np.int)
host_id6 = np.array(infile6).astype(np.int)

sat_num_1 = np.where(host_id1 != -1)[0].size + 0.
sat_num_2 = np.where(host_id2 != -1)[0].size + 0.
sat_num_3 = np.where(host_id3 != -1)[0].size + 0.
sat_num_4 = np.where(host_id4 != -1)[0].size + 0.
sat_num_5 = np.where(host_id5 != -1)[0].size + 0.
sat_num_6 = np.where(host_id6 != -1)[0].size + 0.

print np.where(host_id4 != -1)[0]


print "Final sat # / Initial sat # = ", 1 - sat_num_1 / sat_num_all, sat_num_1, sat_num_all
print "Final sat # / All gal # = ", sat_num_1 / infile1.size
print "\n"
print "Final sat # / Initial sat # = ", 1 - sat_num_2 / sat_num_all, sat_num_2, sat_num_all
print "Final sat # / All gal # = ", sat_num_2 / infile2.size
print "\n"
print "Final sat # / Initial sat # = ", 1 - sat_num_3 / sat_num_all, sat_num_3, sat_num_all
print "Final sat # / All gal # = ", sat_num_3 / infile3.size
print "\n"
print "Final sat # / Initial sat # = ", 1 - sat_num_4 / sat_num_all, sat_num_4, sat_num_all
print "Final sat # / All gal # = ", sat_num_4 / infile4.size
print "\n"
print "Final sat # / Initial sat # = ", 1 - sat_num_5 / sat_num_all, sat_num_5, sat_num_all
print "Final sat # / All gal # = ", sat_num_5 / infile5.size
print "\n"
print "Final sat # / Initial sat # = ", 1 - sat_num_6 / sat_num_all, sat_num_6, sat_num_all
print "Final sat # / All gal # = ", sat_num_6 / infile6.size
print "\n"
