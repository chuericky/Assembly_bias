### Combine the jackknife results into a single file.
### Updated: Nov 8, 2015

import numpy as np
import os, re, glob, sys

f_int = sys.argv[1]
l_int = int(f_int) - 1

b = sys.argv[2]
clus_num = sys.argv[3]
color = sys.argv[4]
str_name = sys.argv[5]   #color and nature (DR, RR)
nat = sys.argv[6]           ### auto or cross
P_th = sys.argv[7]

lenslink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/FOF/W" + f_int + "/"
comlink = lenslink + "combined/"

#fname = glob.glob(lenslink + color + "_" + str_name + "_*_pdf_cross_jack.dat")
fname = glob.glob(lenslink + b + "_" + P_th + "_" + clus_num + "_" + color + "_" + str_name + "*_" + nat + "_jack_wider.dat")
#fname = glob.glob(lenslink + str_name + "_*_pdf_cross_jack.dat")

jack_index = np.loadtxt(fname[0], unpack = True, usecols = (0,))

print jack_index

infile = [0] * len(fname)

for i in range(len(fname)):
    infile[i] = np.loadtxt(fname[i], unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))

matrix = np.vstack((jack_index, np.sum(infile, axis = 0))).T

np.savetxt(comlink + b + "_" + P_th + "_" + clus_num + "_" + color + "_" + str_name + "_" + nat + "_jack_wider.dat", matrix, fmt = '%u', delimiter='\t', newline='\n')
#np.savetxt(comlink + str_name + "_pdf_cross_jack.dat", matrix, fmt = '%u', delimiter='\t', newline='\n')