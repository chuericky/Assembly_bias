### Combine the files into a single file.
### Updated: Nov 8, 2015

import numpy as np
import os, re, glob, sys

f_int = sys.argv[1]
l_int = int(f_int) - 1

b = sys.argv[2]
clus_num = sys.argv[3]
color = sys.argv[4]
str_name = sys.argv[5]
nat = sys.argv[6]
P_th = sys.argv[7]

lenslink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/FOF/W" + f_int + "/"
comlink = lenslink + "combined/"

#fname = glob.glob(lenslink + color + "_" + str_name + "*_pdf_cross_all.dat")
fname = glob.glob(lenslink + b + "_" + P_th + "_" + clus_num + "_" + color + "_" + str_name + "*_" + nat + "_all_wider.dat")
#fname = glob.glob(lenslink + str_name + "*_pdf_cross_all.dat")

infile = [0] * len(fname)

for i in range(len(fname)):
    infile[i] = np.loadtxt(fname[i], unpack = True)

matrix = np.sum(infile, axis = 0)
print matrix

np.savetxt(comlink + b + "_" + P_th + "_" + clus_num + "_" + color + "_" + str_name + "_" + nat + "_all_wider.dat", matrix[None], fmt = '%u', delimiter='\t', newline='\n')
#np.savetxt(comlink + str_name + "_pdf_cross_all.dat", matrix[None], fmt = '%u', delimiter='\t', newline='\n')