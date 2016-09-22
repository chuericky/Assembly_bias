### This program select the shear signals of the concerned galaxies.
### Input: A (argv[1]), rad_fac (argv[2]), mag_b (argv[3]), color (argv[4]), lum (argv[5])
### Updated: Nov 16, 2015


import numpy as np
import os, sys, re, glob

A = sys.argv[1]
r_fac = sys.argv[2]
mag_b = sys.argv[3]
color = sys.argv[4]
lum = sys.argv[5]

shearlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal_2/"
endlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/FOF/"
sortlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_FOF_iso/"
remlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_FOF_iso/"

gal = [0] * 4
gal_new = [0] * 4

for i in range(1,5):
    print i
    
    ffname = "W" + str(i) + "_" + A + "_" + r_fac + "_" + mag_b + "_lens" + color + "_S" + lum + ".tsv"
    
    ingal = np.loadtxt(remlink + ffname, unpack = True, usecols = (0,1,), dtype = {'names': ('id','nam',), 'formats': ('S15','S15',)})

    #bc = np.sort(ingal, order = 'id', kind = 'mergesort')
    #print bc

    gal[i - 1] = np.loadtxt(shearlink + "W" + str(i) + "shear_wider.dat", dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
    #gal[i - 1] = np.loadtxt(shearlink + "W" + str(i) + "shear.dat", dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})

    bc = np.sort(gal[i - 1], order = 'id', kind = 'mergesort')
    
    gal_name = [0] * gal[i - 1].size
    
    for j in range(gal[i - 1].size):
        gal_name[j] = bc[j][0]
    
    c = np.searchsorted(gal_name, ingal[0])

    gal_new[i - 1] = bc[c]

    del bc, gal_name, ingal, c


new_mat = np.hstack((gal_new[0], gal_new[1]))
new_mat = np.hstack((new_mat, gal_new[2]))
new_mat = np.hstack((new_mat, gal_new[3]))


#np.savetxt(endlink + str(b) + "_" + P_th + "_" + clus_num + "_" + color + "_shear.dat", new_mat, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")
np.savetxt(endlink + A + "_" + r_fac + "_" + mag_b + "_" + color + "_S" + lum + "_shear_wider.dat", new_mat, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")

del new_mat

