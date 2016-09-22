### This program selects the Sigma profiles for the designated profiles.
### Updated: catalog name, argv[1], color, argv[2], A (argv[3]), rad_fac (argv[4]), mag_bright (argv[5])
### Updated: May 26, 2016

import numpy as np
import glob, os, sys
import matplotlib.pyplot as plt

catalog = sys.argv[1]
color = sys.argv[2]
fac = sys.argv[3]
rad_fac = sys.argv[4]
mag_b = sys.argv[5]

shear_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/shear_profile/"
#mock_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/erased_assembly_bias_mocks/"
mock_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/"
mock_ex_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/mock/"

color_list = np.loadtxt(mock_link + color + "_" + fac + "_" + rad_fac + "_" + mag_b + ".dat", unpack = True, usecols = (0,1), dtype = {'names': ('id','x'), 'formats': ('int64','float64')})


### Radius of the profiles
radius = np.loadtxt(shear_link + "r.dat")

infile2 = np.loadtxt(shear_link + "DeltaSig_" + catalog + ".txt", dtype = {'names': ('id','sig1','sig2','sig3','sig4','sig5','sig6','sig7','sig8','sig9','sig10','sig11','sig12','sig13','sig14','sig15','sig16','sig17','sig18','sig19','sig20','sig21','sig22','sig23','sig24','sig25','sig26','sig27','sig28','sig29','sig30','sig31','sig32','sig33','sig34','sig35','sig36','sig37','sig38','sig39','sig40','sig41','sig42','sig43','sig44','sig45','sig46','sig47','sig48','sig49','sig50','sig51','sig52','sig53','sig54','sig55','sig56','sig57','sig58','sig59','sig60','sig61','sig62','sig63','sig64','sig65','sig66','sig67','sig68','sig69','sig70','sig71','sig72','sig73','sig74','sig75','sig76','sig77','sig78','sig79','sig80','sig81','sig82','sig83','sig84','sig85','sig86','sig87','sig88','sig89','sig90','sig91','sig92','sig93','sig94','sig95','sig96','sig97','sig98','sig99','sig100'), 'formats':('int64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})

#halo_id = infile[0]

count = 0

b = np.sort(infile2, order = 'id', kind = 'mergesort')

gal_name = [0] * infile2.size

for i in range(infile2.size):
    gal_name[i] = b[i][0]

c = np.searchsorted(gal_name, color_list[0])

gal = b[c]

#a = np.searchsorted(halo_id, color_list[0])

np.savetxt(mock_ex_link + color + "_" + fac + "_" + rad_fac + "_" + mag_b + ".dat", gal, fmt = '%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")



"""
for i in range(color_list.size):
    if (i % 500 == 0): print str(i) + " out of " + str(color_list.size) + " finished!"
    for j in range(len(halo_id)):
        if (color_list[i][0] == halo_id[j]):
            count += 1
            gal.append(infile2[j])
            break

np.savetxt(shear_link + "withoutab_" + color + ".dat", gal, fmt = '%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")


del infile, gal, infile2
"""
