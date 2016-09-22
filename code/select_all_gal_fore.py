### This program adjusts the z_max for the 'all-galaxies' group.
### Updated: Sept 25, 2015

import numpy as np
import os, re, sys, glob

other_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/other_lens/"

os.chdir(other_link)


for i in range(4):
    print i
    all_gal = np.loadtxt(other_link + "W" + str(i + 1) + "_gal_fore2.dat", unpack = True, dtype = {'names': ('id','field', 'ra', 'dec', 'z', 'T_b'), 'formats':('S15','S15','float64','float64','float64','float64')})
    all_gal2 = np.loadtxt(other_link + "W" + str(i + 1) + "_gal_zew.dat", unpack = True, dtype = {'names': ('id','field', 'ra', 'dec', 'z', 'T_b'), 'formats':('S15','S15','float64','float64','float64','float64')})
    
    a = np.where((all_gal[4] <= 0.55))[0]
    b = np.where((all_gal2[4] >= 0.35))[0]
    
    outfile = open("W" + str(i + 1) + "_gal_ccew.dat", "w")
    
    """
    for k in range(all_gal2[0].size):
        outfile.write(str(all_gal2[0][k]) + "\t" + str(all_gal2[1][k]) + "\t" + str(all_gal2[2][k]) + "\t" + str(all_gal2[3][k]) + "\t" + str(all_gal2[4][k]) + "\t" + str(all_gal2[5][k]) + "\n")
    """
    for k in range(b.size):
        outfile.write(str(all_gal2[0][b[k]]) + "\t" + str(all_gal2[1][b[k]]) + "\t" + str(all_gal2[2][b[k]]) + "\t" + str(all_gal2[3][b[k]]) + "\t" + str(all_gal2[4][b[k]]) + "\t" + str(all_gal2[5][b[k]]) + "\n")
    
    for k in range(a.size):
        outfile.write(str(all_gal[0][a[k]]) + "\t" + str(all_gal[1][a[k]]) + "\t" + str(all_gal[2][a[k]]) + "\t" + str(all_gal[3][a[k]]) + "\t" + str(all_gal[4][a[k]]) + "\t" + str(all_gal[5][a[k]]) + "\n")

    outfile.close()

    del a, b

