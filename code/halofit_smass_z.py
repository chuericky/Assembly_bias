### This program selects the lens galaxies according to the luminosities.
### Input: A (argv[1]), r_fac (argv[2]), mag_b (argv[3]), color = argv[4], lum (argv[5])
### Updated: Nov 18, 2015

import glob, os, sys
import numpy as np

A = sys.argv[1]
r_fac = sys.argv[2]
mag_b = sys.argv[3]
color = sys.argv[4]
lum = sys.argv[5]

sortlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_FOF_iso/"
#srcslink = "/Users/rickyccy/Documents/Research_assemblybias/Data/source2/combined2/"
gallink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/FOF/"

os.chdir (sortlink)

fname = glob.glob("W*_" + A + "_" + r_fac + "_" + mag_b + "_lens" + color + "_S" + lum + ".tsv")


for i in range(len(fname)):
    print fname[i]
    ### Load in the lens catalog.
    lensfile = np.loadtxt(sortlink + fname[i], usecols = (0, 1, 4, 5, 24, 27, 32, 36,), dtype = {'names': ('id', 'field', 'alpha', 'delta', 'gal_z', 'color', 'Mr', 'lp_med'), 'formats':('S15', 'S15', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
    
    #lensfile2 = np.loadtxt(sortlink + fname[i], usecols = (0, 1, 4, 5, 24, 27, 32, 36,), dtype = {'names': ('id', 'field', 'alpha', 'delta', 'gal_z', 'color', 'Mr', 'lp_med'), 'formats':('S15', 'S15', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
    
    np.savetxt(gallink + fname[i], lensfile, delimiter = "\t", fmt = "%s")
    
    """
    ra = np.array(lensfile[2])
    dec = np.array(lensfile[3])
    z = np.array(lensfile[4])
    color = np.array(lensfile[5])
    Mr_new = np.array(lensfile[6])
    mass = np.array(lensfile[7])

    ### Luminosity bands
    L1 = np.where((Mr_new <= -20.75) & (Mr_new >= -21.2) & (color <= 1.5))
    L2 = np.where((Mr_new <= -22.0) & (Mr_new >= -22.8) & (color <= 4.0) & (color >= 2.0))
    
    #source = np.where(((Mr_new > -20.0) & (mr <= 24.7)))
    
    list1 = lensfile2[L1]
    list2 = lensfile2[L2]
    #src = lensfile2[source]
    
    np.savetxt(gallink + fname[i][:-4] + 'red.tsv', list1, delimiter = "\t", fmt = "%s")
    np.savetxt(gallink + fname[i][:-4] + 'blu.tsv', list2, delimiter = "\t", fmt = "%s")    #np.savetxt(srcslink + fname[i][:2] + 'src_' + fname[i][2] + '2.tsv', src, delimiter = "\t", fmt = "%s")
    """
