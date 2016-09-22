    ### Divide the lens galaxies into different files according to their pointings.
### Input: Field #: argv[1]
### Updated: May 9, 2016

import numpy as np
import os, glob, re, sys

gal_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/"

f_num = sys.argv[1]

os.chdir (gal_link)

lensgal = np.loadtxt("W{}lens.tsv".format(f_num), unpack = True, usecols = (0,1,4,5,24), dtype = {'names': ('id','field','alpha','delta','z'), 'formats': ('S15','S15','float64','float64','float64')})
lgal = np.loadtxt("W{}lens.tsv".format(f_num), usecols = (0,1,4,5,24), dtype = {'names': ('id','field','alpha','delta','z'), 'formats': ('S15','S15','float64','float64','float64')})

pts = np.unique(lensgal[1])

pts_num = pts.size

gal_info = [0] * pts_num

count = 0

for i in range(pts_num):
    gal_info[i] = lgal[np.where((lensgal[1] == pts[i]))[0]]
    count += gal_info[i].size
    np.savetxt("{}/W{}/{}_lens.tsv".format(gal_link, f_num, pts[i]), gal_info[i], delimiter = "\t", newline = "\n", fmt = "%s\t%s\t%.10g\t%.10g\t%.10g")


if (count != lensgal[0].size):
    print "Something wrong with division!"
    exit(1)
