### Combined shear signals of individual galaxies into a single file.
### Updated: Jul 6, 2015

import numpy as np
import glob, os, sys
import re
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

shearlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal_3/W" + sys.argv[1]
savelink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal_2/"

x = []

os.chdir(shearlink)

fname0 = sorted(glob.glob("*.dat"), key=numericalSort)


for i in range(len(fname0)):
    if (i % 10000 == 0): print i
    x.append(np.loadtxt(fname0[i], dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')}))
    #x.append(np.loadtxt(fname0[i], dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')}))

del fname0

os.chdir(savelink)
np.savetxt("W" + sys.argv[1] + "shear.dat", x, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")
#np.savetxt("W" + sys.argv[1] + "shear_ext.dat", x, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")

del x
