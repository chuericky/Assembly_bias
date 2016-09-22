### Input: A (argv[1]), r_fac (argv[2]), mag_b (argv[3]), lum (argv[4]), color (argv[5])
### Updated: Nov 30, 2015

import numpy as np
import glob, os, sys

A = sys.argv[1]
r_fac = sys.argv[2]
mag_b = sys.argv[3]
lum = sys.argv[4]
color = sys.argv[5]

fstr = A + "_" + r_fac + "_" + mag_b + "_lens" + color + "_S" + lum
fstr2 = A + "_" + r_fac + "_" + mag_b + "_" + color + "_S" + lum

gallink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/FOF/"
shearlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/"

rad = np.loadtxt(shearlink + "dist_bin2.dat", unpack = True, usecols = (1,)) * 1e6
rad2 = np.multiply(np.pi, np.power(rad, 2))

os.chdir(gallink)

red_mass = []
red_z = []

redname = glob.glob("W*_" + fstr + ".tsv")

for i in range(len(redname)):
    redfile = np.loadtxt(redname[i], unpack = True, usecols = (4,7,))
    red_mass.append(np.power(10, redfile[1]))
    red_z.append(redfile[0])

red_mass2 = np.concatenate((red_mass[0], red_mass[1]))
red_mass2 = np.concatenate((red_mass2, red_mass[2]))
red_mass2 = np.concatenate((red_mass2, red_mass[3]))

red_z2 = np.concatenate((red_z[0], red_z[1]))
red_z2 = np.concatenate((red_z2, red_z[2]))
red_z2 = np.concatenate((red_z2, red_z[3]))

del red_mass, red_z

mea = np.mean(red_mass2)

print color, np.log10(mea), mea
star_shear = np.divide(mea, rad2)

os.chdir(gallink)

gg = np.loadtxt(fstr2 + "_gg_wider.dat", unpack = True, usecols = (0,1))

g_shear = np.vstack((gg, star_shear)).T

np.savetxt(fstr2 + "_gg_wider.dat", g_shear, delimiter = "\t", fmt = "%.6E\t%.6E\t%.6E")
"""
outfile = open("stellar_mass.txt", "w")

outfile.write(A + "_" + r_fac + "_" + mag_b + "_blu_S" + lum + "_gg.dat\t")
outfile.write(str(np.log10(np.mean(blu_mass2))))
outfile.write("\t")
outfile.write(str(np.median(blu_z2)))
outfile.write("\n" + A + "_" + r_fac + "_" + mag_b + "_red_S" + lum + "_gg.dat\t")
outfile.write(str(np.log10(np.mean(red_mass2))))
outfile.write("\t")
outfile.write(str(np.median(red_z2)))
outfile.write("\n")


outfile.close()

outfile2 = open("prelim_mass.txt", "w")

outfile2.write(A + "_" + r_fac + "_" + mag_b+ "_blu_S" + lum + "_gg.dat\t150.00\n" + A + "_" + r_fac + "_" + mag_b + "_red_S" + lum + "_gg.dat\t150.00\n")

outfile2.close()

"""


