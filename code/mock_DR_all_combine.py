### This program combines the RR files.
### Updated: Dec 15, 2015

import numpy as np
import glob, os, re

f_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/2pt/DR/projected/"
inter_link = f_link + "inter/"
com_link = f_link + "combined/"

cat_name = ["blu_0.0_0.0_0.0", "red_0.0_0.0_0.0","blu_0.0_0.0_9.0", "red_0.0_0.0_9.0"]
#, "blu_0.0_0.2", "blu_0.0_0.0", "blu_1.0_0.0", "blu_1.2_0.0", "blu_1.5_0.0", "red_0.0_0.2", "red_0.0_0.0", "red_1.0_0.0", "red_1.2_0.0", "red_1.5_0.0"]

os.chdir(f_link)

R_num = 2000000.


for i in range(len(cat_name)):   #20
    fname = glob.glob(cat_name[i] + "*_auto_all.dat")
    infile = [0] * len(fname)

    for j in range(len(fname)):
        infile[j] = np.loadtxt(fname[j], unpack = True, dtype={'names': ('gal_num', 'g0', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12'), 'formats': ('i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15')})

    #print fname[j], str(j), infile[j]

    in_sum = np.sum(infile, axis = 0)

    pair_num = np.multiply(in_sum[0], R_num)

    RR = np.hstack((in_sum[0], np.divide(in_sum[1:], pair_num))).T

    del fname, infile, in_sum, pair_num

    np.savetxt(com_link + cat_name[i] + "_DR_auto_all.dat", RR[None], delimiter = "\t", newline = "\n", fmt = "%d\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E")

    del RR

"""
fname = glob.glob("*_all.dat")
infile = [0] * len(fname)
    
for j in range(len(fname)):
    infile[j] = np.loadtxt(fname[j], unpack = True, dtype={'names': ('gal_num', 'g0', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'g15', 'g16', 'g17', 'g18', 'g19', 'g20', 'g21', 'g22', 'g23', 'g24', 'g25', 'g26', 'g27', 'g28', 'g29'), 'formats': ('i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15')})
    #print fname[j], str(j), infile[j]
    
in_sum = np.sum(infile, axis = 0)

pair_num = in_sum[0] * (in_sum[0] - 1.)
RR = np.hstack((in_sum[0], np.divide(in_sum[1:], pair_num)))

del fname, infile, in_sum, pair_num

np.savetxt(com_link + "RR_all.dat", RR, delimiter = "\t", newline = "\n", fmt = "%d\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E")

del RR
"""