### This program combined the DD files.
### Updated: Dec 15, 2015

import numpy as np
import glob, os, re

f_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/2pt/DD/projected/"
com_link = f_link + "combined/"

os.chdir(f_link)

cat_name = ["blu_0.0_0.0_0.0", "red_0.0_0.0_0.0","blu_0.0_0.0_9.0", "red_0.0_0.0_9.0"]
#, "blu_0.0_0.2", "blu_0.0_0.0", "blu_1.0_0.0", "blu_1.2_0.0", "blu_1.5_0.0", "red_0.0_0.2", "red_0.0_0.0", "red_1.0_0.0", "red_1.2_0.0", "red_1.5_0.0"]

for i in range(len(cat_name)):  #len(cat_name)
    print i

    fname = glob.glob(cat_name[i] + "_DD*auto_all.dat")
    infile = [0] * len(fname)

    for k in range(len(fname)):
        infile[k] = np.loadtxt(fname[k], unpack = True, dtype={'names': ('gal_num', 'g0', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12'), 'formats': ('i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15', 'i15')})
        
    in_sum = np.sum(infile, axis = 0)

    pair_num = np.multiply(infile[0][0], (infile[0][0] - 1.))

    print in_sum[1:]

    DD = np.hstack((infile[0][0], np.divide(in_sum[1:], pair_num))).T


    """
    pair_num = in_sum[1] * (in_sum[1] - 1.)

    DD = np.vstack((id, in_sum[1]))
    DD = np.vstack((DD, np.divide(in_sum[2:], pair_num))).T
    """
    
    del fname, infile, in_sum, pair_num

    np.savetxt(com_link + cat_name[i] + "_DD_auto_all.dat", DD[None], delimiter = "\t", newline = "\n", fmt = "%d\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E")
    del DD

