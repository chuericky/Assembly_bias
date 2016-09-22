### This program saves the shear profile. Inverse covariance matrix is calculated as well.
### Updated: Dec 15, 2015
### Input: argv[1] file.

import numpy as np
import matplotlib.pyplot as plt
import sys
from numpy import linalg

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=6)

file = sys.argv[1]

shear_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/shear_profile/"

r_bin = np.loadtxt(shear_link + "r.dat")
N_ran = 1000

infile = np.loadtxt(shear_link + file, unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100))

bin_num = r_bin.size
f_size = len(infile[0])


gg = [0] * bin_num
avg_gg = [0] * bin_num

for i in range(bin_num):
    gg[i] = infile[i]
    avg_gg[i] = np.average(gg[i])

### The array to store up the bootstrap samples.
DeltaSig = [0] * bin_num

for i in range(bin_num):
    DeltaSig[i] = [0] * N_ran

### Bootstrap
boot_num = [0] * N_ran

for i in range(N_ran):
    boot_num[i] = np.random.randint(0, f_size, size = f_size)

for j in range(bin_num):
    if (j % 10 == 0): print str(j) + " out of " + str(bin_num)
    for i in range(N_ran):
        DeltaSig[j][i] = np.average(gg[j][boot_num[i]])

cov_g_t = np.cov(DeltaSig)
sd_g_t = np.power(np.diag(cov_g_t), 0.5)


cov_fit_g_t = cov_g_t[0:20,0:20]
inv_cov_fit_g_t = linalg.inv(cov_fit_g_t)

print np.diag(np.dot(cov_fit_g_t, inv_cov_fit_g_t))

gg_mat = np.vstack((avg_gg, sd_g_t)).T


np.savetxt(shear_link + "cov_" + file[:-4] + "_gg.dat", cov_fit_g_t, fmt="%.9E", delimiter = "\t", newline = "\n")
np.savetxt(shear_link + file[:-4] + "_gg.dat", gg_mat, fmt="%.6E", delimiter = "\t", newline = "\n")

"""
plt.errorbar(r_bin, avg_gg, yerr = [sd_g_t, sd_g_t], fmt = 'o')
plt.xscale('log')
plt.yscale('log')
plt.show()
plt.close()
"""
del r_bin, infile, gg, avg_gg, DeltaSig, cov_g_t, sd_g_t


