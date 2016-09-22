### Calculate (red - blue)^T * (C_red + C_blue)^-1 * (red - blue)
### Input: A (argv[1]), r_fac (argv[2]), mag_b (argv[3])
### Updated: Jan 22, 2016

import glob, os, sys, re
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg
from scipy.stats import chisqprob
import scipy.stats as st

nr = 4
nr_s = 12
A = sys.argv[1]
r_fac = sys.argv[2]
mag_b = sys.argv[3]

fstr = A + "_" + r_fac + "_" + mag_b

shearlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/FOF/"

os.chdir(shearlink)

b_2pt = np.loadtxt(fstr + "_lensblu_2pt_auto_fx_wider.dat", unpack = True, usecols = (0,))[nr_s:nr_s+nr]
r_2pt = np.loadtxt(fstr + "_lensred_2pt_auto_fx_wider.dat", unpack = True, usecols = (0,))[nr_s:nr_s+nr]

d_2pt = r_2pt - b_2pt

cov_b = np.loadtxt(fstr + "_lensblu_cov.dat")
cov_r = np.loadtxt(fstr + "_lensred_cov.dat")

C = np.add(cov_b, cov_r)
"""
for i in range(4):
    for j in range(4):
        if (i != j):
            C[i][j] = 0
"""
inv_C = linalg.inv(C)

inv_C = np.matrix(inv_C)

print C
print inv_C

print np.diag(np.dot(inv_C, C))
print b_2pt
print r_2pt



dot = np.array(np.dot(np.dot(d_2pt,inv_C),d_2pt))[0][0]

## p-value
p = chisqprob(dot, nr)

print "chi2 = " + str(dot)
print "Deg. of freedom = " + str(nr)
print "p-value = " + str(p)
print "Two-tailed Z-score = " + str(st.norm.ppf(1 - p/2))
