### Calculate (red - blue)^T * (C_red + C_blue)^-1 * (red - blue)
### Input: A (argv[1]), r_fac (argv[2]), mag_b (argv[3]), lum (argv[4])
### Updated: Dec 11, 2015

import glob, os, sys, re
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg

nr = 7
A = sys.argv[1]
r_fac = sys.argv[2]
mag_b = sys.argv[3]
color = sys.argv[4]
color2 = sys.argv[5]
lum = sys.argv[6]
lum2 = sys.argv[7]


fstr = A + "_" + r_fac + "_" + mag_b

shearlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/FOF/"
gallink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/"
savelink = "/Users/rickyccy/Documents/Research_assemblybias/result/Dec_14_2015/"

os.chdir(shearlink)

b_gg = np.loadtxt(fstr + "_" + color + "_S" + lum + "_gg_wider.dat", unpack = True, usecols = (0,))[:nr]
r_gg = np.loadtxt(fstr + "_" + color2 + "_S" + lum2 + "_gg_wider.dat", unpack = True, usecols = (0,))[:nr]

cov_b = np.loadtxt("cov_fit_" + fstr + "_" + color + "_S" + lum + "_gg_wider.dat")
cov_r = np.loadtxt("cov_fit_" + fstr + "_" + color2 + "_S" + lum2 + "_gg_wider.dat")

C = np.add(cov_b, cov_r)

d_gg = r_gg - b_gg

inv_C = linalg.inv(C)

print np.diag(np.dot(C, inv_C))

inv_C = np.matrix(inv_C)

dot = np.array(np.dot(np.dot(d_gg,inv_C),d_gg))[0][0]

print "chi2 = " + str(dot)


r_bin = np.loadtxt(gallink + "dist_bin2.dat", unpack = True, usecols = (1,)) * 1000
b_shear, b_err = np.loadtxt(fstr + "_" + color + "_S" + lum + "_gg_wider.dat", unpack = True, usecols = (0,1))
r_shear, r_err = np.loadtxt(fstr + "_" + color2 + "_S" + lum2 + "_gg_wider.dat", unpack = True, usecols = (0,1))

plt.errorbar(r_bin, b_shear, yerr = [b_err, b_err], color = 'b', fmt = 'o', alpha = 0.7)
plt.errorbar(r_bin, r_shear, yerr = [r_err, r_err], color = 'r', fmt = 'o', alpha = 0.7)
plt.ylabel(r'$\Delta\Sigma [h_{70} M_\odot pc^{-2}]$', fontsize = 20)
plt.xlabel(r'r $[h_{70}^{-1} kpc]$', fontsize = 20)
plt.legend(["gal # = 17312\n9.35 <= " + r"$\log(M_\ast/M_\odot)$" + " <= 10.35", "gal # = 6763\n9.05 <= " + r"$\log(M_\ast/M_\odot)$" + " <= 10.05"], fontsize = 12)
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-1, 2e2])
plt.xlim([20, 1.6e3])
plt.show()
#plt.savefig(savelink + fstr + "_dSig.pdf")
plt.close()
