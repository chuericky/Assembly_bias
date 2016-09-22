### Plots the dSig profiles.
### Updated: Jun 13, 2016

import numpy as np
import glob, os, re, sys
import matplotlib.pyplot as plt

dsig_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/Delta_Sig/"
shearlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/"
save_link = "/Users/rickyccy/Documents/Research_assemblybias/result/Jun_7_2016/density_contrast/"

os.chdir(dsig_link)

quar_num = 4


r_bin = np.loadtxt(shearlink + "dist_bin2.dat", unpack = True, usecols = (1,)) * 1000

fstr = "smass_quartile"

shear = [0] * quar_num
g_err = [0] * quar_num
sm_shear = [0] * quar_num

color = ["b", "g", "c", "r"]
mark = ["^", "s", "o", "p"]

for i in range(quar_num):
    shear[i], g_err[i], sm_shear[i] = np.loadtxt("{}{}_gg.dat".format(fstr,i+1), unpack = True)
    plt.errorbar(r_bin, np.add(shear[i], sm_shear[i]), yerr = [g_err[i], g_err[i]], color = color[i], fmt = mark[i], alpha = 0.6)

plt.ylabel(r'$\Delta\Sigma [h_{70} M_\odot pc^{-2}]$', fontsize = 20)
plt.xlabel(r'r $[h_{70}^{-1} kpc]$', fontsize = 20)
plt.legend(["Lowest mass, 0%-25%","25%-50%","50%-75%","Highest mass, 75%-100%"], fontsize = 10, loc = "upper right")
plt.yscale('log')
plt.xscale('log')
plt.ylim([5e-2, 3e2])
plt.xlim([20, 1.6e3])
plt.savefig(save_link + "dSig_smass.png",bbox_inches='tight')
plt.close()