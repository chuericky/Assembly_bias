### This program plots the 2pt functions.
### Input: parameter name (argv[1]), b (argv[2]), clus_num (argv[3]), parameter1 (argv[4]), parameter2 (argv[5]), P_th (argv[6])
### Updated: Dec 8, 2015

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import os, glob, sys

parameter_name = sys.argv[1]
b = sys.argv[2]
clus_num = sys.argv[3]
#color1 = sys.argv[4]
#color2 = sys.argv[5]
P_th = sys.argv[4]
parameter1 = b + "_" + P_th + "_" + clus_num + "_lensblu"
parameter2 = b + "_" + P_th + "_" + clus_num + "_lensred"

### Number of all lenses, randoms
lens_all_num1 = np.sum(np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/" + parameter_name + "/" + parameter1 + "_lensname.txt", usecols = (0,), unpack = True))
lens_all_num2 = np.sum(np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/" + parameter_name + "/" + parameter2 + "_lensname.txt", usecols = (0,), unpack = True))
#lens_all_num3 = np.sum(np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/onlybright/blu_bew_lensname.txt", usecols = (0,), unpack = True))
#lens_all_num4 = np.sum(np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/onlybright/red_bew_lensname.txt", usecols = (0,), unpack = True))
#lens_all_num5 = np.sum(np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/onlybright/blu_cew_lensname.txt", usecols = (0,), unpack = True))
#lens_all_num6 = np.sum(np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/onlybright/red_cew_lensname.txt", usecols = (0,), unpack = True))
#auto_gal = np.sum(np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/" + parameter_name + "/gal_aew_lensname.txt", usecols = (0,), unpack = True))

#print auto_gal

ranpairlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/"
piclink = "/Users/rickyccy/Documents/Research_assemblybias/result/Feb_20_2016/"

theta_range = np.loadtxt(ranpairlink + "radial_bin2.dat", unpack = True, usecols = (1,))

os.chdir(ranpairlink)

fx1 = np.loadtxt(parameter1 + "_2pt_auto_fx_wider.dat", unpack = True)
fx2 = np.loadtxt(parameter2 + "_2pt_auto_fx_wider.dat", unpack = True)
#fx1 = np.loadtxt(parameter1 + "_2pt_cross_fx_wider.dat", unpack = True)
#fx2 = np.loadtxt(parameter2 + "_2pt_cross_fx_wider.dat", unpack = True)
#fx3 = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/onlybright/blu_bew_2pt_auto_fx.dat", unpack = True)
#fx4 = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/onlybright/red_bew_2pt_auto_fx.dat", unpack = True)
#fx5 = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/onlybright/blu_cew_2pt_auto_fx.dat", unpack = True)
#fx6 = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/onlybright/red_cew_2pt_auto_fx.dat", unpack = True)
"""
fx1 = np.loadtxt(parameter1 + "_2pt_auto_fx.dat", unpack = True)
fx2 = np.loadtxt(parameter2 + "_2pt_auto_fx.dat", unpack = True)
"""
fig = plt.figure()
ax = host_subplot(111, axes_class=AA.Axes)


ax.errorbar(theta_range, fx1[0], yerr = [fx1[1], fx1[1]], color = 'b', fmt = 'o', alpha = 0.7)
ax.errorbar(theta_range, fx2[0], yerr = [fx2[1], fx2[1]], color = 'r', fmt = 'o', alpha = 0.7)
#plt.errorbar(theta_range, fx3[0], yerr = [fx3[1], fx3[1]], color = 'c', fmt = '^', alpha = 0.7)
#plt.errorbar(theta_range, fx4[0], yerr = [fx4[1], fx4[1]], color = 'm', fmt = '^', alpha = 0.7)
#plt.errorbar(theta_range, fx5[0], yerr = [fx5[1], fx5[1]], color = 'g', fmt = 's', alpha = 0.7)
#plt.errorbar(theta_range, fx6[0], yerr = [fx6[1], fx6[1]], color = 'k', fmt = 's', alpha = 0.7)
#ax.legend(["Original red","New red"])
#ax.legend(["gal # = " + str(int(lens_all_num1)) + "\n9.27 <= " + r"$\log(M_\ast/M_\odot)$" + " <= 10.27", "gal # = " + str(int(lens_all_num2)) + "\n10.6 <= " + r"$\log(M_\ast/M_\odot)$" + " <= 10.6"], fontsize = 12)
#plt.legend(["blue, gal # = " + str(int(lens_all_num1)) + r", 0.5 $Mpc h^{-1}$", "red, gal # = " + str(int(lens_all_num2)) + r", 0.5 $Mpc h^{-1}$", "blue, gal # = " + str(int(lens_all_num3)) + r", 0.75 $Mpc h^{-1}$", "red, gal # = " + str(int(lens_all_num4)) + r", 0.75 $Mpc h^{-1}$", "blue, gal # = " + str(int(lens_all_num5)) + r", 1.0 $Mpc h^{-1}$", "red, gal # = " + str(int(lens_all_num6)) + r", 1.0 $Mpc h^{-1}$"], fontsize = 10)
#plt.legend(["gal # = " + str(int(lens_all_num1)) + "\nwhole 'all-galaxies'", "gal # = " + str(int(lens_all_num2)) + "\ngood 'all-galaxies'"], fontsize = 10)
ax.set_ylabel(r'$w(\theta)$')
ax.set_xlabel(r'$\theta (deg)$')
#ax.set_title("auto correlation, A = 1.5", y = 1.1)
ax.set_title(r"Mass $\chi^2$ = 3.425, 2pt fx significance = $0.89\sigma$" , y = 1.1)
#plt.title("Dots - without satellites, triangles - without satellites, auto correlation")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim((0.003, 3))
ax.set_ylim((1e-3, 10))
#plt.vlines(ymin = 1e-5, ymax = 10, x = 0.19, linestyle = '-.')
#plt.vlines(ymin = 1e-5, ymax = 10, x = 0.7, linestyle = '-.')
ax2 = ax.twin()
ax2.set_xticks([0.006195129546466065, 0.06195129546466065, 0.6195129546466065])
ax2.set_xticklabels([r"$10^{-1}$", r"$10^{0}$",
                     r"$10^{1}$"], 'serif')
ax2.set_xlabel(r'$h_{70}^{-1}\rm{Mpc}$')
ax2.axis["right"].major_ticklabels.set_visible(False)
plt.show()
#plt.savefig(piclink + b + "_" + P_th + "_" + clus_num + "_2pt_red_auto.pdf",bbox_inches='tight')
#plt.savefig(piclink + "tew_2pt_all.pdf")


plt.close()


