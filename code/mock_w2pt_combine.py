### This program calculates the projected 2pt function.
### Updated: Dec 15, 2015

import numpy as np
import matplotlib.pyplot as plt
import glob, os, re, sys
import matplotlib.gridspec as gridspec
from pylab import *

dir = "2pt"

rad_fac = sys.argv[1]

#cat_name = ["blu_1.0_0.3_" + rad_fac, "red_1.0_0.3_" + rad_fac, "blu_1.2_0.3_" + rad_fac, "red_1.2_0.3_" + rad_fac, "blu_1.5_0.3_" + rad_fac, "red_1.5_0.3_" + rad_fac]
cat_name = ["blu_1.5_0.5_1.0", "red_1.5_0.5_1.0", "blu_1.5_0.5_0.5", "red_1.5_0.5_0.5", "blu_1.5_0.5_0.1", "red_1.5_0.5_0.1"]
color_opt = ["b", "r"]

DD_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/" + dir + "/DD/projected/combined/"
DR_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/" + dir + "/DR/projected/combined/"
RR_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/" + dir + "/RR/combined/"
mock_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/"
SDSS_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/2pt/SDSS_data/"
twopt_link = mock_link + "2pt/"
result_link = "/Users/rickyccy/Documents/Research_assemblybias/result/Dec_14_2015/"

ran_all_num = 2000000
jack_num = 100

radius = np.loadtxt(mock_link + "radial_bin.dat", unpack = True, usecols = (1,))

RR_file = np.loadtxt(RR_link + "ran_other_RR_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
RR_file_jack = np.loadtxt(RR_link + "ran_other_RR_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))


blu_DD = np.loadtxt(DD_link + cat_name[0] + "_DD_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
red_DD = np.loadtxt(DD_link + cat_name[1] + "_DD_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
blu_DR = np.loadtxt(DR_link + cat_name[0] + "_DR_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
red_DR = np.loadtxt(DR_link + cat_name[1] + "_DR_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
blu2_DD = np.loadtxt(DD_link + cat_name[2] + "_DD_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
red2_DD = np.loadtxt(DD_link + cat_name[3] + "_DD_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
blu2_DR = np.loadtxt(DR_link + cat_name[2] + "_DR_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
red2_DR = np.loadtxt(DR_link + cat_name[3] + "_DR_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
blu3_DD = np.loadtxt(DD_link + cat_name[4] + "_DD_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
red3_DD = np.loadtxt(DD_link + cat_name[5] + "_DD_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
blu3_DR = np.loadtxt(DR_link + cat_name[4] + "_DR_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))
red3_DR = np.loadtxt(DR_link + cat_name[5] + "_DR_auto_all.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13))



blu_DD_jack = np.loadtxt(DD_link + cat_name[0] + "_DD_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
red_DD_jack = np.loadtxt(DD_link + cat_name[1] + "_DD_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
blu_DR_jack = np.loadtxt(DR_link + cat_name[0] + "_DR_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
red_DR_jack = np.loadtxt(DR_link + cat_name[1] + "_DR_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
blu2_DD_jack = np.loadtxt(DD_link + cat_name[2] + "_DD_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
red2_DD_jack = np.loadtxt(DD_link + cat_name[3] + "_DD_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
blu2_DR_jack = np.loadtxt(DR_link + cat_name[2] + "_DR_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
red2_DR_jack = np.loadtxt(DR_link + cat_name[3] + "_DR_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
blu3_DD_jack = np.loadtxt(DD_link + cat_name[4] + "_DD_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
red3_DD_jack = np.loadtxt(DD_link + cat_name[5] + "_DD_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
blu3_DR_jack = np.loadtxt(DR_link + cat_name[4] + "_DR_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))
red3_DR_jack = np.loadtxt(DR_link + cat_name[5] + "_DR_auto_jack.dat", usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14))

blu_twopt = np.divide((blu_DD - 2 * blu_DR), RR_file) + 1
red_twopt = np.divide((red_DD - 2 * red_DR), RR_file) + 1
blu_twopt_jack = np.divide((blu_DD_jack - 2 * blu_DR_jack), RR_file_jack) + 1
red_twopt_jack = np.divide((red_DD_jack - 2 * red_DR_jack), RR_file_jack) + 1
blu2_twopt = np.divide((blu2_DD - 2 * blu2_DR), RR_file) + 1
red2_twopt = np.divide((red2_DD - 2 * red2_DR), RR_file) + 1
blu2_twopt_jack = np.divide((blu2_DD_jack - 2 * blu2_DR_jack), RR_file_jack) + 1
red2_twopt_jack = np.divide((red2_DD_jack - 2 * red2_DR_jack), RR_file_jack) + 1
blu3_twopt = np.divide((blu3_DD - 2 * blu3_DR), RR_file) + 1
red3_twopt = np.divide((red3_DD - 2 * red3_DR), RR_file) + 1
blu3_twopt_jack = np.divide((blu3_DD_jack - 2 * blu3_DR_jack), RR_file_jack) + 1
red3_twopt_jack = np.divide((red3_DD_jack - 2 * red3_DR_jack), RR_file_jack) + 1


print blu_twopt
print blu_twopt_jack[10]

blu_twopt_diff = blu_twopt - blu_twopt_jack
red_twopt_diff = red_twopt - red_twopt_jack
blu2_twopt_diff = blu2_twopt - blu2_twopt_jack
red2_twopt_diff = red2_twopt - red2_twopt_jack
blu3_twopt_diff = blu3_twopt - blu3_twopt_jack
red3_twopt_diff = red3_twopt - red3_twopt_jack

blu_sqdiff = [0] * jack_num
red_sqdiff = [0] * jack_num
blu2_sqdiff = [0] * jack_num
red2_sqdiff = [0] * jack_num
blu3_sqdiff = [0] * jack_num
red3_sqdiff = [0] * jack_num


for i in range(jack_num):
    blu_sqdiff[i] = np.outer(blu_twopt_diff[i], blu_twopt_diff[i])
    red_sqdiff[i] = np.outer(red_twopt_diff[i], red_twopt_diff[i])
    blu2_sqdiff[i] = np.outer(blu2_twopt_diff[i], blu2_twopt_diff[i])
    red2_sqdiff[i] = np.outer(red2_twopt_diff[i], red2_twopt_diff[i])
    blu3_sqdiff[i] = np.outer(blu3_twopt_diff[i], blu3_twopt_diff[i])
    red3_sqdiff[i] = np.outer(red3_twopt_diff[i], red3_twopt_diff[i])

blu_cov_matrix = np.matrix((np.sum(blu_sqdiff, axis = 0) * (jack_num - 1.) / jack_num))
red_cov_matrix = np.matrix((np.sum(red_sqdiff, axis = 0) * (jack_num - 1.) / jack_num))
blu2_cov_matrix = np.matrix((np.sum(blu2_sqdiff, axis = 0) * (jack_num - 1.) / jack_num))
red2_cov_matrix = np.matrix((np.sum(red2_sqdiff, axis = 0) * (jack_num - 1.) / jack_num))
blu3_cov_matrix = np.matrix((np.sum(blu3_sqdiff, axis = 0) * (jack_num - 1.) / jack_num))
red3_cov_matrix = np.matrix((np.sum(red3_sqdiff, axis = 0) * (jack_num - 1.) / jack_num))

#print blu_cov_matrix[0].size

blu_err_bar = np.power(np.diag(blu_cov_matrix), 0.5)
red_err_bar = np.power(np.diag(red_cov_matrix), 0.5)
blu2_err_bar = np.power(np.diag(blu2_cov_matrix), 0.5)
red2_err_bar = np.power(np.diag(red2_cov_matrix), 0.5)
blu3_err_bar = np.power(np.diag(blu3_cov_matrix), 0.5)
red3_err_bar = np.power(np.diag(red3_cov_matrix), 0.5)

blu_twopt_write = np.vstack((blu_twopt, blu_err_bar)).T
red_twopt_write = np.vstack((red_twopt, red_err_bar)).T
blu2_twopt_write = np.vstack((blu2_twopt, blu2_err_bar)).T
red2_twopt_write = np.vstack((red2_twopt, red2_err_bar)).T
blu3_twopt_write = np.vstack((blu3_twopt, blu3_err_bar)).T
red3_twopt_write = np.vstack((red3_twopt, red3_err_bar)).T

np.savetxt(mock_link + dir + "/" + cat_name[0] + "_2pt_auto.dat", blu_twopt_write, fmt = "%.7E", delimiter = "\t")
np.savetxt(mock_link + dir + "/" + cat_name[1] + "_2pt_auto.dat", red_twopt_write, fmt = "%.7E", delimiter = "\t")
np.savetxt(mock_link + dir + "/" + cat_name[2] + "_2pt_auto.dat", blu2_twopt_write, fmt = "%.7E", delimiter = "\t")
np.savetxt(mock_link + dir + "/" + cat_name[3] + "_2pt_auto.dat", red2_twopt_write, fmt = "%.7E", delimiter = "\t")
np.savetxt(mock_link + dir + "/" + cat_name[4] + "_2pt_auto.dat", blu3_twopt_write, fmt = "%.7E", delimiter = "\t")
np.savetxt(mock_link + dir + "/" + cat_name[5] + "_2pt_auto.dat", red3_twopt_write, fmt = "%.7E", delimiter = "\t")

print radius

fig = plt.figure()

gs1 = gridspec.GridSpec(1, 3)
ax1 = fig.add_subplot(gs1[0])
ax2 = fig.add_subplot(gs1[1])
ax3 = fig.add_subplot(gs1[2])

(_, caps1, _) = ax1.errorbar(radius,  blu_twopt, yerr = [ blu_err_bar,  blu_err_bar], color = 'b', fmt = 'o')
(_, caps2, _) = ax1.errorbar(radius,  red_twopt, yerr = [ red_err_bar,  red_err_bar], color = 'r', fmt = 'o')
(_, caps3, _) = ax2.errorbar(radius,  blu2_twopt, yerr = [ blu2_err_bar,  blu2_err_bar], color = 'b', fmt = 'o')
(_, caps4, _) = ax2.errorbar(radius,  red2_twopt, yerr = [ red2_err_bar,  red2_err_bar], color = 'r', fmt = 'o')
(_, caps5, _) = ax3.errorbar(radius,  blu3_twopt, yerr = [ blu3_err_bar,  blu3_err_bar], color = 'b', fmt = 'o')
(_, caps6, _) = ax3.errorbar(radius,  red3_twopt, yerr = [ red3_err_bar,  red3_err_bar], color = 'r', fmt = 'o')
for cap in caps1:
    cap.set_markeredgewidth(1.1)
for cap in caps2:
    cap.set_markeredgewidth(1.1)
for cap in caps3:
    cap.set_markeredgewidth(1.1)
for cap in caps4:
    cap.set_markeredgewidth(1.1)
for cap in caps5:
    cap.set_markeredgewidth(1.1)
for cap in caps6:
    cap.set_markeredgewidth(1.1)

ax1.set_ylabel(r'$w(r_p)$', fontsize = 12)
ax1.set_xlabel(r'$r_p [h^{-1} \rm{Mpc}]$', fontsize = 12)
ax2.set_xlabel(r'$r_p [h^{-1} \rm{Mpc}]$', fontsize = 12)
ax3.set_xlabel(r'$r_p [h^{-1} \rm{Mpc}]$', fontsize = 12)
#plt.title("Dots - without satellites, triangles - without satellites, auto correlation")
setp( ax2.get_yticklabels(), visible=False)
setp( ax3.get_yticklabels(), visible=False)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax1.set_xlim((1.3e-1, 60))
ax1.set_ylim((5e-5, 10))
ax2.set_xlim((1.3e-1, 60))
ax2.set_ylim((5e-5, 10))
ax3.set_xlim((1.3e-1, 60))
ax3.set_ylim((5e-5, 10))
#ax1.text(0.17, 1.3, r'Excluding $1.0r_{200c}$ and $200 h^{-1}$kpc' + '\n' + r'$f_{\rm{excluded}} = 0.47$', fontsize=10)
#ax2.text(0.17, 1.3, r'Excluding $1.2r_{200c}$ and $200 h^{-1}$kpc' + '\n' + r'$f_{\rm{excluded}} = 0.53$', fontsize=10)
#ax3.text(0.17, 1.3, r'Excluding $1.5r_{200c}$ and $200 h^{-1}$kpc' + '\n' + r'$f_{\rm{excluded}} = 0.61$', fontsize=10)
ax1.text(45, 0.6, r'$f_{\rm{excluded}} = 0.984$, $f_{\rm{sat}} = 0.014$' + '\n' + r'$M_{\ast,back}/M_{\ast,lens} \geq 1$', horizontalalignment='right', fontsize=10)
ax2.text(45, 0.6, r'$f_{\rm{excluded}} = 0.996$, $f_{\rm{sat}} = 0.0057$' + '\n' + r'$M_{\ast,back}/M_{\ast,lens} \geq 0.5$', horizontalalignment='right', fontsize=10)
ax3.text(45, 0.6, r'$f_{\rm{excluded}} = 0.998$, $f_{\rm{sat}} = 0.0044$' + '\n' + r'$M_{\ast,back}/M_{\ast,lens} \geq 0.1$', horizontalalignment='right', fontsize=10)
#ax1.set_title('Mock data', fontsize = 12)
#ax1.set_title('Mock data', fontsize = 12)
ax2.set_title(r'Mock Data, excluding $1.5\rm{r_{200c}}$ + $500 h^{-1}$kpc', fontsize = 12)
#ax3.set_title('Mock data', fontsize = 12)
gs1.tight_layout(fig, w_pad = -4.9, rect=[0, 0.4, 1, 1])
#plt.show()
plt.savefig(result_link + "mock_2pt_auto_1.0.eps",bbox_inches='tight')
#plt.savefig(piclink + "tew_2pt_all.pdf")

plt.close()

