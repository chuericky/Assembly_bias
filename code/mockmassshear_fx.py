### Calculates the mass function of the mock halos.  d(gamma)/dM
### Updated: Jun 7, 2016

import numpy as np
import os, re, sys, glob
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

mock_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/mock/"
result_link = "/Users/rickyccy/Documents/Research_assemblybias/result/Jul_27_2016/mock_Sig/"

h = 0.7

os.chdir(mock_link)

hname, hmass = np.loadtxt("sm_1.5_0.5_3.0_catalog.dat", usecols = (0,14), unpack = True, dtype = {'names':('id','m200c'), 'formats':('i8','float64')})

### Convert the halo mass to h_70^-1M_sun instead of h^-1M_sun.
hmass = hmass / h

hmass = np.log10(hmass)

### From M < 10.5, 10.5 <= M < 10.7, ..., M >= 13.7
mass_bin = np.arange(10.5,13.8,0.2)
mass_pltbin = np.arange(10.4,13.9,0.2)

"""
mass_pbin = np.arange(10.6,13.7,0.2)

hist_res, bin_edges_res = np.histogram(hmass, bins = mass_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(mass_pbin, hist_res, 0.2, color = 'm', alpha = 0.3)
plt.show()
plt.close()
"""

inds = np.digitize(hmass, mass_bin)

#print mass_bin

#print inds

#print hmass

#print np.bincount(inds)


shear_bin = 35

### Shear profiles.
hshear = np.loadtxt("sm_1.5_0.5_3.0_shear.dat", usecols = np.arange(1,shear_bin + 1,1))

### Convert the shear profiles from h M_sun pc^-2 to h_70 M_sun pc^-2
hshear = hshear * h

mass_binnum = 18

avg_shear = [0] * mass_binnum

### Find the average shear at each mass bin.
for i in range(mass_binnum):  ### mass_binnum
    pos = np.where((inds == i))[0]
    avg_shear[i] = np.average(hshear[pos], axis = 0)

avg_shear_mat = np.array(avg_shear).T

r = np.arange(0.01,0.36,0.01) / h

m_var = np.var(hmass, ddof = 1)

mass_derv = np.vstack((mass_pltbin, avg_shear_mat)).T

mass_s = np.arange(10.3,13.9,0.02)

derv_weighted_avg = [0] * shear_bin
derv_avg = [0] * shear_bin
spd = [0] * shear_bin


for j in range(shear_bin):
    f, (ax1,ax2) = plt.subplots(1,2, figsize = (10,5))
    ### Deg. = 2
    """
    spl2 = IUS(mass_pltbin, avg_shear_mat[j], k = 2)
    spld2 = spl2.derivative()
    """
    spl3 = IUS(mass_pltbin, avg_shear_mat[j], k = 3)
    spld3 = spl3.derivative()
    """
    spl4 = IUS(mass_pltbin, avg_shear_mat[j], k = 4)
    spld4 = spl4.derivative()
    spl5 = IUS(mass_pltbin, avg_shear_mat[j], k = 5)
    spld5 = spl5.derivative()
    """
    #ax1.plot(mass_pltbin, avg_shear_mat[j], color = 'r')
    #ax1.plot(mass_s, spl2(mass_s), color = 'b')
    ax1.plot(mass_s, spl3(mass_s), color = 'r')
    #ax1.plot(mass_s, spl4(mass_s), color = 'c')
    #ax1.plot(mass_s, spl5(mass_s), color = 'r')
    #ax1.legend(["k = 2", "k = 3", "k = 4", "k = 5"], loc = 'upper left')
    ax1.scatter(mass_pltbin, avg_shear_mat[j], marker = 'x', color = 'r')
    ax1.set_xlabel(r"$\log(M_h) (M_{\odot})$")
    ax1.set_ylabel(r"$\langle\Delta\Sigma(r)\rangle [M_\odot pc^{-2}]$")
    ax1.set_title(r"Bin #{}, r = {} Mpc".format(j+1, r[j]))
    ax1.hlines(xmin = 10.2, xmax = 14, y = 0, linestyle = '-.')
    ax1.set_xlim((10.2,14))
    #ax1.yscale('log')
    #ax1.ylim((1,1000))
    ax1.set_ylim((-1,300))
    #ax1.savefig(result_link + "Sig_{}_linear_interpolated.png".format(j+101))
    #ax1.close()
    #ax2.plot(mass_s, spld2(mass_s), color = 'b')
    #ax2.plot(mass_s, spld4(mass_s), color = 'c')
    #ax2.plot(mass_s, spld5(mass_s), color = 'r')
    #ax2.scatter(mass_pltbin, spld2(mass_pltbin), color = 'b', marker = 'x')
    spd[j] = spld3(mass_pltbin)
    derv_weighted_avg[j] = np.average(spd[j], weights = np.bincount(inds))
    derv_avg[j] = np.mean(spd[j])
    #ax2.scatter(mass_pltbin, spld4(mass_pltbin), color = 'c', marker = 'x')
    #ax2.scatter(mass_pltbin, spld5(mass_pltbin), color = 'r', marker = 'x')
    ax2.hlines(xmin = 10.2, xmax = 14, y = 0, linestyle = '-.')
    ax2.hlines(xmin = 10.2, xmax = 14, y = derv_weighted_avg[j], linestyle = '-', color = 'g')
    ax2.hlines(xmin = 10.2, xmax = 14, y = derv_avg[j], linestyle = '--', color = 'r')
    ax2.legend(["y = 0", "Weighted average", "Simple mean"], loc = 'upper left')
    ax2.scatter(mass_pltbin, spd[j], color = 'b', marker = 'x')
    ax2.plot(mass_s, spld3(mass_s), color = 'b')
    ax2.set_xlabel(r"$\log(M_h) (M_{\odot})$")
    ax2.set_ylabel(r"$\frac{d \langle\Delta\Sigma(r)\rangle}{d\log M_h}$")
    ax2.set_title(r"Bin #{}, r = {} Mpc".format(j+1, r[j]))
    ax2.set_xlim((10.2,14))
    ax2.set_ylim((-50,500))
    plt.tight_layout(w_pad = 0.5)
    plt.savefig(result_link + "MassShear_{}.png".format(j+101), bbox_inches = 'tight')
    plt.close()

spd = np.matrix(spd).T

mass_derv = np.hstack((mass_derv, spd))

np.savetxt(mock_link + "mass_shear_new.dat", mass_derv, fmt = "%.10g", delimiter = "\t")


m_avg_vec = [np.mean(hmass)] * shear_bin

m_var_vec = [m_var] * shear_bin

vec = np.vstack((m_avg_vec, derv_weighted_avg))
vec = np.vstack((vec, derv_avg))
vec = np.vstack((vec, m_var_vec)).T

np.savetxt(mock_link + "shear_mass_derivative_new.dat", vec, fmt = "%.10g", delimiter = "\t")


