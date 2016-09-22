### Check the mass range of the mock galaxies given in the catalog.
### Updated: Oct 10, 2015
### color argv[1]

import numpy as np
import sys, os, glob
import matplotlib.pyplot as plt

savelink = "/Users/rickyccy/Documents/Research_assemblybias/result/Dec_14_2015/"

color = sys.argv[1]
color2 = sys.argv[2]

bins2 = np.arange(9.5, 14.6, 0.5)
bins_mid = np.arange(9.75, 14.26, 0.5)

#mock_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/erased_assembly_bias_mocks/"
mock_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/"

os.chdir(mock_link)

gal = np.loadtxt(color + "_iso.dat", unpack = True, usecols = (4,6,7))
gal2 = np.loadtxt(color2 + "_iso.dat", unpack = True, usecols = (4,6,7))
"""
M_r = gal[0]
print np.average(gal[1]), np.power(np.mean(np.power(gal[1],2)), 0.5)
"""
mass = np.log10(gal[1])
mass2 = np.log10(gal2[1])

print np.where(mass < 10.0)[0]

a = np.histogram(mass, bins = bins2)[0] / float(gal[0].size)
a2 = np.histogram(mass2, bins = bins2)[0] / float(gal2[0].size)

plt.plot(bins_mid, a, color = color[0])
plt.plot(bins_mid, a2, color = color2[0])
plt.scatter(bins_mid, a, color = color[0])
plt.scatter(bins_mid, a2, color = color2[0])

#plt.bar(bins_mid, a, width = 0.5, align = 'center', color = color[0], alpha = 0.7)
#plt.bar(bins_mid, a2, width = 0.5, align = 'center', color = color2[0], alpha = 0.7)
#plt.yscale('log')
plt.xlabel(r'$\log_{10}(M / h^{-1} \rm{M_{\odot}})$', fontsize=12)
plt.ylabel(r'$N(M)$', fontsize=12)
#plt.vlines(x = fmass, ymin = 1, ymax = 1e5, linestyles = 'dashed')
#plt.title(color + ', ' + withornot + ' assembly bias, 1 Mpc separation')
plt.xlim([9.0,14.25])
plt.ylim([0.0,1.0])
plt.xticks(bins2)
plt.show()
#plt.savefig(savelink + "iso_mock_line.eps")
plt.close()


