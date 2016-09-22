### Find the M_200c and calculates the concentrations of the subhalos and host halos.
### Updated: May 26, 2016

import numpy as np
import os, sys, re
import lensing_cosmology as lenco

wlens = lenco.weak_lensing()

mock_list = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/mock/"

os.chdir (mock_list)

### Critical density at z = 0
den_c = wlens.rho_crit(0)

Mr_cons = 800 * np.pi / 3 * den_c

sm_hfile = np.loadtxt("sm_halo_catalog.dat", unpack = True, usecols = (0,2,3), dtype = {'names': ('id','rs','m200c'), 'formats': ('i8','float64','float64')})

### Convert to kpc/h
r_200c = np.power(sm_hfile[2] / Mr_cons, 1./3) / 1e3

### Concentrations
conc = np.divide(r_200c, sm_hfile[1])

sm_hfile = np.vstack((sm_hfile, conc)).T

np.savetxt("sm_halo_catalog_withconc.dat", sm_hfile, delimiter = "\t", fmt = "%lu\t%g\t%g\t%g")



