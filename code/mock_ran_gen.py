### This program generates randoms for the mock catalog.
### Updated: Aug 5, 2015

import numpy as np
import matplotlib.pyplot as plt
import os, glob, sys

x = np.arange(0, 251, 25)

N_size = 200000

ran_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/ran_position/"

os.chdir(ran_link)

x_coord = np.random.uniform(low = 0.0, high = 250.0, size = N_size)
y_coord = np.random.uniform(low = 0.0, high = 250.0, size = N_size)
z_coord = np.random.uniform(low = 0.0, high = 250.0, size = N_size)

x_inds = np.digitize(x_coord, x, right = True) - 1
y_inds = np.digitize(y_coord, x, right = True) - 1
inds = x_inds + 10 * y_inds

coord = np.vstack((x_coord, y_coord, z_coord, inds)).T

np.savetxt("ran_color.dat", coord, fmt = "%.5f\t%.5f\t%.5f\t%d", newline = "\n")

del coord, x, x_coord, y_coord, z_coord, x_inds, y_inds, inds
