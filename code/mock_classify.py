### This program classifies the galaxies into blue and red galaxies, may separate by luminosities as well.
### This is for the catalog with assembly bias.
### Updated: A (argv[1]), rad (argv[2]), star mass ratio (argv[3]), blue upper mass (argv[4]), blue low mass (argv[5]), red up mass (argv[6]), red low mass (argv[7])
### Updated: Dec 15, 2015

import numpy as np
import os, glob, re, itertools, sys

flink = "/Users/rickyccy/Documents/Research_assemblybias/mock/"
os.chdir(flink)

fac = sys.argv[1]
rad_rep = sys.argv[2]
sm = sys.argv[3]
b_up = float(sys.argv[4])
b_low = float(sys.argv[5])
r_up = float(sys.argv[6])
r_low = float(sys.argv[7])

grid = np.arange(0, 251, 25)

#infile = np.loadtxt("Mr19_age_distribution_matching_mock.dat", usecols =(0,1,2,3,7,9,10,11,12), unpack = True, dtype = {'names': ('halo_id', 'x', 'y', 'z', 'm_vir', 'M_r', 'g_r', 'm_host', 'host_id'), 'formats':('S10', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'S10')})
infile = np.loadtxt("sm_probab_" + fac + "_" + rad_rep + "_" + sm + ".dat", usecols =(0,1,2,3,7,9,10,11,12), unpack = True, dtype = {'names': ('halo_id', 'x', 'y', 'z', 'm_vir', 'star_m', 'g_r', 'm_host', 'host_id'), 'formats':('S10', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'S10')})
#infile = np.loadtxt("Mr19_screen.dat", usecols =(0,1,2,3,7,9,10,11,12), unpack = True, dtype = {'names': ('halo_id', 'x', 'y', 'z', 'm_vir', 'M_r', 'g_r', 'm_host', 'host_id'), 'formats':('S10', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'S10')})

halo_id = infile[0].astype(np.int)
x = infile[1]
y = infile[2]
z = infile[3]
m_vir = infile[4]
star_m = infile[5]
g_r = infile[6]
m_host = infile[7]
host_id = infile[8].astype(np.int)

"""
M19_arr = np.where((M_r <= -19) & (M_r > -20))[0]
M20_arr = np.where((M_r <= -20) & (M_r > -21))[0]
M21_arr = np.where((M_r <= -21) & (M_r >= -22))[0]
M19_x_inds = np.digitize(x[M19_arr], grid, right = True) - 1
M19_y_inds = np.digitize(y[M19_arr], grid, right = True) - 1
M19_inds = M19_x_inds + 10 * M19_y_inds
M20_x_inds = np.digitize(x[M20_arr], grid, right = True) - 1
M20_y_inds = np.digitize(y[M20_arr], grid, right = True) - 1
M20_inds = M20_x_inds + 10 * M20_y_inds
M21_x_inds = np.digitize(x[M21_arr], grid, right = True) - 1
M21_y_inds = np.digitize(y[M21_arr], grid, right = True) - 1
M21_inds = M21_x_inds + 10 * M21_y_inds

"""
color_line = 0.76 + 0.15 * (star_m - 10.0)

color = g_r - color_line


#blu_arr = np.where((color < 0) & (host_id == -1) & (m_host >= 3.67e11) & (m_host <= 5.51e11))[0]
#red_arr = np.where((color >= 0) & (host_id == -1) & (m_host >= 3.67e11) & (m_host <= 5.51e11))[0]
#other_arr = np.where((host_id != -1) | (host_id == -1) & (m_host < 3.67e11) | (host_id == -1) & (m_host > 5.51e11))[0]

#blu_arr = np.where((color < 0) & (m_vir >= 3.67e11) & (m_vir <= 5.51e11))[0]
#red_arr = np.where((color >= 0) & (m_vir >= 3.67e11) & (m_vir <= 5.51e11))[0]
#other_arr = np.where((m_vir < 3.67e11) | (m_vir > 5.51e11))[0]
blu_arr = np.where((color < 0) & (star_m >= b_low) & (star_m <= b_up))[0]
red_arr = np.where((color >= 0) & (star_m >= r_low) & (star_m <= r_up))[0]
#other_arr = np.where((m_vir < 3.67e11) | (m_vir > 5.51e11))[0]

print blu_arr
print red_arr

blu_x_inds = np.digitize(x[blu_arr], grid, right = True) - 1
blu_y_inds = np.digitize(y[blu_arr], grid, right = True) - 1
blu_inds = blu_x_inds + 10 * blu_y_inds

red_x_inds = np.digitize(x[red_arr], grid, right = True) - 1
red_y_inds = np.digitize(y[red_arr], grid, right = True) - 1
red_inds = red_x_inds + 10 * red_y_inds

#oth_x_inds = np.digitize(x[other_arr], grid, right = True) - 1
#oth_y_inds = np.digitize(y[other_arr], grid, right = True) - 1
#oth_inds = oth_x_inds + 10 * oth_y_inds

blu = np.vstack((halo_id[blu_arr], x[blu_arr], y[blu_arr], z[blu_arr], star_m[blu_arr], g_r[blu_arr], m_vir[blu_arr], host_id[blu_arr], blu_inds, m_host[blu_arr])).T
red = np.vstack((halo_id[red_arr], x[red_arr], y[red_arr], z[red_arr], star_m[red_arr], g_r[red_arr], m_vir[red_arr], host_id[red_arr], red_inds, m_host[red_arr])).T
#oth = np.vstack((halo_id[other_arr], x[other_arr], y[other_arr], z[other_arr], M_r[other_arr], g_r[other_arr], m_vir[other_arr], host_id[other_arr], oth_inds)).T

np.savetxt("blu_" + fac + "_" + rad_rep + "_" + sm + ".dat", blu, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3e\t%d\t%d\t%.3e', delimiter = '\t', newline = '\n')
np.savetxt("red_" + fac + "_" + rad_rep + "_" + sm + ".dat", red, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3e\t%d\t%d\t%.3e', delimiter = '\t', newline = '\n')
#np.savetxt("other_tobe.dat", oth, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')

#del color_line, color, blu_arr, halo_id, x, y, z, m_vir, M_r, m_host, g_r, host_id, blu, blu_x_inds, blu_y_inds, blu_inds#, red, red_arr, red_x_inds, red_y_inds, red_inds#, other_arr, oth, oth_x_inds, oth_y_inds, oth_inds

"""
M19_blu_arr = np.where((color < 0) & (M_r <= -19) & (M_r > -20))[0]
M19_red_arr = np.where((color > 0) & (M_r <= -19) & (M_r > -20))[0]
M20_blu_arr = np.where((color < 0) & (M_r <= -20) & (M_r > -21))[0]
M20_red_arr = np.where((color > 0) & (M_r <= -20) & (M_r > -21))[0]
M21_blu_arr = np.where((color < 0) & (M_r <= -21) & (M_r >= -22))[0]
M21_red_arr = np.where((color > 0) & (M_r <= -21) & (M_r >= -22))[0]

### Jackknife
M19_blu_x_inds = np.digitize(x[M19_blu_arr], grid, right = True) - 1
M19_blu_y_inds = np.digitize(y[M19_blu_arr], grid, right = True) - 1
M19_blu_inds = M19_blu_x_inds + 10 * M19_blu_y_inds
M19_red_x_inds = np.digitize(x[M19_red_arr], grid, right = True) - 1
M19_red_y_inds = np.digitize(y[M19_red_arr], grid, right = True) - 1
M19_red_inds = M19_red_x_inds + 10 * M19_red_y_inds
M20_blu_x_inds = np.digitize(x[M20_blu_arr], grid, right = True) - 1
M20_blu_y_inds = np.digitize(y[M20_blu_arr], grid, right = True) - 1
M20_blu_inds = M20_blu_x_inds + 10 * M20_blu_y_inds
M20_red_x_inds = np.digitize(x[M20_red_arr], grid, right = True) - 1
M20_red_y_inds = np.digitize(y[M20_red_arr], grid, right = True) - 1
M20_red_inds = M20_red_x_inds + 10 * M20_red_y_inds
M21_blu_x_inds = np.digitize(x[M21_blu_arr], grid, right = True) - 1
M21_blu_y_inds = np.digitize(y[M21_blu_arr], grid, right = True) - 1
M21_blu_inds = M21_blu_x_inds + 10 * M21_blu_y_inds
M21_red_x_inds = np.digitize(x[M21_red_arr], grid, right = True) - 1
M21_red_y_inds = np.digitize(y[M21_red_arr], grid, right = True) - 1
M21_red_inds = M21_red_x_inds + 10 * M21_red_y_inds


M19_blu = np.vstack((halo_id[M19_blu_arr], x[M19_blu_arr], y[M19_blu_arr], z[M19_blu_arr], m_vir[M19_blu_arr], M_r[M19_blu_arr], g_r[M19_blu_arr], m_host[M19_blu_arr], host_id[M19_blu_arr], M19_blu_inds)).T
M20_blu = np.vstack((halo_id[M20_blu_arr], x[M20_blu_arr], y[M20_blu_arr], z[M20_blu_arr], m_vir[M20_blu_arr], M_r[M20_blu_arr], g_r[M20_blu_arr], m_host[M20_blu_arr], host_id[M20_blu_arr], M20_blu_inds)).T
M21_blu = np.vstack((halo_id[M21_blu_arr], x[M21_blu_arr], y[M21_blu_arr], z[M21_blu_arr], m_vir[M21_blu_arr], M_r[M21_blu_arr], g_r[M21_blu_arr], m_host[M21_blu_arr], host_id[M21_blu_arr], M21_blu_inds)).T
M19_red = np.vstack((halo_id[M19_red_arr], x[M19_red_arr], y[M19_red_arr], z[M19_red_arr], m_vir[M19_red_arr], M_r[M19_red_arr], g_r[M19_red_arr], m_host[M19_red_arr], host_id[M19_red_arr], M19_red_inds)).T
M20_red = np.vstack((halo_id[M20_red_arr], x[M20_red_arr], y[M20_red_arr], z[M20_red_arr], m_vir[M20_red_arr], M_r[M20_red_arr], g_r[M20_red_arr], m_host[M20_red_arr], host_id[M20_red_arr], M20_red_inds)).T
M21_red = np.vstack((halo_id[M21_red_arr], x[M21_red_arr], y[M21_red_arr], z[M21_red_arr], m_vir[M21_red_arr], M_r[M21_red_arr], g_r[M21_red_arr], m_host[M21_red_arr], host_id[M21_red_arr], M21_red_inds)).T

#del halo_id, x, y, z, m_vir, M_r, g_r, m_host, host_id, M19_arr, M20_arr, M21_arr, M19_color_line, M19_color, M19_blu_arr, M19_red_arr, M20_color_line, M20_color, M20_blu_arr, M20_red_arr, M21_color_line, M21_color, M21_blu_arr, M21_red_arr

np.savetxt("M19_blu.dat", M19_blu, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')
np.savetxt("M20_blu.dat", M20_blu, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')
np.savetxt("M21_blu.dat", M21_blu, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')
np.savetxt("M19_red.dat", M19_red, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')
np.savetxt("M20_red.dat", M20_red, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')
np.savetxt("M21_red.dat", M21_red, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')

del M19_blu, M19_red, M20_blu, M20_red, M21_blu, M21_red


M19 = np.vstack((halo_id[M19_arr], x[M19_arr], y[M19_arr], z[M19_arr], m_vir[M19_arr], M_r[M19_arr], g_r[M19_arr], m_host[M19_arr], host_id[M19_arr], M19_inds)).T
M20 = np.vstack((halo_id[M20_arr], x[M20_arr], y[M20_arr], z[M20_arr], m_vir[M20_arr], M_r[M20_arr], g_r[M20_arr], m_host[M20_arr], host_id[M20_arr], M20_inds)).T
M21 = np.vstack((halo_id[M21_arr], x[M21_arr], y[M21_arr], z[M21_arr], m_vir[M21_arr], M_r[M21_arr], g_r[M21_arr], m_host[M21_arr], host_id[M21_arr], M21_inds)).T

#del halo_id, x, y, z, m_vir, M_r, g_r, m_host, host_id, M19_arr, M20_arr, M21_arr, M19_color_line, M19_color, M19_blu_arr, M19_red_arr, M20_color_line, M20_color, M20_blu_arr, M20_red_arr, M21_color_line, M21_color, M21_blu_arr, M21_red_arr

np.savetxt("M19.dat", M19, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')
np.savetxt("M20.dat", M20, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')
np.savetxt("M21.dat", M21, fmt = '%d\t%.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%.3f\t%.3e\t%d\t%d', delimiter = '\t', newline = '\n')

del M19, M20, M21
"""