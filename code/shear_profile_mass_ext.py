### This program saves the shear profile. Inverse covariance matrix is calculated as well.
### Input: quartile bin (argv[1])
### Updated: Jun 13, 2016

import numpy as np
import matplotlib.pyplot as plt
import sys
from numpy import linalg

#np.set_printoptions(threshold='nan')
#np.set_printoptions(precision=6)

quar_ind = sys.argv[1]


quar_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/quartiles/"
shearlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/"
dsig_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/Delta_Sig/"

r_bin = np.loadtxt(shearlink + "dist_bin2.dat", unpack = True, usecols = (1,)) * 1000

f_num = 4       ### Number of fields.

fstr = "mass_quartile{}".format(quar_ind)

fin = [0] * f_num
smfile = [0] * f_num


for i in range(f_num):
    fin[i] = np.loadtxt(quar_link + "W{}_{}_shear.dat".format(i+1,fstr), unpack = True, usecols = range(1,43,1))
    smfile[i] = np.loadtxt(quar_link + "W{}_{}_pos_temp.dat".format(i+1,fstr), unpack = True, usecols = (5,))

infile3 = np.hstack((fin[0],fin[1]))
infile3 = np.hstack((infile3, fin[2]))
infile3 = np.hstack((infile3, fin[3]))
sm = np.hstack((smfile[0],smfile[1]))
sm = np.hstack((sm,smfile[2]))
sm = np.hstack((sm,smfile[3]))

print np.power(10,np.median(sm))

sm = np.mean(np.power(10, sm))
print sm

exit(1)

#print np.log10(sm)


rad2 = np.multiply(np.pi, np.power(r_bin * 1000, 2))

### Stellar mass shear
sm_shear = np.divide(sm, rad2)


gg1 = []
gg2 = []
gg3 = []
gg4 = []
gg5 = []
gg6 = []
gg7 = []
gg8 = []
gg9 = []
gg10 = []
gg11 = []
gg12 = []
gg13 = []
gg14 = []
n1 = []
n2 = []
n3 = []
n4 = []
n5 = []
n6 = []
n7 = []
n8 = []
n9 = []
n10 = []
n11 = []
n12 = []
n13 = []
n14 = []

### Number of radial bins and bootstraps
bin_num = 14     ### 14
N_ran = 1000
N_size = len(infile3[0])


weighted_gg = [0] * bin_num
weighted_err_n = [0] * bin_num

print "g1"
for i in range(len(infile3[0])):
    if (np.isnan(infile3[0][i]) == True):
        gg1.append(0)
        n1.append(infile3[0 + bin_num][i])
    else:
        gg1.append(infile3[0][i])
        n1.append(infile3[0 + bin_num][i])

if (len(gg1) != len(n1)): print "g1 does not match n1!"

weighted_gg[0] = np.average(gg1, weights = n1)
weighted_err_n[0] = np.power((np.average(np.power(gg1,2), weights = n1) - np.power(weighted_gg[0], 2)) / (len(gg1) - 1), 0.5)




print "g2"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[1][i]) == True):
        gg2.append(0)
        n2.append(infile3[1 + bin_num][i])
    else:
        gg2.append(infile3[1][i])
        n2.append(infile3[1 + bin_num][i])

if (len(gg2) != len(n2)): print "g2 does not match e2!"

weighted_gg[1] = np.average(gg2, weights = n2)
weighted_err_n[1] = np.power((np.average(np.power(gg2,2), weights = n2) - np.power(weighted_gg[1], 2)) / (len(gg2) - 1), 0.5)



print "g3"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[2][i]) == True):
        gg3.append(0)
        n3.append(infile3[2 + bin_num][i])
    else:
        gg3.append(infile3[2][i])
        n3.append(infile3[2 + bin_num][i])

if (len(gg3) != len(n3)): print "g3 does not match e3!"

weighted_gg[2] = np.average(gg3, weights = n3)
weighted_err_n[2] = np.power((np.average(np.power(gg3,2), weights = n3) - np.power(weighted_gg[2], 2)) / (len(gg3) - 1), 0.5)



print "g4"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[3][i]) == True):
        gg4.append(0)
        n4.append(infile3[3 + bin_num][i])
    else:
        gg4.append(infile3[3][i])
        n4.append(infile3[3 + bin_num][i])

if (len(gg4) != len(n4)): print "g4 does not match e4!"

weighted_gg[3] = np.average(gg4, weights = n4)
weighted_err_n[3] = np.power((np.average(np.power(gg4,2), weights = n4) - np.power(weighted_gg[3], 2)) / (len(gg4) - 1), 0.5)




print "g5"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[4][i]) == True):
        gg5.append(0)
        n5.append(infile3[4 + bin_num][i])
    else:
        gg5.append(infile3[4][i])
        n5.append(infile3[4 + bin_num][i])

if (len(gg5) != len(n5)): print "g5 does not match e5!"

weighted_gg[4] = np.average(gg5, weights = n5)
weighted_err_n[4] = np.power((np.average(np.power(gg5,2), weights = n5) - np.power(weighted_gg[4], 2)) / (len(gg5) - 1), 0.5)



print "g6"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[5][i]) == True):
        gg6.append(0)
        n6.append(infile3[5 + bin_num][i])
    else:
        gg6.append(infile3[5][i])
        n6.append(infile3[5 + bin_num][i])

if (len(gg6) != len(n6)): print "g6 does not match e6!"

weighted_gg[5] = np.average(gg6, weights = n6)
weighted_err_n[5] = np.power((np.average(np.power(gg6,2), weights = n6) - np.power(weighted_gg[5], 2)) / (len(gg6) - 1), 0.5)



print "g7"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[6][i]) == True):
        gg7.append(0)
        n7.append(infile3[6 + bin_num][i])
    else:
        gg7.append(infile3[6][i])
        n7.append(infile3[6 + bin_num][i])

if (len(gg7) != len(n7)): print "g7 does not match e7!"

weighted_gg[6] = np.average(gg7, weights = n7)
weighted_err_n[6] = np.power((np.average(np.power(gg7,2), weights = n7) - np.power(weighted_gg[6], 2)) / (len(gg7) - 1), 0.5)



print "g8"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[7][i]) == True):
        gg8.append(0)
        n8.append(infile3[7 + bin_num][i])
    else:
        gg8.append(infile3[7][i])
        n8.append(infile3[7 + bin_num][i])

if (len(gg8) != len(n8)): print "g8 does not match e8!"

weighted_gg[7] = np.average(gg8, weights = n8)
weighted_err_n[7] = np.power((np.average(np.power(gg8,2), weights = n8) - np.power(weighted_gg[7], 2)) / (len(gg8) - 1), 0.5)



print "g9"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[8][i]) == True):
        gg9.append(0)
        n9.append(infile3[22][i])
    else:
        gg9.append(infile3[8][i])
        n9.append(infile3[22][i])

if (len(gg9) != len(n9)): print "g9 does not match e9!"

weighted_gg[8] = np.average(gg9, weights = n9)
weighted_err_n[8] = np.power((np.average(np.power(gg9,2), weights = n9) - np.power(weighted_gg[8], 2)) / (len(gg9) - 1), 0.5)



print "g10"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[9][i]) == True):
        gg10.append(0)
        n10.append(infile3[23][i])
    else:
        gg10.append(infile3[9][i])
        n10.append(infile3[23][i])

if (len(gg10) != len(n10)): print "g10 does not match e10!"

weighted_gg[9] = np.average(gg10, weights = n10)
weighted_err_n[9] = np.power((np.average(np.power(gg10,2), weights = n10) - np.power(weighted_gg[9], 2)) / (len(gg10) - 1), 0.5)


print "g11"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[10][i]) == True):
        gg11.append(0)
        n11.append(infile3[24][i])
    else:
        gg11.append(infile3[10][i])
        n11.append(infile3[24][i])

if (len(gg11) != len(n11)): print "g11 does not match e11!"

weighted_gg[10] = np.average(gg11, weights = n11)
weighted_err_n[10] = np.power((np.average(np.power(gg11,2), weights = n11) - np.power(weighted_gg[10], 2)) / (len(gg11) - 1), 0.5)


print "g12"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[11][i]) == True):
        gg12.append(0)
        n12.append(infile3[25][i])
    else:
        gg12.append(infile3[11][i])
        n12.append(infile3[25][i])

if (len(gg12) != len(n12)): print "g12 does not match e12!"

weighted_gg[11] = np.average(gg12, weights = n12)
weighted_err_n[11] = np.power((np.average(np.power(gg12,2), weights = n12) - np.power(weighted_gg[11], 2)) / (len(gg12) - 1), 0.5)



print "g13"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[12][i]) == True):
        gg13.append(0)
        n13.append(infile3[26][i])
    else:
        gg13.append(infile3[12][i])
        n13.append(infile3[26][i])

if (len(gg13) != len(n13)): print "g13 does not match e13!"

weighted_gg[12] = np.average(gg13, weights = n13)
weighted_err_n[12] = np.power((np.average(np.power(gg13,2), weights = n13) - np.power(weighted_gg[12], 2)) / (len(gg13) - 1), 0.5)



print "g14"
for i in range(len(infile3[0])):
    if (i % 300000 == 0): print i
    if (np.isnan(infile3[13][i]) == True):
        gg14.append(0)
        n14.append(infile3[27][i])
    else:
        gg14.append(infile3[13][i])
        n14.append(infile3[27][i])

if (len(gg14) != len(n14)): print "g14 does not match e14!"

weighted_gg[13] = np.average(gg14, weights = n14)
weighted_err_n[13] = np.power((np.average(np.power(gg14,2), weights = n14) - np.power(weighted_gg[13], 2)) / (len(gg14) - 1), 0.5)




mat = [0] * bin_num

for i in range(bin_num):
    mat[i] = [0] * 2


mat[0][0] = gg1
mat[1][0] = gg2
mat[2][0] = gg3
mat[3][0] = gg4
mat[4][0] = gg5
mat[5][0] = gg6
mat[6][0] = gg7
mat[7][0] = gg8

mat[8][0] = gg9
mat[9][0] = gg10
mat[10][0] = gg11
mat[11][0] = gg12
mat[12][0] = gg13
mat[13][0] = gg14


mat[0][1] = n1
mat[1][1] = n2
mat[2][1] = n3
mat[3][1] = n4
mat[4][1] = n5
mat[5][1] = n6
mat[6][1] = n7
mat[7][1] = n8

mat[8][1] = n9
mat[9][1] = n10
mat[10][1] = n11
mat[11][1] = n12
mat[12][1] = n13
mat[13][1] = n14


mat = np.array(mat)


#del gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, gg10, gg11, gg12, gg13, gg14, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14




### The array to store up the bootstrap.
g_t = [0] * bin_num

for i in range(bin_num):
    g_t[i] = [0] * N_ran

### Bootstrap
boot_num = [0] * N_ran

for i in range(N_ran):
    boot_num[i] = np.random.randint(0, N_size, size = N_size)

for j in range(bin_num):
    print j
    for i in range(N_ran):
        g_t[j][i] = np.average(np.array(mat[j][0])[boot_num[i]], weights = np.array(mat[j][1])[boot_num[i]])
        #print i,j,boot_num[i],mat[j][0][boot_num[i]],mat[j][1][boot_num[i]],g_t[j][i]

cov_g_t = np.cov(g_t)
sd_g_t = np.power(np.diag(cov_g_t), 0.5)

#inv_cov_g_t = linalg.inv(cov_g_t)

#print np.diag(np.dot(cov_g_t, inv_cov_g_t))

#np.savetxt(galshearlink + "invcov_" + fstr + "_gg_wider.dat", inv_cov_g_t, fmt="%.9E", delimiter = "\t", newline = "\n")

cov_fit_g_t = cov_g_t[0:7,0:7]
inv_cov_fit_g_t = linalg.inv(cov_fit_g_t)

print np.diag(np.dot(cov_fit_g_t, inv_cov_fit_g_t))

np.savetxt(dsig_link + "cov_fit_" + fstr + "_gg.dat", cov_fit_g_t, fmt="%.9E", delimiter = "\t", newline = "\n")
np.savetxt(dsig_link + "invcov_fit_" + fstr + "_gg.dat", inv_cov_fit_g_t, fmt="%.9E", delimiter = "\t", newline = "\n")

del N_size, g_t, boot_num, mat, cov_g_t#, inv_cov_g_t

gg_matrix = np.vstack((np.vstack((weighted_gg, sd_g_t)),sm_shear)).T

np.savetxt(dsig_link + fstr + "_gg.dat", gg_matrix, fmt="%.6E", delimiter = "\t", newline = "\n")

del gg_matrix

"""
plt.errorbar(r_bin, weighted_gg, yerr = sd_g_t)
plt.ylabel(r'$\Delta\Sigma [h_{70} M_\odot pc^{-2}]$', fontsize = 20)
plt.xlabel(r'r $[h_{70}^{-1} kpc]$', fontsize = 20)
plt.yscale('log')
plt.xscale('log')
plt.ylim([1e-1, 4e3])
plt.xlim([10, 6e2])
#plt.ylim([1e-1, 1e2])
#plt.xlim([10, 2e2])
plt.title("W" + sys.argv[1])
plt.plot()
plt.savefig(savelink + "W" + sys.argv[1] + "dSig_weightederr.pdf")

plt.close()
"""
