### This program calculates the 2pt function with jackknife errors of the combined field.
### Input: parameter name (with/without) (argv[1]), b (argv[2]), clus_num (argv[3]), P_th (argv[4])
### Updated: Jan 27, 2016

import numpy as np
import glob, os, re, sys
import matplotlib.pyplot as plt
from numpy import linalg
from scipy.stats import chisqprob
import scipy.stats as st

parameter_name = sys.argv[1]
b = sys.argv[2]
clus_num = sys.argv[3]
P_th = sys.argv[4]

parameter = b + "_" + P_th + "_" + clus_num + "_lensred"
parameter2 = b + "_" + P_th + "_" + clus_num + "_lensblu"

para_tot = b + "_" + P_th + "_" + clus_num

### Number of all lenses, randoms
lens_all_num = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/" + parameter_name + "/" + parameter + "_lensname.txt", usecols = (0,), unpack = True)
lens_all_num2 = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/" + parameter_name + "/" + parameter2 + "_lensname.txt", usecols = (0,), unpack = True)
ran_all_num = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/random_map/ran_num_iso.dat", usecols = (1,), unpack = True)

ranpairlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/"
piclink = "/Users/rickyccy/Documents/Research_assemblybias/result/Sept_29_2015/"

theta_range = np.loadtxt(ranpairlink + "radial_bin2.dat", unpack = True, usecols = (1,))

theta_power = np.power(theta_range, -0.8)

DD_pair = [0] * 4
DR_pair = [0] * 4
RR_pair = [0] * 4
ang_cor = [0] * 4
DD_data = [0] * 4
DR_data = [0] * 4
RR_data = [0] * 4
twopt = [0] * 4
DD_jackgalcount = [0] * 4
DR_jackgalcount = [0] * 4
RR_jackgalcount = [0] * 4
DD_jackgalpair = [0] * 4
DR_jackgalpair = [0] * 4
RR_jackgalpair = [0] * 4
DD_jack = [0] * 4
DR_jack = [0] * 4
RR_jack = [0] * 4
twopt_jack = [0] * 4

DD_pair2 = [0] * 4
DR_pair2 = [0] * 4
ang_cor2 = [0] * 4
DD_data2 = [0] * 4
DR_data2 = [0] * 4
twopt2 = [0] * 4
DD_jackgalcount2 = [0] * 4
DR_jackgalcount2 = [0] * 4
DD_jackgalpair2 = [0] * 4
DR_jackgalpair2 = [0] * 4
DD_jack2 = [0] * 4
DR_jack2 = [0] * 4
twopt_jack2 = [0] * 4

for i in range(4):
    DD_pair[i] = lens_all_num[i] * (lens_all_num[i] - 1)
    DR_pair[i] = lens_all_num[i] * ran_all_num[i]
    RR_pair[i] = ran_all_num[i] * (ran_all_num[i] - 1)
    
    DD_pair2[i] = lens_all_num2[i] * (lens_all_num2[i] - 1)
    DR_pair2[i] = lens_all_num2[i] * ran_all_num[i]

field_num = 4

for i in range(field_num):
    os.chdir("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/W" + str(i + 1) + "/combined/")

    """
    DD_galcount = np.loadtxt(parameter + "_DD_all.dat", unpack = True, usecols = (0,))
    DR_galcount = np.loadtxt(parameter + "_DR_all.dat", unpack = True, usecols = (0,))
    RR_galcount = np.loadtxt("RR_all.dat", unpack = True, usecols = (0,))
    """
    DD_galcount = np.loadtxt(parameter + "_DD_auto_all_wider.dat", unpack = True, usecols = (0,))
    DR_galcount = np.loadtxt(parameter + "_DR_auto_all_wider.dat", unpack = True, usecols = (0,))
    RR_galcount = np.loadtxt("iso_RR_auto_all_wider.dat", unpack = True, usecols = (0,))

    DD_galcount2 = np.loadtxt(parameter2 + "_DD_auto_all_wider.dat", unpack = True, usecols = (0,))
    DR_galcount2 = np.loadtxt(parameter2 + "_DR_auto_all_wider.dat", unpack = True, usecols = (0,))
    
    if (DD_galcount != lens_all_num[i]):
        print "Data-Data numbers don't match!"
        break
    if (DR_galcount != lens_all_num[i]):
        print "Data-Random numbers don't match!"
        break
    if (RR_galcount != ran_all_num[i]):
        print "Random-Random numbers don't match!"
        break

    DD_data[i] = np.loadtxt(parameter + "_DD_auto_all_wider.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,)) / DD_pair[i]
    DR_data[i] = np.loadtxt(parameter + "_DR_auto_all_wider.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,)) / DR_pair[i]
    DD_data2[i] = np.loadtxt(parameter2 + "_DD_auto_all_wider.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,)) / DD_pair2[i]
    DR_data2[i] = np.loadtxt(parameter2 + "_DR_auto_all_wider.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,)) / DR_pair2[i]
    RR_data[i] = np.loadtxt("iso_RR_auto_all_wider.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,)) / RR_pair[i]
    twopt[i] = np.divide((DD_data[i] - 2 * DR_data[i]), RR_data[i]) + 1
    twopt2[i] = np.divide((DD_data2[i] - 2 * DR_data2[i]), RR_data[i]) + 1

    ### Jackknife
    DD_jackgalcount[i] = np.loadtxt(parameter + "_DD_auto_jack_wider.dat", unpack = True, usecols = (1,))
    DR_jackgalcount[i] = np.loadtxt(parameter + "_DR_auto_jack_wider.dat", unpack = True, usecols = (1,))
    DD_jackgalcount2[i] = np.loadtxt(parameter2 + "_DD_auto_jack_wider.dat", unpack = True, usecols = (1,))
    DR_jackgalcount2[i] = np.loadtxt(parameter2 + "_DR_auto_jack_wider.dat", unpack = True, usecols = (1,))
    RR_jackgalcount[i] = np.loadtxt("iso_RR_auto_jack_wider.dat", unpack = True, usecols = (1,))
   
    ### Make sure the total number of galaxies is right in all jackknife regions.
    DD_jackcheck = len(DD_jackgalcount[i]) * lens_all_num[i] - np.sum(DD_jackgalcount[i])
    DR_jackcheck = len(DR_jackgalcount[i]) * lens_all_num[i] - np.sum(DR_jackgalcount[i])
    DD_jackcheck2 = len(DD_jackgalcount2[i]) * lens_all_num2[i] - np.sum(DD_jackgalcount2[i])
    DR_jackcheck2 = len(DR_jackgalcount2[i]) * lens_all_num2[i] - np.sum(DR_jackgalcount2[i])
    RR_jackcheck = len(RR_jackgalcount[i]) * ran_all_num[i] - np.sum(RR_jackgalcount[i])
    
    if (DD_jackcheck != lens_all_num[i]):
        print "Jackknife Data-Data numbers don't match!"
        break
    if (DR_jackcheck != lens_all_num[i]):
        print "Jackknife Data-Random numbers don't match!"
        break
    if (RR_jackcheck != ran_all_num[i]):
        print "Jackknife Random-Random numbers don't match!"
        break

    
    DD_jackgalpair[i] = np.multiply(DD_jackgalcount[i], (DD_jackgalcount[i] - 1))
    DR_jackgalpair[i] = np.multiply(DD_jackgalcount[i], RR_jackgalcount[i])
    DD_jackgalpair2[i] = np.multiply(DD_jackgalcount2[i], (DD_jackgalcount2[i] - 1))
    DR_jackgalpair2[i] = np.multiply(DD_jackgalcount2[i], RR_jackgalcount[i])
    RR_jackgalpair[i] = np.multiply(RR_jackgalcount[i], (RR_jackgalcount[i] - 1))

    DD_jack[i] = np.divide(np.loadtxt(parameter + "_DD_auto_jack_wider.dat", unpack = True, usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,)), DD_jackgalpair[i])
    DR_jack[i] = np.divide(np.loadtxt(parameter + "_DR_auto_jack_wider.dat", unpack = True, usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,)), DR_jackgalpair[i])
    DD_jack2[i] = np.divide(np.loadtxt(parameter2 + "_DD_auto_jack_wider.dat", unpack = True, usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,)), DD_jackgalpair2[i])
    DR_jack2[i] = np.divide(np.loadtxt(parameter2 + "_DR_auto_jack_wider.dat", unpack = True, usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,)), DR_jackgalpair2[i])
    RR_jack[i] = np.divide(np.loadtxt("iso_RR_auto_jack_wider.dat", unpack = True, usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,)), RR_jackgalpair[i])

    twopt_jack[i] = (np.divide((DD_jack[i] - 2 * DR_jack[i]), RR_jack[i]) + 1).T
    twopt_jack2[i] = (np.divide((DD_jack2[i] - 2 * DR_jack2[i]), RR_jack[i]) + 1).T

    del DD_galcount, DR_galcount, RR_galcount, DD_jackcheck, DR_jackcheck, RR_jackcheck

### The combined two pt function.
com_twopt = np.average(twopt, weights = DD_pair, axis = 0)
com_twopt2 = np.average(twopt2, weights = DD_pair2, axis = 0)

com_diff = com_twopt - com_twopt2

DD_jackpair = [0] * 54
jack_twopt = [0] * 54
DD_jackpair2 = [0] * 54
jack_twopt2 = [0] * 54

for i in range(54):
    DD_jackpair[i] = [0] * 4
    jack_twopt[i] = [0] * 4
    DD_jackpair2[i] = [0] * 4
    jack_twopt2[i] = [0] * 4

for i in range(18):
    DD_jackpair[i][0] = DD_jackgalpair[0][i]
    jack_twopt[i][0] = twopt_jack[0][i]
    DD_jackpair2[i][0] = DD_jackgalpair2[0][i]
    jack_twopt2[i][0] = twopt_jack2[0][i]
    for j in range(1,4):
        DD_jackpair[i][j] = DD_pair[j]
        jack_twopt[i][j] = twopt[j]
        DD_jackpair2[i][j] = DD_pair2[j]
        jack_twopt2[i][j] = twopt2[j]

for i in range(18, 24):
    DD_jackpair[i][0] = DD_pair[0]
    jack_twopt[i][0] = twopt[0]
    DD_jackpair2[i][0] = DD_pair2[0]
    jack_twopt2[i][0] = twopt2[0]
    for j in range(2,4):
        DD_jackpair[i][j] = DD_pair[j]
        jack_twopt[i][j] = twopt[j]
        DD_jackpair2[i][j] = DD_pair2[j]
        jack_twopt2[i][j] = twopt2[j]
    DD_jackpair[i][1] = DD_jackgalpair[1][i - 18]
    jack_twopt[i][1] = twopt_jack[1][i - 18]
    DD_jackpair2[i][1] = DD_jackgalpair2[1][i - 18]
    jack_twopt2[i][1] = twopt_jack2[1][i - 18]

for i in range(24, 48):
    for j in range(2):
        DD_jackpair[i][j] = DD_pair[j]
        jack_twopt[i][j] = twopt[j]
        DD_jackpair2[i][j] = DD_pair2[j]
        jack_twopt2[i][j] = twopt2[j]
    DD_jackpair[i][3] = DD_pair[3]
    jack_twopt[i][3] = twopt[3]
    DD_jackpair2[i][3] = DD_pair2[3]
    jack_twopt2[i][3] = twopt2[3]
    DD_jackpair[i][2] = DD_jackgalpair[2][i - 24]
    jack_twopt[i][2] = twopt_jack[2][i - 24]
    DD_jackpair2[i][2] = DD_jackgalpair2[2][i - 24]
    jack_twopt2[i][2] = twopt_jack2[2][i - 24]

for i in range(48, 54):
    for j in range(3):
        DD_jackpair[i][j] = DD_pair[j]
        jack_twopt[i][j] = twopt[j]
        DD_jackpair2[i][j] = DD_pair2[j]
        jack_twopt2[i][j] = twopt2[j]
    DD_jackpair[i][3] = DD_jackgalpair[3][i - 48]
    jack_twopt[i][3] = twopt_jack[3][i - 48]
    DD_jackpair2[i][3] = DD_jackgalpair2[3][i - 48]
    jack_twopt2[i][3] = twopt_jack2[3][i - 48]

jackcom_twopt = [0] * 54
sq_diff = [0] * 54
jackcom_twopt2 = [0] * 54
sq_diff2 = [0] * 54

jackcom_twopt_diff = [0] * 54

for i in range(54):
    ### Combined jackknife two-pt function.
    jackcom_twopt[i] = np.average(jack_twopt[i], weights = DD_jackpair[i], axis = 0)
    jackcom_twopt2[i] = np.average(jack_twopt2[i], weights = DD_jackpair2[i], axis = 0)

for i in range(54):
    jackcom_twopt_diff[i] = jackcom_twopt[i] - jackcom_twopt2[i]



twopt_diff = com_twopt - jackcom_twopt
twopt_diff2 = com_twopt2 - jackcom_twopt2

twopt_tot_diff = com_diff - jackcom_twopt_diff

sq_tot_diff = [0] * 54

for i in range(54):
    sq_diff[i] = np.outer(twopt_diff[i], twopt_diff[i])
    sq_diff2[i] = np.outer(twopt_diff2[i], twopt_diff2[i])
    sq_tot_diff[i] = np.outer(twopt_tot_diff[i], twopt_tot_diff[i])

cov_matrix = np.matrix((np.sum(sq_diff, axis = 0) * 53 / 54.))
cov_matrix2 = np.matrix((np.sum(sq_diff2, axis = 0) * 53 / 54.))
cov_matrix_tot = np.matrix((np.sum(sq_tot_diff, axis = 0) * 53 / 54.))

err_bar = np.power(np.diag(cov_matrix), 0.5)
err_bar2 = np.power(np.diag(cov_matrix2), 0.5)

twopt_write = np.vstack((com_twopt, err_bar)).T
twopt_write2 = np.vstack((com_twopt2, err_bar2)).T

inv_cov_matrix = linalg.inv(cov_matrix)
inv_cov_matrix2 = linalg.inv(cov_matrix2)
inv_cov_matrix_tot = linalg.inv(cov_matrix_tot)

#print np.diag(np.dot(cov_matrix, inv_cov_matrix))
#print np.diag(np.dot(cov_matrix2, inv_cov_matrix2))
#print np.diag(np.dot(cov_matrix_tot, inv_cov_matrix_tot))

cov_mat = cov_matrix[12:16,12:16]
inv_cov_mat = linalg.inv(cov_mat)
cov_mat2 = cov_matrix2[12:16,12:16]
inv_cov_mat2 = linalg.inv(cov_mat2)
cov_mat_tot = cov_matrix_tot[12:16,12:16]
inv_cov_mat_tot = linalg.inv(cov_mat_tot)

err_bar_tot = np.power(np.diag(cov_matrix_tot), 0.5)


twopt_write_tot = np.vstack((com_diff, err_bar_tot)).T

#print np.diag(np.dot(cov_mat, inv_cov_mat))
#print np.diag(np.dot(cov_mat2, inv_cov_mat2))
#print np.diag(np.dot(cov_mat_tot, inv_cov_mat_tot))

cov_two_add = np.add(cov_mat, cov_mat2)
inv_cov_two_add = linalg.inv(cov_two_add)


com_two_cov = com_diff[12:16]

dot = np.array(np.dot(np.dot(com_two_cov,inv_cov_mat_tot),com_two_cov))[0][0]
dot2 = np.array(np.dot(np.dot(com_two_cov,inv_cov_two_add),com_two_cov))[0][0]


p = chisqprob(dot, 4)
p2 = chisqprob(dot2, 4)

s1 = st.norm.ppf(1 - p/2)
s2 = st.norm.ppf(1 - p2/2)

print "New chi^2 = " + str(dot)
print "New p-value = " + str(p)
print "New sigma-value = " + str(s1)

print "Old chi^2 = " + str(dot2)
print "Old p-value = " + str(p2)
print "Old sigma-value = " + str(s2)

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + parameter + "_cov_wider_d.dat", cov_matrix, fmt="%.9E", delimiter = "\t", newline = "\n")

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + parameter + "_cov_d.dat", cov_mat, fmt="%.9E", delimiter = "\t", newline = "\n")

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + parameter + "_2pt_auto_fx_wider_d.dat", twopt_write, fmt = "%.7E", delimiter = "\t")


np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + parameter2 + "_cov_wider_d.dat", cov_matrix2, fmt="%.9E", delimiter = "\t", newline = "\n")

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + parameter2 + "_cov_d.dat", cov_mat2, fmt="%.9E", delimiter = "\t", newline = "\n")

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + parameter2 + "_2pt_auto_fx_wider_d.dat", twopt_write2, fmt = "%.7E", delimiter = "\t")

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + para_tot + "_cov_wider_diff.dat", cov_matrix_tot, fmt="%.9E", delimiter = "\t", newline = "\n")

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + para_tot + "_cov_diff.dat", cov_mat_tot, fmt="%.9E", delimiter = "\t", newline = "\n")

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + para_tot + "_2pt_auto_fx_wider_diff.dat", twopt_write_tot, fmt = "%.7E", delimiter = "\t")





del DD_pair, DR_pair, RR_pair, lens_all_num, ran_all_num, DD_data, DR_data, RR_data, twopt, com_twopt, DD_jackgalcount, DR_jackgalcount, RR_jackgalcount, DD_jackgalpair, DR_jackgalpair, RR_jackgalpair, DD_jack, DR_jack, RR_jack, twopt_jack, DD_jackpair, jack_twopt, twopt_diff, sq_diff, cov_matrix, err_bar, twopt_write