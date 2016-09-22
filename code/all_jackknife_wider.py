### This program calculates the 2pt function with jackknife errors of the combined field.
### Input: parameter name (with/without) (argv[1]), b (argv[2]), clus_num (argv[3]), color (argv[4]), P_th (argv[5])
### Updated: Jan 29, 2016

import numpy as np
import glob, os, re, sys
import matplotlib.pyplot as plt
from numpy import linalg

parameter_name = sys.argv[1]
b = sys.argv[2]
clus_num = sys.argv[3]
color = sys.argv[4]
P_th = sys.argv[5]

parameter = b + "_" + P_th + "_" + clus_num + "_lens" + color

### Number of all lenses, randoms
lens_all_num = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/" + parameter_name + "/" + parameter + "_lensname.txt", usecols = (0,), unpack = True)
ran_all_num = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/random_map/ran_num_iso.dat", usecols = (1,), unpack = True)

ranpairlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/"
piclink = "/Users/rickyccy/Documents/Research_assemblybias/result/Sept_29_2015/"

theta_range = np.loadtxt(ranpairlink + "radial_bin2.dat", unpack = True, usecols = (1,))

theta_power = np.power(theta_range, -0.8)

field_num = 4
jack_num = 54

DD_pair = [0] * field_num
DR_pair = [0] * field_num
RR_pair = [0] * field_num
ang_cor = [0] * field_num
DD_data = [0] * field_num
DR_data = [0] * field_num
RR_data = [0] * field_num
twopt = [0] * field_num
DD_jackgalcount = [0] * field_num
DR_jackgalcount = [0] * field_num
RR_jackgalcount = [0] * field_num
DD_jackgalpair = [0] * field_num
DR_jackgalpair = [0] * field_num
RR_jackgalpair = [0] * field_num
DD_jack = [0] * field_num
DR_jack = [0] * field_num
RR_jack = [0] * field_num
twopt_jack = [0] * field_num

for i in range(field_num):
    DD_pair[i] = lens_all_num[i] * (lens_all_num[i] - 1)
    DR_pair[i] = lens_all_num[i] * ran_all_num[i]
    RR_pair[i] = ran_all_num[i] * (ran_all_num[i] - 1)


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
    RR_data[i] = np.loadtxt("iso_RR_auto_all_wider.dat", unpack = True, usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,)) / RR_pair[i]
    twopt[i] = np.divide((DD_data[i] - 2 * DR_data[i]), RR_data[i]) + 1

    ### Jackknife
    DD_jackgalcount[i] = np.loadtxt(parameter + "_DD_auto_jack_wider.dat", unpack = True, usecols = (1,))
    DR_jackgalcount[i] = np.loadtxt(parameter + "_DR_auto_jack_wider.dat", unpack = True, usecols = (1,))
    RR_jackgalcount[i] = np.loadtxt("iso_RR_auto_jack_wider.dat", unpack = True, usecols = (1,))
   
    ### Make sure the total number of galaxies is right in all jackknife regions.
    DD_jackcheck = len(DD_jackgalcount[i]) * lens_all_num[i] - np.sum(DD_jackgalcount[i])
    DR_jackcheck = len(DR_jackgalcount[i]) * lens_all_num[i] - np.sum(DR_jackgalcount[i])
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
    RR_jackgalpair[i] = np.multiply(RR_jackgalcount[i], (RR_jackgalcount[i] - 1))

    DD_jack[i] = np.divide(np.loadtxt(parameter + "_DD_auto_jack_wider.dat", unpack = True, usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,)), DD_jackgalpair[i])
    DR_jack[i] = np.divide(np.loadtxt(parameter + "_DR_auto_jack_wider.dat", unpack = True, usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,)), DR_jackgalpair[i])
    RR_jack[i] = np.divide(np.loadtxt("iso_RR_auto_jack_wider.dat", unpack = True, usecols = (2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,)), RR_jackgalpair[i])

    twopt_jack[i] = (np.divide((DD_jack[i] - 2 * DR_jack[i]), RR_jack[i]) + 1).T

    del DD_galcount, DR_galcount, RR_galcount, DD_jackcheck, DR_jackcheck, RR_jackcheck

### The combined two pt function.
com_twopt = np.average(twopt, weights = DD_pair, axis = 0)

DD_jackpair = [0] * jack_num
jack_twopt = [0] * jack_num

for i in range(jack_num):
    DD_jackpair[i] = [0] * field_num
    jack_twopt[i] = [0] * field_num

for i in range(18):
    DD_jackpair[i][0] = DD_jackgalpair[0][i]
    jack_twopt[i][0] = twopt_jack[0][i]
    for j in range(1,field_num):
        DD_jackpair[i][j] = DD_pair[j]
        jack_twopt[i][j] = twopt[j]

for i in range(18, 24):
    DD_jackpair[i][0] = DD_pair[0]
    jack_twopt[i][0] = twopt[0]
    for j in range(2,field_num):
        DD_jackpair[i][j] = DD_pair[j]
        jack_twopt[i][j] = twopt[j]
    DD_jackpair[i][1] = DD_jackgalpair[1][i - 18]
    jack_twopt[i][1] = twopt_jack[1][i - 18]

for i in range(24, 48):
    for j in range(2):
        DD_jackpair[i][j] = DD_pair[j]
        jack_twopt[i][j] = twopt[j]
    DD_jackpair[i][3] = DD_pair[3]
    jack_twopt[i][3] = twopt[3]
    DD_jackpair[i][2] = DD_jackgalpair[2][i - 24]
    jack_twopt[i][2] = twopt_jack[2][i - 24]

for i in range(48, 54):
    for j in range(3):
        DD_jackpair[i][j] = DD_pair[j]
        jack_twopt[i][j] = twopt[j]
    DD_jackpair[i][3] = DD_jackgalpair[3][i - 48]
    jack_twopt[i][3] = twopt_jack[3][i - 48]

jackcom_twopt = [0] * jack_num
sq_diff = [0] * jack_num

for i in range(jack_num):
    ### Combined jackknife two-pt function.
    jackcom_twopt[i] = np.average(jack_twopt[i], weights = DD_jackpair[i], axis = 0)

twopt_diff = com_twopt - jackcom_twopt

for i in range(jack_num):
    sq_diff[i] = np.outer(twopt_diff[i], twopt_diff[i])

cov_matrix = np.matrix((np.sum(sq_diff, axis = 0) * (jack_num - 1.) / jack_num))

err_bar = np.power(np.diag(cov_matrix), 0.5)

twopt_write = np.vstack((com_twopt, err_bar)).T

inv_cov_matrix = linalg.inv(cov_matrix)

print np.diag(np.dot(cov_matrix, inv_cov_matrix))

cov_mat = cov_matrix[12:16,12:16]
inv_cov_mat = linalg.inv(cov_mat)

print np.diag(np.dot(cov_mat, inv_cov_mat))

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + parameter + "_cov_wider.dat", cov_matrix, fmt="%.9E", delimiter = "\t", newline = "\n")

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + parameter + "_cov.dat", cov_mat, fmt="%.9E", delimiter = "\t", newline = "\n")

np.savetxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/" + parameter_name + "/" + parameter + "_2pt_auto_fx_wider.dat", twopt_write, fmt = "%.7E", delimiter = "\t")






del DD_pair, DR_pair, RR_pair, lens_all_num, ran_all_num, DD_data, DR_data, RR_data, twopt, com_twopt, DD_jackgalcount, DR_jackgalcount, RR_jackgalcount, DD_jackgalpair, DR_jackgalpair, RR_jackgalpair, DD_jack, DR_jack, RR_jack, twopt_jack, DD_jackpair, jack_twopt, twopt_diff, sq_diff, cov_matrix, err_bar, twopt_write