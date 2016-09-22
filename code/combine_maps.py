### This program combines the random maps to a bigger maps.
### Updated : Sep 13, 2015

import os, glob, sys
import numpy as np

field_num = sys.argv[1]
gal_cat = sys.argv[2]

field_link = "/Users/rickyccy/Documents/Research_assemblybias/random_map/W" + field_num
field_link_com = "/Users/rickyccy/Documents/Research_assemblybias/random_map/W" + field_num + "/combined/"

os.chdir(field_link)

fname = glob.glob("*ran_" + gal_cat + ".dat")

com_num = len(fname)

print com_num


big_map = []
for j in range(com_num):   ### com_num
    infile = np.loadtxt(fname[j])
    stri = fname[j][:6]
    
    field_ind = 0
    
    if (stri == "W1p4p3" or stri == "W1p3p3" or stri == "W1p2p3" or stri == "W1p1p3"):
        field_ind = 0
    if (stri == "W1m0p3" or stri == "W1m1p3" or stri == "W1m2p3" or stri == "W1m3p3"):
        field_ind = 1
    if (stri == "W1p4p2" or stri == "W1p3p2" or stri == "W1p2p2" or stri == "W1p1p2"):
        field_ind = 2
    if (stri == "W1m0p2" or stri == "W1m1p2" or stri == "W1m2p2" or stri == "W1m3p2"):
        field_ind = 3
    if (stri == "W1p4p1" or stri == "W1p3p1" or stri == "W1p2p1" or stri == "W1p1p1"):
        field_ind = 4
    if (stri == "W1m0p1" or stri == "W1m1p1" or stri == "W1m2p1" or stri == "W1m3p1"):
        field_ind = 5
    if (stri == "W1p4m0" or stri == "W1p3m0" or stri == "W1p2m0" or stri == "W1p1m0"):
        field_ind = 6
    if (stri == "W1m0m0" or stri == "W1m1m0" or stri == "W1m2m0" or stri == "W1m3m0"):
        field_ind = 7
    if (stri == "W1p4m1" or stri == "W1p3m1" or stri == "W1p2m1" or stri == "W1p1m1"):
        field_ind = 8
    if (stri == "W1m0m1" or stri == "W1m1m1" or stri == "W1m2m1" or stri == "W1m3m1"):
        field_ind = 9
    if (stri == "W1p4m2" or stri == "W1p3m2" or stri == "W1p2m2" or stri == "W1p1m2"):
        field_ind = 10
    if (stri == "W1m0m2" or stri == "W1m1m2" or stri == "W1m2m2" or stri == "W1m3m2"):
        field_ind = 11
    if (stri == "W1p4m3" or stri == "W1p3m3" or stri == "W1p2m3" or stri == "W1p1m3"):
        field_ind = 12
    if (stri == "W1m0m3" or stri == "W1m1m3" or stri == "W1m2m3" or stri == "W1m3m3"):
        field_ind = 13
    if (stri == "W1p4m4" or stri == "W1p3m4" or stri == "W1p2m4" or stri == "W1p1m4"):
        field_ind = 14
    if (stri == "W1m0m4" or stri == "W1m1m4" or stri == "W1m2m4" or stri == "W1m3m4"):
        field_ind = 15
    if (stri == "W1m4p3" or stri == "W1m4p2" or stri == "W1m4p1" or stri == "W1m4m0"):
        field_ind = 16
    if (stri == "W1m4m1" or stri == "W1m4m2" or stri == "W1m4m3" or stri == "W1m4m4"):
        field_ind = 17

    if (stri == "W2p3p2" or stri == "W2p3p1" or stri == "W2p3m0" or stri == "W2p3m1"):
        field_ind = 18
    if (stri == "W2p2p2" or stri == "W2p2p1" or stri == "W2p2m0" or stri == "W2p2m1"):
        field_ind = 19
    if (stri == "W2p1p2" or stri == "W2p1p1" or stri == "W2p1m0" or stri == "W2p1m1"):
        field_ind = 20
    if (stri == "W2m0p2" or stri == "W2m0p1" or stri == "W2m0m0" or stri == "W2m0m1"):
        field_ind = 21
    if (stri == "W2m1p2" or stri == "W2m1p1" or stri == "W2m1m0" or stri == "W2m1m1"):
        field_ind = 22
    if (stri == "W2p3p3" or stri == "W2p2p3" or stri == "W2p1p3" or stri == "W2m0p3" or stri == "W2m1p3"):
        field_ind = 23

    if (stri == "W3p3p3" or stri == "W3p2p3"):
        field_ind = 24
    if (stri == "W3p1p3" or stri == "W3m0p3"):
        field_ind = 25
    if (stri == "W3m1p3" or stri == "W3m2p3"):
        field_ind = 26
    if (stri == "W3p3p2" or stri == "W3p2p2"):
        field_ind = 27
    if (stri == "W3p1p2" or stri == "W3m0p2"):
        field_ind = 28
    if (stri == "W3m1p2" or stri == "W3m2p2"):
        field_ind = 29
    if (stri == "W3p3p1" or stri == "W3p2p1"):
        field_ind = 30
    if (stri == "W3p1p1" or stri == "W3m0p1"):
        field_ind = 31
    if (stri == "W3m1p1" or stri == "W3m2p1"):
        field_ind = 32
    if (stri == "W3p3m0" or stri == "W3p2m0"):
        field_ind = 33
    if (stri == "W3p1m0" or stri == "W3m0m0"):
        field_ind = 34
    if (stri == "W3m1m0" or stri == "W3m2m0"):
        field_ind = 35
    if (stri == "W3p3m1" or stri == "W3p2m1"):
        field_ind = 36
    if (stri == "W3p1m1" or stri == "W3m0m1"):
        field_ind = 37
    if (stri == "W3m1m1" or stri == "W3m2m1"):
        field_ind = 38
    if (stri == "W3p3m2" or stri == "W3p2m2"):
        field_ind = 39
    if (stri == "W3p1m2" or stri == "W3m0m2"):
        field_ind = 40
    if (stri == "W3m1m2" or stri == "W3m2m2"):
        field_ind = 41
    if (stri == "W3p3m3" or stri == "W3p2m3"):
        field_ind = 42
    if (stri == "W3p1m3" or stri == "W3m0m3"):
        field_ind = 43
    if (stri == "W3m1m3" or stri == "W3m2m3" or stri == "W3m3m3"):
        field_ind = 44
    if (stri == "W3m3m2" or stri == "W3m3m1"):
        field_ind = 45
    if (stri == "W3m3m0" or stri == "W3m3p1"):
        field_ind = 46
    if (stri == "W3m3p2" or stri == "W3m3p3"):
        field_ind = 47

    if (stri == "W4p2m2" or stri == "W4p1m2" or stri == "W4m0m2" or stri == "W4m1m2"):
        field_ind = 48
    if (stri == "W4p2m1" or stri == "W4p1m1" or stri == "W4m0m1" or stri == "W4m1m1"):
        field_ind = 49
    if (stri == "W4p2m0" or stri == "W4p1m0" or stri == "W4m0m0" or stri == "W4m1m0"):
        field_ind = 50
    if (stri == "W4p1p1" or stri == "W4m0p1" or stri == "W4m1p1" or stri == "W4m1p2" or stri == "W4m1p3"):
        field_ind = 51
    if (stri == "W4m2m0" or stri == "W4m2p1" or stri == "W4m2p2" or stri == "W4m2p3"):
        field_ind = 52
    if (stri == "W4m3m0" or stri == "W4m3p1" or stri == "W4m3p2" or stri == "W4m3p3"):
        field_ind = 53

    field_arr = [field_ind] * len(infile)

    print stri, len(infile), field_ind

    infile = np.vstack((field_arr, infile.T)).T

    if (j == 0):
        big_map = infile
    else:
        big_map = np.append(big_map, infile, axis = 0)

    del field_arr

np.savetxt(field_link_com + "ran_" + gal_cat + ".dat", big_map, fmt='%d\t%.6f\t%.6f')

del big_map
