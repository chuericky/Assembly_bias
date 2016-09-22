### This program combines the pixel count columns and the neighbors into a single file.
### Input: Field # (argv[1])
### Updated: May 22, 2016

import numpy as np
import os, sys, re, glob

f_num = sys.argv[1]

pix_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/W{}/area_count/".format(f_num)
gal_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/"

os.chdir(pix_link)

area_file = glob.glob("*area_count.dat")

lenarea = len(area_file)

area = []

for i in range(lenarea):  #len(area_file)
    area.append(np.loadtxt(area_file[i], dtype = {'names': ('id','p1','p2','p3','p4','p5','n1','n2','n3','n4','n5','f1','f2','f3','f4','f5'), 'formats': ('S15','i20','i20','i20','i20','i20','i20','i20','i20','i20','i20','float64','float64','float64','float64','float64')}))


area_stack = np.hstack((area[0],area[1]))

for i in range(2,lenarea):
    area_stack = np.hstack((area_stack,area[i]))

del area

bc = np.sort(area_stack, order = 'id', kind = 'mergesort')

os.chdir(gal_link)

ingal = np.loadtxt("W{}_over.dat".format(f_num), unpack = True, dtype = {'names': ('id','n1','n2','n3','n4','n5'), 'formats': ('S15','int64','int64','int64','int64','int64')})
ingal2 = np.loadtxt("W{}_over.dat".format(f_num), dtype = {'names': ('id','n1','n2','n3','n4','n5'), 'formats': ('S15','int64','int64','int64','int64','int64')})

lengal = bc.size

gal_name = [0] * lengal

for j in range(lengal):
    gal_name[j] = bc[j][0]

c = np.searchsorted(gal_name, ingal[0])

del ingal, area_stack, gal_name

gal_new = bc[c]

### Make sure the two catalogs have the same sequence.
for i in range(lengal):
    if (ingal2[i][0] != gal_new[i][0]):
        print "Catalogs are not in the same sequence! {}, {}".format(ingal2[i][0],gal_new[i][0])
        exit(1)


gal_mat = [0] * lengal

for i in range(lengal):
    gal_mat[i] = [0] * 25

for i in range(lengal):
    for j in range(1,6):
        gal_mat[i][j - 1] = ingal2[i][j]
    for j in range(1,16):
        gal_mat[i][j + 9] = gal_new[i][j]

#gal_mat = np.matrix(gal_mat)

for i in range(lengal):
    for j in range(5,10):
        gal_mat[i][j] = (gal_mat[i][j - 5] * gal_mat[i][j + 10] + 0.) / gal_mat[i][j + 5]


outfile = open("W{}_overdensity.txt".format(f_num), "w")

for i in range(lengal):
    outfile.write(str(ingal2[i][0]) + "\t")
    for j in range(5):
        outfile.write(str(int(gal_mat[i][j])) + "\t")
    for j in range(5,10):
        outfile.write(str('%.10E' % gal_mat[i][j]) + "\t")
    for j in range(10,20):
        outfile.write(str(int(gal_mat[i][j])) + "\t")
    for j in range(20,25):
        outfile.write(str('%.7E' % gal_mat[i][j]) + "\t")
    outfile.write("\n")

outfile.close()



del bc, ingal2, c, gal_new, gal_mat
