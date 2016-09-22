### This program plots the halo coordinates to check their positions.
### Updated: Aug 4, 2015

import numpy as np
import matplotlib.pyplot as plt
import os, glob, re, itertools

flink = "/Users/rickyccy/Documents/Research_assemblybias/mock/"
figlink = "/Users/rickyccy/Documents/Research_assemblybias/result/Aug_4_2015/mock_distribution/"
os.chdir(flink)

infile = open("M19_blu.dat", "r")

x = []
y = []
z = []

for line in itertools.islice(infile, 10, None):
    sline = line.split()
    x.append(float(sline[1]))
    y.append(float(sline[2]))
    z.append(float(sline[3]))

del sline
infile.close()

plt.scatter(x, y, alpha = 0.1)
#plt.xlim((0,100))
#plt.ylim((0,100))
plt.plot()
plt.show()
#plt.savefig(figlink + "M19_xy.png")
plt.close()

del x, y, z
