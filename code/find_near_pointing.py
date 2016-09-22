### Find the neighboring 8 pointings.
### Field # (argv[1])
### Updated: May 12, 2016

import os, re, sys, glob
import numpy as np
from astropy.io import fits

def ang_sep(a1,d1,a2,d2):   ### Calculates the angular separation between (a1,d1) and (a2,d2) in celestial coordinates.  Ans in radian
    return np.arccos(np.cos(a1 - a2) * np.cos(d1) * np.cos(d2) + np.sin(d1) * np.sin(d2))

d2r = np.pi / 180

f_man = sys.argv[1]

masklink = "/Volumes/RICKYCHUE/Masks/W{}/".format(f_man)

os.chdir (masklink)

### FITS files
fitsname = np.loadtxt("fitsname.txt", dtype = {'names': ('id',), 'formats': ('S50',)})

fits_num = fitsname.size

### Coordinates of the center of the pointings.
fits_cen = [0] * fits_num

for i in range(fits_num):
    hdu = fits.open(fitsname[i][0])
    prihdr = hdu[0].header
    ### Read in RA and dec
    fits_cen[i] = [prihdr['CRVAL1'] * d2r, prihdr['CRVAL2'] * d2r]
    hdu.close(fitsname[i][0])
    #print fitsname[i][0], fits_cen[i]

### Padding files
padsname = np.loadtxt("padsname.txt", dtype = {'names': ('id',), 'formats': ('S50',)})

pads_num = padsname.size

### Coordinates of the center of the padding pointings.
pads_cen = [0] * pads_num

for i in range(pads_num):
    hdup = fits.open(padsname[i][0])
    prihdrp = hdup[0].header
    ### Read in RA and dec
    pads_cen[i] = [prihdrp['CRVAL1'] * d2r, prihdrp['CRVAL2'] * d2r]
    hdup.close(padsname[i][0])
    #print padsname[i][0], pads_cen[i]

tag = np.array(["m5", "m4", "m3", "m2", "m1", "m0", "p1", "p2", "p3", "p4", "p5"])
tag_num = np.arange(-5,6,1)


### Loop over all the pointings.
for i in range(fits_num):
    consider_pt = [0] * 9   ### 9 pointings should be considered, including itself.
    count = 0
    
    for j in range(fits_num):
        if (fitsname[i][0] == fitsname[j][0]):
            consider_pt[count] = fitsname[j][0]
            count += 1
        else:
            a_sep = ang_sep(fits_cen[i][0], fits_cen[i][1], fits_cen[j][0], fits_cen[j][1])
            if (a_sep <= 0.026179938779914945):     ### Radian for 1.5 deg
                consider_pt[count] = fitsname[j][0]
                count += 1

    for j in range(pads_num):
        a_sep = ang_sep(fits_cen[i][0], fits_cen[i][1], pads_cen[j][0], pads_cen[j][1])
        if (a_sep <= 0.026179938779914945):     ### Radian for 1.5 deg
            consider_pt[count] = padsname[j][0]
            count += 1

    if (count != 9):
        print "Something wrong with counting!"
        exit(1)

    y = [0] * count

    print "Consider: {}".format(fitsname[i][0][:6])
    for j in range(count):
        y[j] = consider_pt[j]

    ### Sort the files according to dec, then RA.
    x = [0] * count

    for j in range(count):
        x[j] = [y[j], tag_num[np.where(y[j][2:4] == tag)[0][0]], tag_num[np.where(y[j][4:6] == tag)[0][0]]]

    b = sorted(x, key = lambda z: (z[2], z[1]), reverse = True)

    c = np.array(b).T[0]

    np.savetxt("{}_fitsname.txt".format(fitsname[i][0][:6]), c, fmt = "%s", delimiter = "\t", newline = "\n")

    del consider_pt, y, x, c, b



