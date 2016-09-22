### This program checks whether the lens galaxies lie inside the unmasked areas.
### Updated: Jun 8, 2015

import numpy as np
import os, glob
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib as mpl


mask_dir = "/Users/rickyccy/Documents/Research_assemblybias/Masks/W"
mask_map_dir = "/Users/rickyccy/Documents/Research_assemblybias/result/Jul_6_2015/"
lens_dir = "/Users/rickyccy/Documents/Research_assemblybias/random_map/W"
#lens_dir = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_isolated/"


for i in range(3,4):                            ### Different fields.
    print i
    os.chdir(lens_dir + str(i) + "/combined/")
    #os.chdir(lens_dir)
    #lens_file = glob.glob("W" + str(i) + "lens.tsv")   ### Lens files.
    lens_file = glob.glob("ran*")   ### Lens files.
    
    lens_coord = [0] * len(lens_file)
    
    for j in range(len(lens_file)):
        #lens_coord[j] = np.loadtxt(lens_file[j], unpack = True, usecols = (4,5,))
        lens_coord[j] = np.loadtxt(lens_file[j], unpack = True, usecols = (1,2,))

    os.chdir(mask_dir + str(i))

    border_masks_file = "border.reg.txt"

    border_mask = []

    with open(border_masks_file) as myfile2:      ### Open each mask files.
        data2=myfile2.readlines()
        for k in range(len(data2)):        ### len(data)
            x2 = np.array(data2[k][12:-2].split(",")) ### Split the string by commas.
            x2 = x2.astype(float)
            y2 = x2.reshape(x2.size/2, 2)
            border_mask.append(y2)                 ### Masks

    z_bor = np.array([0] * len(border_mask))

    mask_file = glob.glob("C*.txt")
    
    masks = []                                  ### An array to store up the polygon masks.
    
    for j in range(len(mask_file)):                          ### len(mask_file)
        with open(mask_file[j]) as myfile:      ### Open each mask files.
            data=myfile.readlines()
            for k in range(8,len(data)):        ### len(data)
                x = np.array(data[k][12:-2].split(",")) ### Split the string by commas.
                x = x.astype(float)
                y = x.reshape(x.size/2, 2)
                masks.append(y)                 ### Masks

    z = np.array([0] * len(masks))

    """
    os.chdir(mask_dir + str(i) + "/other/")                 ### Load in masks for other colors.
    mask_file_2 = glob.glob("*.txt")
    
    masks_2 = []                                  ### An array to store up the polygon masks.
    
    for j in range(len(mask_file_2)):                          ### len(mask_file)
        with open(mask_file_2[j]) as myfile_2:      ### Open each mask files.
            data_2=myfile_2.readlines()
            for k in range(4,len(data_2)):        ### len(data)
                x_2 = np.array(data_2[k][12:-2].split(",")) ### Split the string by commas.
                x_2 = x_2.astype(float)
                y_2 = x_2.reshape(x_2.size/2, 2)
                masks_2.append(y_2)                 ### Masks

    z_2 = np.array([0] * len(masks_2))
    """
    for j in range(len(lens_file)):
        print j
        fig, ax = plt.subplots()

        # Make the collection and add it to the plot.
        coll = PolyCollection(masks, array=z, cmap=mpl.cm.cool, edgecolor = 'None')
        ax.add_collection(coll)

        coll2 = PolyCollection(border_mask, array=z_bor, cmap=mpl.cm.cool, edgecolor = 'None')
        ax.add_collection(coll2)
        """
            coll_2 = PolyCollection(masks_2, array=z_2, cmap=mpl.cm.cool, edgecolor = 'None', alpha = 0.3)
            ax.add_collection(coll_2)
        """
        ax.invert_xaxis()
        ax.autoscale_view()
        plt.scatter(lens_coord[j][0], lens_coord[j][1], color = 'r', marker = '.', alpha = 0.05, s = 1)
        plt.xlabel('RA, (deg)')
        plt.ylabel('dec, (deg)')

        plt.title("W" + str(i) + "_ran_" + str(j + 1))

        # Add a colorbar for the PolyCollection
        fig.savefig(mask_map_dir + "W" + str(i) + "_ran" + str(j + 1) + ".png", dpi = 3000)
        plt.close()

    del masks, mask_file, lens_file, lens_coord, border_masks_file, border_mask, data, data2#, mask_file_2, masks_2
