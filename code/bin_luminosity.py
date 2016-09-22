### This program selects the lens galaxies according to the luminosities.
### Updated: Oct 2, 2015

import glob, os, sys
import numpy as np

red_low = float(sys.argv[1])
red_up = float(sys.argv[2])
red_name = sys.argv[3]
blu_low = float(sys.argv[4])
blu_up = float(sys.argv[5])
blu_name = sys.argv[6]

print red_low, red_up, red_name
print blu_low, blu_up, blu_name

"""
lenslink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_isolated/"
sortlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_with/"
"""

lenslink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens/"
sortlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_onlybright/"
#srcslink = "/Users/rickyccy/Documents/Research_assemblybias/Data/source2/combined2/"

os.chdir (lenslink)

fname = glob.glob("*.tsv")

for i in range(len(fname)):
    print fname[i]
    ### Load in the lens catalog.
    lensfile = np.loadtxt(lenslink + fname[i], unpack = True, dtype = {'names': ('id', 'field', 'x', 'y', 'alpha', 'delta', 'bg', 'lv', 'mu_max', 'mu_thres', 'flag', 'a', 'b', 'theta', 'a_err', 'b_err', 'theta_err', 'class', 'e1', 'e2', 'weight', 'fitclass', 'SNratio', 'mask', 'gal_z', 'gal_z_min', 'gal_z_max', 'T_b', 'modify', 'c2', 'Mu', 'Mg', 'Mr', 'Mi', 'Mz', 'star_flag', 'lp_med', 'lp_inf', 'lp_sup', 'mu', 'mu_err', 'ext_u', 'mi', 'mi_err', 'ext_i', 'mr', 'mr_err', 'ext_r', 'mg', 'mg_err', 'ext_g', 'my', 'my_err', 'ext_y', 'mz', 'mz_err', 'ext_z', 'size', 'FWHM_image', 'FWHM_world', 'Kron_rad', 'Flux_rad', 'bulge_frac', 'model_flux', 'isoarea', 'PDF_1' , 'PDF_2' , 'PDF_3' , 'PDF_4' , 'PDF_5' , 'PDF_6' , 'PDF_7' , 'PDF_8' , 'PDF_9' , 'PDF_10' , 'PDF_11' , 'PDF_12' , 'PDF_13' , 'PDF_14' , 'PDF_15' , 'PDF_16' , 'PDF_17' , 'PDF_18' , 'PDF_19' , 'PDF_20' , 'PDF_21' , 'PDF_22' , 'PDF_23' , 'PDF_24' , 'PDF_25' , 'PDF_26' , 'PDF_27' , 'PDF_28' , 'PDF_29' , 'PDF_30' , 'PDF_31' , 'PDF_32' , 'PDF_33' , 'PDF_34' , 'PDF_35' , 'PDF_36' , 'PDF_37' , 'PDF_38' , 'PDF_39' , 'PDF_40' , 'PDF_41' , 'PDF_42' , 'PDF_43' , 'PDF_44' , 'PDF_45' , 'PDF_46' , 'PDF_47' , 'PDF_48' , 'PDF_49' , 'PDF_50' , 'PDF_51' , 'PDF_52' , 'PDF_53' , 'PDF_54' , 'PDF_55' , 'PDF_56' , 'PDF_57' , 'PDF_58' , 'PDF_59' , 'PDF_60' , 'PDF_61' , 'PDF_62' , 'PDF_63' , 'PDF_64' , 'PDF_65' , 'PDF_66' , 'PDF_67' , 'PDF_68' , 'PDF_69' , 'PDF_70'), 'formats':('S15', 'S6', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64')})
    
    lensfile2 = np.loadtxt(lenslink + fname[i], dtype = {'names': ('id', 'field', 'x', 'y', 'alpha', 'delta', 'bg', 'lv', 'mu_max', 'mu_thres', 'flag', 'a', 'b', 'theta', 'a_err', 'b_err', 'theta_err', 'class', 'e1', 'e2', 'weight', 'fitclass', 'SNratio', 'mask', 'gal_z', 'gal_z_min', 'gal_z_max', 'T_b', 'modify', 'c2', 'Mu', 'Mg', 'Mr', 'Mi', 'Mz', 'star_flag', 'lp_med', 'lp_inf', 'lp_sup', 'mu', 'mu_err', 'ext_u', 'mi', 'mi_err', 'ext_i', 'mr', 'mr_err', 'ext_r', 'mg', 'mg_err', 'ext_g', 'my', 'my_err', 'ext_y', 'mz', 'mz_err', 'ext_z', 'size', 'FWHM_image', 'FWHM_world', 'Kron_rad', 'Flux_rad', 'bulge_frac', 'model_flux', 'isoarea', 'PDF_1' , 'PDF_2' , 'PDF_3' , 'PDF_4' , 'PDF_5' , 'PDF_6' , 'PDF_7' , 'PDF_8' , 'PDF_9' , 'PDF_10' , 'PDF_11' , 'PDF_12' , 'PDF_13' , 'PDF_14' , 'PDF_15' , 'PDF_16' , 'PDF_17' , 'PDF_18' , 'PDF_19' , 'PDF_20' , 'PDF_21' , 'PDF_22' , 'PDF_23' , 'PDF_24' , 'PDF_25' , 'PDF_26' , 'PDF_27' , 'PDF_28' , 'PDF_29' , 'PDF_30' , 'PDF_31' , 'PDF_32' , 'PDF_33' , 'PDF_34' , 'PDF_35' , 'PDF_36' , 'PDF_37' , 'PDF_38' , 'PDF_39' , 'PDF_40' , 'PDF_41' , 'PDF_42' , 'PDF_43' , 'PDF_44' , 'PDF_45' , 'PDF_46' , 'PDF_47' , 'PDF_48' , 'PDF_49' , 'PDF_50' , 'PDF_51' , 'PDF_52' , 'PDF_53' , 'PDF_54' , 'PDF_55' , 'PDF_56' , 'PDF_57' , 'PDF_58' , 'PDF_59' , 'PDF_60' , 'PDF_61' , 'PDF_62' , 'PDF_63' , 'PDF_64' , 'PDF_65' , 'PDF_66' , 'PDF_67' , 'PDF_68' , 'PDF_69' , 'PDF_70'), 'formats':('S15', 'S6', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64')})
    
    T_b = np.array(lensfile[27])
    Mr_new = np.array(lensfile[32])
    mr = np.array(lensfile[45])

    ### Luminosity bands
    """
    L1 = np.where((Mr_new <= -20.0) & (Mr_new > -21.0))
    L2 = np.where((Mr_new <= -21.0) & (Mr_new > -21.5))
    L3 = np.where((Mr_new <= -21.5) & (Mr_new > -22.0))
    L4 = np.where((Mr_new <= -22.0) & (Mr_new > -22.5))
    L5 = np.where((Mr_new <= -22.5) & (Mr_new > -23.0))
    L6 = np.where((Mr_new <= -23.0) & (Mr_new > -23.5))
    L7 = np.where((Mr_new <= -23.5) & (Mr_new > -24.0))
    L8 = np.where((Mr_new <= -24.0) & (Mr_new > -24.5))
    """
    
    L1 = np.where((Mr_new <= red_low) & (Mr_new >= red_up) & (T_b <= 1.5))
    L2 = np.where((Mr_new <= blu_low) & (Mr_new >= blu_up) & (T_b <= 4.0) & (T_b >= 2.0))
    #L1 = np.where((Mr_new <= -19.0) & (Mr_new >= -24.5) & ((T_b <= 1.5) | ((T_b <= 4.0) & (T_b >= 2.0))))
    #L2 = np.where((Mr_new <= -19.0) & (Mr_new >= -24.5) & (T_b <= 4.0) & (T_b >= 2.0))
    
    #source = np.where(((Mr_new > -20.0) & (mr <= 24.7)))
    
    list1 = lensfile2[L1]
    list2 = lensfile2[L2]
    
    """
    list2 = lensfile2[L2]
    list3 = lensfile2[L3]
    list4 = lensfile2[L4]
    list5 = lensfile2[L5]
    list6 = lensfile2[L6]
    list7 = lensfile2[L7]
    list8 = lensfile2[L8]
    """
    #src = lensfile2[source]
    
    
    np.savetxt(sortlink + fname[i][:-4] + 'red_' + red_name + '.tsv', list1, delimiter = "\t", fmt = "%s")
    np.savetxt(sortlink + fname[i][:-4] + 'blu_' + blu_name + '.tsv', list2, delimiter = "\t", fmt = "%s")
    
    """
    np.savetxt(sortlink + fname[i][:-4] + '_L2.tsv', list2, delimiter = "\t", fmt = "%s")
    np.savetxt(sortlink + fname[i][:-4] + '_L3.tsv', list3, delimiter = "\t", fmt = "%s")
    np.savetxt(sortlink + fname[i][:-4] + '_L4.tsv', list4, delimiter = "\t", fmt = "%s")
    np.savetxt(sortlink + fname[i][:-4] + '_L5.tsv', list5, delimiter = "\t", fmt = "%s")
    np.savetxt(sortlink + fname[i][:-4] + '_L6.tsv', list6, delimiter = "\t", fmt = "%s")
    np.savetxt(sortlink + fname[i][:-4] + '_L7.tsv', list7, delimiter = "\t", fmt = "%s")
    np.savetxt(sortlink + fname[i][:-4] + '_L8.tsv', list8, delimiter = "\t", fmt = "%s")
    """
    #np.savetxt(srcslink + fname[i][:2] + 'src_' + fname[i][2] + '2.tsv', src, delimiter = "\t", fmt = "%s")
