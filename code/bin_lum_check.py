### This program selects the lens galaxies according to the luminosities.
### Input: clus_num (argv[1]), red_low (argv[2]), red_up (argv[3]), red_name (argv[4]), blu_low (argv[5]), blu_up (argv[6]), blu_name (argv[7]), A (argv[8]), mag_b (argv[9]), lum class (argv[10])
### Updated: Dec 6, 2015

import glob, os, sys
import numpy as np

clus_num = int(sys.argv[1])
red_low = float(sys.argv[2])
red_up = float(sys.argv[3])
red_name = sys.argv[4]
blu_low = float(sys.argv[5])
blu_up = float(sys.argv[6])
blu_name = sys.argv[7]
P_th = sys.argv[8]
mag_b = sys.argv[9]
lum = sys.argv[10]

print red_low, red_up, red_name
print blu_low, blu_up, blu_name

lenslink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_further/"
sortlink = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_FOF_iso/"
#srcslink = "/Users/rickyccy/Documents/Research_assemblybias/Data/source2/combined2/"

os.chdir (lenslink)

fname = glob.glob(P_th + "_0." + str(clus_num) + "_" + mag_b + "_W*lens.tsv")


for i in range(len(fname)):
    print fname[i]
    ### Load in the lens catalog.
    lensfile = np.loadtxt(lenslink + fname[i], unpack = True, dtype = {'names': ('id', 'field', 'x', 'y', 'alpha', 'delta', 'bg', 'lv', 'mu_max', 'mu_thres', 'flag', 'a', 'b', 'theta', 'a_err', 'b_err', 'theta_err', 'class', 'e1', 'e2', 'weight', 'fitclass', 'SNratio', 'mask', 'gal_z', 'gal_z_min', 'gal_z_max', 'T_b', 'modify', 'c2', 'Mu', 'Mg', 'Mr', 'Mi', 'Mz', 'star_flag', 'lp_med', 'lp_inf', 'lp_sup', 'mu', 'mu_err', 'ext_u', 'mi', 'mi_err', 'ext_i', 'mr', 'mr_err', 'ext_r', 'mg', 'mg_err', 'ext_g', 'my', 'my_err', 'ext_y', 'mz', 'mz_err', 'ext_z', 'size', 'FWHM_image', 'FWHM_world', 'Kron_rad', 'Flux_rad', 'bulge_frac', 'model_flux', 'isoarea', 'PDF_1' , 'PDF_2' , 'PDF_3' , 'PDF_4' , 'PDF_5' , 'PDF_6' , 'PDF_7' , 'PDF_8' , 'PDF_9' , 'PDF_10' , 'PDF_11' , 'PDF_12' , 'PDF_13' , 'PDF_14' , 'PDF_15' , 'PDF_16' , 'PDF_17' , 'PDF_18' , 'PDF_19' , 'PDF_20' , 'PDF_21' , 'PDF_22' , 'PDF_23' , 'PDF_24' , 'PDF_25' , 'PDF_26' , 'PDF_27' , 'PDF_28' , 'PDF_29' , 'PDF_30' , 'PDF_31' , 'PDF_32' , 'PDF_33' , 'PDF_34' , 'PDF_35' , 'PDF_36' , 'PDF_37' , 'PDF_38' , 'PDF_39' , 'PDF_40' , 'PDF_41' , 'PDF_42' , 'PDF_43' , 'PDF_44' , 'PDF_45' , 'PDF_46' , 'PDF_47' , 'PDF_48' , 'PDF_49' , 'PDF_50' , 'PDF_51' , 'PDF_52' , 'PDF_53' , 'PDF_54' , 'PDF_55' , 'PDF_56' , 'PDF_57' , 'PDF_58' , 'PDF_59' , 'PDF_60' , 'PDF_61' , 'PDF_62' , 'PDF_63' , 'PDF_64' , 'PDF_65' , 'PDF_66' , 'PDF_67' , 'PDF_68' , 'PDF_69' , 'PDF_70'), 'formats':('S15', 'S6', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64')})
    
    lensfile2 = np.loadtxt(lenslink + fname[i], dtype = {'names': ('id', 'field', 'x', 'y', 'alpha', 'delta', 'bg', 'lv', 'mu_max', 'mu_thres', 'flag', 'a', 'b', 'theta', 'a_err', 'b_err', 'theta_err', 'class', 'e1', 'e2', 'weight', 'fitclass', 'SNratio', 'mask', 'gal_z', 'gal_z_min', 'gal_z_max', 'T_b', 'modify', 'c2', 'Mu', 'Mg', 'Mr', 'Mi', 'Mz', 'star_flag', 'lp_med', 'lp_inf', 'lp_sup', 'mu', 'mu_err', 'ext_u', 'mi', 'mi_err', 'ext_i', 'mr', 'mr_err', 'ext_r', 'mg', 'mg_err', 'ext_g', 'my', 'my_err', 'ext_y', 'mz', 'mz_err', 'ext_z', 'size', 'FWHM_image', 'FWHM_world', 'Kron_rad', 'Flux_rad', 'bulge_frac', 'model_flux', 'isoarea', 'PDF_1' , 'PDF_2' , 'PDF_3' , 'PDF_4' , 'PDF_5' , 'PDF_6' , 'PDF_7' , 'PDF_8' , 'PDF_9' , 'PDF_10' , 'PDF_11' , 'PDF_12' , 'PDF_13' , 'PDF_14' , 'PDF_15' , 'PDF_16' , 'PDF_17' , 'PDF_18' , 'PDF_19' , 'PDF_20' , 'PDF_21' , 'PDF_22' , 'PDF_23' , 'PDF_24' , 'PDF_25' , 'PDF_26' , 'PDF_27' , 'PDF_28' , 'PDF_29' , 'PDF_30' , 'PDF_31' , 'PDF_32' , 'PDF_33' , 'PDF_34' , 'PDF_35' , 'PDF_36' , 'PDF_37' , 'PDF_38' , 'PDF_39' , 'PDF_40' , 'PDF_41' , 'PDF_42' , 'PDF_43' , 'PDF_44' , 'PDF_45' , 'PDF_46' , 'PDF_47' , 'PDF_48' , 'PDF_49' , 'PDF_50' , 'PDF_51' , 'PDF_52' , 'PDF_53' , 'PDF_54' , 'PDF_55' , 'PDF_56' , 'PDF_57' , 'PDF_58' , 'PDF_59' , 'PDF_60' , 'PDF_61' , 'PDF_62' , 'PDF_63' , 'PDF_64' , 'PDF_65' , 'PDF_66' , 'PDF_67' , 'PDF_68' , 'PDF_69' , 'PDF_70'), 'formats':('S15', 'S6', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'i2', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64' , 'float64')})
    
    T_b = np.array(lensfile[27])
    Mr_new = np.array(lensfile[32])
    smass = np.array(lensfile[36])
    #print np.min(smass), np.max(smass)
    #mr = np.array(lensfile[45])

    
    #L1 = np.where((Mr_new >= red_low) & (Mr_new <= red_up) & (T_b <= 1.5))
    #L2 = np.where((Mr_new >= blu_low) & (Mr_new <= blu_up) & (T_b <= 4.0) & (T_b >= 2.0))
    ### Avoid the problematic blue samples.
    L1 = np.where((smass >= red_low) & (smass <= red_up) & (T_b <= 1.5))
    L2 = np.where((smass >= blu_low) & (smass <= blu_up) & (T_b <= 4.0) & (T_b >= 2.0) & (Mr_new >= -21.0))
    
    list1 = lensfile2[L1]
    list2 = lensfile2[L2]
    
    np.savetxt(sortlink + "W" + str(i + 1) + "_" + str(P_th) + "_0." + str(clus_num) + "_" + mag_b + "_lens" + red_name + '_S' + lum + '.tsv', list1, delimiter = "\t", fmt = "%s")
    np.savetxt(sortlink + "W" + str(i + 1) + "_" + str(P_th) + "_0." + str(clus_num) +  "_" + mag_b +"_lens" + blu_name + '_S' + lum + '.tsv', list2, delimiter = "\t", fmt = "%s")
    
    #np.savetxt(sortlink + fname[i][4:6] + '_' + str(b) + '_' + str(P_th) + '_' + str(clus_num) + '_lens' + red_name + '.tsv', list1, delimiter = "\t", fmt = "%s")
    #np.savetxt(sortlink + fname[i][4:6] + '_' + str(b) + '_' + str(P_th) + '_' + str(clus_num) + '_lens' + blu_name + '.tsv', list2, delimiter = "\t", fmt = "%s")
    
