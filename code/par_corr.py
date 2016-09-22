### It computes the correlations among parameters.  Calculates the predicted concentration difference and mass difference by Gaussian statistics, then calculates the density of the environments splitting by concentration quartiles, and then further separate the galaxies within the same mass bin w.r.t. concentrations, using Gaussian statistics.
### Matrix: a_world, b_world, a/b (elongation), 1 - b/a (ellpticity), Mu, Mg, Mr, Mi, Mz, flux_rad, isoarea_world, overdensity rad 1 - 5, shear 1 - 14.
### Updated: Aug 19, 2016

import numpy as np
import os, sys, re, glob
from scipy.stats import norm
import matplotlib.pyplot as plt
from numpy import linalg
from scipy.optimize import curve_fit
from scipy.stats import norm

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

### Weighted average, variance and skewness.
def w_avg_var_skew(values, weights):
    average = np.average(values, weights=weights)
    w_sum = np.sum(weights)
    biased_var = np.average((values-average)**2, weights=weights)
    variance = biased_var * (w_sum / (w_sum - 1.))
    biased_skew = np.average((values-average)**3, weights=weights) / biased_var**1.5
    skew = biased_skew * (w_sum / (w_sum - 1.))**0.5
    return (average, variance**0.5, skew)

gal_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/"
mass_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/mass/"
conc_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/concentration/"
shear_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/shear/"
mock_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/mock/"
shear_dist_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/"
result_link = "/Users/rickyccy/Documents/Research_assemblybias/result/Jul_27_2016/"
result_shear_link = "/Users/rickyccy/Documents/Research_assemblybias/result/Jul_27_2016/shear_compare/"
quartile_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/quartiles/"
quarvec_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/quartile_resvec/"

os.chdir(gal_link)

### Dimensionless Hubble constant
h = 0.7

file_num = 4

lens_delta = [0] * file_num
lens_shear = [0] * file_num
lens_par = [0] * file_num
lens_name = [0] * file_num

for i in range(file_num):
    lens_delta[i] = np.loadtxt("W{}_overdensity.dat".format(i+1), usecols = (6,7,8,9,10))
    lens_shear[i] = np.loadtxt(shear_link + "W{}shear.dat".format(i+1), usecols = (1,2,3,4,5,6,7,8,9,10,11,12,13,14))
    ### a_world, b_world, a/b, 1 - b/a, flux_rad, isoarea_world, other than magnitude, everything is in log scale.
    lens_par[i] = np.loadtxt("W{}lens.tsv".format(i+1), usecols = (11,12,61,64,30,31,32,33,34,36))
    lens_name[i] = np.loadtxt("W{}_overdensity.dat".format(i+1), unpack = True, usecols = (0,11), dtype = {'names':('id','x'), 'formats':('S20','float64')})[0]

delta_mat = [0]
shear_mat = [0]
par_mat = [0]
name_mat = [0]


### Combining overdensities, shear and parameter matrices.
delta_mat = np.vstack((lens_delta[0], lens_delta[1]))
shear_mat = np.vstack((lens_shear[0], lens_shear[1]))
par_mat = np.vstack((lens_par[0], lens_par[1]))
name_mat = np.hstack((lens_name[0], lens_name[1]))

for i in range(2,file_num):
    delta_mat = np.vstack((delta_mat, lens_delta[i]))
    shear_mat = np.vstack((shear_mat, lens_shear[i]))
    par_mat = np.vstack((par_mat, lens_par[i]))
    name_mat = np.hstack((name_mat, lens_name[i]))

delta_mat = np.log10(delta_mat).T
shear_mat = shear_mat.T
par_mat = par_mat.T
name_mat = name_mat.T


"""
a = par_mat[0]
b = par_mat[1]

## Ellipticity (a^2 - b^2)/(a^2 + b^2)
a_sq = np.power(a,2)
b_sq = np.power(b,2)
ellp2 = np.divide((a_sq - b_sq),(a_sq + b_sq))

## Elongation a/b
elong = np.divide(a,b)

## Ellipticity 1 - b/a
ellp = 1 - np.divide(b,a)
"""
#par_mat = np.insert(par_mat, 2, elong, axis = 0)
#par_mat = np.insert(par_mat, 3, ellp, axis = 0)

del lens_delta, lens_shear, lens_par, lens_name #, a, b, a_sq, b_sq, elong, ellp

### Put a_world, b_world, flux_radius and isoarea_world into log scale.
for i in range(4):
    par_mat[i] = np.log10(par_mat[i])

mat = np.vstack((par_mat, delta_mat))
mat = np.vstack((mat, shear_mat))

false_pos = np.where((np.isnan(mat[15]) == False))[0]

### Some shear profiles have NAN, avoid them.
gdmat = mat.T[false_pos].T

gdname = name_mat.T[false_pos].T

"""
for i in range(5):
    step = 0.02
    res_bin = np.arange(2.5,5.0,step)
    res_pltbin = np.arange(2.51,4.989,step)

    hist_res, bin_edges_res = np.histogram(gdshear[i], bins = res_bin)
    hist_res = hist_res / (np.sum(hist_res) + 0.)
    plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
    #plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
    plt.ylabel("Probability")
    #plt.title(r"$\log(\bar{M_h}/h_{70}^{-1}M_\odot) = $" + str('%.5g' % mass_avg))
    #plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
    plt.xlabel(r"$\delta(r)$")
    plt.xlim((2.5,5.0))
    plt.show()
    plt.close()

"""


del mat, delta_mat, shear_mat, par_mat, name_mat


### Covariance matrices, for all the galaxies.
covmat = np.cov(gdmat)
corrmat = np.corrcoef(gdmat)

gd_std = np.power(np.diag(covmat),0.5)

### Mean of the parameters, for all the galaxies.
gd_mean = np.mean(gdmat,axis = 1)

gd_meanstd = np.vstack((gd_mean,gd_std)).T

### Mean of the parameters, overdensity and shear
par_mean = np.matrix(gd_mean[:10])
den_mean = np.matrix(gd_mean[10:15])
shear_mean = np.matrix(gd_mean[-14:-5])        ### If all shear bins, use -14:, first 9 bins.

np.savetxt(gal_link + "par_correlation_matrix.dat", corrmat)
np.savetxt(gal_link + "par_mean.dat", gd_meanstd)

### Compute the predicted shear from the correlation functions.
### gamma_p (vec p) = <gamma p> * <pp>^-1 vec p

### Covariance and inverse covariance matrices of the galaxy parameters, <pp> and <pp>^-1.
cov_par = covmat[:10,:10]
inv_cov_par = linalg.inv(cov_par)

print np.dot(cov_par, inv_cov_par)
print "\n\n"

### Cross-Covariance matrix of the shear with galaxy parameters, and the matrix of overdensity with galaxy parameters, for all the galaxies.
cov_gamma_par = covmat[-14:-5,:10]     ### If all shear bins, use -14:
cov_den_par = covmat[10:15,:10]         ### <np> correlation between number density and parameter
cov_gamma_par_tsp = cov_gamma_par.T




### Linear operator, <gamma * p> * <p * p>^-1, for all the galaxies, and < n * p> * <p * p>^-1;
f_gpar = np.dot(cov_gamma_par, inv_cov_par)
f_dpar = np.dot(cov_den_par, inv_cov_par)



### Operator of the parameters, delta p of CFHTLenS galaxies, all the galaxies.
par_opr = gdmat[:10,:] - par_mean.T

### Standard deviations of the parameters
par_std = gd_std[:10]

### Standardized residuals.
std_opr = np.divide(par_opr.T, par_std).T

"""
aind = np.arange(0,par_opr[0].size,1)

for i in range(10):
    nor_score = np.divide(np.array(par_opr[i])[0], par_std[i])
    plt.scatter(aind, nor_score, alpha = 0.3)
    plt.show()
    plt.close()
"""



### Delta n = <n * p> * <p * p>^(-1) * dp
delta_numden = np.dot(f_dpar, par_opr)
delta_gamma = np.dot(f_gpar, par_opr)



shear_mean = np.array(shear_mean)

shear_cpr = gdmat[-14:-5,:] - shear_mean.T


### Radial bins in CFHTLenS.
rr = obs_r_bin = np.loadtxt(shear_dist_link + "dist_bin2.dat", unpack = True, usecols = (1,))[:9]



### Covariance matrix of shears, <gamma * gamma>, and overdensities, <logn * logn>
cov_gamma = covmat[-14:-5,-14:-5]
inv_cov_gamma = linalg.inv(cov_gamma)
#cov_den = covmat[10:15,10:15]
#inv_cov_den = linalg.inv(cov_den)

print np.dot(cov_gamma, inv_cov_gamma)
print "\n\n"

#print np.dot(cov_den,inv_cov_den)
#print "\n\n"

### Kernel matrix, <gamma * gamma>^-1 * <gamma * p> * <p * p>^-1, <logn * logn>^-1 * <logn * p> * <p * p>^-1
gau_ker = np.dot(inv_cov_gamma, f_gpar)
#den_gau_ker = np.dot(inv_cov_den, f_dpar)

### Matrix for Residual of the parameter, 10 columns
np.savetxt(gal_link + "res_par.dat", par_opr.T, fmt = "%.10E")

### <gamma * gamma>^-1 * <gamma * p> * <p * p>^-1, <logn * p> * <p * p>^-1
np.savetxt(gal_link + "shear_gaussian_res_kernel.dat", gau_ker, delimiter = "\t", fmt = "%.10E")
np.savetxt(gal_link + "overdensity_gaussian_res_kernel.dat", f_dpar, delimiter = "\t", fmt = "%.10E")

"""
aaaaa = np.power(np.diag(cov_par), 0.5)

print np.multiply(f_dpar, aaaaa)

#print cov_par

exit()
"""

den_num = 5

den_name = ["3Mpc", "5Mpc", "7Mpc", "9Mpc", "11Mpc"]
par_bin = 4

delta_numden = np.array(delta_numden)

for m in range(1,1):
    
    print m
    
    ### 25%, 50% and 75% quartiles
    quartile1 = np.percentile(delta_numden[m],25)
    quartile2 = np.percentile(delta_numden[m],50)
    quartile3 = np.percentile(delta_numden[m],75)
    
    quar_ind = [0] * par_bin
    
    ### The indices of the galaxies which belong to the corresponding density quartiles.
    quar_ind[0] = np.where((delta_numden[m] <= quartile1))[0]
    quar_ind[1] = np.where(((delta_numden[m] > quartile1) & (delta_numden[m] <= quartile2)))[0]
    quar_ind[2] = np.where((delta_numden[m] > quartile2) & (delta_numden[m] <= quartile3))[0]
    quar_ind[3] = np.where((delta_numden[m] > quartile3))[0]
    
    quar = [0] * par_bin
    quar_name = [0] * par_bin
    
    for i in range(par_bin):
        quar[i] = delta_numden[m][quar_ind[i]]
        quar_name[i] = gdname[quar_ind[i]]

    
    quar_num = 4
    
    quar_tit = [0] * quar_num
    for k in range(quar_num):
        quar_tit[k] = [0] * file_num
        for j in range(file_num):
            aaa = np.core.defchararray.find(quar_name[k], "W{}".format(j+1))
            quar_tit[k][j] = quar_name[k][np.where((aaa == 0))[0]]
            #print m,k,j,quar_tit[k][j].size


    lens_numden = [0] * file_num
    lens_pos = [0] * file_num
    lens_gamma = [0] * file_num
    
    for k in range(file_num):
        lens_numden[k] = np.loadtxt("W{}_overdensity.dat".format(k+1), usecols = (0,6,7,8,9,10),dtype = {'names':('id','n1','n2','n3','n4','n5'), 'formats':('S20','float64','float64','float64','float64','float64')})
        lens_pos[k] = np.loadtxt("W{}lens.tsv".format(k+1), usecols = (0,1,4,5,24,36),dtype = {'names':('id','field','ra','dec','z','sm'), 'formats':('S20','S20','float64','float64','float64','float64')})
        lens_gamma[k] = np.loadtxt(shear_link + "W{}shear.dat".format(k+1), dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
        
        bc = np.sort(lens_numden[k], order = 'id', kind = 'mergesort')
        cd = np.sort(lens_pos[k], order = 'id', kind = 'mergesort')
        ce = np.sort(lens_gamma[k], order = 'id', kind = 'mergesort')
        
        g_name = [0] * lens_numden[k].size
        
        for j in range(lens_numden[k].size):
            g_name[j] = bc[j][0]
    
    
        for j in range(quar_num):
            c = np.searchsorted(g_name, quar_tit[j][k])
            gal_new = bc[c]
            gal_pos = cd[c]
            gal_gamma = ce[c]
            
            np.savetxt(quartile_link + "W{}_den{}bin{}_overdensity.dat".format(k+1,den_name[m], j+1), gal_new, fmt = "%s\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E", delimiter = "\t", newline = "\n")
            np.savetxt(quartile_link + "W{}_den{}bin{}_pos_temp.dat".format(k+1,den_name[m],j+1), gal_pos, fmt = "%s\t%s\t%.9g\t%.9g\t%.9g\t%.9g", delimiter = "\t", newline = "\n")
            np.savetxt(quartile_link + "W{}_den{}bin{}_shear.dat".format(k+1,den_name[m],j+1), gal_gamma, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")
        
        del bc, cd, ce, g_name, gal_new, gal_pos, gal_gamma

    ### Gaussian residuals of the parameters for all CFHTLenS galaxies.
    #np.savetxt(mass_link + "shear_den{}_res.dat".format(den_name[m]), par_opr[m].T, fmt = "%.10E")
    

    del quar, quar_tit, quar_ind, quar_name, lens_numden, lens_pos, lens_gamma, quartile1, quartile2, quartile3#, quar_stdscore







### CFHTLenS distance bin.
obs_r_bin = np.loadtxt(shear_dist_link + "dist_bin2.dat", unpack = True, usecols = (1,))[:9]

mock_r_bin = np.arange(0.01,0.36,0.01) / h

### Mass derivatives and mass variance, convert the hubble factor back.
g_mass_derv_load = np.loadtxt(mock_link + "shear_mass_derivative_new.dat", unpack = True, usecols = (1,))
mass_var = np.loadtxt(mock_link + "shear_mass_derivative_new.dat", unpack = True, usecols = (3,))[0]
mass_avg = np.loadtxt(mock_link + "shear_mass_derivative_new.dat", unpack = True, usecols = (0,))[0]
#g_mass_den_derv_load = np.loadtxt(mock_link + "overden_mass_derivative.dat", unpack = True, usecols = (1,))

### Interpolate the derivatives times the variance of mass to the r bins for the observational data.
g_mass_derv = np.interp(obs_r_bin, mock_r_bin, g_mass_derv_load) * mass_var

#g_mass_den_derv = g_mass_den_derv_load * mass_var

### Mass residual kernel.
mass_ker = np.dot(g_mass_derv, gau_ker)
#mass_den_ker = np.dot(g_mass_den_derv, den_gau_ker)

### d(gamma)/d(logM) * <(logM)^2> * <gamma * gamma>^-1 * <gamma * p> * <p * p>^-1; d(logn)/d(logM) * <(logM)^2> * <logn * logn>^-1 * <logn * p> * <p * p>^-1
np.savetxt(mass_link + "shear_mass_res_kernel.dat", mass_ker, fmt = "%.10E")
#np.savetxt(mass_link + "overden_mass_res_kernel.dat", mass_den_ker, fmt = "%.10E")

### Number of parameters
par_num = cov_gamma_par[0].size

par_ker = [0] * par_num

for i in range(par_num):
    par_ker[i] = np.dot(cov_gamma_par_tsp[i], gau_ker)


par_name = ["a_world", "b_world", "flux_rad", "iso_area", "Mag_u", "Mag_g", "Mag_r", "Mag_i", "Mag_z", "starmass"]

### <p_i * gamma> * <gamma * gamma>^-1 * <gamma * p> * <p * p>^-1, where p_i is any parameter.
for i in range(par_num):
    np.savetxt(mass_link + "shear_{}_res_kernel.dat".format(par_name[i]), par_ker[i], fmt = "%.10E")


#########################################

temp_vec = np.dot(cov_par, mass_ker)
temp_vec_denom = np.dot(mass_ker, temp_vec)

prep_vec = [0] * den_num

par_bin = 4


for m in range(den_num):  ### den_num
    
    print m
    
    ### Calculate the orthogonal vector to the vector pointing to the direction for largest increase in M_halo.
    temp_vec_numt = np.dot(f_dpar[m], temp_vec)
    temp_fac = temp_vec_numt / temp_vec_denom
    
    prep_vec[m] = f_dpar[m] - np.multiply(temp_fac, mass_ker)
    
    np.savetxt(mass_link + "prep_vec_{}.dat".format(den_name[m]), prep_vec[m], fmt = "%.10E")
    
    ### residuals of the prep vector
    prep_res = np.array(np.dot(prep_vec[m], par_opr))[0]

    ### 25%, 50% and 75% quartiles
    quartile1 = np.percentile(prep_res,25)
    quartile2 = np.percentile(prep_res,50)
    quartile3 = np.percentile(prep_res,75)
    
    quar_ind = [0] * par_bin
    
    ### The indices of the galaxies which belong to the corresponding mass quartiles.
    quar_ind[0] = np.where((prep_res <= quartile1))[0]
    quar_ind[1] = np.where(((prep_res > quartile1) & (prep_res <= quartile2)))[0]
    quar_ind[2] = np.where((prep_res > quartile2) & (prep_res <= quartile3))[0]
    quar_ind[3] = np.where((prep_res > quartile3))[0]
    
    quar = [0] * par_bin
    quar_name = [0] * par_bin
    
    
    for i in range(par_bin):
        #quar[i] = par_res[m][quar_ind[i]]
        quar[i] = prep_res[quar_ind[i]]
        quar_name[i] = gdname[quar_ind[i]]
    
    quar_num = 4
    
    quar_tit = [0] * quar_num
    for k in range(quar_num):
        quar_tit[k] = [0] * file_num
        for j in range(file_num):
            aaa = np.core.defchararray.find(quar_name[k], "W{}".format(j+1))
            quar_tit[k][j] = quar_name[k][np.where((aaa == 0))[0]]
            #print m,k,j,quar_tit[k][j].size

    lens_numden = [0] * file_num
    lens_pos = [0] * file_num
    lens_gamma = [0] * file_num
    
    for k in range(file_num):
        lens_numden[k] = np.loadtxt("W{}_overdensity.dat".format(k+1), usecols = (0,6,7,8,9,10),dtype = {'names':('id','n1','n2','n3','n4','n5'), 'formats':('S20','float64','float64','float64','float64','float64')})
        lens_pos[k] = np.loadtxt("W{}lens.tsv".format(k+1), usecols = (0,1,4,5,24,36),dtype = {'names':('id','field','ra','dec','z','sm'), 'formats':('S20','S20','float64','float64','float64','float64')})
        lens_gamma[k] = np.loadtxt(shear_link + "W{}shear.dat".format(k+1), dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
        
        bc = np.sort(lens_numden[k], order = 'id', kind = 'mergesort')
        cd = np.sort(lens_pos[k], order = 'id', kind = 'mergesort')
        ce = np.sort(lens_gamma[k], order = 'id', kind = 'mergesort')
        
        g_name = [0] * lens_numden[k].size
        
        for j in range(lens_numden[k].size):
            g_name[j] = bc[j][0]
    
    
        for j in range(quar_num):
            c = np.searchsorted(g_name, quar_tit[j][k])
            gal_new = bc[c]
            gal_pos = cd[c]
            gal_gamma = ce[c]
            
            np.savetxt(quartile_link + "W{}_{}prepdir_bin{}_overdensity.dat".format(k+1,den_name[m], j+1), gal_new, fmt = "%s\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E", delimiter = "\t", newline = "\n")
            np.savetxt(quartile_link + "W{}_{}prepdir_bin{}_pos_temp.dat".format(k+1,den_name[m],j+1), gal_pos, fmt = "%s\t%s\t%.9g\t%.9g\t%.9g\t%.9g", delimiter = "\t", newline = "\n")
            np.savetxt(quartile_link + "W{}_{}prepdir_bin{}_shear.dat".format(k+1,den_name[m],j+1), gal_gamma, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")
        
        del bc, cd, ce, g_name, gal_new, gal_pos, gal_gamma

    ### Gaussian residuals of the parameters for all CFHTLenS galaxies.
    np.savetxt(mass_link + "shear_{}prepdir_res.dat".format(den_name[m]), prep_res.T, fmt = "%.10E")
    

    del quar, quar_tit, quar_ind, quar_name, lens_numden, lens_pos, lens_gamma, quartile1, quartile2, quartile3#, quar_stdscore










#########################################
"""
### Check
mass_res_check = np.dot(g_mass_derv, np.dot(inv_cov_gamma, shear_cpr))

print np.min(mass_res_check), np.max(mass_res_check)

step = 5e-3
res_bin = np.arange(-0.2,0.2,step)
res_pltbin = np.arange(-0.1975,0.1975,step)

hist_res, bin_edges_res = np.histogram(mass_res_check, bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
#plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Shear estimate, " + r"$\log(\bar{M_h}/h_{70}^{-1}M_\odot) = $" + str('%.5g' % mass_avg))
#plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\log M_h$")
plt.xlim((-0.13,0.13))
plt.savefig(result_link + "actualshear_mass_res.png")
plt.close()
######################
"""
"""
smass = gdmat[9,:]

print smass

### 25%, 50% and 75% quartiles
quartile1 = np.percentile(smass,25)
quartile2 = np.percentile(smass,50)
quartile3 = np.percentile(smass,75)

quar1_ind = np.where((smass <= quartile1))[0]
quar2_ind = np.where(((smass > quartile1) & (smass <= quartile2)))[0]
quar3_ind = np.where(((smass > quartile2) & (smass <= quartile3)))[0]
quar4_ind = np.where((smass > quartile3))[0]

quar = [0] * 4
quar_name = [0] * 4

quar[0] = smass[quar1_ind]
quar[1] = smass[quar2_ind]
quar[2] = smass[quar3_ind]
quar[3] = smass[quar4_ind]

quar_name[0] = gdname[quar1_ind]
quar_name[1] = gdname[quar2_ind]
quar_name[2] = gdname[quar3_ind]
quar_name[3] = gdname[quar4_ind]


quar_num = 4

quar_tit = [0] * quar_num
for k in range(quar_num):
    quar_tit[k] = [0] * file_num
    for j in range(file_num):
        aaa = np.core.defchararray.find(quar_name[k], "W{}".format(j+1))
        quar_tit[k][j] = quar_name[k][np.where((aaa == 0))[0]]
        print k,j,quar_tit[k][j].size

lens_numden = [0] * file_num
lens_pos = [0] * file_num
lens_gamma = [0] * file_num

for k in range(file_num):
    lens_numden[k] = np.loadtxt("W{}_overdensity.dat".format(k+1), usecols = (0,6,7,8,9,10),dtype = {'names':('id','n1','n2','n3','n4','n5'), 'formats':('S20','float64','float64','float64','float64','float64')})
    lens_pos[k] = np.loadtxt("W{}lens.tsv".format(k+1), usecols = (0,1,4,5,24,36),dtype = {'names':('id','field','ra','dec','z','sm'), 'formats':('S20','S20','float64','float64','float64','float64')})
    lens_gamma[k] = np.loadtxt(shear_link + "W{}shear.dat".format(k+1), dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
    
    bc = np.sort(lens_numden[k], order = 'id', kind = 'mergesort')
    cd = np.sort(lens_pos[k], order = 'id', kind = 'mergesort')
    ce = np.sort(lens_gamma[k], order = 'id', kind = 'mergesort')
    
    g_name = [0] * lens_numden[k].size
    
    for j in range(lens_numden[k].size):
        g_name[j] = bc[j][0]
    
    
    for j in range(quar_num):
        c = np.searchsorted(g_name, quar_tit[j][k])
        gal_new = bc[c]
        gal_pos = cd[c]
        gal_gamma = ce[c]
        np.savetxt(quartile_link + "W{}_smass_quartile{}_overdensity.dat".format(k+1,j+1), gal_new, fmt = "%s\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E", delimiter = "\t", newline = "\n")
        np.savetxt(quartile_link + "W{}_smass_quartile{}_pos_temp.dat".format(k+1,j+1), gal_pos, fmt = "%s\t%s\t%.9g\t%.9g\t%.9g\t%.9g", delimiter = "\t", newline = "\n")
        np.savetxt(quartile_link + "W{}_smass_quartile{}_shear.dat".format(k+1,j+1), gal_gamma, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")
    
    del bc, cd, ce, g_name, gal_new, gal_pos, gal_gamma

del quar, quar_tit

exit(1)
"""


par_res = [0] * par_num

### d (parameter)
for i in range(par_num):
    par_res[i] = np.array(np.dot(par_ker[i], par_opr))[0]
    #print i, par_res[i], par_res[i].size, np.min(par_res[i]), np.max(par_res[i])



### dlogM
mass_res = np.array(np.dot(mass_ker, par_opr))[0]
#mass_den_res = np.dot(mass_den_ker, par_opr)

print np.min(mass_res), np.max(mass_res)
#print np.min(mass_den_res), np.max(mass_den_res)


par_opr = np.array(par_opr)



par_bin = 4

for m in range(par_num):
    
    print m, par_name[m]
    
    """
    ### 25%, 50% and 75% quartiles
    quartile1 = np.percentile(par_res[m],25)
    quartile2 = np.percentile(par_res[m],50)
    quartile3 = np.percentile(par_res[m],75)
    
    quar_ind = [0] * par_bin

    ### The indices of the galaxies which belong to the corresponding mass quartiles.
    quar_ind[0] = np.where((par_res[m] <= quartile1))[0]
    quar_ind[1] = np.where(((par_res[m] > quartile1) & (par_res[m] <= quartile2)))[0]
    quar_ind[2] = np.where((par_res[m] > quartile2) & (par_res[m] <= quartile3))[0]
    quar_ind[3] = np.where((par_res[m] > quartile3))[0]
    """
    ### 25%, 50% and 75% quartiles
    quartile1 = np.percentile(par_opr[m],25)
    quartile2 = np.percentile(par_opr[m],50)
    quartile3 = np.percentile(par_opr[m],75)
    
    quar_ind = [0] * par_bin
    
    ### The indices of the galaxies which belong to the corresponding mass quartiles.
    quar_ind[0] = np.where((par_opr[m] <= quartile1))[0]
    quar_ind[1] = np.where(((par_opr[m] > quartile1) & (par_opr[m] <= quartile2)))[0]
    quar_ind[2] = np.where((par_opr[m] > quartile2) & (par_opr[m] <= quartile3))[0]
    quar_ind[3] = np.where((par_opr[m] > quartile3))[0]
    
    """
    ### Standardized scores of the parameters.
    quar_stdscore = [0] * par_bin
    
    for p in range(par_bin):
        quar_stdscore[p] = std_opr.T[quar_ind[p]]
        quar_stdmean = np.mean(quar_stdscore[p], axis = 0)
        quar_stdsdv = np.power(np.diag(np.cov(quar_stdscore[p].T)), 0.5)
        
        quar_stat = np.vstack((quar_stdmean, quar_stdsdv)).T

        np.savetxt(quarvec_link + "{}bin{}_res_par_standardized.dat".format(par_name[m], p+1), quar_stdscore[p], fmt='%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E', delimiter = "\t", newline = "\n")
        np.savetxt(quarvec_link + "{}bin{}_res_par_statistics.dat".format(par_name[m], p+1), quar_stat, fmt='%.10E\t%.10E', newline = "\n")
    """
    
    quar = [0] * par_bin
    quar_name = [0] * par_bin
    

    for i in range(par_bin):
        #quar[i] = par_res[m][quar_ind[i]]
        quar[i] = par_opr[m][quar_ind[i]]
        quar_name[i] = gdname[quar_ind[i]]

    quar_num = 4

    quar_tit = [0] * quar_num
    for k in range(quar_num):
        quar_tit[k] = [0] * file_num
        for j in range(file_num):
            aaa = np.core.defchararray.find(quar_name[k], "W{}".format(j+1))
            quar_tit[k][j] = quar_name[k][np.where((aaa == 0))[0]]
            #print m,k,j,quar_tit[k][j].size


    lens_numden = [0] * file_num
    lens_pos = [0] * file_num
    lens_gamma = [0] * file_num

    for k in range(file_num):
        lens_numden[k] = np.loadtxt("W{}_overdensity.dat".format(k+1), usecols = (0,6,7,8,9,10),dtype = {'names':('id','n1','n2','n3','n4','n5'), 'formats':('S20','float64','float64','float64','float64','float64')})
        lens_pos[k] = np.loadtxt("W{}lens.tsv".format(k+1), usecols = (0,1,4,5,24,36),dtype = {'names':('id','field','ra','dec','z','sm'), 'formats':('S20','S20','float64','float64','float64','float64')})
        lens_gamma[k] = np.loadtxt(shear_link + "W{}shear.dat".format(k+1), dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
    
        bc = np.sort(lens_numden[k], order = 'id', kind = 'mergesort')
        cd = np.sort(lens_pos[k], order = 'id', kind = 'mergesort')
        ce = np.sort(lens_gamma[k], order = 'id', kind = 'mergesort')
    
        g_name = [0] * lens_numden[k].size
    
        for j in range(lens_numden[k].size):
            g_name[j] = bc[j][0]


        for j in range(quar_num):
            c = np.searchsorted(g_name, quar_tit[j][k])
            gal_new = bc[c]
            gal_pos = cd[c]
            gal_gamma = ce[c]
        
            np.savetxt(quartile_link + "W{}_{}bin{}_overdensity.dat".format(k+1,par_name[m], j+1), gal_new, fmt = "%s\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E", delimiter = "\t", newline = "\n")
            np.savetxt(quartile_link + "W{}_{}bin{}_pos_temp.dat".format(k+1,par_name[m],j+1), gal_pos, fmt = "%s\t%s\t%.9g\t%.9g\t%.9g\t%.9g", delimiter = "\t", newline = "\n")
            np.savetxt(quartile_link + "W{}_{}bin{}_shear.dat".format(k+1,par_name[m],j+1), gal_gamma, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")

        del bc, cd, ce, g_name, gal_new, gal_pos, gal_gamma

    ### Gaussian residuals of the parameters for all CFHTLenS galaxies.
    np.savetxt(mass_link + "shear_{}_res.dat".format(par_name[m]), par_opr[m].T, fmt = "%.10E")

    """
    par_min = np.min(par_res[m])
    par_max = np.max(par_res[m])

    print par_min, par_max


    step = 4e-4
    res_bin = np.arange(-0.01,0.015,step)
    res_pltbin = np.arange(-0.0098,0.0148,step)
    
    hist_res, bin_edges_res = np.histogram(par_res[m], bins = res_bin)
    hist_res = hist_res / (np.sum(hist_res) + 0.)
    plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
    #plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
    plt.ylabel("Probability")
    #plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
    plt.xlabel(r"$\delta \log${}".format(par_name[m]))
    plt.xlim((-0.007,0.007))
    plt.savefig(result_link + "shear_{}_res.png".format(par_name[m]))
    #plt.savefig(result_link + "shear_mass_res.png")
    plt.close()
    
    
    del step, res_bin, res_pltbin, hist_res, bin_edges_res, par_min, par_max
    """
    del quar, quar_tit, quar_ind, quar_name, lens_numden, lens_pos, lens_gamma, quartile1, quartile2, quartile3#, quar_stdscore







###############################################

### 25%, 50% and 75% quartiles
quartile1 = np.percentile(mass_res,25)
quartile2 = np.percentile(mass_res,50)
quartile3 = np.percentile(mass_res,75)

mass_bin = 4

quar_ind = [0] * mass_bin

### The indices of the galaxies which belong to the corresponding mass quartiles.
quar_ind[0] = np.where((mass_res <= quartile1))[0]
quar_ind[1] = np.where(((mass_res > quartile1) & (mass_res <= quartile2)))[0]
quar_ind[2] = np.where((mass_res > quartile2) & (mass_res <= quartile3))[0]
quar_ind[3] = np.where((mass_res > quartile3))[0]

### Standardized scores of the parameters.
quar_stdscore = [0] * par_bin
    
for p in range(par_bin):
    quar_stdscore[p] = std_opr.T[quar_ind[p]]
    quar_stdmean = np.mean(quar_stdscore[p], axis = 0)
    quar_stdsdv = np.power(np.diag(np.cov(quar_stdscore[p].T)), 0.5)
    
    quar_stat = np.vstack((quar_stdmean, quar_stdsdv)).T
    
    print p,quar_ind[p],quar_stdscore[p]
        
    np.savetxt(quarvec_link + "halomassbin{}_res_par_statistics.dat".format(p+1), quar_stat, fmt='%.10E\t%.10E', newline = "\n")
    np.savetxt(quarvec_link + "halomassbin{}_res_par_standardized.dat".format(p+1), quar_stdscore[p], fmt='%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E', delimiter = "\t", newline = "\n")

quar = [0] * mass_bin
quar_name = [0] * mass_bin

for i in range(mass_bin):
    quar[i] = mass_res[quar_ind[i]]
    quar_name[i] = gdname[quar_ind[i]]

quar_num = 4
    
quar_tit = [0] * quar_num
for k in range(quar_num):
    quar_tit[k] = [0] * file_num
    for j in range(file_num):
        aaa = np.core.defchararray.find(quar_name[k], "W{}".format(j+1))
        quar_tit[k][j] = quar_name[k][np.where((aaa == 0))[0]]
        print k,j,quar_tit[k][j].size

lens_numden = [0] * file_num
lens_pos = [0] * file_num
lens_gamma = [0] * file_num
    
for k in range(file_num):
    lens_numden[k] = np.loadtxt("W{}_overdensity.dat".format(k+1), usecols = (0,6,7,8,9,10),dtype = {'names':('id','n1','n2','n3','n4','n5'), 'formats':('S20','float64','float64','float64','float64','float64')})
    lens_pos[k] = np.loadtxt("W{}lens.tsv".format(k+1), usecols = (0,1,4,5,24,36),dtype = {'names':('id','field','ra','dec','z','sm'), 'formats':('S20','S20','float64','float64','float64','float64')})
    lens_gamma[k] = np.loadtxt(shear_link + "W{}shear.dat".format(k+1), dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
        
    bc = np.sort(lens_numden[k], order = 'id', kind = 'mergesort')
    cd = np.sort(lens_pos[k], order = 'id', kind = 'mergesort')
    ce = np.sort(lens_gamma[k], order = 'id', kind = 'mergesort')
        
    g_name = [0] * lens_numden[k].size
        
    for j in range(lens_numden[k].size):
        g_name[j] = bc[j][0]
        
        
    for j in range(quar_num):
        c = np.searchsorted(g_name, quar_tit[j][k])
        gal_new = bc[c]
        gal_pos = cd[c]
        gal_gamma = ce[c]
        
        """
        gal_name = [a[0] for a in gal_new]
        gal_name2 = [a[0] for a in gal_pos]
        gal_name3 = [a[0] for a in gal_gamma]
        
        for m in range(len(gal_name)):
            if (gal_name[m] != gal_name2[m]):
                print m, gal_name[m], gal_name2[m]
            if (gal_name[m] != gal_name3[m]):
                print m, gal_name[m], gal_name3[m]

        print k,j,"Check done"
        """
        
        np.savetxt(quartile_link + "W{}_massbin{}_overdensity.dat".format(k+1,j+1), gal_new, fmt = "%s\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E", delimiter = "\t", newline = "\n")
        np.savetxt(quartile_link + "W{}_massbin{}_pos_temp.dat".format(k+1,j+1), gal_pos, fmt = "%s\t%s\t%.9g\t%.9g\t%.9g\t%.9g", delimiter = "\t", newline = "\n")
        np.savetxt(quartile_link + "W{}_massbin{}_shear.dat".format(k+1,j+1), gal_gamma, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")

    del bc, cd, ce, g_name, gal_new, gal_pos, gal_gamma

del quar, quar_tit


np.savetxt(mass_link + "shear_mass_res.dat", mass_res.T, fmt = "%.10E")
#np.savetxt(mass_link + "overden_mass_res.dat", mass_den_res.T, fmt = "%.10E")

exit()

conc_quarbin = 4


print par_mean
print "\n"

c_mbin = [0] * mass_bin

### Mass bins ranges.
c_mbin[0] = ["10.5", "10.8"]
c_mbin[1] = ["11.3", "11.5"]
c_mbin[2] = ["12.0", "12.2"]
c_mbin[3] = ["12.4", "12.7"]


### 4 mass bins.
for i in range(mass_bin):
    
    ### Parameters of the galaxies divided by concentration.
    conc_par = par_opr.T[quar_ind[i]].T
    conc_gdname = gdname[quar_ind[i]]
    
    
    f_str = mock_link + "shear_conc_{}_{}_massbin{}_derivative_intp.dat".format(c_mbin[i][0], c_mbin[i][1], i+1)
    
    print i, quar_ind[i], quar_ind[i].size
    print conc_par[0].size, conc_par.T[0].size
    
    g_conc_derv_load = np.loadtxt(f_str, unpack = True, usecols = (1,))
    conc_var = np.loadtxt(f_str, unpack = True, usecols = (3,))[0]
    conc_avg = np.loadtxt(f_str, unpack = True, usecols = (0,))[0]
    
    ### Interpolate the derivatives times the variance of concentration to the r bins for the observational data.
    g_conc_derv = np.interp(obs_r_bin, mock_r_bin, g_conc_derv_load) * conc_var
    
    ### Concentration residual kernel, d(gamma)/dc * <c^2> * <gamma * gamma>^-1 * <gamma * p> * <p * p>^-1, at fixed mass.
    c_ker = np.dot(g_conc_derv, gau_ker)
    
    np.savetxt(mass_link + "conc/shear_massbin{}_conc_res_kernel.dat".format(i+1), c_ker, fmt = "%.10E")
    
    ### dc
    c_res = np.array(np.dot(c_ker, conc_par))[0]
    
    print np.min(c_res), np.max(c_res)
    
    ### 25%, 50% and 75% quartiles
    c_quartile1 = np.percentile(c_res,25)
    c_quartile2 = np.percentile(c_res,50)
    c_quartile3 = np.percentile(c_res,75)
    
    ### The indices of the galaxies which belong to the corresponding mass quartiles.
    c_quar1_ind = np.where((c_res <= c_quartile1))[0]
    c_quar2_ind = np.where(((c_res > c_quartile1) & (c_res <= c_quartile2)))[0]
    c_quar3_ind = np.where(((c_res > c_quartile2) & (c_res <= c_quartile3)))[0]
    c_quar4_ind = np.where((c_res > c_quartile3))[0]
    
    c_name = [0] * file_num
    
    c_name[0] = conc_gdname[c_quar1_ind]
    c_name[1] = conc_gdname[c_quar2_ind]
    c_name[2] = conc_gdname[c_quar3_ind]
    c_name[3] = conc_gdname[c_quar4_ind]
    
    conc_quar_tit = [0] * quar_num
    
    for k in range(quar_num):
        conc_quar_tit[k] = [0] * file_num
        for j in range(file_num):
            abb = np.core.defchararray.find(c_name[k], "W{}".format(j+1))
            conc_quar_tit[k][j] = c_name[k][np.where((abb == 0))[0]]
            #print k,j
            #print abb
            #print conc_quar_tit[k][j],conc_quar_tit[k][j].size

    #exit()

    conc_numden = [0] * file_num
    conc_pos = [0] * file_num
    conc_gamma = [0] * file_num
    conc_name = [0] * file_num
    
    for k in range(file_num):
        qk_str = quartile_link + "W{}_massbin{}_".format(k+1,i+1)
        
        conc_numden[k] = np.loadtxt("{}overdensity.dat".format(qk_str), usecols =(0,1,2,3,4,5),dtype = {'names':('id','n1','n2','n3','n4','n5'), 'formats':('S20','float64','float64','float64','float64','float64')})
        conc_pos[k] = np.loadtxt("{}pos_temp.dat".format(qk_str), usecols = (0,1,2,3,4,5), dtype = {'names':('id','field','ra','dec','z','sm'), 'formats':('S20','S20','float64','float64','float64','float64')})
        conc_gamma[k] = np.loadtxt("{}shear.dat".format(qk_str), usecols = range(0,43), dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
    
        ab_numden = np.sort(conc_numden[k], order = 'id', kind = 'mergesort')
        ab_pos = np.sort(conc_pos[k], order = 'id', kind = 'mergesort')
        ab_gamma = np.sort(conc_gamma[k], order = 'id', kind = 'mergesort')
        
        cg_name = [0] * conc_numden[k].size
        
        for j in range(conc_numden[k].size):
            cg_name[j] = ab_numden[j][0]

        for j in range(quar_num):
            cc = np.searchsorted(cg_name, conc_quar_tit[j][k])
            
            cgal_new = ab_numden[cc]
            cgal_pos = ab_pos[cc]
            cgal_gamma = ab_gamma[cc]
            
            qkc_str = quartile_link + "conc/W{}_massbin{}_concbin{}_".format(k+1,i+1,j+1)
            
            np.savetxt("{}overdensity.dat".format(qkc_str), cgal_new, fmt = "%s\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E", delimiter = "\t", newline = "\n")
            np.savetxt("{}pos_temp.dat".format(qkc_str), cgal_pos, fmt = "%s\t%s\t%.9g\t%.9g\t%.9g\t%.9g", delimiter = "\t", newline = "\n")
            np.savetxt("{}shear.dat".format(qkc_str), cgal_gamma, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")

            del cc, cgal_pos, cgal_new, cgal_gamma

    
        del qk_str, ab_numden, ab_pos, ab_gamma, cg_name
    
    """
    step = 5e-4
    res_bin = np.arange(-0.02,0.002,step)
    res_pltbin = np.arange(-0.01975,0.00175,step)

    hist_res, bin_edges_res = np.histogram(c_res, bins = res_bin)
    hist_res = hist_res / (np.sum(hist_res) + 0.)
    plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
    #plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
    plt.ylabel("Probability")
    #plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
    plt.title("Mass bin {}".format(i+1))
    plt.xlabel(r"$\delta c$")
    plt.xlim((-0.02,0.002))
    plt.savefig(result_link + "shear_massbin{}_conc_res.png".format(i+1))
    #plt.savefig(result_link + "shear_mass_res.png")
    plt.close()
    """
    
    
    del conc_par, g_conc_derv_load, conc_var, conc_avg, g_conc_derv, c_ker, c_res, c_quartile1, c_quartile2, c_quartile3, c_quar1_ind, c_quar2_ind, c_quar3_ind, c_quar4_ind, conc_numden, conc_pos, conc_gamma, conc_name, conc_gdname, c_name, conc_quar_tit





del gdmat, quar_name, c_mbin



##################################################################




"""
### Standard deviation of dlog(M_h)
mass_std = np.std(mass_res)

### Index of the galaxies within 1-sigma of the mean log(M_h)
mass_ind = np.where(((mass_res <= mass_std) & (mass_res >= -mass_std)))[0]

### Interested galaxies parameters and id
gal_par = par_opr[:,mass_ind]
gal_id = gdname[mass_ind]


"""
"""
step = 2e-4
res_bin = np.arange(-0.005,0.005,step)
res_pltbin = np.arange(-0.0049,0.0049,step)

hist_res, bin_edges_res = np.histogram(mass_res[0], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
#plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Shear estimate, " + r"$\log(\bar{M_h}/h_{70}^{-1}M_\odot) = $" + str('%.5g' % mass_avg))
#plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\log M_h$")
plt.xlim((-0.005,0.005))
plt.savefig(result_link + "shear_mass_res.png")
plt.close()
"""


"""
step = 5e-3
res_bin = np.arange(-0.10,0.10,step)
res_pltbin = np.arange(-0.0975,0.0975,step)

hist_res, bin_edges_res = np.histogram(mass_den_res[0], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
#plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Overdensity estimate, " + r"$\log(\bar{M_h}/h_{70}^{-1}M_\odot) = $" + str('%.5g' % mass_avg))
#plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\log M_h$")
plt.xlim((-0.10,0.10))
plt.savefig(result_link + "overdensity_mass_res.png")
plt.close()
"""
"""
conc_bin = ["11.2", "11.4", "11.6", "11.8", "12.0", "12.2", "12.4", "12.6", "12.8"]

conc_binnum = len(conc_bin) - 1



for i in range(conc_binnum):        ### conc_binnum
    ### Concentration derivatives and concentration variance, convert the hubble factor back.
    g_conc_derv_load = np.loadtxt(mock_link + "shear_conc_{}_{}_derivative_new.dat".format(conc_bin[i], conc_bin[i+1]), unpack = True, usecols = (1,))
    conc_var = np.loadtxt(mock_link + "shear_conc_{}_{}_derivative_new.dat".format(conc_bin[i], conc_bin[i+1]), unpack = True, usecols = (3,))[0]
    conc_avg = np.loadtxt(mock_link + "shear_conc_{}_{}_derivative_new.dat".format(conc_bin[i], conc_bin[i+1]), unpack = True, usecols = (0,))[0]
    #g_conc_den_derv_load = np.loadtxt(mock_link + "overden_conc_{}_{}_derivative.dat".format(conc_bin[i], conc_bin[i+1]), unpack = True, usecols = (1,))
    
    g_conc_derv = np.interp(obs_r_bin, mock_r_bin, g_conc_derv_load) * conc_var
    #g_conc_den_derv = g_conc_den_derv_load * conc_var
    
    #print i, conc_var, g_conc_derv_load, g_conc_derv

    ### Concentration residual kernel.
    conc_ker = np.dot(g_conc_derv, gau_ker)
    #conc_den_ker = np.dot(g_conc_den_derv, den_gau_ker)

    ### d(gamma)/dc * <(dc)^2> * <gamma * gamma>^-1 * <gamma * p> * <p * p>^-1, d(logn)/dc * <(dc)^2> * <logn * logn>^-1 * <logn * p> * <p * p>^-1
    np.savetxt(conc_link + "shear_conc_{}_{}_res_kernel.dat".format(conc_bin[i], conc_bin[i+1]), conc_ker, fmt = "%.10E")
    #np.savetxt(conc_link + "overden_conc_{}_{}_res_kernel.dat".format(conc_bin[i], conc_bin[i+1]), conc_den_ker, fmt = "%.10E")

    ### dc
    conc_res = np.array(np.dot(conc_ker, gal_par))
    #conc_den_res = np.dot(conc_den_ker, par_opr)

    print np.min(conc_res), np.max(conc_res)
    
    ### Rank the galaxies according to their predicted concentration in ascending order.
    #conc_ind = np.argsort(conc_res, kind = 'mergesort')[0]
    #conc_order = conc_res[0][conc_ind]

    ### 25%, 50% and 75% quartiles
    quartile1 = np.percentile(conc_res,25)
    quartile2 = np.percentile(conc_res,50)
    quartile3 = np.percentile(conc_res,75)

    quar1_ind = np.where((conc_res[0] <= quartile1))[0]
    quar2_ind = np.where(((conc_res[0] > quartile1) & (conc_res[0] <= quartile2)))[0]
    quar3_ind = np.where(((conc_res[0] > quartile2) & (conc_res[0] <= quartile3)))[0]
    quar4_ind = np.where((conc_res[0] > quartile3))[0]

    quar = [0] * 4
    quar_name = [0] * 4

    quar[0] = conc_res[0][quar1_ind]
    quar[1] = conc_res[0][quar2_ind]
    quar[2] = conc_res[0][quar3_ind]
    quar[3] = conc_res[0][quar4_ind]

    quar_name[0] = gal_id[quar1_ind]
    quar_name[1] = gal_id[quar2_ind]
    quar_name[2] = gal_id[quar3_ind]
    quar_name[3] = gal_id[quar4_ind]
    
    np.savetxt(conc_link + "shear_conc_{}_{}_res_1sigma.dat".format(conc_bin[i], conc_bin[i+1]), conc_res.T, fmt = "%.10E")
    #np.savetxt(conc_link + "overden_conc_{}_{}_res.dat".format(conc_bin[i], conc_bin[i+1]), conc_den_res.T, fmt = "%.10E")
    
    
    quar_num = 4

    quar_tit = [0] * quar_num
    for k in range(quar_num):
        quar_tit[k] = [0] * file_num
        for j in range(file_num):
            aaa = np.core.defchararray.find(quar_name[k], "W{}".format(j+1))
            quar_tit[k][j] = quar_name[k][np.where((aaa == 0))[0]]
            print k,j,quar_tit[k][j].size

    lens_numden = [0] * file_num
    lens_pos = [0] * file_num
    lens_gamma = [0] * file_num

    for k in range(file_num):
        lens_numden[k] = np.loadtxt("W{}_overdensity.dat".format(k+1), usecols = (0,6,7,8,9,10),dtype = {'names':('id','n1','n2','n3','n4','n5'), 'formats':('S20','float64','float64','float64','float64','float64')})
        lens_pos[k] = np.loadtxt("W{}lens.tsv".format(k+1), usecols = (0,1,4,5,24,36),dtype = {'names':('id','field','ra','dec','z','sm'), 'formats':('S20','S20','float64','float64','float64','float64')})
        lens_gamma[k] = np.loadtxt(shear_link + "W{}shear.dat".format(k+1), dtype = {'names': ('id', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'e8', 'e9', 'e10', 'e11', 'e12', 'e13', 'e14', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n11', 'n12', 'n13', 'n14'), 'formats':('S20', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64')})
    
        bc = np.sort(lens_numden[k], order = 'id', kind = 'mergesort')
        cd = np.sort(lens_pos[k], order = 'id', kind = 'mergesort')
        ce = np.sort(lens_gamma[k], order = 'id', kind = 'mergesort')
    
        g_name = [0] * lens_numden[k].size
    
        for j in range(lens_numden[k].size):
            g_name[j] = bc[j][0]
    
    
        for j in range(quar_num):
            c = np.searchsorted(g_name, quar_tit[j][k])
            gal_new = bc[c]
            gal_pos = cd[c]
            gal_gamma = ce[c]
            np.savetxt(quartile_link + "W{}_{}_{}_quartile{}_overdensity.dat".format(k+1,conc_bin[i], conc_bin[i+1],j+1), gal_new, fmt = "%s\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E", delimiter = "\t", newline = "\n")
            np.savetxt(quartile_link + "W{}_{}_{}_quartile{}_pos_temp.dat".format(k+1,conc_bin[i], conc_bin[i+1],j+1), gal_pos, fmt = "%s\t%s\t%.9g\t%.9g\t%.9g\t%.9g", delimiter = "\t", newline = "\n")
            np.savetxt(quartile_link + "W{}_{}_{}_quartile{}_shear.dat".format(k+1,conc_bin[i], conc_bin[i+1],j+1), gal_gamma, fmt='%s\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E\t%.6E', delimiter = "\t", newline = "\n")
    
    del bc, g_name, gal_new, gal_pos, gal_gamma, quar_tit, quar_name, c, lens_pos, lens_numden
"""
"""
    if (i == 0):
        step = 1e-2
        res_bin = np.arange(-0.3,0.3,step)
        res_pltbin = np.arange(-0.295,0.295,step)
    else:
        step = 5e-4
        res_bin = np.arange(-0.015,0.015,step)
        res_pltbin = np.arange(-0.01475,0.01475,step)
    

    hist_res, bin_edges_res = np.histogram(conc_res[0], bins = res_bin)
    hist_res = hist_res / (np.sum(hist_res) + 0.)
    plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
    #plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
    plt.ylabel("Probability")
    plt.title(r"$\bar{c} = $" + str('%.5g' % conc_avg))
    #plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
    plt.xlabel(r"$\delta c$")
    if (i == 0):
        plt.xlim((-0.3,0.3))
    else:
        plt.xlim((-0.015,0.015))
    plt.savefig(result_link + "shear_conc_{}_{}_res_1sigma.png".format(conc_bin[i], conc_bin[i+1]))
    plt.close()
"""


"""
    conc_res_check = np.dot(g_conc_derv, np.dot(inv_cov_gamma, shear_cpr))

    if (i == 0):
        step = 0.25
        res_bin = np.arange(-6,6,step)
        res_pltbin = np.arange(-5.875,5.875,step)
    else:
        step = 1.25e-2
        res_bin = np.arange(-0.4,0.4,step)
        res_pltbin = np.arange(-0.39375,0.39374,step)

    hist_res, bin_edges_res = np.histogram(conc_res_check, bins = res_bin)
    hist_res = hist_res / (np.sum(hist_res) + 0.)
    plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
    #plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
    plt.ylabel("Probability")
    plt.title(r"$\bar{c} = $" + str('%.5g' % conc_avg))
    #plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
    plt.xlabel(r"$\delta c$")
    if (i == 0):
        plt.xlim((-6,6))
    else:
        plt.xlim((-0.4,0.4))
    plt.savefig(result_link + "actualshear_conc_{}_{}_res_1sigma.png".format(conc_bin[i], conc_bin[i+1]))
    plt.close()
"""
"""
    if (i == 0):
        step = 0.1
        res_bin = np.arange(-6,6,step)
        res_pltbin = np.arange(-5.95, 5.95,step)
    else:
        step = 0.05
        res_bin = np.arange(-6,6,step)
        res_pltbin = np.arange(-5.975,5.975,step)


    hist_res, bin_edges_res = np.histogram(conc_den_res[0], bins = res_bin)
    hist_res = hist_res / (np.sum(hist_res) + 0.)
    plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
    #plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
    plt.ylabel("Probability")
    plt.title(r"$\bar{c} = $" + str('%.5g' % conc_avg))
    #plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
    plt.xlabel(r"$\delta c$")
    plt.xlim((-6,6))
    plt.savefig(result_link + "overden_conc_{}_{}_res.png".format(conc_bin[i], conc_bin[i+1]))
    plt.close()
"""

#del g_conc_derv_load, conc_var, conc_avg, g_conc_derv, conc_ker, quar#,  g_conc_den_derv_load, g_conc_den_derv, conc_den_ker


"""
print quar_name[0]
print "\n\n"

quar_num = 4

quar_tit = [0] * quar_num
for i in range(quar_num):
    quar_tit[i] = [0] * file_num
    for j in range(file_num):
        aaa = np.core.defchararray.find(quar_name[i], "W{}".format(j+1))
        quar_tit[i][j] = quar_name[i][np.where((aaa == 0))[0]]
        print i,j,quar_tit[i][j].size

lens_numden = [0] * file_num
lens_pos = [0] * file_num

for i in range(file_num):
    lens_numden[i] = np.loadtxt("W{}_overdensity.dat".format(i+1), usecols = (0,6,7,8,9,10),dtype = {'names':('id','n1','n2','n3','n4','n5'), 'formats':('S20','float64','float64','float64','float64','float64')})
    lens_pos[i] = np.loadtxt("W{}lens.tsv".format(i+1), usecols = (0,4,5,24),dtype = {'names':('id','ra','dec','z'), 'formats':('S20','float64','float64','float64')})
    
    bc = np.sort(lens_numden[i], order = 'id', kind = 'mergesort')
    cd = np.sort(lens_pos[i], order = 'id', kind = 'mergesort')
    
    g_name = [0] * lens_numden[i].size
    
    for j in range(lens_numden[i].size):
        g_name[j] = bc[j][0]
    
    
    for j in range(quar_num):
        c = np.searchsorted(g_name, quar_tit[j][i])
        gal_new = bc[c]
        gal_pos = cd[c]
        np.savetxt(quartile_link + "W{}_quartile{}_overdensity.dat".format(i+1,j+1), gal_new, fmt = "%s\t%.10E\t%.10E\t%.10E\t%.10E\t%.10E", delimiter = "\t", newline = "\n")
        np.savetxt(quartile_link + "W{}_quartile{}_pos.dat".format(i+1,j+1), gal_pos, fmt = "%s\t%.9g\t%.9g\t%g", delimiter = "\t", newline = "\n")

    del bc, g_name, gal_new
"""















"""
### Predicted shear, delta gamma_p (gau_shear)
gau_shear = np.dot(f_gpar, par_opr).T
pre_shear = np.add(gau_shear, shear_mean)

np.savetxt(gal_link + "gaussian_res_shear.dat", gau_shear, delimiter = "\t", fmt = "%.10E")


### Actual shear
act_shear = gdmat[-14:-5,:].T       ### If all shear bins, use -14:

### <gamma * gamma>^-1
shear_pre_inv = linalg.inv(np.cov(gau_shear.T))

imd_mat = np.dot(shear_pre_inv, gau_shear.T).T

np.savetxt(gal_link + "imd_res_shear.dat", imd_mat, delimiter = "\t", fmt = "%.10E")
"""
"""

### Residual shear
res_shear = np.subtract(act_shear, pre_shear)

resres_mean = np.array(np.mean(res_shear,axis = 0))
resres_var = np.array(np.var(res_shear,axis = 0))


### Check the mean of the scatter, should be zero.
res_mean = np.array(np.mean(gau_shear,axis = 0))
res_var = np.array(np.var(gau_shear,axis = 0))

rad_bin = np.loadtxt("/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/dist_bin2.dat", usecols = (1,), unpack = True)


step = 5


res_bin = np.arange(-150,150,step)
res_pltbin = np.arange(-147.5,147.5,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[0], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][0]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 1a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[0]))
plt.savefig(result_link + "shear1.png")
plt.close()

step = 200

res_bin = np.arange(-5000,5000,step)
res_pltbin = np.arange(-4900,4900,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[0], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(3000,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][0]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][0]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 1b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[0]))
plt.savefig(result_link + "res_shear1.png")
plt.close()

step = 5

res_bin = np.arange(-150,150,step)
res_pltbin = np.arange(-147.5,147.5,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[1], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][1]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][1]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 2a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[1]))
plt.savefig(result_link + "shear2.png")
plt.close()

step = 200

res_bin = np.arange(-5000,5000,step)
res_pltbin = np.arange(-4900,4900,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[1], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(3000,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][1]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][1]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 2b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[1]))
plt.savefig(result_link + "res_shear2.png")
plt.close()

step = 5

res_bin = np.arange(-150,150,step)
res_pltbin = np.arange(-147.5,147.5,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[2], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][2]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][2]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 3a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[2]))
plt.savefig(result_link + "shear3.png")
plt.close()

step = 200

res_bin = np.arange(-4000,4000,step)
res_pltbin = np.arange(-3900,3900,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[2], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(1500,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][2]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][2]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 3b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[2]))
plt.savefig(result_link + "res_shear3.png")
plt.close()


step = 5

res_bin = np.arange(-100,100,step)
res_pltbin = np.arange(-97.5,97.5,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[3], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][3]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][3]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 4a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[3]))
plt.savefig(result_link + "shear4.png")
plt.close()

step = 200

res_bin = np.arange(-4000,4000,step)
res_pltbin = np.arange(-3900,3900,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[3], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(1500,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][3]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][3]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 4b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[3]))
plt.savefig(result_link + "res_shear4.png")
plt.close()

step = 5

res_bin = np.arange(-100,100,step)
res_pltbin = np.arange(-97.5,97.5,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[4], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(50,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][4]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][4]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 5a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[4]))
plt.savefig(result_link + "shear5.png")
plt.close()

step = 150

res_bin = np.arange(-3000,3000,step)
res_pltbin = np.arange(-2925,2925,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[4], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(1000,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][4]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][4]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 5b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[4]))
plt.savefig(result_link + "res_shear5.png")
plt.close()


step = 2

res_bin = np.arange(-50,50,step)
res_pltbin = np.arange(-49,49,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[5], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(20,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][5]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][5]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 6a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[5]))
plt.savefig(result_link + "shear6.png")
plt.close()

step = 100

res_bin = np.arange(-2000,2000,step)
res_pltbin = np.arange(-1950,1950,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[5], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(750,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][5]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][5]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 6b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[5]))
plt.savefig(result_link + "res_shear6.png")
plt.close()

step = 2

res_bin = np.arange(-40,40,step)
res_pltbin = np.arange(-39,39,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[6], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(20,0.06,r'$\mu = $' + str('%.3f' % res_mean[0][6]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][6]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 7a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[6]))
plt.savefig(result_link + "shear7.png")
plt.close()

step = 100

res_bin = np.arange(-2000,2000,step)
res_pltbin = np.arange(-1950,1950,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[6], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(750,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][6]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][6]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 7b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[6]))
plt.savefig(result_link + "res_shear7.png")
plt.close()

step = 1

res_bin = np.arange(-30,30,step)
res_pltbin = np.arange(-29.5,29.5,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[7], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(15,0.05,r'$\mu = $' + str('%.3f' % res_mean[0][7]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][7]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 8a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[7]))
plt.savefig(result_link + "shear8.png")
plt.close()

step = 50

res_bin = np.arange(-1500,1500,step)
res_pltbin = np.arange(-1475,1475,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[7], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(750,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][7]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][7]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 8b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[7]))
plt.savefig(result_link + "res_shear8.png")
plt.close()

step = 1

res_bin = np.arange(-20,20,step)
res_pltbin = np.arange(-19.5,19.5,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[8], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(10,0.05,r'$\mu = $' + str('%.3f' % res_mean[0][8]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][8]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 9a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[8]))
plt.savefig(result_link + "shear9.png")
plt.close()

step = 50

res_bin = np.arange(-1000,1000,step)
res_pltbin = np.arange(-975,975,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[8], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(500,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][8]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][8]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 9b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[8]))
plt.savefig(result_link + "res_shear9.png")
plt.close()


step = 1

res_bin = np.arange(-20,20,step)
res_pltbin = np.arange(-19.5,19.5,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[9], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(10,0.05,r'$\mu = $' + str('%.3f' % res_mean[0][9]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][9]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 10a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[9]))
plt.savefig(result_link + "shear10.png")
plt.close()

step = 50

res_bin = np.arange(-750,750,step)
res_pltbin = np.arange(-725,725,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[9], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(300,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][9]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][9]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 10b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[9]))
plt.savefig(result_link + "res_shear10.png")
plt.close()


step = 0.5

res_bin = np.arange(-15,15,step)
res_pltbin = np.arange(-14.75,14.75,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[10], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(7.5,0.05,r'$\mu = $' + str('%.3f' % res_mean[0][10]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][10]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 11a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[10]))
plt.savefig(result_link + "shear11.png")
plt.close()

step = 25

res_bin = np.arange(-400,400,step)
res_pltbin = np.arange(-387.5,387.5,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[10], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(150,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][10]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][10]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 11b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[10]))
plt.savefig(result_link + "res_shear11.png")
plt.close()

step = 0.5

res_bin = np.arange(-15,15,step)
res_pltbin = np.arange(-14.75,14.75,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[11], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(7.5,0.05,r'$\mu = $' + str('%.3f' % res_mean[0][11]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][11]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 12a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[11]))
plt.savefig(result_link + "shear12.png")
plt.close()

step = 25

res_bin = np.arange(-400,400,step)
res_pltbin = np.arange(-387.5,387.5,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[11], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(150,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][11]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][11]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 12b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[11]))
plt.savefig(result_link + "res_shear12.png")
plt.close()

step = 0.2

res_bin = np.arange(-5,5,step)
res_pltbin = np.arange(-4.9,4.9,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[12], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(3.5,0.05,r'$\mu = $' + str('%.3f' % res_mean[0][12]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][12]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 13a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[12]))
plt.savefig(result_link + "shear13.png")
plt.close()

step = 10

res_bin = np.arange(-200,200,step)
res_pltbin = np.arange(-195,195,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[12], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(50,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][12]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][12]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 13b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[12]))
plt.savefig(result_link + "res_shear13.png")
plt.close()

step = 0.2

res_bin = np.arange(-5,5,step)
res_pltbin = np.arange(-4.9,4.9,step)

hist_res, bin_edges_res = np.histogram(gau_shear.T[13], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(3.5,0.05,r'$\mu = $' + str('%.3f' % res_mean[0][13]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % res_var[0][13]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 14a " + r"$\langle\gamma p\rangle\langle pp\rangle^{-1}\vec{p}$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[13]))
plt.savefig(result_link + "shear14.png")
plt.close()

step = 10

res_bin = np.arange(-200,200,step)
res_pltbin = np.arange(-195,195,step)

hist_res, bin_edges_res = np.histogram(res_shear.T[13], bins = res_bin)
hist_res = hist_res / (np.sum(hist_res) + 0.)
plt.bar(res_pltbin, hist_res, step, color = 'm', alpha = 0.3)
plt.text(50,0.06,r'$\mu = $' + str('%.3f' % resres_mean[0][13]) + '\n' + r'$\sigma^2 = $' + str('%.3f' % resres_var[0][13]),fontsize = 12)
plt.ylabel("Probability")
plt.title("Fig. 14b " + r"$\gamma_{actual}(r) - \gamma_{predicted}(r)$")
plt.xlabel(r"$\delta\gamma(r = {} Mpc/h)$".format('%.4f' % rad_bin[13]))
plt.savefig(result_link + "res_shear14.png")
plt.close()
"""

"""
y_la = [0] * 29

y_la[0] = "log(a_world), 0"
y_la[1] = "log(b_world), 1"
y_la[2] = "log(flux radius), 2"
y_la[3] = "log(iso area), 3"
y_la[4] = r"$M_u$, 4"
y_la[5] = r"$M_g$, 5"
y_la[6] = r"$M_r$, 6"
y_la[7] = r"$M_i$, 7"
y_la[8] = r"$M_z$, 8"
y_la[9] = r"$\log(M_\star)$, 9"
y_la[10] = r"$\delta(3 Mpc)$, 10"
y_la[11] = r"$\delta(5 Mpc)$, 11"
y_la[12] = r"$\delta(7 Mpc)$, 12"
y_la[13] = r"$\delta(9 Mpc)$, 13"
y_la[14] = r"$\delta(11 Mpc)$, 14"
y_la[15] = r"$\gamma$, innermost, 15"
y_la[28] = r"$\gamma$, outermost, 28"

for i in range(16,28):
    y_la[i] = r"$\gamma$"

plt.pcolor(corrmat)
plt.colorbar()
plt.yticks(np.arange(0.5,29.5),y_la)
plt.xticks(np.arange(0.5,29.5),range(0,29))
plt.xlim((0,29))
plt.ylim((0,29))
plt.title("Correlation matrix")
plt.savefig(result_link + "corrmat.png", bbox_inches='tight')
plt.close()


x_bin = np.arange(7,12,0.1)
x_pltbin = np.arange(7.05,11.95,0.1)
#x_bin = np.arange(-4000,4000,100)
#x_pltbin = np.arange(-3950,3950,100)

print x_bin.size, x_pltbin.size

#shear_mat[13] = np.log10(shear_mat[13])

#x_bin = np.arange(10000,60000,2000)
#x_pltbin = np.arange(11000,59000,2000)

bar_width = 0.1

hist, bin_edge = np.histogram(gdmat[9], bins = x_bin)

plt.bar(x_pltbin, hist, bar_width, color = 'm', alpha = 0.3)
plt.xlabel(r'$\log(M_\star)$')
#plt.ylabel('Probability')
#plt.show()
plt.savefig(result_link + "starmass.png")
plt.close()
"""

