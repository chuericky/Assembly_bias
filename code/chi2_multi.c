/* This program fits the shear profile with the central NFW profile, fixing baryonic contribution to obtain the halo mass and concentration of the halos. M_h and conc. are two free parameters.
    Generalized chi-square fitting by covariance matrix, fitting using brute force method.
    Input: with / without (argv[1])
    Updated: Aug 15, 2015
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tool.h"
#include "cosmopar.h"
#include "halo_model.h"

typedef struct {
    double g_t, d_sig;
} shear_profile;

typedef struct {
    char name[40];
    double star_mass_log, star_mass_linear, R_200, R_200_Mpc, z, ang_dis;
} file_name;

int main(int argc, char *argv[]) {
    char shear_cat[150], pre_mass[150], shear_name[150], para_name[150], plot_name[150], mass_name[150], invcov_name[150], parameter_name[40];
    FILE *shear_cat_file, *pre_mass_file, *shear_file, *para_file, *plot_file, *mass_file, *invcov_file;
    int i = 0, j = 0, k = 0, s = 0, t = 0, p = 0, r_bin = 0, r_plot_bin = 0, r_plot_bin2 = 0, shear_cat_num = 0, pre_mass_num = 0, spot = 0, file_count = 0, s_file_count = 0, sr_bin = 0, fitnum = 0, m_num = 0, fitmin_mass = 0, fitmin_conc = 0, fitmin_mass2 = 0, fitmin_conc2 = 0, c_num = 0, fitnum_count = 0, cov_num = 0;
    double zl_ref = 0, zs_ref = 0, ang_dis_ref = 0, eff_ref = 0, Sig_cr_ref = 0, r_min = 0, r_max = 0, r_fitmin = 0, dlogr = 0, r_plot_min = 0, r_plot_max = 0, dlogr_plot = 0, d2r = M_PI / 180., arcs2d = 1 / 3600., m_min = 0, m_max = 0, dlogm = 0, c_min = 0, c_max = 0, c_diff = 0, min_chisq = 0, min_red_chisq = 0, m_opt = 0, m_opt2 = 0, m_opt3 = 0, c_opt = 0, c_opt2 = 0, c_opt3 = 0, r_new = 0, eps_c = 0, eps_m = 0, eps_compare = 0, eps2_c = 0, eps2_m = 0, eps2_compare = 1, eps3_m = 0, eps3_c = 0, eps3_compare = 1, onesigma_dchi = 0, onesigma_dredchi = 0, onesigma = 0, m_up = 0, m_low = 0, c_up = 0, c_low = 0, fit_r = 0;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 2) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    zl_ref = 0.27;
    ang_dis_ref = Da(zl_ref);
    printf("%lE\n", Da(0.3));
    
    getchar();
    
    // Number of rows (columns) in the inverse covariance matrix.
    cov_num = 11;
    
    
    // Read in stellar masses in each bin.
    strcpy (shear_cat, shear_gal_link);
    strcpy (parameter_name, argv[1]);
    strcat (shear_cat, parameter_name);
    strcat (shear_cat, "/stellar_mass.txt");
    
    shear_cat_file = fopen (shear_cat, "r");
    
    // Make sure the source name file is present.
    if (!shear_cat_file) {
        printf("Cannot find the shear name file!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count number of lines in the shear name file.
    shear_cat_num = count_line (spot, shear_cat_file, file_count);
    
    // An array to store the radius bin and gamma_tan.
    file_name *lens_name = malloc (shear_cat_num * sizeof *lens_name);
    
    file_count = 0;
    
    while (fscanf(shear_cat_file, "%s %lf %lf", &lens_name[file_count].name, &lens_name[file_count].star_mass_log, &lens_name[file_count].z) == 3) {
        lens_name[file_count].star_mass_linear = pow(10, lens_name[file_count].star_mass_log);
        lens_name[file_count].ang_dis = Da(lens_name[file_count].z);
        file_count++;
    }
    
    strcpy (pre_mass, shear_gal_link);
    strcat (pre_mass, parameter_name);
    strcat (pre_mass, "/prelim_mass.txt");
    
    pre_mass_file = fopen (pre_mass, "r");
    
    pre_mass_num = 0;
    
    while (fscanf(pre_mass_file, "%*s %lf", &lens_name[pre_mass_num].R_200) == 1) {
        lens_name[pre_mass_num].R_200_Mpc = lens_name[pre_mass_num].R_200 / kpc_2_pc;
        pre_mass_num++;
    }
    
    if (file_count != pre_mass_num) {
        printf("Line count in %s and %s don't match!!!\n", shear_cat, pre_mass);
        exit(1);
    }
    
    // Log radius bins, in arcsecs.
    r_bin = 30;
    r_min = 5;
    r_max = 2500;
    r_fitmin = 0;   // Range of r to be fit
    dlogr = (log10(r_max) - log10(r_min)) / r_bin;
    
    fit_r = 0.2;    // Fitted radius of the halos, in Mpc h^-1.
    
    // Radius bins for plotting, in arcsecs.
    r_plot_bin = 230;
    r_plot_bin2 = 200;
    r_plot_min = 2;
    r_plot_max = 200;
    dlogr_plot = (log10(r_plot_max) - log10(r_plot_min)) / r_plot_bin2;
    
    double *r_phy = malloc (r_bin * sizeof (double));
    
    // Centers of the annulus width (in arcsec), and convert to the physical distance away from the lens center (in h_70^-1 Mpc).
    for (i = 0; i < r_bin; i++) {
        r_phy[i] = ((r_min * pow(10, (i + 0.5) * dlogr)) * d2r * arcs2d) * ang_dis_ref;
    }
    
    // The radius bin for plotting, in pc.
    double *r_plot_phy = malloc (r_plot_bin * sizeof (double));
    
    for (i = 0; i < r_plot_bin; i++) {
        r_plot_phy[i] = ((r_plot_min * pow(10, i * dlogr_plot)) * d2r * arcs2d) * ang_dis_ref * Mpc_2_pc;
    }
    
    double *r_fitmax = malloc (file_count * sizeof (double));
    int *r_fitnum = malloc (file_count * sizeof (int));
    
    
    // Open each shear file.
    for (i = 0; i < file_count; i++) { //file_count
        
        printf("\nChi-squared fitting of %s started!\n", lens_name[i].name);
        
        r_fitnum[i] = 0;
        r_fitmax[i] = fit_r; // lens_name[i].R_200_Mpc;
        
        for (j = 0; j < r_bin; j++) {
            if (r_phy[j] > r_fitmin && r_phy[j] < r_fitmax[i])
                r_fitnum[i]++;
        }
        
        // Read in the shear files for each bin.
        strcpy (shear_name, shear_gal_link);
        strcat (shear_name, parameter_name);
        strcat (shear_name, "/");
        strcat (shear_name, lens_name[i].name);
        
        shear_file = fopen (shear_name, "r");
        
        spot = 0;
        s_file_count = -1;
        
        sr_bin = count_line (spot, shear_file, s_file_count);
        
        shear_profile *s_info = malloc (sr_bin * sizeof *s_info);
        
        s_file_count = 0;
        
        while (fscanf(shear_file, "%lf %*lf", &s_info[s_file_count].d_sig) == 1) {
            s_file_count++;
        }
        
        // Input the inverse covariance matrix.
        strcpy (invcov_name, shear_gal_link);
        strcat (invcov_name, parameter_name);
        strcat (invcov_name, "/invcov_");
        strcat (invcov_name, lens_name[i].name);
        
        invcov_file = fopen (invcov_name, "r");
        
        double **invcov = malloc (cov_num * sizeof (double *));
        
        for (j = 0; j < cov_num; j++) {
            invcov[j] = malloc (cov_num * sizeof (double));
            for (k = 0; k < cov_num; k++) {
                fscanf(invcov_file, "%lE", &invcov[j][k]);
            }
        }
        
        // Setting the inverse covariance matrix diagonal.
        /*for (j = 0; j < cov_num; j++) {
            for (k = 0; k < cov_num; k++) {
                if (j != k) {
                    invcov[j][k] = 0;
                }
            }
        }*/
        
        
        // Iterate until the M_200 converges
        p = 0;
        eps2_compare = 1;
        
        do {
            // Iteration after the first time.
            if (p > 0) {
                r_fitnum[i] = 0;
                r_fitmax[i] = fit_r; // r_new
                
                for (j = 0; j < r_bin; j++) {
                    if (r_phy[j] > r_fitmin && r_phy[j] < r_fitmax[i])
                        r_fitnum[i]++;
                }
            }
            
            double *r_fit = malloc (r_fitnum[i] * sizeof (double));
            double *dsig_data = malloc (r_fitnum[i] * sizeof (double));
        
            fitnum = 0;
        
            // Input the region where fitting takes place.
            // r_fit in pc h_70^-1, dsig_data in M_sun pc^-2 h_70.
            for (j = 0; j < s_file_count; j++) {
                if (r_phy[j] > r_fitmin && r_phy[j] < r_fitmax[i]) {
                    r_fit[fitnum] = r_phy[j] * Mpc_2_pc;
                    dsig_data[fitnum] = s_info[j].d_sig;
                    fitnum++;
                }
            }
            
            // Double check if the number of data points being fit is correct.
            if (r_fitnum[i] != fitnum) {
                printf("Number of fitting distance bin is not correct! %s \n", lens_name[i].name);
                exit(1);
            }
        
            m_num = 200;
            c_num = 200;
        
            double *m_bin = malloc ((m_num + 1) * sizeof (double));
            double *c_bin = malloc ((c_num + 1) * sizeof (double));
        
            // Chi-squared fitting.
            double *dsig_bar_fit = malloc (fitnum * sizeof (double));
            
            // d_sigma due to baryon, this is not to be fit.
            dsig_baryon(lens_name[i].star_mass_linear, r_fit, dsig_bar_fit, fitnum);
    
            // Log mass bins, for chi-square fitting.
            m_min = 1e10;
            m_max = 1e14;
            dlogm = (log10(m_max) - log10(m_min)) / m_num;
        
            // Bins for concentration, for chi-square fitting.
            c_min = 0.2;
            c_max = lens_name[i].R_200_Mpc / r_phy[0];  // Make sure the maximum concentration aligns with the innermost bin.
            if (p > 0) {
                c_max = r_new / r_phy[0];
            }
            
            c_diff = (c_max - c_min) / c_num;
        
        
            for (j = 0; j <= m_num; j++) {
                m_bin[j] = m_min * pow(10, j * dlogm);
            }
        
            for (j = 0; j <= c_num; j++) {
                c_bin[j] = c_min + j * c_diff;
            }
        
            // Making a coarse map to locate the local minima.
            double **chi_sq_map = malloc ((m_num + 1) * sizeof (double *));
            double **redchi_sq_map = malloc ((m_num + 1) * sizeof (double *));
        
            for (j = 0; j <= m_num; j++) {
                // j = 0 for mass fitting, j = 1 for concentration fitting.
                chi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                redchi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
            
                for (k = 0; k <= c_num; k++) {
                    chi_sq_map[j][k] = 0;
                    redchi_sq_map[j][k] = 0;
                }
            }
            
            double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
            double *diff_vector = malloc (fitnum * sizeof (double));
            double *intermediate_vector = malloc (fitnum * sizeof (double));
        
            // Find chi^2
            for (j = 0; j <= m_num; j++) { // m_num
                for (k = 0; k <= c_num; k++) { // c_num
            
                    dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                    
                    fitnum_count = 0;
                    
                    // Vector storing the difference between the data and model.
                    for (s = 0; s < fitnum; s++) {
                        diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                        intermediate_vector[s] = 0;
                    }
                    
                    // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                    for (s = 0; s < fitnum; s++) {
                        for (t = 0; t < fitnum; t++) {
                            intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                        }
                    }
            
                    for (s = 0; s < fitnum; s++) {
                        // Assuming all central galaxies. Skipping over all the NAN data points.
                        /*if (dsig_data[s] != dsig_data[s]) {
                            continue;
                        }*/
                        // Chi-square value, chi^2 = c_1*N * a_N*1
                        chi_sq_map[j][k] += intermediate_vector[s] * diff_vector[s];
                        fitnum_count++;
                    }
                    redchi_sq_map[j][k] = chi_sq_map[j][k] / (fitnum_count - 3);
                }
            }
            
            free (diff_vector);
            free (intermediate_vector);
            free (dsig_cen1h_fit);
        
            // Find the minimum index for mass and concentration in brute force fitting.
            fitmin_mass = 0;
            fitmin_conc = 0;
            min_matrix((m_num + 1), (c_num + 1), chi_sq_map, &fitmin_mass, &fitmin_conc, 0);
        
            // Assigning the values of the arrays for the next iteration.
            if (fitmin_mass != 0 && fitmin_mass != m_num) {
                m_min = m_bin[fitmin_mass - 1];
                m_max = m_bin[fitmin_mass + 1];
                dlogm = (log10(m_max) - log10(m_min)) / m_num;
                for (j = 0; j <= m_num; j++) {
                    m_bin[j] = m_min * pow(10, j * dlogm);
                }
            }
        
            if (fitmin_conc != 0 && fitmin_conc != c_num) {
                c_min = c_bin[fitmin_conc - 1];
                c_max = c_bin[fitmin_conc + 1];
                c_diff = (c_max - c_min) / c_num;
                for (j = 0; j <= c_num; j++) {
                    c_bin[j] = c_min + j * c_diff;
                }
            }
        
            if (fitmin_mass == 0) {
                m_min = m_bin[fitmin_mass];
                m_max = m_bin[fitmin_mass + 1];
                dlogm = (log10(m_max) - log10(m_min)) / m_num;
                for (j = 0; j <= m_num; j++) {
                    m_bin[j] = m_min * pow(10, j * dlogm);
                }
            }
        
            if (fitmin_mass == m_num) {
                m_min = m_bin[fitmin_mass - 1];
                m_max = m_bin[fitmin_mass];
                dlogm = (log10(m_max) - log10(m_min)) / m_num;
                for (j = 0; j <= m_num; j++) {
                    m_bin[j] = m_min * pow(10, j * dlogm);
                }
            }
        
            if (fitmin_conc == 0) {
                c_min = c_bin[fitmin_conc];
                c_max = c_bin[fitmin_conc + 1];
                c_diff = (c_max - c_min) / c_num;
                for (j = 0; j <= c_num; j++) {
                    c_bin[j] = c_min + j * c_diff;
                }
            }
        
            if (fitmin_conc == c_num) {
                c_min = c_bin[fitmin_conc - 1];
                c_max = c_bin[fitmin_conc];
                c_diff = (c_max - c_min) / c_num;
                for (j = 0; j <= c_num; j++) {
                    c_bin[j] = c_min + j * c_diff;
                }
            }
        
            free (chi_sq_map);
            free (redchi_sq_map);
        
    
            // Find the minimum chi-squared value by grid-searching method.
            t = 0;
        
            do {
                if (t > 0) {
                    if (fitmin_mass != 0 && fitmin_mass != m_num) {
                        m_min = m_bin[fitmin_mass - 1];
                        m_max = m_bin[fitmin_mass + 1];
                        dlogm = (log10(m_max) - log10(m_min)) / m_num;
                        for (j = 0; j <= m_num; j++) {
                            m_bin[j] = m_min * pow(10, j * dlogm);
                        }
                    }
                
                    if (fitmin_conc != 0 && fitmin_conc != c_num) {
                        c_min = c_bin[fitmin_conc - 1];
                        c_max = c_bin[fitmin_conc + 1];
                        c_diff = (c_max - c_min) / c_num;
                        for (j = 0; j <= c_num; j++) {
                            c_bin[j] = c_min + j * c_diff;
                        }
                    }
                
                    if (fitmin_mass == 0) {
                        m_min = m_bin[fitmin_mass];
                        m_max = m_bin[fitmin_mass + 1];
                        dlogm = (log10(m_max) - log10(m_min)) / m_num;
                        for (j = 0; j <= m_num; j++) {
                            m_bin[j] = m_min * pow(10, j * dlogm);
                        }
                    }
                
                    if (fitmin_mass == m_num) {
                        m_min = m_bin[fitmin_mass - 1];
                        m_max = m_bin[fitmin_mass];
                        dlogm = (log10(m_max) - log10(m_min)) / m_num;
                        for (j = 0; j <= m_num; j++) {
                            m_bin[j] = m_min * pow(10, j * dlogm);
                        }
                    }
                
                    if (fitmin_conc == 0) {
                        c_min = c_bin[fitmin_conc];
                        c_max = c_bin[fitmin_conc + 1];
                        c_diff = (c_max - c_min) / c_num;
                        for (j = 0; j <= c_num; j++) {
                            c_bin[j] = c_min + j * c_diff;
                        }
                    }
                
                    if (fitmin_conc == c_num) {
                        c_min = c_bin[fitmin_conc - 1];
                        c_max = c_bin[fitmin_conc];
                        c_diff = (c_max - c_min) / c_num;
                        for (j = 0; j <= c_num; j++) {
                            c_bin[j] = c_min + j * c_diff;
                        }
                    }
                }

                // Making a coarse map to locate the local minima.
                double **chi_sq_map = malloc ((m_num + 1) * sizeof (double *));
                double **redchi_sq_map = malloc ((m_num + 1) * sizeof (double *));
            
                for (j = 0; j <= m_num; j++) {
                    // j = 0 for mass fitting, j = 1 for concentration fitting.
                    chi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                    redchi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                
                    for (k = 0; k <= c_num; k++) {
                        chi_sq_map[j][k] = 0;
                        redchi_sq_map[j][k] = 0;
                    }
                }
                
                double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
                double *diff_vector = malloc (fitnum * sizeof (double));
                double *intermediate_vector = malloc (fitnum * sizeof (double));
            
                // Find chi^2
                for (j = 0; j <= m_num; j++) { // m_num
                    for (k = 0; k <= c_num; k++) { // c_num
                    
                        dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                        
                        fitnum_count = 0;
                        
                        // Vector storing the difference between the data and model.
                        for (s = 0; s < fitnum; s++) {
                            diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                            intermediate_vector[s] = 0;
                        }
                        
                        // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                        for (s = 0; s < fitnum; s++) {
                            for (t = 0; t < fitnum; t++) {
                                intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                            }
                        }
                        
                        for (s = 0; s < fitnum; s++) {
                            // Assuming all central galaxies. Skipping over all the NAN data points.
                            /*if (dsig_data[s] != dsig_data[s]) {
                             continue;
                             }*/
                            // Chi-square value, chi^2 = c_1*N * a_N*1
                            chi_sq_map[j][k] += intermediate_vector[s] * diff_vector[s];
                            fitnum_count++;
                        }
                        
                        redchi_sq_map[j][k] = chi_sq_map[j][k] / (fitnum_count - 3);
                        // First loop over conc., then mass.
                        
                    }
                }
                
                free (diff_vector);
                free (intermediate_vector);
                free (dsig_cen1h_fit);
                
            
                // Find the minimum index for mass and concentration in brute force fitting.
                fitmin_mass = 0;
                fitmin_conc = 0;
                min_matrix((m_num + 1), (c_num + 1), chi_sq_map, &fitmin_mass, &fitmin_conc, 0);
                
                min_red_chisq = redchi_sq_map[fitmin_mass][fitmin_conc];
            
                // do-while loop stopping criteria.
                if (fitmin_mass != 0 && fitmin_mass != m_num) {
                    eps_m = fabs((m_bin[fitmin_mass + 1] - m_bin[fitmin_mass - 1]) / m_bin[fitmin_mass]);
                }
            
                if (fitmin_conc != 0 && fitmin_conc != c_num) {
                    eps_c = fabs((c_bin[fitmin_conc + 1] - c_bin[fitmin_conc - 1]) / c_bin[fitmin_conc]);
                }
            
                if (fitmin_mass == 0) {
                    eps_m = fabs((m_bin[fitmin_mass + 1] - m_bin[fitmin_mass]) / m_bin[fitmin_mass]);
                }
            
                if (fitmin_mass == m_num) {
                    eps_m = fabs((m_bin[fitmin_mass] - m_bin[fitmin_mass - 1]) / m_bin[fitmin_mass]);
                }
            
                if (fitmin_conc == 0) {
                    eps_c = fabs((c_bin[fitmin_conc + 1] - c_bin[fitmin_conc]) / c_bin[fitmin_conc]);
                }
            
                if (fitmin_conc == c_num) {
                    eps_c = fabs((c_bin[fitmin_conc] - c_bin[fitmin_conc - 1]) / c_bin[fitmin_conc]);
                }
            
                m_opt = m_bin[fitmin_mass];
                c_opt = c_bin[fitmin_conc];
            
                eps_compare = max_compare(eps_m, eps_c);
                

                t++;
                free (chi_sq_map);
                free (redchi_sq_map);
            
            } while (eps_compare >= 1e-5);
            
            r_new = rad_200(rho_crit(lens_name[i].z), m_opt) / Mpc_2_pc;
            
            if (p > 0) {
                eps2_m = fabs((m_opt2 - m_opt) / m_opt2);
                eps2_c = fabs((c_opt2 - c_opt) / c_opt2);
                eps2_compare = max_compare(eps2_m, eps2_c);
            }
            
            c_opt2 = c_opt;
            m_opt2 = m_opt;
        
            p++;
            
            if (p == 50) {
                printf("# of iterations exceed 50!\n");
                break;
            }
            
            free (m_bin);
            free (c_bin);
            free (r_fit);
            free (dsig_data);
            free (dsig_bar_fit);
            
        } while (eps2_compare >= 1e-3);
        
        
        // Find the 1-sigma level of m_opt and c_opt.
        onesigma = 2.3;
        onesigma_dredchi = (min_red_chisq * (fitnum_count - 3) + onesigma) / (fitnum_count - 3);
        onesigma_dchi = min_red_chisq * (fitnum_count - 3) + onesigma;
        
        m_opt3 = m_opt;
        c_opt3 = c_opt;

        printf("%lf\t%f\t%f\n", min_red_chisq,  min_red_chisq * (fitnum_count - 3), onesigma_dchi);
    
        
        p = 0;
        
        // Upper quartile of 1-sigma mass.
        // Iterate until the M_200 converges.
        
        do {
            // Iteration after the first time.
            /*if (p > 0) {
                r_fitnum[i] = 0;
                r_fitmax[i] = r_new;
                
                for (j = 0; j < r_bin; j++) {
                    if (r_phy[j] > r_fitmin && r_phy[j] < r_fitmax[i])
                    r_fitnum[i]++;
                }
            }*/

            
            double *r_fit = malloc (r_fitnum[i] * sizeof (double));
            double *dsig_data = malloc (r_fitnum[i] * sizeof (double));
            
            fitnum = 0;
            
            // Input the region where fitting takes place.
            // r_fit in pc h_70^-1, dsig_data in M_sun pc^-2 h_70.
            for (j = 0; j < s_file_count; j++) {
                if (r_phy[j] > r_fitmin && r_phy[j] < r_fitmax[i]) {
                    r_fit[fitnum] = r_phy[j] * Mpc_2_pc;
                    dsig_data[fitnum] = s_info[j].d_sig;
                    fitnum++;
                }
            }
            
            // Double check if the number of data points being fit is correct.
            if (r_fitnum[i] != fitnum) {
                printf("Number of fitting distance bin is not correct! %s \n", lens_name[i].name);
                exit(1);
            }
            
            m_num = 200;
            c_num = 200;
            
            double *m_bin = malloc ((m_num + 1) * sizeof (double));
            double *c_bin = malloc ((c_num + 1) * sizeof (double));
            
            // Chi-squared fitting.
            double *dsig_bar_fit = malloc (fitnum * sizeof (double));
            
            // d_sigma due to baryon, this is not to be fit.
            dsig_baryon(lens_name[i].star_mass_linear, r_fit, dsig_bar_fit, fitnum);
            
            // Log mass bins, for chi-square fitting.
            m_min = m_opt3;
            m_max = 1e14;
            dlogm = (log10(m_max) - log10(m_min)) / m_num;
            
            // Bins for concentration, for chi-square fitting.
            c_min = 0.2;
            c_max = lens_name[i].R_200_Mpc / r_phy[0];  // Make sure the maximum concentration aligns with the innermost bin.
            /*if (p > 0) {
                c_max = r_new / r_phy[0];
            }*/
            
            c_diff = (c_max - c_min) / c_num;
            
            for (j = 0; j <= m_num; j++) {
                m_bin[j] = m_min * pow(10, j * dlogm);
            }
            
            for (j = 0; j <= c_num; j++) {
                c_bin[j] = c_min + j * c_diff;
            }
            
            // Making a coarse map to locate the local minima.
            double **chi_sq_map2 = malloc ((m_num + 1) * sizeof (double *));
            double **redchi_sq_map2 = malloc ((m_num + 1) * sizeof (double *));
            
            for (j = 0; j <= m_num; j++) {
                // j = 0 for mass fitting, j = 1 for concentration fitting.
                chi_sq_map2[j] = malloc ((c_num + 1) * sizeof (double));
                redchi_sq_map2[j] = malloc ((c_num + 1) * sizeof (double));
                
                for (k = 0; k <= c_num; k++) {
                    chi_sq_map2[j][k] = 0;
                    redchi_sq_map2[j][k] = 0;
                }
            }
            
            double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
            double *diff_vector = malloc (fitnum * sizeof (double));
            double *intermediate_vector = malloc (fitnum * sizeof (double));
            
            // Find chi^2
            for (j = 0; j <= m_num; j++) { // m_num
                for (k = 0; k <= c_num; k++) { // c_num
                    
                    dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                    
                    fitnum_count = 0;
                    
                    // Vector storing the difference between the data and model.
                    for (s = 0; s < fitnum; s++) {
                        diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                        intermediate_vector[s] = 0;
                    }
                    
                    // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                    for (s = 0; s < fitnum; s++) {
                        for (t = 0; t < fitnum; t++) {
                            intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                        }
                    }
                    
                    for (s = 0; s < fitnum; s++) {
                        // Assuming all central galaxies. Skipping over all the NAN data points.
                        /*if (dsig_data[s] != dsig_data[s]) {
                         continue;
                         }*/
                        // Chi-square value, chi^2 = c_1*N * a_N*1
                        chi_sq_map2[j][k] += intermediate_vector[s] * diff_vector[s];
                        fitnum_count++;
                    }
                    redchi_sq_map2[j][k] = chi_sq_map2[j][k] / (fitnum_count - 3);
                }
            }
            
            free (diff_vector);
            free (intermediate_vector);
            free (dsig_cen1h_fit);
            
            int *mass_s_count = malloc ((m_num + 1) * sizeof (int));
            
            for (j = 0; j <= m_num; j++) {
                mass_s_count[j] = 0;
            }
            
            for (j = 0; j <= m_num; j++) {
                for (k = 0; k <= c_num; k++) {
                    //printf("%d\t%d\t%lE\t%lf\t%lf\t%lf\n", j, k, m_bin[j], c_bin[k], chi_sq_map2[j][k], onesigma_dchi);
                    if (chi_sq_map2[j][k] < onesigma_dchi) {
                        mass_s_count[j]++;
                    }
                }
            }
            
            for (j = 0; j <= m_num; j++) {
                if (mass_s_count[j] == 0) {
                    break;
                }
            }
            
            
            if (j == m_num + 1) {
                m_up = m_bin[m_num];
                printf("Lower bound of the upper 1-sigma mass: %lf\n", m_up);
                break;
            }
            
            
            // Assigning the values of the arrays for the next iteration.
            for (j = 0; j < m_num; j++) {
                if (mass_s_count[j] != 0 && mass_s_count[j + 1] == 0) {
                    m_min = m_bin[j];
                    m_max = m_bin[j + 1];
                    break;
                }
                if (mass_s_count[j] == 0 && mass_s_count[j + 1] != 0) {
                    m_min = m_bin[j];
                    m_max = m_bin[j + 1];
                    break;
                }
            }
            
            /*for (j = 0; j <= m_num; j++) {
                printf("%d\t%lE\t%d\n", j, m_bin[j], mass_s_count[j]);
            }*/
            
            dlogm = (log10(m_max) - log10(m_min)) / m_num;
            for (j = 0; j <= m_num; j++) {
                m_bin[j] = m_min * pow(10, j * dlogm);
            }
        
            
            free (mass_s_count);
            
            
            for (j = 0; j <= m_num; j++) {
                free (chi_sq_map2[j]);
                free (redchi_sq_map2[j]);
            }
            
            free (chi_sq_map2);
            free (redchi_sq_map2);
            
            // Find the minimum chi-squared value by grid-searching method.
            t = 0;
            
            do {
                // Making a coarse map to locate the local minima.
                double **chi_sq_map = malloc ((m_num + 1) * sizeof (double *));
                double **redchi_sq_map = malloc ((m_num + 1) * sizeof (double *));
                
                for (j = 0; j <= m_num; j++) {
                    // j = 0 for mass fitting, j = 1 for concentration fitting.
                    chi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                    redchi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                    
                    for (k = 0; k <= c_num; k++) {
                        chi_sq_map[j][k] = 0;
                        redchi_sq_map[j][k] = 0;
                    }
                }
                
                double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
                double *diff_vector = malloc (fitnum * sizeof (double));
                double *intermediate_vector = malloc (fitnum * sizeof (double));
                
                // Find chi^2
                for (j = 0; j <= m_num; j++) { // m_num
                    for (k = 0; k <= c_num; k++) { // c_num
                        
                        dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                        
                        fitnum_count = 0;
                        
                        // Vector storing the difference between the data and model.
                        for (s = 0; s < fitnum; s++) {
                            diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                            intermediate_vector[s] = 0;
                        }
                        
                        // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                        for (s = 0; s < fitnum; s++) {
                            for (t = 0; t < fitnum; t++) {
                                intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                            }
                        }
                        
                        for (s = 0; s < fitnum; s++) {
                            // Assuming all central galaxies. Skipping over all the NAN data points.
                            /*if (dsig_data[s] != dsig_data[s]) {
                             continue;
                             }*/
                            // Chi-square value, chi^2 = c_1*N * a_N*1
                            chi_sq_map[j][k] += intermediate_vector[s] * diff_vector[s];
                            fitnum_count++;
                        }
                        
                        redchi_sq_map[j][k] = chi_sq_map[j][k] / (fitnum_count - 3);
                        // First loop over conc., then mass.
                    }
                }
                
                free (diff_vector);
                free (intermediate_vector);
                free (dsig_cen1h_fit);
                
                int *mass_s_count = malloc ((m_num + 1) * sizeof (int));
                
                for (j = 0; j <= m_num; j++) {
                    mass_s_count[j] = 0;
                }
                
                for (j = 0; j <= m_num; j++) {
                    for (k = 0; k <= c_num; k++) {
                        //printf("%d\t%d\t%lE\t%lf\t%lf\t%lf\n", j, k, m_bin[j], c_bin[k], chi_sq_map2[j][k], onesigma_dchi);
                        if (chi_sq_map[j][k] < onesigma_dchi) {
                            mass_s_count[j]++;
                        }
                    }
                }
            
                
                
                
                // Assigning the values of the arrays for the next iteration.
                for (j = 0; j < m_num; j++) {
                    if (mass_s_count[j] != 0 && mass_s_count[j + 1] == 0) {
                        eps3_m = fabs((m_bin[j + 1] - m_bin[j]) / m_bin[j]);
                        m_min = m_bin[j];
                        m_max = m_bin[j + 1];
                        dlogm = (log10(m_max) - log10(m_min)) / m_num;
                        m_up = m_bin[j];
                        break;
                    }
                    if (mass_s_count[j] == 0 && mass_s_count[j + 1] != 0) {
                        eps3_m = fabs((m_bin[j + 1] - m_bin[j]) / m_bin[j]);
                        m_min = m_bin[j];
                        m_max = m_bin[j + 1];
                        dlogm = (log10(m_max) - log10(m_min)) / m_num;
                        m_up = m_bin[j];
                        break;
                    }
                }
                
                printf("Fractional error in 1-sigma mass: %lE\n", eps3_m);
                printf("Upper quartile of 1-sig mass fit: %lE\n", m_up);
                
                /*for (j = 0; j <= m_num; j++) {
                    printf("%d\t%lE\t%d\n", j, m_bin[j], mass_s_count[j]);
                }*/
                
                for (j = 0; j <= m_num; j++) {
                    m_bin[j] = m_min * pow(10, j * dlogm);
                }

                free (mass_s_count);
                
                for (j = 0; j <= m_num; j++) {
                    free (chi_sq_map[j]);
                    free (redchi_sq_map[j]);
                }
                
                free (chi_sq_map);
                free (redchi_sq_map);
            } while (eps3_m > 1e-3);
            
            
            
            p++;
            
            free (m_bin);
            free (c_bin);
            free (r_fit);
            free (dsig_data);
            free (dsig_bar_fit);
            
        } while (p < 1); // do loop for 1-sigma
        
        
        p = 0;
        
        // Lower quartile of 1-sigma mass.
        // Iterate until the M_200 converges.
        
        do {
            
            
            double *r_fit = malloc (r_fitnum[i] * sizeof (double));
            double *dsig_data = malloc (r_fitnum[i] * sizeof (double));
            
            fitnum = 0;
            
            // Input the region where fitting takes place.
            // r_fit in pc h_70^-1, dsig_data in M_sun pc^-2 h_70.
            for (j = 0; j < s_file_count; j++) {
                if (r_phy[j] > r_fitmin && r_phy[j] < r_fitmax[i]) {
                    r_fit[fitnum] = r_phy[j] * Mpc_2_pc;
                    dsig_data[fitnum] = s_info[j].d_sig;
                    fitnum++;
                }
            }
            
            // Double check if the number of data points being fit is correct.
            if (r_fitnum[i] != fitnum) {
                printf("Number of fitting distance bin is not correct! %s \n", lens_name[i].name);
                exit(1);
            }
            
            m_num = 200;
            c_num = 200;
            
            double *m_bin = malloc ((m_num + 1) * sizeof (double));
            double *c_bin = malloc ((c_num + 1) * sizeof (double));
            
            // Chi-squared fitting.
            double *dsig_bar_fit = malloc (fitnum * sizeof (double));
            
            // d_sigma due to baryon, this is not to be fit.
            dsig_baryon(lens_name[i].star_mass_linear, r_fit, dsig_bar_fit, fitnum);
            
            // Log mass bins, for chi-square fitting.
            m_min = 1e10;
            m_max = m_opt3;
            dlogm = (log10(m_max) - log10(m_min)) / m_num;
            
            // Bins for concentration, for chi-square fitting.
            c_min = 0.2;
            c_max = lens_name[i].R_200_Mpc / r_phy[0];  // Make sure the maximum concentration aligns with the innermost bin.

            
            c_diff = (c_max - c_min) / c_num;
            
            for (j = 0; j <= m_num; j++) {
                m_bin[j] = m_min * pow(10, j * dlogm);
            }
            
            for (j = 0; j <= c_num; j++) {
                c_bin[j] = c_min + j * c_diff;
            }
            
            // Making a coarse map to locate the local minima.
            double **chi_sq_map2 = malloc ((m_num + 1) * sizeof (double *));
            double **redchi_sq_map2 = malloc ((m_num + 1) * sizeof (double *));
            
            for (j = 0; j <= m_num; j++) {
                // j = 0 for mass fitting, j = 1 for concentration fitting.
                chi_sq_map2[j] = malloc ((c_num + 1) * sizeof (double));
                redchi_sq_map2[j] = malloc ((c_num + 1) * sizeof (double));
                
                for (k = 0; k <= c_num; k++) {
                    chi_sq_map2[j][k] = 0;
                    redchi_sq_map2[j][k] = 0;
                }
            }
            
            double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
            double *diff_vector = malloc (fitnum * sizeof (double));
            double *intermediate_vector = malloc (fitnum * sizeof (double));
            
            // Find chi^2
            for (j = 0; j <= m_num; j++) { // m_num
                for (k = 0; k <= c_num; k++) { // c_num
                    
                    dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                    
                    fitnum_count = 0;
                    
                    // Vector storing the difference between the data and model.
                    for (s = 0; s < fitnum; s++) {
                        diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                        intermediate_vector[s] = 0;
                    }
                    
                    // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                    for (s = 0; s < fitnum; s++) {
                        for (t = 0; t < fitnum; t++) {
                            intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                        }
                    }
                    
                    for (s = 0; s < fitnum; s++) {
                        // Assuming all central galaxies. Skipping over all the NAN data points.

                        // Chi-square value, chi^2 = c_1*N * a_N*1
                        chi_sq_map2[j][k] += intermediate_vector[s] * diff_vector[s];
                        fitnum_count++;
                    }
                    redchi_sq_map2[j][k] = chi_sq_map2[j][k] / (fitnum_count - 3);
                }
            }
            
            free (diff_vector);
            free (intermediate_vector);
            free (dsig_cen1h_fit);
            
            int *mass_u_count = malloc ((m_num + 1) * sizeof (int));
            
            for (j = 0; j <= m_num; j++) {
                mass_u_count[j] = 0;
            }
            
            for (j = 0; j <= m_num; j++) {
                for (k = 0; k <= c_num; k++) {
                    //printf("%d\t%d\t%lE\t%lf\t%lf\t%lf\n", j, k, m_bin[j], c_bin[k], chi_sq_map2[j][k], onesigma_dchi);
                    if (chi_sq_map2[j][k] < onesigma_dchi) {
                        mass_u_count[j]++;
                    }
                }
            }
            
            for (j = 0; j <= m_num; j++) {
                if (mass_u_count[j] == 0) {
                    break;
                }
            }
            
            
            if (j == m_num + 1) {
                m_low = m_bin[0];
                printf("Upper bound of the lower 1-sigma mass: %lf\n", m_low);
                break;
            }
            
            
            // Assigning the values of the arrays for the next iteration.
            for (j = 1; j < m_num; j++) {
                if (mass_u_count[j] == 0 && mass_u_count[j + 1] != 0) {
                    m_min = m_bin[j];
                    m_max = m_bin[j + 1];
                    break;
                }
                if (mass_u_count[j] != 0 && mass_u_count[j + 1] == 0) {
                    m_min = m_bin[j];
                    m_max = m_bin[j + 1];
                    break;
                }
            }
            
            /*for (j = 0; j <= m_num; j++) {
                printf("%d\t%lE\t%d\n", j, m_bin[j], mass_u_count[j]);
            }*/
            
            dlogm = (log10(m_max) - log10(m_min)) / m_num;
            for (j = 0; j <= m_num; j++) {
                m_bin[j] = m_min * pow(10, j * dlogm);
            }
            
            free (mass_u_count);
            
            
            for (j = 0; j <= m_num; j++) {
                free (chi_sq_map2[j]);
                free (redchi_sq_map2[j]);
            }
            
            free (chi_sq_map2);
            free (redchi_sq_map2);
            
            // Find the minimum chi-squared value by grid-searching method.
            t = 0;
            
            do {
                // Making a coarse map to locate the local minima.
                double **chi_sq_map = malloc ((m_num + 1) * sizeof (double *));
                double **redchi_sq_map = malloc ((m_num + 1) * sizeof (double *));
                
                for (j = 0; j <= m_num; j++) {
                    // j = 0 for mass fitting, j = 1 for concentration fitting.
                    chi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                    redchi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                    
                    for (k = 0; k <= c_num; k++) {
                        chi_sq_map[j][k] = 0;
                        redchi_sq_map[j][k] = 0;
                    }
                }
                
                double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
                double *diff_vector = malloc (fitnum * sizeof (double));
                double *intermediate_vector = malloc (fitnum * sizeof (double));
                
                // Find chi^2
                for (j = 0; j <= m_num; j++) { // m_num
                    for (k = 0; k <= c_num; k++) { // c_num
                        
                        dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                        
                        fitnum_count = 0;
                        
                        // Vector storing the difference between the data and model.
                        for (s = 0; s < fitnum; s++) {
                            diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                            intermediate_vector[s] = 0;
                        }
                        
                        // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                        for (s = 0; s < fitnum; s++) {
                            for (t = 0; t < fitnum; t++) {
                                intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                            }
                        }
                        
                        for (s = 0; s < fitnum; s++) {
                            // Assuming all central galaxies. Skipping over all the NAN data points.

                            // Chi-square value, chi^2 = c_1*N * a_N*1
                            chi_sq_map[j][k] += intermediate_vector[s] * diff_vector[s];
                            fitnum_count++;
                        }
                        
                        redchi_sq_map[j][k] = chi_sq_map[j][k] / (fitnum_count - 3);
                        // First loop over conc., then mass.
                    }
                }
                
                free (diff_vector);
                free (intermediate_vector);
                free (dsig_cen1h_fit);
                
                int *mass_u_count = malloc ((m_num + 1) * sizeof (int));
                
                for (j = 0; j <= m_num; j++) {
                    mass_u_count[j] = 0;
                }
                
                for (j = 0; j <= m_num; j++) {
                    for (k = 0; k <= c_num; k++) {
                        //printf("%d\t%d\t%lE\t%lf\t%lf\t%lf\n", j, k, m_bin[j], c_bin[k], chi_sq_map2[j][k], onesigma_dchi);
                        if (chi_sq_map[j][k] < onesigma_dchi) {
                            mass_u_count[j]++;
                        }
                    }
                }
                
                
                
                
                // Assigning the values of the arrays for the next iteration.
                for (j = 1; j < m_num; j++) {
                    if (mass_u_count[j] == 0 && mass_u_count[j + 1] != 0) {
                        eps3_m = fabs((m_bin[j + 1] - m_bin[j]) / m_bin[j]);
                        m_min = m_bin[j];
                        m_max = m_bin[j + 1];
                        dlogm = (log10(m_max) - log10(m_min)) / m_num;
                        m_low = m_bin[j];
                        break;
                    }
                    if (mass_u_count[j] != 0 && mass_u_count[j + 1] == 0) {
                        eps3_m = fabs((m_bin[j + 1] - m_bin[j]) / m_bin[j]);
                        m_min = m_bin[j];
                        m_max = m_bin[j + 1];
                        dlogm = (log10(m_max) - log10(m_min)) / m_num;
                        m_low = m_bin[j];
                        break;
                    }
                }
                
                printf("Fractional error in 1-sigma mass: %lE\n", eps3_m);
                printf("Lower quartile of 1-sig mass fit: %lE\n", m_low);
                
                
                /*for (j = 0; j <= m_num; j++) {
                    printf("%d\t%lE\t%d\n", j, m_bin[j], mass_u_count[j]);
                }*/
                
                for (j = 0; j <= m_num; j++) {
                    m_bin[j] = m_min * pow(10, j * dlogm);
                }
                
                free (mass_u_count);
                
                for (j = 0; j <= m_num; j++) {
                    free (chi_sq_map[j]);
                    free (redchi_sq_map[j]);
                }
                
                free (chi_sq_map);
                free (redchi_sq_map);
            } while (eps3_m > 1e-3);
            
            
            
            p++;
            
            free (m_bin);
            free (c_bin);
            free (r_fit);
            free (dsig_data);
            free (dsig_bar_fit);
            
        } while (p < 1); // do loop for 1-sigma
        
        
        
        
        p = 0;
        
        // Upper quartile of 1-sigma conc.
        // Iterate until the M_200 converges.
        
        do {
            
            double *r_fit = malloc (r_fitnum[i] * sizeof (double));
            double *dsig_data = malloc (r_fitnum[i] * sizeof (double));
            
            fitnum = 0;
            
            // Input the region where fitting takes place.
            // r_fit in pc h_70^-1, dsig_data in M_sun pc^-2 h_70.
            for (j = 0; j < s_file_count; j++) {
                if (r_phy[j] > r_fitmin && r_phy[j] < r_fitmax[i]) {
                    r_fit[fitnum] = r_phy[j] * Mpc_2_pc;
                    dsig_data[fitnum] = s_info[j].d_sig;
                    fitnum++;
                }
            }
            
            // Double check if the number of data points being fit is correct.
            if (r_fitnum[i] != fitnum) {
                printf("Number of fitting distance bin is not correct! %s \n", lens_name[i].name);
                exit(1);
            }
            
            m_num = 200;
            c_num = 200;
            
            double *m_bin = malloc ((m_num + 1) * sizeof (double));
            double *c_bin = malloc ((c_num + 1) * sizeof (double));
            
            // Chi-squared fitting.
            double *dsig_bar_fit = malloc (fitnum * sizeof (double));
            
            // d_sigma due to baryon, this is not to be fit.
            dsig_baryon(lens_name[i].star_mass_linear, r_fit, dsig_bar_fit, fitnum);
            
            // Log mass bins, for chi-square fitting.
            m_min = 1e10;
            m_max = 1e14;
            dlogm = (log10(m_max) - log10(m_min)) / m_num;
            
            // Bins for concentration, for chi-square fitting.
            c_min = c_opt3;
            c_max = lens_name[i].R_200_Mpc / r_phy[0];  // Make sure the maximum concentration aligns with the innermost bin.
            /*if (p > 0) {
             c_max = r_new / r_phy[0];
             }*/
            
            c_diff = (c_max - c_min) / c_num;
            
            for (j = 0; j <= m_num; j++) {
                m_bin[j] = m_min * pow(10, j * dlogm);
            }
            
            for (j = 0; j <= c_num; j++) {
                c_bin[j] = c_min + j * c_diff;
            }
            
            // Making a coarse map to locate the local minima.
            double **chi_sq_map2 = malloc ((m_num + 1) * sizeof (double *));
            double **redchi_sq_map2 = malloc ((m_num + 1) * sizeof (double *));
            
            for (j = 0; j <= m_num; j++) {
                // j = 0 for mass fitting, j = 1 for concentration fitting.
                chi_sq_map2[j] = malloc ((c_num + 1) * sizeof (double));
                redchi_sq_map2[j] = malloc ((c_num + 1) * sizeof (double));
                
                for (k = 0; k <= c_num; k++) {
                    chi_sq_map2[j][k] = 0;
                    redchi_sq_map2[j][k] = 0;
                }
            }
            
            double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
            double *diff_vector = malloc (fitnum * sizeof (double));
            double *intermediate_vector = malloc (fitnum * sizeof (double));
            
            // Find chi^2
            for (j = 0; j <= m_num; j++) { // m_num
                for (k = 0; k <= c_num; k++) { // c_num
                    
                    dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                    
                    fitnum_count = 0;
                    
                    // Vector storing the difference between the data and model.
                    for (s = 0; s < fitnum; s++) {
                        diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                        intermediate_vector[s] = 0;
                    }
                    
                    // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                    for (s = 0; s < fitnum; s++) {
                        for (t = 0; t < fitnum; t++) {
                            intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                        }
                    }
                    
                    for (s = 0; s < fitnum; s++) {
                        // Assuming all central galaxies. Skipping over all the NAN data points.
                        /*if (dsig_data[s] != dsig_data[s]) {
                         continue;
                         }*/
                        // Chi-square value, chi^2 = c_1*N * a_N*1
                        chi_sq_map2[j][k] += intermediate_vector[s] * diff_vector[s];
                        fitnum_count++;
                    }
                    redchi_sq_map2[j][k] = chi_sq_map2[j][k] / (fitnum_count - 3);
                }
            }
            
            free (diff_vector);
            free (intermediate_vector);
            free (dsig_cen1h_fit);
            
            int *conc_s_count = malloc ((c_num + 1) * sizeof (int));
            
            for (k = 0; k <= c_num; k++) {
                conc_s_count[k] = 0;
            }
            
            for (j = 0; j <= m_num; j++) {
                for (k = 0; k <= c_num; k++) {
                    //printf("%d\t%d\t%lE\t%lf\t%lf\t%lf\n", j, k, m_bin[j], c_bin[k], chi_sq_map2[j][k], onesigma_dchi);
                    if (chi_sq_map2[j][k] < onesigma_dchi) {
                        conc_s_count[k]++;
                    }
                }
            }
            
            for (k = 0; k <= c_num; k++) {
                if (conc_s_count[k] == 0) {
                    break;
                }
            }
            
            
            if (k == c_num + 1) {
                c_up = c_bin[c_num];
                printf("Lower bound of the upper 1-sigma conc: %lf\n", c_up);
                break;
            }
            
            
            // Assigning the values of the arrays for the next iteration.
            for (k = 1; k < c_num; k++) {
                if (conc_s_count[k] != 0 && conc_s_count[k + 1] == 0) {
                    c_min = c_bin[k];
                    c_max = c_bin[k + 1];
                    break;
                }
                if (conc_s_count[k] == 0 && conc_s_count[k + 1] != 0) {
                    c_min = c_bin[k];
                    c_max = c_bin[k + 1];
                    break;
                }
            }
            
            /*for (k = 0; k <= c_num; k++) {
             printf("%d\t%lE\t%d\n", k, c_bin[k], conc_s_count[k]);
             }*/
            
            
            c_diff = (c_max - c_min) / c_num;
            
            for (j = 0; j <= c_num; j++) {
                c_bin[j] = c_min + j * c_diff;
            }
            
            
            free (conc_s_count);
            
            
            for (j = 0; j <= m_num; j++) {
                free (chi_sq_map2[j]);
                free (redchi_sq_map2[j]);
            }
            
            free (chi_sq_map2);
            free (redchi_sq_map2);
            
            // Find the minimum chi-squared value by grid-searching method.
            t = 0;
            
            do {
                // Making a coarse map to locate the local minima.
                double **chi_sq_map = malloc ((m_num + 1) * sizeof (double *));
                double **redchi_sq_map = malloc ((m_num + 1) * sizeof (double *));
                
                for (j = 0; j <= m_num; j++) {
                    // j = 0 for mass fitting, j = 1 for concentration fitting.
                    chi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                    redchi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                    
                    for (k = 0; k <= c_num; k++) {
                        chi_sq_map[j][k] = 0;
                        redchi_sq_map[j][k] = 0;
                    }
                }
                
                double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
                double *diff_vector = malloc (fitnum * sizeof (double));
                double *intermediate_vector = malloc (fitnum * sizeof (double));
                
                // Find chi^2
                for (j = 0; j <= m_num; j++) { // m_num
                    for (k = 0; k <= c_num; k++) { // c_num
                        
                        dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                        
                        fitnum_count = 0;
                        
                        // Vector storing the difference between the data and model.
                        for (s = 0; s < fitnum; s++) {
                            diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                            intermediate_vector[s] = 0;
                        }
                        
                        // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                        for (s = 0; s < fitnum; s++) {
                            for (t = 0; t < fitnum; t++) {
                                intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                            }
                        }
                        
                        for (s = 0; s < fitnum; s++) {
                            // Assuming all central galaxies. Skipping over all the NAN data points.
                            /*if (dsig_data[s] != dsig_data[s]) {
                             continue;
                             }*/
                            // Chi-square value, chi^2 = c_1*N * a_N*1
                            chi_sq_map[j][k] += intermediate_vector[s] * diff_vector[s];
                            fitnum_count++;
                        }
                        
                        redchi_sq_map[j][k] = chi_sq_map[j][k] / (fitnum_count - 3);
                        // First loop over conc., then mass.
                    }
                }
                
                free (diff_vector);
                free (intermediate_vector);
                free (dsig_cen1h_fit);
                
                int *conc_s_count = malloc ((c_num + 1) * sizeof (int));
                
                for (k = 0; k <= c_num; k++) {
                    conc_s_count[k] = 0;
                }
                
                for (j = 0; j <= m_num; j++) {
                    for (k = 0; k <= c_num; k++) {
                        //printf("%d\t%d\t%lE\t%lf\t%lf\t%lf\n", j, k, m_bin[j], c_bin[k], chi_sq_map2[j][k], onesigma_dchi);
                        if (chi_sq_map[j][k] < onesigma_dchi) {
                            conc_s_count[k]++;
                        }
                    }
                }
                
                
                
                
                // Assigning the values of the arrays for the next iteration.
                for (k = 0; k < c_num; k++) {
                    if (conc_s_count[k] != 0 && conc_s_count[k + 1] == 0) {
                        eps3_c = fabs((c_bin[k + 1] - c_bin[k]) / c_bin[k]);
                        c_min = c_bin[k];
                        c_max = c_bin[k + 1];
                        c_diff = (c_max - c_min) / c_num;
                        c_up = c_bin[k];
                        break;
                    }
                }
                
                
                for (j = 0; j <= c_num; j++) {
                    c_bin[j] = c_min + j * c_diff;
                }
                
                printf("Fractional error in 1-sigma conc: %lE\n", eps3_c);
                printf("Upper quartile of 1-sig conc fit: %lE\n", c_up);
                
                /*for (j = 0; j <= m_num; j++) {
                 printf("%d\t%lE\t%d\n", j, m_bin[j], mass_s_count[j]);
                 }*/
                
                
                free (conc_s_count);
                
                for (j = 0; j <= m_num; j++) {
                    free (chi_sq_map[j]);
                    free (redchi_sq_map[j]);
                }
                
                free (chi_sq_map);
                free (redchi_sq_map);
            } while (eps3_c > 1e-3);
            
            
            
            p++;
            
            free (m_bin);
            free (c_bin);
            free (r_fit);
            free (dsig_data);
            free (dsig_bar_fit);
            
        } while (p < 1); // do loop for 1-sigma

        
        
        p = 0;
        
        // Lower quartile of 1-sigma conc.
        // Iterate until the M_200 converges.
        
        do {
            
            double *r_fit = malloc (r_fitnum[i] * sizeof (double));
            double *dsig_data = malloc (r_fitnum[i] * sizeof (double));
            
            fitnum = 0;
            
            // Input the region where fitting takes place.
            // r_fit in pc h_70^-1, dsig_data in M_sun pc^-2 h_70.
            for (j = 0; j < s_file_count; j++) {
                if (r_phy[j] > r_fitmin && r_phy[j] < r_fitmax[i]) {
                    r_fit[fitnum] = r_phy[j] * Mpc_2_pc;
                    dsig_data[fitnum] = s_info[j].d_sig;
                    fitnum++;
                }
            }
            
            // Double check if the number of data points being fit is correct.
            if (r_fitnum[i] != fitnum) {
                printf("Number of fitting distance bin is not correct! %s \n", lens_name[i].name);
                exit(1);
            }
            
            m_num = 200;
            c_num = 200;
            
            double *m_bin = malloc ((m_num + 1) * sizeof (double));
            double *c_bin = malloc ((c_num + 1) * sizeof (double));
            
            // Chi-squared fitting.
            double *dsig_bar_fit = malloc (fitnum * sizeof (double));
            
            // d_sigma due to baryon, this is not to be fit.
            dsig_baryon(lens_name[i].star_mass_linear, r_fit, dsig_bar_fit, fitnum);
            
            // Log mass bins, for chi-square fitting.
            m_min = 1e10;
            m_max = 1e14;
            dlogm = (log10(m_max) - log10(m_min)) / m_num;
            
            // Bins for concentration, for chi-square fitting.
            c_min = 0.2;
            c_max = c_opt3;  // Make sure the maximum concentration aligns with the innermost bin.
            /*if (p > 0) {
             c_max = r_new / r_phy[0];
             }*/
            
            c_diff = (c_max - c_min) / c_num;
            
            for (j = 0; j <= m_num; j++) {
                m_bin[j] = m_min * pow(10, j * dlogm);
            }
            
            for (j = 0; j <= c_num; j++) {
                c_bin[j] = c_min + j * c_diff;
            }
            
            // Making a coarse map to locate the local minima.
            double **chi_sq_map2 = malloc ((m_num + 1) * sizeof (double *));
            double **redchi_sq_map2 = malloc ((m_num + 1) * sizeof (double *));
            
            for (j = 0; j <= m_num; j++) {
                // j = 0 for mass fitting, j = 1 for concentration fitting.
                chi_sq_map2[j] = malloc ((c_num + 1) * sizeof (double));
                redchi_sq_map2[j] = malloc ((c_num + 1) * sizeof (double));
                
                for (k = 0; k <= c_num; k++) {
                    chi_sq_map2[j][k] = 0;
                    redchi_sq_map2[j][k] = 0;
                }
            }
            
            double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
            double *diff_vector = malloc (fitnum * sizeof (double));
            double *intermediate_vector = malloc (fitnum * sizeof (double));
            
            // Find chi^2
            for (j = 0; j <= m_num; j++) { // m_num
                for (k = 0; k <= c_num; k++) { // c_num
                    
                    dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                    
                    fitnum_count = 0;
                    
                    // Vector storing the difference between the data and model.
                    for (s = 0; s < fitnum; s++) {
                        diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                        intermediate_vector[s] = 0;
                    }
                    
                    // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                    for (s = 0; s < fitnum; s++) {
                        for (t = 0; t < fitnum; t++) {
                            intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                        }
                    }
                    
                    for (s = 0; s < fitnum; s++) {
                        // Assuming all central galaxies. Skipping over all the NAN data points.
                        /*if (dsig_data[s] != dsig_data[s]) {
                         continue;
                         }*/
                        // Chi-square value, chi^2 = c_1*N * a_N*1
                        chi_sq_map2[j][k] += intermediate_vector[s] * diff_vector[s];
                        fitnum_count++;
                    }
                    redchi_sq_map2[j][k] = chi_sq_map2[j][k] / (fitnum_count - 3);
                }
            }
            
            free (diff_vector);
            free (intermediate_vector);
            free (dsig_cen1h_fit);
            
            int *conc_u_count = malloc ((c_num + 1) * sizeof (int));
            
            for (k = 0; k <= c_num; k++) {
                conc_u_count[k] = 0;
            }
            
            for (j = 0; j <= m_num; j++) {
                for (k = 0; k <= c_num; k++) {
                    //printf("%d\t%d\t%lE\t%lf\t%lf\t%lf\n", j, k, m_bin[j], c_bin[k], chi_sq_map2[j][k], onesigma_dchi);
                    if (chi_sq_map2[j][k] < onesigma_dchi) {
                        conc_u_count[k]++;
                    }
                }
            }
            
            for (k = 0; k <= c_num; k++) {
                if (conc_u_count[k] == 0) {
                    break;
                }
            }
            
            
            if (k == c_num + 1) {
                c_low = c_bin[0];
                printf("Upper bound of the lower 1-sigma conc: %lf\n", c_low);
                break;
            }
            
            
            // Assigning the values of the arrays for the next iteration.
            for (k = 1; k < c_num; k++) {
                if (conc_u_count[k] != 0 && conc_u_count[k + 1] == 0) {
                    c_min = c_bin[k];
                    c_max = c_bin[k + 1];
                    break;
                }
                if (conc_u_count[k] == 0 && conc_u_count[k + 1] != 0) {
                    c_min = c_bin[k];
                    c_max = c_bin[k + 1];
                    break;
                }
            }
            
            /*for (k = 0; k <= c_num; k++) {
             printf("%d\t%lE\t%d\n", k, c_bin[k], conc_u_count[k]);
             }*/
            
            
            c_diff = (c_max - c_min) / c_num;
            
            for (j = 0; j <= c_num; j++) {
                c_bin[j] = c_min + j * c_diff;
            }
            
            
            free (conc_u_count);
            
            
            for (j = 0; j <= m_num; j++) {
                free (chi_sq_map2[j]);
                free (redchi_sq_map2[j]);
            }
            
            free (chi_sq_map2);
            free (redchi_sq_map2);
            
            // Find the minimum chi-squared value by grid-searching method.
            t = 0;
            
            do {
                // Making a coarse map to locate the local minima.
                double **chi_sq_map = malloc ((m_num + 1) * sizeof (double *));
                double **redchi_sq_map = malloc ((m_num + 1) * sizeof (double *));
                
                for (j = 0; j <= m_num; j++) {
                    // j = 0 for mass fitting, j = 1 for concentration fitting.
                    chi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                    redchi_sq_map[j] = malloc ((c_num + 1) * sizeof (double));
                    
                    for (k = 0; k <= c_num; k++) {
                        chi_sq_map[j][k] = 0;
                        redchi_sq_map[j][k] = 0;
                    }
                }
                
                double *dsig_cen1h_fit = malloc (fitnum * sizeof (double));
                double *diff_vector = malloc (fitnum * sizeof (double));
                double *intermediate_vector = malloc (fitnum * sizeof (double));
                
                // Find chi^2
                for (j = 0; j <= m_num; j++) { // m_num
                    for (k = 0; k <= c_num; k++) { // c_num
                        
                        dsig_cen1h_twopara(m_bin[j], lens_name[i].z, c_bin[k], r_fit, dsig_cen1h_fit, fitnum);
                        
                        fitnum_count = 0;
                        
                        // Vector storing the difference between the data and model.
                        for (s = 0; s < fitnum; s++) {
                            diff_vector[s] = dsig_data[s] - (dsig_cen1h_fit[s] + dsig_bar_fit[s]);
                            intermediate_vector[s] = 0;
                        }
                        
                        // Intermediate vector: c_1*N = a_1*N * b_N*N, where b is the inverse covariance matrix.
                        for (s = 0; s < fitnum; s++) {
                            for (t = 0; t < fitnum; t++) {
                                intermediate_vector[s] += diff_vector[t] * invcov[t][s];
                            }
                        }
                        
                        for (s = 0; s < fitnum; s++) {
                            // Assuming all central galaxies. Skipping over all the NAN data points.
                            /*if (dsig_data[s] != dsig_data[s]) {
                             continue;
                             }*/
                            // Chi-square value, chi^2 = c_1*N * a_N*1
                            chi_sq_map[j][k] += intermediate_vector[s] * diff_vector[s];
                            fitnum_count++;
                        }
                        
                        redchi_sq_map[j][k] = chi_sq_map[j][k] / (fitnum_count - 3);
                        // First loop over conc., then mass.
                    }
                }
                
                free (diff_vector);
                free (intermediate_vector);
                free (dsig_cen1h_fit);
                
                int *conc_u_count = malloc ((c_num + 1) * sizeof (int));
                
                for (k = 0; k <= c_num; k++) {
                    conc_u_count[k] = 0;
                }
                
                for (j = 0; j <= m_num; j++) {
                    for (k = 0; k <= c_num; k++) {
                        //printf("%d\t%d\t%lE\t%lf\t%lf\t%lf\n", j, k, m_bin[j], c_bin[k], chi_sq_map2[j][k], onesigma_dchi);
                        if (chi_sq_map[j][k] < onesigma_dchi) {
                            conc_u_count[k]++;
                        }
                    }
                }
                
                
                
                
                // Assigning the values of the arrays for the next iteration.
                for (k = 0; k < c_num; k++) {
                    if (conc_u_count[k] == 0 && conc_u_count[k + 1] != 0) {
                        eps3_c = fabs((c_bin[k + 1] - c_bin[k]) / c_bin[k]);
                        c_min = c_bin[k];
                        c_max = c_bin[k + 1];
                        c_diff = (c_max - c_min) / c_num;
                        c_low = c_bin[k];
                        break;
                    }
                }
                
                
                printf("Fractional error in 1-sigma conc: %lE\n", eps3_c);
                printf("Lower quartile of 1-sig conc fit: %lf\n", c_low);
                
                /*for (j = 0; j <= c_num; j++) {
                 printf("%d\t%lE\t%d\n", j, c_bin[j], conc_u_count[j]);
                 }*/
                
                for (j = 0; j <= c_num; j++) {
                    c_bin[j] = c_min + j * c_diff;
                }
                
                
                free (conc_u_count);
                
                for (j = 0; j <= m_num; j++) {
                    free (chi_sq_map[j]);
                    free (redchi_sq_map[j]);
                }
                
                free (chi_sq_map);
                free (redchi_sq_map);
            } while (eps3_c > 1e-3);
            
            
            
            p++;
            
            free (m_bin);
            free (c_bin);
            free (r_fit);
            free (dsig_data);
            free (dsig_bar_fit);
            
        } while (p < 1); // do loop for 1-sigma
        
        
        
        

        
        
        strcpy (plot_name, shear_gal_link);
        strcat (plot_name, parameter_name);
        strcat (plot_name, "/plot_");
        strcat (plot_name, lens_name[i].name);
        
        plot_file = fopen (plot_name, "w");
        
        strcpy (mass_name, shear_gal_link);
        strcat (mass_name, parameter_name);
        strcat (mass_name, "/prelim_mass2.txt");
        
        mass_file = fopen (mass_name, "w");
        
        // Make sure the plot file is opened successfully.
        if (!mass_file) {
            printf("Cannot open the prelim_mass file!\n");
            exit(1);
        }
        
        /*r_plot_bin = 100;
        lens_name[i].star_mass_linear = -99;
        lens_name[i].z = 0;
        
        // The radius bin for plotting, in pc.
        double *r_plot_phy2 = malloc (r_plot_bin * sizeof (double));
        
        for (i = 0; i < r_plot_bin; i++) {
            r_plot_phy2[i] = (0.01 + i * 0.01) * 1000000;
        }
        
        m_opt = 7.472196e+11;
        c_opt = 11.1938;*/
        
        // Array to store up the dsigmas for plotting.
        double *dsig_bar_plot = malloc (r_plot_bin * sizeof (double));
        double *dsig_cen1h_plot = malloc (r_plot_bin * sizeof (double));
        double *dsig_cen1h_plot_low = malloc (r_plot_bin * sizeof (double));
        double *dsig_cen1h_plot_up = malloc (r_plot_bin * sizeof (double));
        
        dsig_baryon(lens_name[i].star_mass_linear, r_plot_phy, dsig_bar_plot, r_plot_bin);
        dsig_cen1h_twopara(m_opt, lens_name[i].z, c_opt, r_plot_phy, dsig_cen1h_plot, r_plot_bin);
        dsig_cen1h_twopara(m_low, lens_name[i].z, c_opt, r_plot_phy, dsig_cen1h_plot_low, r_plot_bin);
        dsig_cen1h_twopara(m_up, lens_name[i].z, c_opt, r_plot_phy, dsig_cen1h_plot_up, r_plot_bin);
        
        for (j = 0; j < r_plot_bin; j++) {
            // Correct the r_plot_phy back to kpc for easier plotting.
            fprintf(plot_file, "%E\t%f\t%E\t%E\t%E\t%E\t%E\t%E\t%f\t%f\t%f\t%E\t%E\t%E\n", r_plot_phy[j] / kpc_2_pc, lens_name[i].z, m_opt, m_low, m_up, min_red_chisq, dsig_bar_plot[j], dsig_cen1h_plot[j], c_opt, c_low, c_up, rad_200(rho_crit(lens_name[i].z), m_opt) / kpc_2_pc, dsig_cen1h_plot_low[j], dsig_cen1h_plot_up[j]);
            //fprintf(plot_file, "%E\n", dsig_cen1h_plot[j]);
        }
        
        // free (r_plot_phy2);
        
        fprintf(mass_file, "%s\t%f\t%f\n", lens_name[i].name, lens_name[i].z, rad_200(rho_crit(lens_name[i].z), m_opt) / kpc_2_pc);

        printf("%s\nBest-fit mass\t=\t%E\t%E\t%E\nBest-fit conc.\t=\t%f\t%f\t%f\nRed chi-squared\t=\t%f\nBest-fit R_200\t=\t%f\n", lens_name[i].name, m_low, m_opt, m_up, c_low, c_opt, c_up, min_red_chisq, r_new * kpc_2_pc);
        printf("Percentage err. in m\t=\t%E\nPercentage err. in c\t=\t%E\n", eps2_m, eps2_c);
        
        fclose (plot_file);
        fclose (mass_file);
        fclose (shear_file);
        fclose (invcov_file);
        
        free (s_info);
        free (dsig_bar_plot);
        free (dsig_cen1h_plot);
        free (dsig_cen1h_plot_low);
        free (dsig_cen1h_plot_up);
        free (invcov);
    }
    
    free (r_phy);
    free (r_plot_phy);
    free (r_fitnum);
    free (r_fitmax);
    free (lens_name);
    fclose (shear_cat_file);
    fclose (pre_mass_file);
    
    
    return 0;
}