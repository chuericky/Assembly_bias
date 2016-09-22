 /* This program selects the corresponding source galaxies for a lens. And calculates the lensing signal by correlation function measurement.
    The searching of source galaxies is done by chaining mesh technique.
    Error bars of the shear profiles are determined by bootstrapping 10^4 times the entire CFHTLens area.
    Only central galaxies are involved.
    // Multiplicative biases are applied individually. <Sigma_crit> for each lens-source pair is calculated as <Sigma_cr^-1>/<Sigma_cr^-2>.
    Input: Field # (argv[1]), multiple of n lens (argv[2]), value of n (argv[3]).

    Updated: Nov 7, 2015
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "tool.h"
#include "cosmopar.h"

// Struct for the galaxy catalog.
typedef struct {
    char name[20];
    int num;
} namelens;

typedef struct {
    char name[15];
    int num;
} namesrcs;

typedef struct {
    double ra, dec, ra_low, ra_up, dec_low, dec_up;
} fieldpos;

typedef struct {
    char id[15], field[6];
    double x, y, alpha, delta, bg, lv, mu_max, mu_thres, a, b, theta, a_err, b_err, theta_err, class, e1, e2, weight, SNratio, gal_z, gal_z_min, gal_z_max, T_BPZ, modify, c2, lp_mu, lp_mg, lp_mr, lp_mi, lp_mz, lp_med, lp_inf, lp_sup, mu, mu_err, ext_u, mi, mi_err, ext_i, mr, mr_err, ext_r, mg, mg_err, ext_g, my, my_err, ext_y, mz, mz_err, ext_z, size, FWHM_image, FWHM_world, kron_rad, flux_rad, bulge_frac, model_flux, isoarea;
    double PDF[70];
    int flag, fitclass, mask, star_flag;
} gallens;

typedef struct {
    char id[15], field[6];
    double x, y, alpha, delta, bg, lv, mu_max, mu_thres, a, b, theta, a_err, b_err, theta_err, class, e1, e2, weight, SNratio, gal_z, gal_z_min, gal_z_max, T_BPZ, modify, c2, lp_mu, lp_mg, lp_mr, lp_mi, lp_mz, lp_med, lp_inf, lp_sup, mu, mu_err, ext_u, mi, mi_err, ext_i, mr, mr_err, ext_r, mg, mg_err, ext_g, my, my_err, ext_y, mz, mz_err, ext_z, size, FWHM_image, FWHM_world, kron_rad, flux_rad, bulge_frac, model_flux, isoarea;
    double PDF[70];
    int flag, fitclass, mask, star_flag;
} galsrcs;




int main(int argc, char *argv[]) {
    char lensname_file[150], srcsname_file[150], field_file[150], field_manual[3], lgal_manum[8], lgal_manual[8], distname_file[150], color[20];// parameter_name[20];
    FILE *lensname, *srcsname, *fieldname, *distname;
    int i = 0, j = 0, k = 0, m = 0, p = 0, q = 0, s = 0, t = 0, u = 0, v = 0, lensfile_num = 0, spot = 0, file_count = -1, lensindex = 0, srcsfile_num = 0, srcsindex = 0, sgal_index = 0, lgal_index = 0, field_num = 0, fieldindex = 0, r_bincount = 0, ramesh_num = 0, decmesh_num = 0, r_ind = 0, d_ind = 0, src_cou = 0, lens_class = 0, lens2_class = 0, low_pc = 0, high_pc = 0, N_ran = 0, boot_count = 0, PDF_bin = 70, field_int = 0, lgal_num = 0, lgal_end = 0, lgal_eva = 0, r_concern = 15, r_count = 0;
    double r_min = 0, r_max = 0, r_intermediate = 0, dlogr = 0, ramap_size = 0, decmap_size = 0, ramesh_size = 0, decmesh_size = 0, decmax = 0, gl_a_sep = 0, ls_a_sep = 0, ls_a_sep_rad = 0, gcir_ang = 0, s2ang_sgal = 0, c2ang_sgal = 0, e_tan = 0, e_cross = 0, zl_ref = 0, zs_ref = 0, eff_ref = 0, Sig_cr_ref = 0, eff_ls = 0, eff_weight = 0, percentile = 0, d2r = M_PI / 180., r2d = 180. / M_PI, sig_norm = 0, gtn = 0, gxn = 0, gd = 0, ann_in = 0, ann_out = 0;

    // Make sure the exact number of arguments are passed.
    if (argc != 4) {
        printf("Invalid passing arguments!\n");
        return 1;
    }

    srand (time(NULL));
    
    // Reference lens and source redshifts, lensing efficiency and Sigma_crit. See Valender et. al. 2011
    zl_ref = 0.27;
    zs_ref = 0.98;
    eff_ref = Da(zl_ref) * Da12(zl_ref, zs_ref) / Da(zs_ref);
    Sig_cr_ref = Sig_cri(zl_ref, zs_ref);
    
    // Redshift bins, Da bins and Da12 bins.
    double *z_hist = malloc (PDF_bin * sizeof (double));
    double *Da_hist = malloc (PDF_bin * sizeof (double));
    double **Da12_hist = malloc (PDF_bin * sizeof (double *));
    
    // k-loop is for lens PDF, m-loop is for source PDF.
    for (k = 0; k < PDF_bin; k++) {
        z_hist[k] = 0.025 + k * 0.05;
        Da_hist[k] = Da(z_hist[k]);
        Da12_hist[k] = malloc (PDF_bin * sizeof (double));
        for (m = 0; m < PDF_bin; m++) {
            Da12_hist[k][m] = 0;
        }
    }
    
    for (k = 0; k < PDF_bin; k++) {
        for (m = k + 1; m < PDF_bin; m++) {
            Da12_hist[k][m] = Da12(z_hist[k], z_hist[m]);
        }
    }
    
    // Parameter to be looked for.
    // strcpy (parameter_name, argv[4]);
    
    // Log radius bins, in arcsecs.
    r_bincount = 20;
    r_min = 5;
    r_intermediate = 2500;
    r_count = 21;
    dlogr = (log10(r_intermediate) - log10(r_min)) / r_bincount;
    
    double *r_width = malloc ((r_count + 1) * sizeof (double));
    double *phyr_cen = malloc ((r_count + 1) * sizeof (double)); // Center of the annulus
    double *phyr_width = malloc ((r_count + 1) * sizeof (double));
    
    // Write in the file to store the location of the bins.
    strcpy (distname_file, shear_multi_link);
    strcat (distname_file, "dist_bin2.dat");
    
    distname = fopen (distname_file, "w");
    
    // "Widths" of the annulus width.
    for (k = 0; k < (r_count + 1); k++) {
        r_width[k] = r_min * pow(10, k * dlogr);
        phyr_width[k] = r_min * pow(10, k * dlogr) * Da(zl_ref) / 3600 * d2r;
        phyr_cen[k] = r_min * pow(10, (k - 0.5) * dlogr) * Da(zl_ref) / 3600 * d2r;   // In Mpc
        //printf("%d\t%lf\t%lf\t%lf\t%lf\n", k, phyr_width[k], phyr_cen[k], r_min * pow(10, (k - 0.5) * dlogr) * Da(0.025) / 3600 * d2r, r_min * pow(10, k * dlogr) * Da(0.025) / 3600 * d2r);
    }
    //getchar();
    for (k = 1; k < r_concern; k++) {
        //printf("%d\t%lE\n", k, phyr_cen[k]);
        fprintf(distname, "%d\t%lE\n", k, phyr_cen[k]);
    }
    
    fclose (distname);
    

    // Write in the source files.
    strcpy (srcsname_file, srcs_multi_link);
    strcat (srcsname_file, "srcsname.txt");
    
    srcsname = fopen (srcsname_file, "r");
    
    // Make sure the source name file is present.
    if (!srcsname) {
        printf("Cannot find the source name file!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count number of lines in the source name file.
    srcsfile_num = count_line (spot, srcsname, file_count);
    
    // An array storing all file names.
    namesrcs *srcs = malloc (srcsfile_num * sizeof *srcs);
    
    srcsindex = 0;
    
    while (fscanf(srcsname, "%d %s", &srcs[srcsindex].num, &srcs[srcsindex].name) == 2) {
        srcsindex++;
    }
    
    fclose (srcsname);
    
    
    // Write in the lens files.
    strcpy (lensname_file, lens_multi_link);   // lensmultisortlink
    //strcat (lensname_file, parameter_name);
    strcat (lensname_file, "lensname.txt");
    
    lensname = fopen (lensname_file, "r");
    
    // Make sure the lens name file is present.
    if (!lensname) {
        printf("Cannot find the lens name file!\n");
        exit(1);
    }
    
    // Count numer of lines in the lens name file.
    lensfile_num = count_line (spot, lensname, file_count);
    
    // An array storing all file names.
    namelens *lens = malloc (lensfile_num * sizeof *lens);
    
    lensindex = 0;
    
    while (fscanf(lensname, "%d %s", &lens[lensindex].num, &lens[lensindex].name) == 2) {
        lensindex++;
    }
    
    spot = 0;
    file_count = -1;
    
    // Write in the position file of the wide fields.
    strcpy (field_file, field_link);
    strcat (field_file, "position.dat");
    
    fieldname = fopen (field_file, "r");
    
    // Make sure the field position file is present.
    if (!fieldname) {
        printf("Cannot find the file for position of the fields!\n");
        exit(1);
    }
    
    // Count number of lines in the field position file.
    field_num = count_line (spot, fieldname, file_count);
    
    // An array storing all file names.
    fieldpos *fpos = malloc (field_num * sizeof *fpos);
    
    fieldindex = 0;
    
    while (fscanf(fieldname, "%*s %lf %lf %lf %lf %lf %lf", &fpos[fieldindex].ra, &fpos[fieldindex].dec, &fpos[fieldindex].ra_low, &fpos[fieldindex].ra_up, &fpos[fieldindex].dec_low, &fpos[fieldindex].dec_up) == 6) {
        fieldindex++;
    }
    
    fclose (fieldname);
    
    N_ran = 500;  // Number of bootstraps.
    
    // Sigma_crit normalization factor, unit in M_sun pc^-1.
    sig_norm = (pow(c, 2) / (4 * M_PI * G_N)) * (pc_2_cm / (Msun_2_g * Mpc_2_pc));
    
    
    // For Field W1 to W4.
    strcpy (field_manual, argv[1]);
    field_int = atoi(field_manual);
    
    for (i = field_int - 1; i < field_int; i++) {
        
        printf("\nField W%d\n", i + 1);
        
        FILE *srcsfile;
        char W_fieldsrcs[2], srcsfile_string[150];
        
        strcpy (srcsfile_string, srcs_multi_link);
        strcat (srcsfile_string, srcs[i].name);
        
        strncpy (W_fieldsrcs, srcs[i].name, 2);
        W_fieldsrcs[2] = 0;
        
        srcsfile = fopen (srcsfile_string, "r");
        
        // Make sure the source name file is present.
        if (!srcsfile) {
            printf("Cannot find the source file %s!\n", srcs[i].name);
            exit(1);
        }
        
        
        // Open up the information of the source galaxies.
        galsrcs *sgal = malloc (srcs[i].num * sizeof *sgal);
        
        sgal_index = 0;
        
        // Total 135 columns of the source catalog.
        while (fscanf(srcsfile, "%s %s %lf %lf %lf %lf %lE %lE %lf %lf %d %lE %lE %lf %lE %lE %lf %lf %lf %lf %lE %d %lf %d %lf %lf %lf %lf %lE %lE %lE %lE %lE %lE %lE %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lE %lf %lE %lf %lf %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE", &sgal[sgal_index].id, &sgal[sgal_index].field, &sgal[sgal_index].x, &sgal[sgal_index].y, &sgal[sgal_index].alpha, &sgal[sgal_index].delta, &sgal[sgal_index].bg, &sgal[sgal_index].lv, &sgal[sgal_index].mu_max, &sgal[sgal_index].mu_thres, &sgal[sgal_index].flag, &sgal[sgal_index].a, &sgal[sgal_index].b, &sgal[sgal_index].theta, &sgal[sgal_index].a_err, &sgal[sgal_index].b_err, &sgal[sgal_index].theta_err, &sgal[sgal_index].class, &sgal[sgal_index].e1, &sgal[sgal_index].e2, &sgal[sgal_index].weight, &sgal[sgal_index].fitclass, &sgal[sgal_index].SNratio, &sgal[sgal_index].mask, &sgal[sgal_index].gal_z, &sgal[sgal_index].gal_z_min, &sgal[sgal_index].gal_z_max, &sgal[sgal_index].T_BPZ, &sgal[sgal_index].modify, &sgal[sgal_index].c2, &sgal[sgal_index].lp_mu, &sgal[sgal_index].lp_mg, &sgal[sgal_index].lp_mr, &sgal[sgal_index].lp_mi, &sgal[sgal_index].lp_mz, &sgal[sgal_index].star_flag, &sgal[sgal_index].lp_med, &sgal[sgal_index].lp_inf, &sgal[sgal_index].lp_sup, &sgal[sgal_index].mu, &sgal[sgal_index].mu_err, &sgal[sgal_index].ext_u, &sgal[sgal_index].mi, &sgal[sgal_index].mi_err, &sgal[sgal_index].ext_i, &sgal[sgal_index].mr, &sgal[sgal_index].mr_err, &sgal[sgal_index].ext_r, &sgal[sgal_index].mg, &sgal[sgal_index].mg_err, &sgal[sgal_index].ext_g, &sgal[sgal_index].my, &sgal[sgal_index].my_err, &sgal[sgal_index].ext_y, &sgal[sgal_index].mz, &sgal[sgal_index].mz_err, &sgal[sgal_index].ext_z, &sgal[sgal_index].size, &sgal[sgal_index].FWHM_image, &sgal[sgal_index].FWHM_world, &sgal[sgal_index].kron_rad, &sgal[sgal_index].flux_rad, &sgal[sgal_index].bulge_frac, &sgal[sgal_index].model_flux, &sgal[sgal_index].isoarea, &sgal[sgal_index].PDF[0], &sgal[sgal_index].PDF[1], &sgal[sgal_index].PDF[2], &sgal[sgal_index].PDF[3], &sgal[sgal_index].PDF[4], &sgal[sgal_index].PDF[5], &sgal[sgal_index].PDF[6], &sgal[sgal_index].PDF[7], &sgal[sgal_index].PDF[8], &sgal[sgal_index].PDF[9], &sgal[sgal_index].PDF[10], &sgal[sgal_index].PDF[11], &sgal[sgal_index].PDF[12], &sgal[sgal_index].PDF[13], &sgal[sgal_index].PDF[14], &sgal[sgal_index].PDF[15], &sgal[sgal_index].PDF[16], &sgal[sgal_index].PDF[17], &sgal[sgal_index].PDF[18], &sgal[sgal_index].PDF[19], &sgal[sgal_index].PDF[20], &sgal[sgal_index].PDF[21], &sgal[sgal_index].PDF[22], &sgal[sgal_index].PDF[23], &sgal[sgal_index].PDF[24], &sgal[sgal_index].PDF[25], &sgal[sgal_index].PDF[26], &sgal[sgal_index].PDF[27], &sgal[sgal_index].PDF[28], &sgal[sgal_index].PDF[29], &sgal[sgal_index].PDF[30], &sgal[sgal_index].PDF[31], &sgal[sgal_index].PDF[32], &sgal[sgal_index].PDF[33], &sgal[sgal_index].PDF[34], &sgal[sgal_index].PDF[35], &sgal[sgal_index].PDF[36], &sgal[sgal_index].PDF[37], &sgal[sgal_index].PDF[38], &sgal[sgal_index].PDF[39], &sgal[sgal_index].PDF[40], &sgal[sgal_index].PDF[41], &sgal[sgal_index].PDF[42], &sgal[sgal_index].PDF[43], &sgal[sgal_index].PDF[44], &sgal[sgal_index].PDF[45], &sgal[sgal_index].PDF[46], &sgal[sgal_index].PDF[47], &sgal[sgal_index].PDF[48], &sgal[sgal_index].PDF[49], &sgal[sgal_index].PDF[50], &sgal[sgal_index].PDF[51], &sgal[sgal_index].PDF[52], &sgal[sgal_index].PDF[53], &sgal[sgal_index].PDF[54], &sgal[sgal_index].PDF[55], &sgal[sgal_index].PDF[56], &sgal[sgal_index].PDF[57], &sgal[sgal_index].PDF[58], &sgal[sgal_index].PDF[59], &sgal[sgal_index].PDF[60], &sgal[sgal_index].PDF[61], &sgal[sgal_index].PDF[62], &sgal[sgal_index].PDF[63], &sgal[sgal_index].PDF[64], &sgal[sgal_index].PDF[65], &sgal[sgal_index].PDF[66], &sgal[sgal_index].PDF[67], &sgal[sgal_index].PDF[68], &sgal[sgal_index].PDF[69]) == 135) {
            sgal_index++;
            
            //if (sgal_index == 10000) break;
        }
        
        // Construct the chaining mesh grids, in arcsecs. Grid size along dec is defined as usual. Grid width along RA is defined according to the dec of the extremes of the fields.
        ramap_size = fabs(fpos[i].ra_up - fpos[i].ra_low) * 3600;
        decmap_size = fabs(fpos[i].dec_up - fpos[i].dec_low) * 3600;
        
        decmax = max_compare(fabs(fpos[i].dec_up), fabs(fpos[i].dec_low));
        
        // Number of meshes along RA and DEC, followed by the size of each grid along RA and DEC (convert to arcsecs). // Sizes in deg.
        ramesh_size = 0.1;
        decmesh_size = 0.1;
        
        ramesh_num = ramap_size / (3600 * ramesh_size);
        decmesh_num = decmap_size / (3600 * decmesh_size);
        
        // A map to store the particles in the corresponding mesh. Loop over RA first, then DEC.
        double *ra_center = malloc (ramesh_num * sizeof (double));
        double *dec_center = malloc (decmesh_num * sizeof (double));
        
        for (k = 0; k < ramesh_num; k++) {
            ra_center[k] = fpos[i].ra_low + (k + 0.5) * ramesh_size;
        }
        for (k = 0; k < decmesh_num; k++) {
            dec_center[k] = fpos[i].dec_low + (k + 0.5) * decmesh_size;
        }
        
        // Create a map to store up which mesh do the source galaxies belong to, first loop over RA, then DEC. Another array is to store the index of the galaxies in the source list.
        int **sgal_count = malloc (decmesh_num * sizeof (int *));
        int **sgal_cum = malloc (decmesh_num * sizeof (int *));     // Cumulative number of source galaxies.
        double ***sgal_coord = malloc (decmesh_num * sizeof (double **));     // Coordinates of the center of the grids.
        int **sgal_meshsort = malloc (sgal_index * sizeof (int *));
        
        for (k = 0; k < sgal_index; k++) {
            sgal_meshsort[k] = malloc (3 * sizeof (int));
            for (p = 0; p < 3; p++) {
                sgal_meshsort[k][p] = 0;    // [k][0] is mesh num in along RA, [k][1] is mesh num in along DEC, [k][2] is index of the galaxy.
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            sgal_count[k] = malloc (ramesh_num * sizeof (int));
            sgal_cum[k] = malloc (ramesh_num * sizeof (int));
            sgal_coord[k] = malloc (ramesh_num * sizeof (double *));
            for (p = 0; p < ramesh_num; p++) {
                sgal_count[k][p] = 0;
                sgal_cum[k][p] = 0;
                sgal_coord[k][p] = malloc (2 * sizeof (double));
                for (q = 0; q < 2; q++) {
                    sgal_coord[k][p][0] = 0;
                    sgal_coord[k][p][1] = 0;
                }
            }
        }
        
        // Store up the number of counts in the mesh grid
        for (k = 0; k < sgal_index; k++) {
            r_ind = min(ra_center, ramesh_num, sgal[k].alpha);
            d_ind = min(dec_center, decmesh_num, sgal[k].delta);
            sgal_meshsort[k][0] = r_ind;
            sgal_meshsort[k][1] = d_ind;
            sgal_meshsort[k][2] = k;
            sgal_count[d_ind][r_ind]++;
        }
        
        // Cumulative number of galaxies in each grid.
        for (k = 0; k < 1; k++) {
            sgal_cum[k][0] = sgal_count[0][0];
            for (p = 1; p < ramesh_num; p++) {
                sgal_cum[k][p] = sgal_cum[k][p - 1] + sgal_count[k][p];
            }
        }
        
        for (k = 1; k < decmesh_num; k++) {
            sgal_cum[k][0] = sgal_cum[k - 1][ramesh_num - 1] + sgal_count[k][0];
            for (p = 1; p < ramesh_num; p++) {
                sgal_cum[k][p] = sgal_cum[k][p - 1] + sgal_count[k][p];
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            for (p = 0; p < ramesh_num; p++) {
                for (q = 0; q < 2; q++) {
                    sgal_coord[k][p][0] = ra_center[p];
                    sgal_coord[k][p][1] = dec_center[k];
                }
            }
        }
        
        // Sort according to grid of RA; if the same, then sort according to grid of DEC; if still the same, then sort by galaxy index in the source galaxy catalog.
        qsort(sgal_meshsort, sgal_index, sizeof sgal_meshsort[0], compare);
        
        FILE *lensfile;
        char lensfile_string[150];
        
        strcpy (lensfile_string, lens_multi_link);   //lensmultisortlink
        //strcat (lensfile_string, parameter_name);
        //strcat (lensfile_string, "/");
        strcat (lensfile_string, lens[i].name);
        
        lensfile = fopen (lensfile_string, "r");
        
        // Make sure the source name file is present.
        if (!lensfile) {
            printf("Cannot find the lens file %s!\n", lens[i].name);
            exit(1);
        }
        
        // Open up the information of the source galaxies.
        gallens *lgal = malloc (lens[i].num * sizeof *lgal);
        
        lgal_index = 0;
        
        // Total 135 columns of the source catalog.
        while (fscanf(lensfile, "%s %s %lf %lf %lf %lf %lE %lE %lf %lf %d %lE %lE %lf %lE %lE %lf %lf %lf %lf %lE %d %lf %d %lf %lf %lf %lf %lE %lE %lE %lE %lE %lE %lE %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lE %lf %lE %lf %lf %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE", &lgal[lgal_index].id, &lgal[lgal_index].field, &lgal[lgal_index].x, &lgal[lgal_index].y, &lgal[lgal_index].alpha, &lgal[lgal_index].delta, &lgal[lgal_index].bg, &lgal[lgal_index].lv, &lgal[lgal_index].mu_max, &lgal[lgal_index].mu_thres, &lgal[lgal_index].flag, &lgal[lgal_index].a, &lgal[lgal_index].b, &lgal[lgal_index].theta, &lgal[lgal_index].a_err, &lgal[lgal_index].b_err, &lgal[lgal_index].theta_err, &lgal[lgal_index].class, &lgal[lgal_index].e1, &lgal[lgal_index].e2, &lgal[lgal_index].weight, &lgal[lgal_index].fitclass, &lgal[lgal_index].SNratio, &lgal[lgal_index].mask, &lgal[lgal_index].gal_z, &lgal[lgal_index].gal_z_min, &lgal[lgal_index].gal_z_max, &lgal[lgal_index].T_BPZ, &lgal[lgal_index].modify, &lgal[lgal_index].c2, &lgal[lgal_index].lp_mu, &lgal[lgal_index].lp_mg, &lgal[lgal_index].lp_mr, &lgal[lgal_index].lp_mi, &lgal[lgal_index].lp_mz, &lgal[lgal_index].star_flag, &lgal[lgal_index].lp_med, &lgal[lgal_index].lp_inf, &lgal[lgal_index].lp_sup, &lgal[lgal_index].mu, &lgal[lgal_index].mu_err, &lgal[lgal_index].ext_u, &lgal[lgal_index].mi, &lgal[lgal_index].mi_err, &lgal[lgal_index].ext_i, &lgal[lgal_index].mr, &lgal[lgal_index].mr_err, &lgal[lgal_index].ext_r, &lgal[lgal_index].mg, &lgal[lgal_index].mg_err, &lgal[lgal_index].ext_g, &lgal[lgal_index].my, &lgal[lgal_index].my_err, &lgal[lgal_index].ext_y, &lgal[lgal_index].mz, &lgal[lgal_index].mz_err, &lgal[lgal_index].ext_z, &lgal[lgal_index].size, &lgal[lgal_index].FWHM_image, &lgal[lgal_index].FWHM_world, &lgal[lgal_index].kron_rad, &lgal[lgal_index].flux_rad, &lgal[lgal_index].bulge_frac, &lgal[lgal_index].model_flux, &lgal[lgal_index].isoarea, &lgal[lgal_index].PDF[0], &lgal[lgal_index].PDF[1], &lgal[lgal_index].PDF[2], &lgal[lgal_index].PDF[3], &lgal[lgal_index].PDF[4], &lgal[lgal_index].PDF[5], &lgal[lgal_index].PDF[6], &lgal[lgal_index].PDF[7], &lgal[lgal_index].PDF[8], &lgal[lgal_index].PDF[9], &lgal[lgal_index].PDF[10], &lgal[lgal_index].PDF[11], &lgal[lgal_index].PDF[12], &lgal[lgal_index].PDF[13], &lgal[lgal_index].PDF[14], &lgal[lgal_index].PDF[15], &lgal[lgal_index].PDF[16], &lgal[lgal_index].PDF[17], &lgal[lgal_index].PDF[18], &lgal[lgal_index].PDF[19], &lgal[lgal_index].PDF[20], &lgal[lgal_index].PDF[21], &lgal[lgal_index].PDF[22], &lgal[lgal_index].PDF[23], &lgal[lgal_index].PDF[24], &lgal[lgal_index].PDF[25], &lgal[lgal_index].PDF[26], &lgal[lgal_index].PDF[27], &lgal[lgal_index].PDF[28], &lgal[lgal_index].PDF[29], &lgal[lgal_index].PDF[30], &lgal[lgal_index].PDF[31], &lgal[lgal_index].PDF[32], &lgal[lgal_index].PDF[33], &lgal[lgal_index].PDF[34], &lgal[lgal_index].PDF[35], &lgal[lgal_index].PDF[36], &lgal[lgal_index].PDF[37], &lgal[lgal_index].PDF[38], &lgal[lgal_index].PDF[39], &lgal[lgal_index].PDF[40], &lgal[lgal_index].PDF[41], &lgal[lgal_index].PDF[42], &lgal[lgal_index].PDF[43], &lgal[lgal_index].PDF[44], &lgal[lgal_index].PDF[45], &lgal[lgal_index].PDF[46], &lgal[lgal_index].PDF[47], &lgal[lgal_index].PDF[48], &lgal[lgal_index].PDF[49], &lgal[lgal_index].PDF[50], &lgal[lgal_index].PDF[51], &lgal[lgal_index].PDF[52], &lgal[lgal_index].PDF[53], &lgal[lgal_index].PDF[54], &lgal[lgal_index].PDF[55], &lgal[lgal_index].PDF[56], &lgal[lgal_index].PDF[57], &lgal[lgal_index].PDF[58], &lgal[lgal_index].PDF[59], &lgal[lgal_index].PDF[60], &lgal[lgal_index].PDF[61], &lgal[lgal_index].PDF[62], &lgal[lgal_index].PDF[63], &lgal[lgal_index].PDF[64], &lgal[lgal_index].PDF[65], &lgal[lgal_index].PDF[66], &lgal[lgal_index].PDF[67], &lgal[lgal_index].PDF[68], &lgal[lgal_index].PDF[69]) == 135) {
            
            lgal_index++;
            //if (lgal_index == 20) break;
        }
        
        // Index of the map. Assign each lens to the nearest grid. And do the lensing analysis on each of the lens.
        
        // Multiple numbers of lens galaxies, starting from 0.
        strcpy (lgal_manual, argv[2]);
        strcpy (lgal_manum, argv[3]);
        lgal_eva = atoi(lgal_manum);
        lgal_num = atoi(lgal_manual) * lgal_eva;
        
        lgal_end = lgal_index;
        
        if ((atoi(lgal_manual) + 1) * lgal_eva < lgal_end) {
            lgal_end = (atoi(lgal_manual) + 1) * lgal_eva;
        }
        
        // For a particular subset of lens
        for (p = lgal_num; p < lgal_end; p++) {
            
            if (p % 10 == 0)   printf("%d done!\n", p);
            
            // Create the arrays to store up the bootstrap sources.
            double **srcs_g_t = malloc (srcs[i].num * sizeof(double *));
            double **srcs_g_d = malloc (srcs[i].num * sizeof(double *));
            for (u = 0; u < srcs[i].num; u++) {
                srcs_g_t[u] = malloc (r_concern * sizeof(double));
                srcs_g_d[u] = malloc (r_concern * sizeof(double));
                for (s = 0; s < r_concern; s++) {
                    srcs_g_t[u][s] = 0;
                    srcs_g_d[u][s] = 0;
                }
            }
            
            double *sum_eff_ls = malloc (r_concern * sizeof (double));
            double *sum_eff_weight = malloc (r_concern * sizeof (double));
            double *gtn_phy = malloc (r_concern * sizeof (double));
            double *gd_phy = malloc (r_concern * sizeof (double));
            double *gtn_phy_acc = malloc (r_concern * sizeof (double));
            double *gd_phy_acc = malloc (r_concern * sizeof (double));
            //double *gxn_phy = malloc ((r_bincount + 1) * sizeof (double));
            
            double *srcs_count = malloc (r_concern * sizeof (double)); // Count the number of sources at each bin.
            
            src_cou = 0;    // A counter to count the total number of sources.
            
            for (u = 0; u < r_concern; u++) {
                gtn_phy_acc[u] = 0;
                gd_phy_acc[u] = 0;
                srcs_count[u] = 0;
            }
            
            // Count
            if ((lgal_end - p) % 100 == 99) printf("%d / %d lens galaxies calculations done!\t#%d\n", (p - lgal_num), (lgal_end - lgal_num), p);
            
            
            // Chaining mesh.
            for (q = 0; q < decmesh_num; q++) { //decmesh_num
                for (s = 0; s < ramesh_num; s++) {  //ramesh_num
                    
                    // Angular separation between the lens and the grid.
                    gl_a_sep = ang_sep (lgal[p].alpha, lgal[p].delta, sgal_coord[q][s][0], sgal_coord[q][s][1]);
                    
                    for (m = 1; m < r_count + 1; m++) {    // 32
                    
                        // Angular radius of the annulus.
                        ann_in = r_width[m - 1];
                        ann_out = r_width[m];
                    
                        // Considering also the edging effect. Make sure there is at least 1 source galaxy in the grid.
                        if ((gl_a_sep - ann_out <= 1.5 * ramesh_size * 3600) && (ann_in - gl_a_sep <= 1.5 * ramesh_size * 3600) && sgal_count[q][s] > 0) {
                        
                            // Create a dynamic array to store the index of the corresponding source galaxies.
                            int *sgrid_sort = malloc (sgal_count[q][s] * sizeof (int));
                            int *s_ind = malloc (sgal_count[q][s] * sizeof (int));
                        
                            for (t = 0; t < sgal_count[q][s]; t++) {    // sgal_count[q][s]
                            
                                sgrid_sort[t] = sgal_cum[q][s] - sgal_count[q][s] + t;  // Number of concerned source galaxy
                                s_ind[t] = sgal_meshsort[sgrid_sort[t]][2];     // The index of the concerned source galaxy in the source catalog.
                                
                            
                                // Make sure the code is reading in the source galaxies in the desired grid.
                                if (sgal_meshsort[sgrid_sort[t]][0] != s || sgal_meshsort[sgrid_sort[t]][1] != q) {
                                    printf("Error! A source galaxy is not in the concerned grid.\n");
                                    exit (1);
                                }
                            
                            
                                // Source galaxy selection criteria. (Redshift requirements)
                                if ((sgal[s_ind[t]].gal_z_min > lgal[p].gal_z) && (sgal[s_ind[t]].gal_z - lgal[p].gal_z >= 0.1)) {
                                
                                    // Lens-source separation.
                                    ls_a_sep = ang_sep (lgal[p].alpha, lgal[p].delta, sgal[s_ind[t]].alpha, sgal[s_ind[t]].delta);
                                    
                                    // For sources inside the annulus.
                                    if (ls_a_sep <= ann_out && ls_a_sep >= ann_in) {
                                        src_cou++;  // A counter to count the total number of sources in angular bins.
                                    
                                        ls_a_sep_rad = ls_a_sep * d2r / 3600;
                                    
                                        // These quantities are defined in physical r unit.
                                        for (u = 0; u < r_concern; u++) {
                                            sum_eff_ls[u] = 0;
                                            sum_eff_weight[u] = 0;
                                            gtn_phy[u] = 0;
                                            gd_phy[u] = 0;
                                            //gxn_phy[u] = 0;
                                        }
                                    
                                        // Lensing efficiency for the lens-source pair, followed by the lensing efficiency weight factor (1/Sig_cr^2).
                                        //  sum_eff_weight = (1/Sig_cr)
                                        // r_concern up to the 16th data point.
                                        exp_invsig2(Da_hist, Da12_hist, lgal[p].PDF, sgal[s_ind[t]].PDF, PDF_bin, phyr_width, r_concern, ls_a_sep_rad, sum_eff_ls);
                                        exp_invsigsq2(Da_hist, Da12_hist, lgal[p].PDF, sgal[s_ind[t]].PDF, PDF_bin, phyr_width, r_concern, ls_a_sep_rad, sum_eff_weight);
                                    
                                        // Angle between e1 and the great circle.
                                        gcir_ang = ori_ang(lgal[p].alpha, lgal[p].delta, sgal[s_ind[t]].alpha, sgal[s_ind[t]].delta, ls_a_sep);
                                    
                                        // Sine and cosine of the rotation matrix.
                                        s2ang_sgal = sin(2 * gcir_ang);
                                        c2ang_sgal = cos(2 * gcir_ang);
                                    
                                        // Ellipticity in tangential and cross components. e2 should be subtracted by c2 first, then flipped by -1.
                                        e_tan = -c2ang_sgal * sgal[s_ind[t]].e1 + s2ang_sgal * (sgal[s_ind[t]].e2 - sgal[s_ind[t]].c2);
                                        //e_cross = s2ang_sgal * sgal[s_ind[t]].e1 + c2ang_sgal * (sgal[s_ind[t]].e2 - sgal[s_ind[t]].c2);
                                        
                                        numsrcs(Da_hist, lgal[p].PDF, PDF_bin, phyr_width, r_concern, ls_a_sep_rad, srcs_count);
            
                                        for (u = 1; u < r_concern; u++) {
                                            gtn_phy[u] = sgal[s_ind[t]].weight * e_tan * sum_eff_ls[u] / (1 + sgal[s_ind[t]].modify);
                                            gd_phy[u] = sgal[s_ind[t]].weight * sum_eff_weight[u];
                                            gtn_phy_acc[u] += gtn_phy[u];
                                            gd_phy_acc[u] += gd_phy[u];
                                            srcs_g_t[src_cou][u] = gtn_phy[u];
                                            srcs_g_d[src_cou][u] = gd_phy[u];
                                            //printf("%d\t%d\t%lE\t%lE\t%lE\t%lE\n", src_cou, u, gtn_phy[u], gd_phy[u], gtn_phy_acc[u], gd_phy_acc[u]);
                                            //gxn_phy[u] = sgal[s_ind[t]].weight * e_cross * sum_eff_ls[u] / (1 + sgal[s_ind[t]].modify);
                                        }
                                    } // Loop for selecting the galaxies inside the annulus.
                                } // Loop for the redshift separation of the source and the lens.
                            }     // t loop, for source galaxies inside the grid.
                                     
                            free (sgrid_sort);
                            free (s_ind);
                        }         // condition if a grid if inside the annulus
                    }             // m loop, annulus bins
                }                 // s loop, dec of the grids
            }                     // q loop, ra of the grids
            
            
            /************* Bootstrap ***************/
            
            int **bootstrap = malloc (N_ran * sizeof (int *));
            double **g_t_boot = malloc (N_ran * sizeof (double *));
            double **g_d_boot = malloc (N_ran * sizeof (double *));
            double **g_boot = malloc (N_ran * sizeof (double *));
            double *g_boot_sum = malloc (r_concern * sizeof (double));
            double *g_boot_mean = malloc (r_concern * sizeof (double));
            double *g_boot_diff_sq = malloc (r_concern * sizeof (double));
            double *g_boot_bin_err = malloc (r_concern * sizeof (double));
            
            for (u = 0; u < r_concern; u++) {
                g_boot_sum[u] = 0;
                g_boot_mean[u] = 0;
                g_boot_diff_sq[u] = 0;
                g_boot_bin_err[u] = 0;
            }
            
            for (q = 0; q < N_ran; q++) {
                bootstrap[q] = malloc (src_cou * sizeof (int));
                g_t_boot[q] = malloc (r_concern * sizeof (double));
                g_d_boot[q] = malloc (r_concern * sizeof (double));
                g_boot[q] = malloc (r_concern * sizeof (double));
                
                for (u = 0; u < r_concern; u++) {
                    g_t_boot[q][u] = 0;
                    g_d_boot[q][u] = 0;
                    g_boot[q][u] = 0;
                }
                
                for (s = 0; s < src_cou; s++) {
                    // Random generator for bootstrapping.
                    bootstrap[q][s] = (int)(rand() / (RAND_MAX / (src_cou + 0.)));
                    
                    for (u = 1; u < r_concern; u++) {
                        g_t_boot[q][u] += srcs_g_t[bootstrap[q][s]][u];
                        g_d_boot[q][u] += srcs_g_d[bootstrap[q][s]][u];
                    }
                }
                
                for (u = 1; u < r_concern; u++) {
                    g_boot[q][u] = g_t_boot[q][u] / g_d_boot[q][u];
                    // Nan value for g_boot[q][u]
                    if (g_boot[q][u] != g_boot[q][u]) {
                        g_boot[q][u] = 0;
                    }
                    g_boot_sum[u] += g_boot[q][u];     // Sum of bootstraps at that r bin.
                }
            }
            
            // Calculate the mean of each bin.
            for (u = 1; u < r_concern; u++) {
                g_boot_mean[u] = g_boot_sum[u] / N_ran;
            }
            
            // Diagonal elements of the covariance matrix.
            for (q = 0; q < N_ran; q++) {
                for (u = 1; u < r_concern; u++) {
                    g_boot_diff_sq[u] += pow((g_boot[q][u] - g_boot_mean[u]), 2);
                }
            }
            
            char lensshear_file[150];
            FILE *lensshear_name;
            
            strcpy (lensshear_file, shear_multi_link);
            // strcat (lensshear_file, parameter_name);
            strcat (lensshear_file, "gal_3/");
            strcat (lensshear_file, lgal[p].id);
            strcat (lensshear_file, "_wider.dat");
            
            lensshear_name = fopen (lensshear_file, "w");
            
            fprintf(lensshear_name, "%s\t", lgal[p].id);
            
            // Error bars of each bin, unbiased variance square root.
            for (u = 1; u < r_concern; u++) {
                g_boot_bin_err[u] = pow((g_boot_diff_sq[u] / (N_ran - 1)), 0.5);
                fprintf(lensshear_name, "%lE\t", sig_norm * gtn_phy_acc[u] / gd_phy_acc[u]);
            }
            
            for (u = 1; u < r_concern; u++) {
                fprintf(lensshear_name, "%lE\t", sig_norm * g_boot_bin_err[u]);
            }
            
            for (u = 1; u < r_concern; u++) {
                //fprintf(lensshear_name, "%lE\t", gd_phy_acc[u]);
                fprintf(lensshear_name, "%lE\t", srcs_count[u]);
            }
            
            fprintf(lensshear_name, "\n");
            
            fclose (lensshear_name);
            
            
            for (q = 0; q < N_ran; q++) {
                free (bootstrap[q]);
                free (g_t_boot[q]);
                free (g_d_boot[q]);
                free (g_boot[q]);
            }
            
            free (bootstrap);
            free (g_t_boot);
            free (g_d_boot);
            free (g_boot);
            free (g_boot_sum);
            free (g_boot_mean);
            free (g_boot_diff_sq);
            free (g_boot_bin_err);
            
            /*
            char lensshear_file[150];
            FILE *lensshear_name;
            
            strcpy (lensshear_file, shear_multi2_link);
            strcat (lensshear_file, lgal[p].id);
            strcat (lensshear_file, ".dat");
            
            lensshear_name = fopen (lensshear_file, "w");
            
            fprintf(lensshear_name, "%s\t", lgal[p].id);
            
            // Error bars of each bin, unbiased variance square root.
            for (u = 1; u < r_concern; u++) {
                fprintf(lensshear_name, "%lE\t", sig_norm * gtn_phy_acc[u]);
            }
            
            for (u = 1; u < r_concern; u++) {
                fprintf(lensshear_name, "%lE\t", gd_phy_acc[u]);
            }
            fprintf(lensshear_name, "\n");
            
            fclose (lensshear_name);*/

            
            for (u = 0; u < srcs[i].num; u++) {
                free(srcs_g_t[u]);
                free(srcs_g_d[u]);
            }
            
            free (srcs_g_t);
            free (srcs_g_d);
            free (sum_eff_ls);
            free (sum_eff_weight);
            free (gtn_phy);
            free (gd_phy);
            free (gtn_phy_acc);
            free (gd_phy_acc);
            free (srcs_count);

        }                         // p loop, lens galaxies
        
        for (k = 0; k < decmesh_num; k++) {
            free(sgal_count[k]);
            free(sgal_cum[k]);
            for (p = 0; p < ramesh_num; p++) {
                free(sgal_coord[k][p]);
            }
            free(sgal_coord[k]);
        }
        
        for (k = 0; k < sgal_index; k++) {
            free (sgal_meshsort[k]);
        }
    
        free (sgal_meshsort);
        free (sgal_cum);
        free (sgal_count);
        free (sgal_coord);
        free (ra_center);
        free (dec_center);
        free (sgal);
        free (lgal);
        
        fclose (srcsfile);
        fclose (lensfile);
        
    }
    
    free(r_width);
    free(phyr_width);
    free(lens);
    free(srcs);
    free(fpos);
    free(Da12_hist);
    free(Da_hist);
    free(z_hist);
    
    return 0;
}

