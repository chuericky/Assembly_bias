/*  This program further screens potential satellites which might reside in groups by finding if there is any galaxies brighter than it by at least 1 mag within 300 kpc in projected plane
    - prependicular to L.O.S.: D_0 = 0.3 Mpc (physical)
    Input: field #: (argv[1]), D_0 (perp to L.O.S.) (argv[2]), fac (argv[3]), stellar mass ratio (argv[4])
 
    Updated: Dec 14, 2015
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
    double ra, dec, ra_low, ra_up, dec_low, dec_up;
} fieldpos;

typedef struct {
    double z, Da, Dc;
    double numden_mod[4];
} numden;

typedef struct {
    char id[15], field[6];
    double x, y, alpha, delta, bg, lv, mu_max, mu_thres, a, b, theta, a_err, b_err, theta_err, class, e1, e2, weight, SNratio, gal_z, gal_z_min, gal_z_max, T_BPZ, modify, c2, lp_mu, lp_mg, lp_mr, lp_mi, lp_mz, lp_med, lp_inf, lp_sup, mu, mu_err, ext_u, mi, mi_err, ext_i, mr, mr_err, ext_r, mg, mg_err, ext_g, my, my_err, ext_y, mz, mz_err, ext_z, size, FWHM_image, FWHM_world, kron_rad, flux_rad, bulge_frac, model_flux, isoarea, app_m_com;
    double PDF[70];
    int flag, fitclass, mask, star_flag, clus_num, clus_doublenum, sat, FOF_layer;
} gallens;

typedef struct {
    char id[15], field[6];
    double alpha, delta, lp_mr, gal_z, mi, my, app_m_com, lp_med;
} galfull;


void writematrix(FILE* infile, gallens *w);
int field_compare(gallens *w);
int field_compare2(galfull *w);

int main(int argc, char *argv[]) {
    int i = 0, j = 0, k = 0, p = 0, q = 0, s = 0, t = 0, u = 0, w = 0, lensindex = 0, spot = 0, file_count = 0, lensfile_num = 0, field_num = 0, fieldindex = 0, field_int = 0, lgal_index = 0, ramesh_num = 0, decmesh_num = 0, r_ind = 0, d_ind = 0, fd_count = 0, fd_cou_diff = 0, fd_array_count = 0, clus_id_num = 0, clus_id_doublenum = 0, fdfd_count = 0, fdfd_count2 = 0, conc_count = 0, fd_cycle = 100, gal_counter = 0, max_clus_id = 0, bad_count = 0, PDF_bin = 70, numden_index = 0, avg_z_ind = 0, partner_count = 0, p_restore = 0, par_counter = 0, max_lum_ind = 0, par_negcount = 0, iso_count = 0, cen_count = 0, fullfile_num = 0, fgal_index = 0, fullindex = 0;
    double ramap_size = 0, decmap_size = 0, ramesh_size = 0, decmesh_size = 0, gl_a_sep = 0, ll_a_sep = 0, D_0 = 0, V_0 = 0, r2d = 180 / M_PI, sm = 0;
    char lensname_file[150], goodname_file[150], othername_file[150], field_file[150], field_manual[3], numdensity_file[150], D_b_man[10], fullname_file[150], fac[6], star_mass[5];
    FILE *lensname, *fieldname, *goodname, *othername, *numdensityname, *fullname;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 5) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    strcpy (fac, argv[3]);
    
    strcpy (star_mass, argv[4]);
    sm = log10(atof (star_mass));
    
    strcpy (numdensity_file, lensmultisortlink);
    strcat (numdensity_file, "redshift_num_density.dat");
    
    numdensityname = fopen (numdensity_file, "r");
    
    // Make sure the redshift file is present.
    if (!numdensityname) {
        printf("Cannot find the redshift file!\n");
        exit(1);
    }
    
    // An array storing the redshift and number density -1/3 power.
    numden *num_den = malloc (21 * sizeof *num_den);
    
    numden_index = 0;
    
    while (fscanf(numdensityname, "%lf %*lf %*lf %*lf %*lf %lf %lf %lf %lf", &num_den[numden_index].z, &num_den[numden_index].numden_mod[0], &num_den[numden_index].numden_mod[1], &num_den[numden_index].numden_mod[2], &num_den[numden_index].numden_mod[3]) == 5) {
        num_den[numden_index].Da = Da(0.3);
        //printf("%d\t%lf\t%lf\n", numden_index, num_den[numden_index].z, num_den[numden_index].Da);
        numden_index++;
    }

    fclose (numdensityname);
    
    // Redshift bins, Da bins and Dc bins.
    double *z_hist = malloc (PDF_bin * sizeof (double));
    double *Da_hist = malloc (PDF_bin * sizeof (double));
    double *Dc_hist = malloc (PDF_bin * sizeof (double));
    double *z_slice = malloc (numden_index * sizeof (double));
    
    for (k = 0; k < numden_index; k++) {
        z_slice[k] = num_den[k].z;
    }
    
    // k-loop is for lens PDF, m-loop is for source PDF.
    for (k = 0; k < PDF_bin; k++) {
        z_hist[k] = 0.025 + k * 0.05;
        Da_hist[k] = Da(z_hist[k]);
        Dc_hist[k] = Dc(z_hist[k]);
    }
    
    // Write in the lens files.
    strcpy (lensname_file, lensisolink);
    strcat (lensname_file, fac);
    strcat (lensname_file, "_lensname.txt");
    
    lensname = fopen (lensname_file, "r");
    
    // Make sure the lens name file is present.
    if (!lensname) {
        printf("Cannot find the lens name file!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count numer of lines in the lens name file.
    lensfile_num = count_line (spot, lensname, file_count);
    
    // An array storing all file names.
    namelens *lens = malloc (lensfile_num * sizeof *lens);
    
    lensindex = 0;
    
    while (fscanf(lensname, "%d %s", &lens[lensindex].num, &lens[lensindex].name) == 2) {
        lensindex++;
    }

    fclose (lensname);
    
    
    // Write in the lens files.
    strcpy (fullname_file, lensisolink);     // lens_with_link
    strcat (fullname_file, fac);
    strcat (fullname_file, "_lensname.txt");
    
    fullname = fopen (fullname_file, "r");
    
    // Make sure the lens name file is present.
    if (!fullname) {
        printf("Cannot find the full catalog name file!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count number of lines in the full name file.
    fullfile_num = count_line (spot, fullname, file_count);
    
    // An array storing all file names.
    namelens *full_cat = malloc (fullfile_num * sizeof *lens);
    
    fullindex = 0;
    
    while (fscanf(fullname, "%d %s", &full_cat[fullindex].num, &full_cat[fullindex].name) == 2) {
        fullindex++;
    }
    
    fclose (fullname);
    

    
    
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

    // For Field W1 to W4.
    strcpy (field_manual, argv[1]);
    field_int = atoi(field_manual);
    
    // Transverse comoving distance in the projected plane (assume all lenses at z = 0.3), in Mpc h^-1. The cutoff in transverse direction and along the l.o.s. (D_0, V_0). Threshold probability.
    strcpy (D_b_man, argv[2]);
    D_0 = atof (D_b_man);


    
    for (i = field_int - 1; i < field_int; i++) {
        printf("\nField W%d\n", i + 1);
        
        FILE *lensfile, *f_file;
        char lensfile_string[150], full_string[150];
        
        strcpy (lensfile_string, lensisolink);   // lens_with_link
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
            // 0 -> i' band, 1 -> y' band
            //if (field_compare(&lgal[lgal_index]) == 0)  lgal[lgal_index].app_m_com = lgal[lgal_index].mi;
            //if (field_compare(&lgal[lgal_index]) == 1)  lgal[lgal_index].app_m_com = lgal[lgal_index].my;
            lgal[lgal_index].sat = 0;
            //printf("%d\t%s\t%lf\t%lf\t%lf\n", lgal_index, lgal[lgal_index].id, lgal[lgal_index].mi, lgal[lgal_index].my, lgal[lgal_index].app_m_com);
            lgal_index++;
            //if (lgal_index == 100) break;
        }
        fclose (lensfile);
        
        
        // Full catalog.
        strcpy (full_string, lensisolink);
        strcat (full_string, full_cat[i].name);
        
        f_file = fopen (full_string, "r");
        
        // Make sure the fuull name file is present.
        if (!f_file) {
            printf("Cannot find the full file %s!\n", full_cat[i].name);
            exit(1);
        }
        
        // Open up the information of the source galaxies.
        galfull *fgal = malloc (full_cat[i].num * sizeof *fgal);
        
        fgal_index = 0;
        
        // Total 135 columns of the source catalog.
        while (fscanf(f_file, "%s %s %*lf %*lf %lf %lf %*lE %*lE %*lf %*lf %*d %*lE %*lE %*lf %*lE %*lE %*lf %*lf %*lf %*lf %*lE %*d %*lf %*d %lf %*lf %*lf %*lf %*lE %*lE %*lE %*lE %lE %*lE %*lE %*d %lf %*lf %*lf %*lf %*lf %*lf %lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %lf %*lf %*lf %*lf %*lf %*lf %*lE %*lf %*lE %*lf %*lf %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE", &fgal[fgal_index].id, &fgal[fgal_index].field, &fgal[fgal_index].alpha, &fgal[fgal_index].delta, &fgal[fgal_index].gal_z, &fgal[fgal_index].lp_mr, &fgal[fgal_index].lp_med, &fgal[fgal_index].mi, &fgal[fgal_index].my) == 9) {
            //if (field_compare2(&fgal[fgal_index]) == 0)  fgal[fgal_index].app_m_com = fgal[fgal_index].mi;
            //if (field_compare2(&fgal[fgal_index]) == 1)  fgal[fgal_index].app_m_com = fgal[fgal_index].my;
            //printf("%d\t%s\t%lf\t%lf\t%lf\t%lf\n", fgal_index, fgal[fgal_index].id, fgal[fgal_index].alpha, fgal[fgal_index].delta, fgal[fgal_index].gal_z, fgal[fgal_index].lp_mr);
            //printf("%d\t%s\t%lf\t%lf\t%lf\n", fgal_index, fgal[fgal_index].id, fgal[fgal_index].mi, fgal[fgal_index].my, fgal[fgal_index].app_m_com);
            
            fgal_index++;
            //if (fgal_index == 100) break;
        }
        fclose (f_file);

        
        // Construct the chaining mesh grids, in arcsecs. Grid size along dec is defined as usual. Grid width along RA is defined according to the dec of the extremes of the fields.
        ramap_size = fabs(fpos[i].ra_up - fpos[i].ra_low) * 3600;
        decmap_size = fabs(fpos[i].dec_up - fpos[i].dec_low) * 3600;
        
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
        
        
        // Create a map to store up which mesh do the lens galaxies belong to, first loop over RA, then DEC. Another array is to store the index of the galaxies in the lens list.
        int **fgal_count = malloc (decmesh_num * sizeof (int *));
        int **fgal_cum = malloc (decmesh_num * sizeof (int *));     // Cumulative number of lens galaxies.
        double ***fgal_coord = malloc (decmesh_num * sizeof (double **));     // Coordinates of the center of the grids.
        int **fgal_meshsort = malloc (fgal_index * sizeof (int *));
        
        for (k = 0; k < fgal_index; k++) {
            fgal_meshsort[k] = malloc (3 * sizeof (int));
            for (p = 0; p < 3; p++) {
                fgal_meshsort[k][p] = 0;    // [k][0] is mesh num in along RA, [k][1] is mesh num in along DEC, [k][2] is index of the galaxy.
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            fgal_count[k] = malloc (ramesh_num * sizeof (int));
            fgal_cum[k] = malloc (ramesh_num * sizeof (int));
            fgal_coord[k] = malloc (ramesh_num * sizeof (double *));
            for (p = 0; p < ramesh_num; p++) {
                fgal_count[k][p] = 0;
                fgal_cum[k][p] = 0;
                fgal_coord[k][p] = malloc (2 * sizeof (double));
                for (q = 0; q < 2; q++) {
                    fgal_coord[k][p][0] = 0;
                    fgal_coord[k][p][1] = 0;
                }
            }
        }
        
        // Store up the number of counts in the mesh grid
        for (k = 0; k < fgal_index; k++) {
            r_ind = min(ra_center, ramesh_num, fgal[k].alpha);
            d_ind = min(dec_center, decmesh_num, fgal[k].delta);
            fgal_meshsort[k][0] = r_ind;
            fgal_meshsort[k][1] = d_ind;
            fgal_meshsort[k][2] = k;
            fgal_count[d_ind][r_ind]++;
        }
        
        // Cumulative number of galaxies in each grid.
        for (k = 0; k < 1; k++) {
            fgal_cum[k][0] = fgal_count[0][0];
            for (p = 1; p < ramesh_num; p++) {
                fgal_cum[k][p] = fgal_cum[k][p - 1] + fgal_count[k][p];
            }
        }
        
        for (k = 1; k < decmesh_num; k++) {
            fgal_cum[k][0] = fgal_cum[k - 1][ramesh_num - 1] + fgal_count[k][0];
            for (p = 1; p < ramesh_num; p++) {
                fgal_cum[k][p] = fgal_cum[k][p - 1] + fgal_count[k][p];
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            for (p = 0; p < ramesh_num; p++) {
                for (q = 0; q < 2; q++) {
                    fgal_coord[k][p][0] = ra_center[p];
                    fgal_coord[k][p][1] = dec_center[k];
                }
            }
        }

        
        // Sort according to grid of RA; if the same, then sort according to grid of DEC; if still the same, then sort by galaxy index in the lens galaxy catalog.
        qsort(fgal_meshsort, fgal_index, sizeof fgal_meshsort[0], compare);
        
        
        // Looping over each seed lens over the grids.
        for (p = 0; p < lgal_index; p++) { // lgal_index
            
            if (p % 10000 == 0) {
                printf("Seed Galaxy, p = %d out of %d\n", p, lens[i].num);
            }

            for (q = 0; q < decmesh_num; q++) { //decmesh_num
                for (s = 0; s < ramesh_num; s++) {  //ramesh_num
                    // Angular separation between the lens and the grid.
                    gl_a_sep = ang_sep (lgal[p].alpha, lgal[p].delta, fgal_coord[q][s][0], fgal_coord[q][s][1]);
                            
                    // See which grid to look at.
                    if (gl_a_sep <= 1.5 * ramesh_size * 3600 && fgal_count[q][s] > 0) {
                        // Create a dynamic array to store the index of the corresponding lens galaxies.
                        int *fgrid_sort = malloc (fgal_count[q][s] * sizeof (int));
                        int *f_ind = malloc (fgal_count[q][s] * sizeof (int));
                                
                        for (t = 0; t < fgal_count[q][s]; t++) {    // fgal_count[q][s]
                                    
                            fgrid_sort[t] = fgal_cum[q][s] - fgal_count[q][s] + t;  // Number of concerned lens galaxy
                            f_ind[t] = fgal_meshsort[fgrid_sort[t]][2];     // The index of the concerned lens galaxy in the lens catalog.
                            
                            ll_a_sep = ang_sep(lgal[p].alpha, lgal[p].delta, fgal[f_ind[t]].alpha, fgal[f_ind[t]].delta);
                                        
                            // Those within transverse physical distance D_0.
                            if (num_den[0].Da * ll_a_sep / r2d / 3600 <= D_0 && num_den[0].Da * ll_a_sep / r2d / 3600 >= 1e-3) {
                                
                                //if (lgal[p].app_m_com - fgal[f_ind[t]].app_m_com >= mag_bright) {
                                if (fgal[f_ind[t]].lp_med - lgal[p].lp_med >= sm) {
                                    
                                //if (lgal[p].lp_mr - fgal[f_ind[t]].lp_mr >= mag_bright) {
                                    //printf("%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", p, f_ind[t], lgal[p].alpha, fgal[f_ind[t]].alpha, lgal[p].delta, fgal[f_ind[t]].delta, num_den[avg_z_ind].Da * ll_a_sep / r2d / 3600, lgal[p].lp_med, fgal[f_ind[t]].lp_med);
                                    lgal[p].sat = 1;
                                }
                                
                                if (lgal[p].sat == 1) break;

                            } // Those lenses within physical distance D_0
                                        
                        } // t loop, for the lenses inside the grid.
                                    
                    } // which grid to look at.
                            
                } // s loop, same dec.
            } // q loop, same ra.
                    
        } // p-loop (lens catalog)
        
        strcpy (goodname_file, lensfurtherlink);   // lens_without_link
        strcat (goodname_file, fac);
        strcat (goodname_file, "_");
        strcat (goodname_file, D_b_man);
        strcat (goodname_file, "_");
        strcat (goodname_file, star_mass);
        strcat (goodname_file, "_W");
        strcat (goodname_file, field_manual);
        strcat (goodname_file, "lens.tsv");
        
        goodname = fopen (goodname_file, "w");
        
        iso_count = 0;
        cen_count = 0;
        
        // See which galaxies are to be determined inside a cluster. And write those 'isolated' and central galaxies down.
        for (k = 0; k < lgal_index; k++) {
            if (lgal[k].sat == 0) {
                iso_count++;
                writematrix(goodname, &lgal[k]);
            }
            if (lgal[k].sat != 0) {
                cen_count++;
                //writematrix(othername, &lgal[k]);
            }
        }
        
        printf("Total # of galaxies = %d\n", lens[i].num);
        printf("# of central galaxies = %d (%.*f %%)\n", iso_count, 3, iso_count * 100. / lens[i].num);
        printf("# of satellite galaxies = %d (%.3f %%)\n", cen_count, 3, cen_count * 100. / lens[i].num);
        
        
        fclose (goodname);
        
        for (k = 0; k < decmesh_num; k++) {
            free(fgal_count[k]);
            free(fgal_cum[k]);
            for (p = 0; p < ramesh_num; p++) {
                free(fgal_coord[k][p]);
            }
            free(fgal_coord[k]);
        }
        
        for (k = 0; k < fgal_index; k++) {
            free (fgal_meshsort[k]);
        }
        
        free (fgal_meshsort);
        free (fgal_cum);
        free (fgal_count);
        free (fgal_coord);
        free (lgal);
        free (fgal);
        free (ra_center);
        free (dec_center);
    }   // i loop, Field loop
    
    free (fpos);
    free (num_den);
    free (z_slice);
    free (lens);
    free (z_hist);
    free (Da_hist);
    free (Dc_hist);
    
    return 0;
}

// The function to write the entire matrix to a txt file.
void writematrix(FILE *infile, gallens *w) {
    fprintf(infile, "%s\t%s\t%lf\t%lf\t%lf\t%lf\t%lE\t%lE\t%lf\t%lf\t%d\t%lE\t%lE\t%lf\t%lE\t%lE\t%lf\t%lf\t%lf\t%lf\t%lE\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lE\t%lf\t%lE\t%lf\t%lf\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", w->id, w->field, w->x, w->y, w->alpha, w->delta, w->bg, w->lv, w->mu_max, w->mu_thres, w->flag, w->a, w->b, w->theta, w->a_err, w->b_err, w->theta_err, w->class, w->e1, w->e2, w->weight, w->fitclass, w->SNratio, w->mask, w->gal_z, w->gal_z_min, w->gal_z_max, w->T_BPZ, w->modify, w->c2, w->lp_mu, w->lp_mg, w->lp_mr, w->lp_mi, w->lp_mz, w->star_flag, w->lp_med, w->lp_inf, w->lp_sup, w->mu, w->mu_err, w->ext_u, w->mi, w->mi_err, w->ext_i, w->mr, w->mr_err, w->ext_r, w->mg, w->mg_err, w->ext_g, w->my, w->my_err, w->ext_y, w->mz, w->mz_err, w->ext_z, w->size, w->FWHM_image, w->FWHM_world, w->kron_rad, w->flux_rad, w->bulge_frac, w->model_flux, w->isoarea, w->PDF[0], w->PDF[1], w->PDF[2], w->PDF[3], w->PDF[4], w->PDF[5], w->PDF[6], w->PDF[7], w->PDF[8], w->PDF[9], w->PDF[10], w->PDF[11], w->PDF[12], w->PDF[13], w->PDF[14], w->PDF[15], w->PDF[16], w->PDF[17], w->PDF[18], w->PDF[19], w->PDF[20], w->PDF[21], w->PDF[22], w->PDF[23], w->PDF[24], w->PDF[25], w->PDF[26], w->PDF[27], w->PDF[28], w->PDF[29], w->PDF[30], w->PDF[31], w->PDF[32], w->PDF[33], w->PDF[34], w->PDF[35], w->PDF[36], w->PDF[37], w->PDF[38], w->PDF[39], w->PDF[40], w->PDF[41], w->PDF[42], w->PDF[43], w->PDF[44], w->PDF[45], w->PDF[46], w->PDF[47], w->PDF[48], w->PDF[49], w->PDF[50], w->PDF[51], w->PDF[52], w->PDF[53], w->PDF[54], w->PDF[55], w->PDF[56], w->PDF[57], w->PDF[58], w->PDF[59], w->PDF[60], w->PDF[61], w->PDF[62], w->PDF[63], w->PDF[64], w->PDF[65], w->PDF[66], w->PDF[67], w->PDF[68], w->PDF[69]);
}

// The function to determine if the field uses the i' (0) or y' (1) band.
int field_compare(gallens *w) {
    if (strcmp("W1m0m4", w->field) == 0 || strcmp("W1m1m4", w->field) == 0 || strcmp("W1m2m4", w->field) == 0 || strcmp("W1m3m4", w->field) == 0 || strcmp("W1m4m4", w->field) == 0 || strcmp("W1p1m4", w->field) == 0 || strcmp("W1p1p1", w->field) == 0 || strcmp("W1p2m4", w->field) == 0 || strcmp("W1p3m4", w->field) == 0 || strcmp("W1p3p1", w->field) == 0 || strcmp("W1p4m4", w->field) == 0 || strcmp("W3m0m1", w->field) == 0 || strcmp("W3m2m1", w->field) == 0 || strcmp("W3m2p1", w->field) == 0 || strcmp("W3p2m3", w->field) == 0 || strcmp("W4m1p1", w->field) == 0 || strcmp("W4m1p2", w->field) == 0 || strcmp("W4m1p3", w->field) == 0 || strcmp("W4m2p2", w->field) == 0 || strcmp("W4m2p3", w->field) == 0 || strcmp("W4m3p3", w->field) == 0)    return 1;
    else    return 0;
}

int field_compare2(galfull *w) {
    if (strcmp("W1m0m4", w->field) == 0 || strcmp("W1m1m4", w->field) == 0 || strcmp("W1m2m4", w->field) == 0 || strcmp("W1m3m4", w->field) == 0 || strcmp("W1m4m4", w->field) == 0 || strcmp("W1p1m4", w->field) == 0 || strcmp("W1p1p1", w->field) == 0 || strcmp("W1p2m4", w->field) == 0 || strcmp("W1p3m4", w->field) == 0 || strcmp("W1p3p1", w->field) == 0 || strcmp("W1p4m4", w->field) == 0 || strcmp("W3m0m1", w->field) == 0 || strcmp("W3m2m1", w->field) == 0 || strcmp("W3m2p1", w->field) == 0 || strcmp("W3p2m3", w->field) == 0 || strcmp("W4m1p1", w->field) == 0 || strcmp("W4m1p2", w->field) == 0 || strcmp("W4m1p3", w->field) == 0 || strcmp("W4m2p2", w->field) == 0 || strcmp("W4m2p3", w->field) == 0 || strcmp("W4m3p3", w->field) == 0)    return 1;
    else    return 0;
}
