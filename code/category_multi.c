/* This program selects the lens galaxies according to the following criteria. Then it selects source galaxies. Lens and source files are stored in binary format.
    a. 0.2 <= z_photo <= 0.4.
    e. MAG_i / y < 23.
 
    Some fields use i' band, some use y' band.
 
    It normalizes the redshift pdf of each galaxy as well.
 
    Updated: Jul 1, 2015
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tool.h"

// Struct for the galaxy catalog.
typedef struct {
    char id[15], field[6];
    double x, y, alpha, delta, bg, lv, mu_max, mu_thres, a, b, theta, a_err, b_err, theta_err, class, e1, e2, weight, SNratio, gal_z, gal_z_min, gal_z_max, T_BPZ, modify, c2, lp_mu, lp_mg, lp_mr, lp_mi, lp_mz, lp_med, lp_inf, lp_sup, mu, mu_err, ext_u, mi, mi_err, ext_i, mr, mr_err, ext_r, mg, mg_err, ext_g, my, my_err, ext_y, mz, mz_err, ext_z, size, FWHM_image, FWHM_world, kron_rad, flux_rad, bulge_frac, model_flux, isoarea;
    double PDF[70];
    int flag, fitclass, mask, star_flag;
} gal;

void writematrix(FILE* infile, gal *w);
int field_compare(gal *w);


int main() {
    int i = 0, j = 0, gal_num = 0, spot = 0, gal_count = -1, gal_index = 0, ry_count = 0, ri_count = 0, by_count = 0, bi_count = 0, lenscount = 0, totcount = 0, totcount2 = 0, sri_count = 0, sry_count = 0, sby_count = 0, sbi_count = 0, si_count = 0, sy_count = 0, totcount3 = 0, s2ri_count = 0, s2ry_count = 0, PDF_bin = 70;
    FILE *W1lens, *W2lens, *W3lens, *W4lens, *W1src, *W2src, *W3src, *W4src;
    char W1lens_file[150], W2lens_file[150], W3lens_file[150], W4lens_file[150], W1src_file[150], W2src_file[150], W3src_file[150], W4src_file[150];
    double PDF_norm = 0;

    // Write the catalog of the lens and sources for W1 to W4. For lens galaxies, red and blue galaxy catalogs are separated.
    strcpy (W1lens_file, lens_multi_link);
    strcpy (W2lens_file, lens_multi_link);
    strcpy (W3lens_file, lens_multi_link);
    strcpy (W4lens_file, lens_multi_link);
    strcpy (W1src_file, srcs_multi_link);
    strcpy (W2src_file, srcs_multi_link);
    strcpy (W3src_file, srcs_multi_link);
    strcpy (W4src_file, srcs_multi_link);
    
    strcat (W1lens_file, "W1lens2.tsv");
    strcat (W2lens_file, "W2lens2.tsv");
    strcat (W3lens_file, "W3lens2.tsv");
    strcat (W4lens_file, "W4lens2.tsv");
    strcat (W1src_file, "W1src2.tsv");
    strcat (W2src_file, "W2src2.tsv");
    strcat (W3src_file, "W3src2.tsv");
    strcat (W4src_file, "W4src2.tsv");
    
    W1lens = fopen (W1lens_file, "w");
    W2lens = fopen (W2lens_file, "w");
    W3lens = fopen (W3lens_file, "w");
    W4lens = fopen (W4lens_file, "w");
    W1src = fopen (W1src_file, "w");
    W2src = fopen (W2src_file, "w");
    W3src = fopen (W3src_file, "w");
    W4src = fopen (W4src_file, "w");
    
    // An array to store up all the normalized PDF's, and redshifts bins.
    double *norm_PDF = malloc (PDF_bin * sizeof (double));
    
    for (j = 0; j < PDF_bin; j++) {
        norm_PDF[j] = 0;
    }
    
    for (i = 0; i < 203; i++) { // 203
        FILE* catalog_all;
        char catadir[150], file_num[4], W_field[2];
        
        printf("%d\n", i);
        
        // Read in the entire catalog.
        strcpy (catadir, catalog_multi_link);
        strcat (catadir, "CFHTLens_");
        sprintf (file_num, "%d", i);
        strcat (catadir, file_num);
        strcat (catadir, ".tsv");
        
        catalog_all = fopen (catadir, "r");
        
        // Make sure the catalog is available.
        if (!catalog_all) {
            printf("Cannot find the CFHTLenS catalog!\n");
            exit(1);
        }
        
        // Number of lines in the catalog.
        gal_num = count_line(spot, catalog_all, gal_count);
       
        // A struct to store all information of the galaxies.
        gal *galx = malloc (gal_num * sizeof *galx);
        
        gal_index = 0;
        
        // Total 135 columns of the catalog.
        while (fscanf(catalog_all, "%s %s %lf %lf %lf %lf %lE %lE %lf %lf %d %lE %lE %lf %lE %lE %lf %lf %lf %lf %lE %d %lf %d %lf %lf %lf %lf %lE %lE %lE %lE %lE %lE %lE %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lE %lf %lE %lf %lf %lE %lE %lE %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE, %lE", &galx[gal_index].id, &galx[gal_index].field, &galx[gal_index].x, &galx[gal_index].y, &galx[gal_index].alpha, &galx[gal_index].delta, &galx[gal_index].bg, &galx[gal_index].lv, &galx[gal_index].mu_max, &galx[gal_index].mu_thres, &galx[gal_index].flag, &galx[gal_index].a, &galx[gal_index].b, &galx[gal_index].theta, &galx[gal_index].a_err, &galx[gal_index].b_err, &galx[gal_index].theta_err, &galx[gal_index].class, &galx[gal_index].e1, &galx[gal_index].e2, &galx[gal_index].weight, &galx[gal_index].fitclass, &galx[gal_index].SNratio, &galx[gal_index].mask, &galx[gal_index].gal_z, &galx[gal_index].gal_z_min, &galx[gal_index].gal_z_max, &galx[gal_index].T_BPZ, &galx[gal_index].modify, &galx[gal_index].c2, &galx[gal_index].lp_mu, &galx[gal_index].lp_mg, &galx[gal_index].lp_mr, &galx[gal_index].lp_mi, &galx[gal_index].lp_mz, &galx[gal_index].star_flag, &galx[gal_index].lp_med, &galx[gal_index].lp_inf, &galx[gal_index].lp_sup, &galx[gal_index].mu, &galx[gal_index].mu_err, &galx[gal_index].ext_u, &galx[gal_index].mi, &galx[gal_index].mi_err, &galx[gal_index].ext_i, &galx[gal_index].mr, &galx[gal_index].mr_err, &galx[gal_index].ext_r, &galx[gal_index].mg, &galx[gal_index].mg_err, &galx[gal_index].ext_g, &galx[gal_index].my, &galx[gal_index].my_err, &galx[gal_index].ext_y, &galx[gal_index].mz, &galx[gal_index].mz_err, &galx[gal_index].ext_z, &galx[gal_index].size, &galx[gal_index].FWHM_image, &galx[gal_index].FWHM_world, &galx[gal_index].kron_rad, &galx[gal_index].flux_rad, &galx[gal_index].bulge_frac, &galx[gal_index].model_flux, &galx[gal_index].isoarea, &galx[gal_index].PDF[0], &galx[gal_index].PDF[1], &galx[gal_index].PDF[2], &galx[gal_index].PDF[3], &galx[gal_index].PDF[4], &galx[gal_index].PDF[5], &galx[gal_index].PDF[6], &galx[gal_index].PDF[7], &galx[gal_index].PDF[8], &galx[gal_index].PDF[9], &galx[gal_index].PDF[10], &galx[gal_index].PDF[11], &galx[gal_index].PDF[12], &galx[gal_index].PDF[13], &galx[gal_index].PDF[14], &galx[gal_index].PDF[15], &galx[gal_index].PDF[16], &galx[gal_index].PDF[17], &galx[gal_index].PDF[18], &galx[gal_index].PDF[19], &galx[gal_index].PDF[20], &galx[gal_index].PDF[21], &galx[gal_index].PDF[22], &galx[gal_index].PDF[23], &galx[gal_index].PDF[24], &galx[gal_index].PDF[25], &galx[gal_index].PDF[26], &galx[gal_index].PDF[27], &galx[gal_index].PDF[28], &galx[gal_index].PDF[29], &galx[gal_index].PDF[30], &galx[gal_index].PDF[31], &galx[gal_index].PDF[32], &galx[gal_index].PDF[33], &galx[gal_index].PDF[34], &galx[gal_index].PDF[35], &galx[gal_index].PDF[36], &galx[gal_index].PDF[37], &galx[gal_index].PDF[38], &galx[gal_index].PDF[39], &galx[gal_index].PDF[40], &galx[gal_index].PDF[41], &galx[gal_index].PDF[42], &galx[gal_index].PDF[43], &galx[gal_index].PDF[44], &galx[gal_index].PDF[45], &galx[gal_index].PDF[46], &galx[gal_index].PDF[47], &galx[gal_index].PDF[48], &galx[gal_index].PDF[49], &galx[gal_index].PDF[50], &galx[gal_index].PDF[51], &galx[gal_index].PDF[52], &galx[gal_index].PDF[53], &galx[gal_index].PDF[54], &galx[gal_index].PDF[55], &galx[gal_index].PDF[56], &galx[gal_index].PDF[57], &galx[gal_index].PDF[58], &galx[gal_index].PDF[59], &galx[gal_index].PDF[60], &galx[gal_index].PDF[61], &galx[gal_index].PDF[62], &galx[gal_index].PDF[63], &galx[gal_index].PDF[64], &galx[gal_index].PDF[65], &galx[gal_index].PDF[66], &galx[gal_index].PDF[67], &galx[gal_index].PDF[68], &galx[gal_index].PDF[69]) == 135) {
            
            PDF_norm = 0;
            
            for (j = 0; j < PDF_bin; j++) {
                PDF_norm += galx[gal_index].PDF[j];
            }
            
            for (j = 0; j < PDF_bin; j++) {
                norm_PDF[j] = galx[gal_index].PDF[j] / PDF_norm;
                galx[gal_index].PDF[j] = norm_PDF[j];
            }

            // Determine which field does the galaxy belong to.
            strncpy(W_field, galx[gal_index].field, 2);
            W_field[2] = 0;
            
            // Lens selection, using photo-z is fine.
            //if (galx[gal_index].gal_z >= 0.2 && galx[gal_index].gal_z <= 0.4 && ((galx[gal_index].T_BPZ >= 1.0 && galx[gal_index].T_BPZ <= 1.5 && galx[gal_index].lp_med >= 9.0 && galx[gal_index].lp_med <= 11.0) || (galx[gal_index].T_BPZ >= 2.0 && galx[gal_index].T_BPZ <= 4.0 && galx[gal_index].lp_med >= 9.0 && galx[gal_index].lp_med <= 11.0))) {
            //if (galx[gal_index].gal_z >= 0.2 && galx[gal_index].gal_z <= 0.4 && ((galx[gal_index].T_BPZ >= 1.0 && galx[gal_index].T_BPZ <= 1.5 && galx[gal_index].lp_mr >= -21.5 && galx[gal_index].lp_mr <= -19.0) || (galx[gal_index].T_BPZ >= 2.0 && galx[gal_index].T_BPZ <= 4.0 && galx[gal_index].lp_mr >= -21.5 && galx[gal_index].lp_mr <= -19.0))) {
            if (galx[gal_index].gal_z >= 0.2 && galx[gal_index].gal_z <= 0.4) { // && ((galx[gal_index].T_BPZ >= 2.0 && galx[gal_index].T_BPZ <= 4.0 && galx[gal_index].lp_med >= 9.0 && galx[gal_index].lp_med <= 9.25 && galx[gal_index].lp_mr >= -21.0))) {
                
                totcount++;
                lenscount++;
                // Boolean test. 1 -> y' band, 0 -> i' band.
                if (field_compare(&galx[gal_index]) == 0 && galx[gal_index].mi <= 23.0 && galx[gal_index].mi > -98.0 && galx[gal_index].lp_mu > -98.0 && galx[gal_index].lp_mg > -98.0 && galx[gal_index].lp_mr > -98.0 && galx[gal_index].lp_mi > -98.0 && galx[gal_index].lp_mz > -98.0 && galx[gal_index].lp_med > -98.0 && galx[gal_index].lp_inf > -98.0 && galx[gal_index].lp_sup > -98.0) {
                    ri_count++;
                    if (strcmp("W1", W_field) == 0)   writematrix(W1lens, &galx[gal_index]);
                    if (strcmp("W2", W_field) == 0)   writematrix(W2lens, &galx[gal_index]);
                    if (strcmp("W3", W_field) == 0)   writematrix(W3lens, &galx[gal_index]);
                    if (strcmp("W4", W_field) == 0)   writematrix(W4lens, &galx[gal_index]);
                }
                if (field_compare(&galx[gal_index]) == 1 && galx[gal_index].my <= 23.0 && galx[gal_index].my > -98.0 && galx[gal_index].lp_mu > -98.0 && galx[gal_index].lp_mg > -98.0 && galx[gal_index].lp_mr > -98.0 && galx[gal_index].lp_mi > -98.0 && galx[gal_index].lp_mz > -98.0 && galx[gal_index].lp_med > -98.0 && galx[gal_index].lp_inf > -98.0 && galx[gal_index].lp_sup > -98.0) {
                    ry_count++;
                    if (strcmp("W1", W_field) == 0)   writematrix(W1lens, &galx[gal_index]);
                    if (strcmp("W3", W_field) == 0)   writematrix(W3lens, &galx[gal_index]);
                    if (strcmp("W4", W_field) == 0)   writematrix(W4lens, &galx[gal_index]);
                }
            }
            
            
            // Source selection. (0.3 is set as the lower bound as lens and source have to be at least 0.1 separated in redshift space.)
            if (galx[gal_index].gal_z >= 0.3 && galx[gal_index].gal_z <= 0.4) {
                totcount2++;
                // Boolean test. 1 -> y' band, 0 -> i' band.
                if (field_compare(&galx[gal_index]) == 0 && galx[gal_index].mi > 23.0 && galx[gal_index].mi <= 24.7 && galx[gal_index].modify >= -0.5 && galx[gal_index].modify <= 0.2) {
                    sri_count++;
                    if (strcmp("W1", W_field) == 0)   writematrix(W1src, &galx[gal_index]);
                    if (strcmp("W2", W_field) == 0)   writematrix(W2src, &galx[gal_index]);
                    if (strcmp("W3", W_field) == 0)   writematrix(W3src, &galx[gal_index]);
                    if (strcmp("W4", W_field) == 0)   writematrix(W4src, &galx[gal_index]);
                }
                if (field_compare(&galx[gal_index]) == 1 && galx[gal_index].my > 23.0 && galx[gal_index].my <= 24.7 && galx[gal_index].modify >= -0.5 && galx[gal_index].modify <= 0.2) {
                    sry_count++;
                    if (strcmp("W1", W_field) == 0)   writematrix(W1src, &galx[gal_index]);
                    if (strcmp("W3", W_field) == 0)   writematrix(W3src, &galx[gal_index]);
                    if (strcmp("W4", W_field) == 0)   writematrix(W4src, &galx[gal_index]);
                }
            }
            
            
            if (galx[gal_index].gal_z > 0.4) {
                totcount3++;
                // "Red" Sources.
                    // Boolean test. 1 -> y' band, 0 -> i' band.
                if (field_compare(&galx[gal_index]) == 0 && galx[gal_index].mi <= 24.7 && galx[gal_index].modify >= -0.5 && galx[gal_index].modify <= 0.2) {
                    s2ri_count++;
                    if (strcmp("W1", W_field) == 0)   writematrix(W1src, &galx[gal_index]);
                    if (strcmp("W2", W_field) == 0)   writematrix(W2src, &galx[gal_index]);
                    if (strcmp("W3", W_field) == 0)   writematrix(W3src, &galx[gal_index]);
                    if (strcmp("W4", W_field) == 0)   writematrix(W4src, &galx[gal_index]);
                }
                if (field_compare(&galx[gal_index]) == 1 && galx[gal_index].my <= 24.7 && galx[gal_index].modify >= -0.5 && galx[gal_index].modify <= 0.2) {
                    s2ry_count++;
                    if (strcmp("W1", W_field) == 0)   writematrix(W1src, &galx[gal_index]);
                    if (strcmp("W3", W_field) == 0)   writematrix(W3src, &galx[gal_index]);
                    if (strcmp("W4", W_field) == 0)   writematrix(W4src, &galx[gal_index]);
                }
            }
            
            gal_index++;
        }
        
        free(galx);
        galx = NULL;
        
        fclose(catalog_all);
    }
    
    printf("%d\t%d\t%d\t%d\t%d\t%d\n", totcount, lenscount, ri_count, ry_count, bi_count, by_count);
    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n", totcount2, sri_count, sry_count, sbi_count, sby_count, si_count, sy_count);
    printf("%d\t%d\t%d\n", totcount3, s2ri_count, s2ry_count);
    
    fclose (W1lens);
    fclose (W2lens);
    fclose (W3lens);
    fclose (W4lens);
    fclose (W1src);
    fclose (W2src);
    fclose (W3src);
    fclose (W4src);
    
    free (norm_PDF);
    
    return 0;
}


// The function to write the entire matrix to a txt file.
void writematrix(FILE *infile, gal *w) {
    //fwrite(w, sizeof(gal), 1, infile);
    fprintf(infile, "%s\t%s\t%lf\t%lf\t%lf\t%lf\t%lE\t%lE\t%lf\t%lf\t%d\t%lE\t%lE\t%lf\t%lE\t%lE\t%lf\t%lf\t%lf\t%lf\t%lE\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lE\t%lf\t%lE\t%lf\t%lf\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", w->id, w->field, w->x, w->y, w->alpha, w->delta, w->bg, w->lv, w->mu_max, w->mu_thres, w->flag, w->a, w->b, w->theta, w->a_err, w->b_err, w->theta_err, w->class, w->e1, w->e2, w->weight, w->fitclass, w->SNratio, w->mask, w->gal_z, w->gal_z_min, w->gal_z_max, w->T_BPZ, w->modify, w->c2, w->lp_mu, w->lp_mg, w->lp_mr, w->lp_mi, w->lp_mz, w->star_flag, w->lp_med, w->lp_inf, w->lp_sup, w->mu, w->mu_err, w->ext_u, w->mi, w->mi_err, w->ext_i, w->mr, w->mr_err, w->ext_r, w->mg, w->mg_err, w->ext_g, w->my, w->my_err, w->ext_y, w->mz, w->mz_err, w->ext_z, w->size, w->FWHM_image, w->FWHM_world * 3600, w->kron_rad, w->flux_rad, w->bulge_frac, w->model_flux, w->isoarea * 3600 * 3600, w->PDF[0], w->PDF[1], w->PDF[2], w->PDF[3], w->PDF[4], w->PDF[5], w->PDF[6], w->PDF[7], w->PDF[8], w->PDF[9], w->PDF[10], w->PDF[11], w->PDF[12], w->PDF[13], w->PDF[14], w->PDF[15], w->PDF[16], w->PDF[17], w->PDF[18], w->PDF[19], w->PDF[20], w->PDF[21], w->PDF[22], w->PDF[23], w->PDF[24], w->PDF[25], w->PDF[26], w->PDF[27], w->PDF[28], w->PDF[29], w->PDF[30], w->PDF[31], w->PDF[32], w->PDF[33], w->PDF[34], w->PDF[35], w->PDF[36], w->PDF[37], w->PDF[38], w->PDF[39], w->PDF[40], w->PDF[41], w->PDF[42], w->PDF[43], w->PDF[44], w->PDF[45], w->PDF[46], w->PDF[47], w->PDF[48], w->PDF[49], w->PDF[50], w->PDF[51], w->PDF[52], w->PDF[53], w->PDF[54], w->PDF[55], w->PDF[56], w->PDF[57], w->PDF[58], w->PDF[59], w->PDF[60], w->PDF[61], w->PDF[62], w->PDF[63], w->PDF[64], w->PDF[65], w->PDF[66], w->PDF[67], w->PDF[68], w->PDF[69]);
}

// The function to determine if the field uses the i' (0) or y' (1) band.
int field_compare(gal *w) {
    if (strcmp("W1m0m4", w->field) == 0 || strcmp("W1m1m4", w->field) == 0 || strcmp("W1m2m4", w->field) == 0 || strcmp("W1m3m4", w->field) == 0 || strcmp("W1m4m4", w->field) == 0 || strcmp("W1p1m4", w->field) == 0 || strcmp("W1p1p1", w->field) == 0 || strcmp("W1p2m4", w->field) == 0 || strcmp("W1p3m4", w->field) == 0 || strcmp("W1p3p1", w->field) == 0 || strcmp("W1p4m4", w->field) == 0 || strcmp("W3m0m1", w->field) == 0 || strcmp("W3m2m1", w->field) == 0 || strcmp("W3m2p1", w->field) == 0 || strcmp("W3p2m3", w->field) == 0 || strcmp("W4m1p1", w->field) == 0 || strcmp("W4m1p2", w->field) == 0 || strcmp("W4m1p3", w->field) == 0 || strcmp("W4m2p2", w->field) == 0 || strcmp("W4m2p3", w->field) == 0 || strcmp("W4m3p3", w->field) == 0)    return 1;
    else    return 0;
}

