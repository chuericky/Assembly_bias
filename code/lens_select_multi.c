/* This program selects the lens galaxies which do not belong to any clusters.
 Input: A, factor on virial radius (argv[1])
 Updated: Oct 27, 2015
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "tool.h"
#include "cosmopar.h"

typedef struct {
    char name[15];
    int num;
} namelens;

// Struct for the cluster catalog.
typedef struct {
    double ra, dec, z, sig, mass, r;
} clus;

// Struct for the lens galaxy catalog.
typedef struct {
    char id[20], field[6];
    double x, y, alpha, delta, bg, lv, mu_max, mu_thres, a, b, theta, a_err, b_err, theta_err, class, e1, e2, weight, SNratio, gal_z, gal_z_min, gal_z_max, T_BPZ, modify, c2, lp_mu, lp_mg, lp_mr, lp_mi, lp_mz, lp_med, lp_inf, lp_sup, mu, mu_err, ext_u, mi, mi_err, ext_i, mr, mr_err, ext_r, mg, mg_err, ext_g, my, my_err, ext_y, mz, mz_err, ext_z, size, FWHM_image, FWHM_world, kron_rad, flux_rad, bulge_frac, model_flux, isoarea;
    double PDF[70];
    int flag, fitclass, mask, star_flag;
} gallens;

void writematrix(FILE *infile, gallens *w);


int main(int argc, char *argv[]) {
    char lens_file[150], clus_file[150], lensnumfile[150], r_fac[5];
    FILE *clusname, *lensnumname;
    int i = 0, j = 0, k = 0, m = 0, spot = 0, file_count = -1, clus_num = 0, clus_count = 0, lensfile_num = 0, lensindex = 0, lgal_index = 0, aux_count = 0;
    double c_l_ang = 0, c_l_dis = 0, R_factor = 0;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 2) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    strcpy (r_fac, argv[1]);
    R_factor = atof(r_fac);
    
    // Read in the lens name file.
    strcpy (lensnumfile, lensisosortlink);
    strcat (lensnumfile, "lensname.txt");
    
    lensnumname = fopen (lensnumfile, "r");

    // Make sure the lens name file is present.
    if (!lensnumname) {
        printf("Cannot find the lens name file!\n");
        exit(1);
    }
    
    // Count numer of lines in the lens name file.
    lensfile_num = count_line (spot, lensnumname, file_count);
    
    // An array to store all file names.
    namelens *lens = malloc (lensfile_num * sizeof *lens);
    
    lensindex = 0;
    
    while (fscanf(lensnumname, "%d %s", &lens[lensindex].num, &lens[lensindex].name) == 2) {
        lensindex++;
    }
    
    
    // 4 wide fields.
    for (i = 1; i < 5; i++) {
        
        char field[2];
        
        // Read in the cluster catalog.
        sprintf(field, "%d", i);
        strcpy (clus_file, cluster_link);
        strcat (clus_file, "clusters_W");
        strcat (clus_file, field);
        strcat (clus_file, ".cat");
        
        clusname = fopen (clus_file, "r");
        
        // Make sure the cluster catalog file is present.
        if (!clusname) {
            printf("Cannot find the cluster name file!\n");
            exit(1);
        }
        
        spot = 0;
        file_count = -1;

        // Count number of lines in the cluster catalog file.
        clus_num = count_line (spot, clusname, file_count);
        
        clus *clus_info = malloc (clus_num * sizeof *clus_info);
        
        clus_count = 0;
        
        while (fscanf (clusname, "%*d %lf %lf %lf %lf %lf %lf", &clus_info[clus_count].ra, &clus_info[clus_count].dec, &clus_info[clus_count].z, &clus_info[clus_count].sig, &clus_info[clus_count].mass, &clus_info[clus_count].r) == 6) {
            clus_count++;
        }
             
        // Firstly calculate the angular diameter distances for the clusters. To prevent repetitive calculation.
        double *Da_clus = malloc (clus_count * sizeof (double));
             
        for (j = 0; j < clus_count; j++) {
            Da_clus[j] = Da(clus_info[j].z);
        }
        
        // Lens catalog.
        for (j = 0; j < lensindex; j++) { // lensindex
            
            char W_fieldclus[2];
            
            strncpy (W_fieldclus, lens[j].name+1, 1);
            W_fieldclus[2] = 0;
            
            
            if (strcmp(W_fieldclus, field) == 0) {
                
                printf("%s search begins!\n", lens[j].name);
                
                FILE *lensname, *lensout;//, *lensout2;
                char lensfile_string[150], lensout_string[150];//, lensout2_string[150];
                
                strcpy (lensfile_string, lensisosortlink);
                strcat (lensfile_string, lens[j].name);
                
                lensname = fopen (lensfile_string, "r");
                
                // Make sure the source name file is present.
                if (!lensname) {
                    printf("Cannot find the lens file %s!\n", lens[j].name);
                    exit(1);
                }
                
                strcpy (lensout_string, lensisolink);
                strcat (lensout_string, r_fac);
                strcat (lensout_string, "_");
                strcat (lensout_string, lens[j].name);
                /*strcpy (lensout2_string, lensisolink);
                strcat (lensout2_string, "2");
                strcat (lensout2_string, lens[j].name);
                 
                lensout2 = fopen (lensout2_string, "w");*/
                
                lensout = fopen (lensout_string, "w");
                
                // Make sure the source name file is present.
                if (!lensout) {
                    printf("Cannot open the lens file %s!\n", lens[j].name);
                    exit(1);
                }
                /*if (!lensout2) {
                    printf("Cannot open the lens file %s!\n", lens[j].name);
                    exit(1);
                }*/
                
                // Open up the information of the source galaxies.
                gallens *lgal = malloc (lens[j].num * sizeof *lgal);
                
                lgal_index = 0;
                
                // Total 135 columns of the source catalog.
                while (fscanf(lensname, "%s %s %lf %lf %lf %lf %lE %lE %lf %lf %d %lE %lE %lf %lE %lE %lf %lf %lf %lf %lE %d %lf %d %lf %lf %lf %lf %lE %lE %lE %lE %lE %lE %lE %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lE %lf %lE %lf %lf %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE", &lgal[lgal_index].id, &lgal[lgal_index].field, &lgal[lgal_index].x, &lgal[lgal_index].y, &lgal[lgal_index].alpha, &lgal[lgal_index].delta, &lgal[lgal_index].bg, &lgal[lgal_index].lv, &lgal[lgal_index].mu_max, &lgal[lgal_index].mu_thres, &lgal[lgal_index].flag, &lgal[lgal_index].a, &lgal[lgal_index].b, &lgal[lgal_index].theta, &lgal[lgal_index].a_err, &lgal[lgal_index].b_err, &lgal[lgal_index].theta_err, &lgal[lgal_index].class, &lgal[lgal_index].e1, &lgal[lgal_index].e2, &lgal[lgal_index].weight, &lgal[lgal_index].fitclass, &lgal[lgal_index].SNratio, &lgal[lgal_index].mask, &lgal[lgal_index].gal_z, &lgal[lgal_index].gal_z_min, &lgal[lgal_index].gal_z_max, &lgal[lgal_index].T_BPZ, &lgal[lgal_index].modify, &lgal[lgal_index].c2, &lgal[lgal_index].lp_mu, &lgal[lgal_index].lp_mg, &lgal[lgal_index].lp_mr, &lgal[lgal_index].lp_mi, &lgal[lgal_index].lp_mz, &lgal[lgal_index].star_flag, &lgal[lgal_index].lp_med, &lgal[lgal_index].lp_inf, &lgal[lgal_index].lp_sup, &lgal[lgal_index].mu, &lgal[lgal_index].mu_err, &lgal[lgal_index].ext_u, &lgal[lgal_index].mi, &lgal[lgal_index].mi_err, &lgal[lgal_index].ext_i, &lgal[lgal_index].mr, &lgal[lgal_index].mr_err, &lgal[lgal_index].ext_r, &lgal[lgal_index].mg, &lgal[lgal_index].mg_err, &lgal[lgal_index].ext_g, &lgal[lgal_index].my, &lgal[lgal_index].my_err, &lgal[lgal_index].ext_y, &lgal[lgal_index].mz, &lgal[lgal_index].mz_err, &lgal[lgal_index].ext_z, &lgal[lgal_index].size, &lgal[lgal_index].FWHM_image, &lgal[lgal_index].FWHM_world, &lgal[lgal_index].kron_rad, &lgal[lgal_index].flux_rad, &lgal[lgal_index].bulge_frac, &lgal[lgal_index].model_flux, &lgal[lgal_index].isoarea, &lgal[lgal_index].PDF[0], &lgal[lgal_index].PDF[1], &lgal[lgal_index].PDF[2], &lgal[lgal_index].PDF[3], &lgal[lgal_index].PDF[4], &lgal[lgal_index].PDF[5], &lgal[lgal_index].PDF[6], &lgal[lgal_index].PDF[7], &lgal[lgal_index].PDF[8], &lgal[lgal_index].PDF[9], &lgal[lgal_index].PDF[10], &lgal[lgal_index].PDF[11], &lgal[lgal_index].PDF[12], &lgal[lgal_index].PDF[13], &lgal[lgal_index].PDF[14], &lgal[lgal_index].PDF[15], &lgal[lgal_index].PDF[16], &lgal[lgal_index].PDF[17], &lgal[lgal_index].PDF[18], &lgal[lgal_index].PDF[19], &lgal[lgal_index].PDF[20], &lgal[lgal_index].PDF[21], &lgal[lgal_index].PDF[22], &lgal[lgal_index].PDF[23], &lgal[lgal_index].PDF[24], &lgal[lgal_index].PDF[25], &lgal[lgal_index].PDF[26], &lgal[lgal_index].PDF[27], &lgal[lgal_index].PDF[28], &lgal[lgal_index].PDF[29], &lgal[lgal_index].PDF[30], &lgal[lgal_index].PDF[31], &lgal[lgal_index].PDF[32], &lgal[lgal_index].PDF[33], &lgal[lgal_index].PDF[34], &lgal[lgal_index].PDF[35], &lgal[lgal_index].PDF[36], &lgal[lgal_index].PDF[37], &lgal[lgal_index].PDF[38], &lgal[lgal_index].PDF[39], &lgal[lgal_index].PDF[40], &lgal[lgal_index].PDF[41], &lgal[lgal_index].PDF[42], &lgal[lgal_index].PDF[43], &lgal[lgal_index].PDF[44], &lgal[lgal_index].PDF[45], &lgal[lgal_index].PDF[46], &lgal[lgal_index].PDF[47], &lgal[lgal_index].PDF[48], &lgal[lgal_index].PDF[49], &lgal[lgal_index].PDF[50], &lgal[lgal_index].PDF[51], &lgal[lgal_index].PDF[52], &lgal[lgal_index].PDF[53], &lgal[lgal_index].PDF[54], &lgal[lgal_index].PDF[55], &lgal[lgal_index].PDF[56], &lgal[lgal_index].PDF[57], &lgal[lgal_index].PDF[58], &lgal[lgal_index].PDF[59], &lgal[lgal_index].PDF[60], &lgal[lgal_index].PDF[61], &lgal[lgal_index].PDF[62], &lgal[lgal_index].PDF[63], &lgal[lgal_index].PDF[64], &lgal[lgal_index].PDF[65], &lgal[lgal_index].PDF[66], &lgal[lgal_index].PDF[67], &lgal[lgal_index].PDF[68], &lgal[lgal_index].PDF[69]) == 135) {
                    lgal_index++;
                }
                
                // For each lens galaxies.
                for (k = 0; k < lgal_index; k++) { // lgal_index
                    
                    // Check progress
                    if (k % 10000 == 9999) {
                        printf("%d / %d lens completed!\n", k, lgal_index);
                    }
                    
                    // Magnitude cutoff by Ford et. al. 2015. These galaxies are not classified to be inside any clusters.
                    /*if (lgal[k].lp_mi > -19.35) {
                        writematrix(lensout, &lgal[k]);
                    }
                    
                    if (lgal[k].lp_mi <= -19.35) {*/
                        
                    aux_count = 0; // Have to change with line 153.
                        
                    // Loop over cluster catalog for each lens galaxies.
                    for (m = 0; m < clus_count; m++) { // clus_count, have to change with line 150.

                        // Redshift cut
                        // Outside many sigma.
                        if (fabs(clus_info[m].z - lgal[k].gal_z) > 0.2 * (1 + clus_info[m].z)) {
                            //c_l_ang = ang_sep (clus_info[m].ra, clus_info[m].dec, lgal[k].alpha, lgal[k].delta) * M_PI / (3600. * 180);
                            // Physical transverse separation of the cluster and the lens galaxy, in pc.
                            //c_l_dis = c_l_ang * Da_clus[m] * Mpc_2_pc;
                                
                            aux_count++;
                            continue;   // Proceed to the next cluster.
                        }
                            
                        if (fabs(clus_info[m].z - lgal[k].gal_z) <= 0.2 * (1 + clus_info[m].z)) {
                                
                            // Angular separation between the cluster and the lens galaxy, in radian.
                            c_l_ang = ang_sep (clus_info[m].ra, clus_info[m].dec, lgal[k].alpha, lgal[k].delta) * M_PI / (3600. * 180);
                            // Physical transverse separation of the cluster and the lens galaxy, in pc.
                            c_l_dis = c_l_ang * Da_clus[m] * Mpc_2_pc;
                            
                            //printf("%d\t%lf\t%d\t%lf\t%lf\t%lf\n", m, clus_info[m].z, k, lgal[k].gal_z, fabs(clus_info[m].z - lgal[k].gal_z), 0.2 * (1 + clus_info[m].z));
                            //printf("%lf\t%lE\t%lE\t%lE\n", c_l_ang, Da_clus[m], c_l_dis, clus_info[m].r);
                            
                            // Outside 2 R_vir of the cluster. Just to be safe.
                            if (c_l_dis > R_factor * clus_info[m].r) {
                                aux_count++;
                                continue;
                            }
                                
                            if (c_l_dis <= R_factor * clus_info[m].r) {
                                
                                //writematrix(lensout2, &lgal[k]);
                                break;
                            }
                        }
                    }
                    
                    if (aux_count == clus_count) {
                        writematrix(lensout, &lgal[k]);
                    }
                }
                
                free (lgal);
                fclose (lensname);
                fclose (lensout);;
                //fclose (lensout2);
                
                
                printf("%s search ends!\n", lens[j].name);
            }
            
        }
        
        
        free (clus_info);
        free (Da_clus);
        fclose (clusname);
        
    }
    
    free (lens);
    fclose (lensnumname);
    
    return 0;
}


// The function to write the entire matrix to a txt file.
void writematrix(FILE *infile, gallens *w) {
    fprintf(infile, "%s\t%s\t%lf\t%lf\t%lf\t%lf\t%lE\t%lE\t%lf\t%lf\t%d\t%lE\t%lE\t%lf\t%lE\t%lE\t%lf\t%lf\t%lf\t%lf\t%lE\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lE\t%lf\t%lE\t%lf\t%lf\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", w->id, w->field, w->x, w->y, w->alpha, w->delta, w->bg, w->lv, w->mu_max, w->mu_thres, w->flag, w->a, w->b, w->theta, w->a_err, w->b_err, w->theta_err, w->class, w->e1, w->e2, w->weight, w->fitclass, w->SNratio, w->mask, w->gal_z, w->gal_z_min, w->gal_z_max, w->T_BPZ, w->modify, w->c2, w->lp_mu, w->lp_mg, w->lp_mr, w->lp_mi, w->lp_mz, w->star_flag, w->lp_med, w->lp_inf, w->lp_sup, w->mu, w->mu_err, w->ext_u, w->mi, w->mi_err, w->ext_i, w->mr, w->mr_err, w->ext_r, w->mg, w->mg_err, w->ext_g, w->my, w->my_err, w->ext_y, w->mz, w->mz_err, w->ext_z, w->size, w->FWHM_image, w->FWHM_world, w->kron_rad, w->flux_rad, w->bulge_frac, w->model_flux, w->isoarea, w->PDF[0], w->PDF[1], w->PDF[2], w->PDF[3], w->PDF[4], w->PDF[5], w->PDF[6], w->PDF[7], w->PDF[8], w->PDF[9], w->PDF[10], w->PDF[11], w->PDF[12], w->PDF[13], w->PDF[14], w->PDF[15], w->PDF[16], w->PDF[17], w->PDF[18], w->PDF[19], w->PDF[20], w->PDF[21], w->PDF[22], w->PDF[23], w->PDF[24], w->PDF[25], w->PDF[26], w->PDF[27], w->PDF[28], w->PDF[29], w->PDF[30], w->PDF[31], w->PDF[32], w->PDF[33], w->PDF[34], w->PDF[35], w->PDF[36], w->PDF[37], w->PDF[38], w->PDF[39], w->PDF[40], w->PDF[41], w->PDF[42], w->PDF[43], w->PDF[44], w->PDF[45], w->PDF[46], w->PDF[47], w->PDF[48], w->PDF[49], w->PDF[50], w->PDF[51], w->PDF[52], w->PDF[53], w->PDF[54], w->PDF[55], w->PDF[56], w->PDF[57], w->PDF[58], w->PDF[59], w->PDF[60], w->PDF[61], w->PDF[62], w->PDF[63], w->PDF[64], w->PDF[65], w->PDF[66], w->PDF[67], w->PDF[68], w->PDF[69]);
}