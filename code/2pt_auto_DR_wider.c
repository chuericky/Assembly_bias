/* This program calculates the 2pt auto-correlation function.
 Searching of random galaxies by chaining mesh technique.
 Use Landy-Szalay estimator. (DD - 2DR + RR) / RR.
 Including jack-knifed error estimation.
 
 Input: Field # (argv[1]), random realization (argv[2]), multiple of n lens (argv[3]), value of n (argv[4]), with / without (argv[5]), b (argv[6]), clus_num (argv[7]), color (argv[8]), P_th (argv[9])
 
 Updated: Nov 8, 2015
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
    char name[50];
    int num;
} namelens;

typedef struct {
    double ra, dec, ra_low, ra_up, dec_low, dec_up;
} fieldpos;

typedef struct {
    double ra, dec;
    int jack;
} rangal;

typedef struct {
    char id[20];
    double alpha, delta;
    int jack;
} gallens;

typedef struct {
    int jack_num, jack_min, jack_max;
} jackreg;


int main(int argc, char *argv[]) {
    int i = 0, k = 0, m = 0, p = 0, q = 0, s = 0, t = 0, u = 0, r_bincount = 0, ranindex = 0, lensfile_num = 0, spot = 0, file_count = 0, field_num = 0, fieldindex = 0, field_int = 0, rgal_index = 0, lensindex = 0, lgal_index = 0, r_ind = 0, d_ind = 0, ramesh_num = 0, decmesh_num = 0, rgal_eva = 0, rgal_end = 0, rgal_num = 0, lgal_num = 0, lgal_end = 0, lgal_eva = 0, gal_count = 0;
    double r_min = 0, r_max = 0, dlogr = 0, d2r = M_PI / 180., ramap_size = 0, decmap_size = 0, ramesh_size = 0, decmesh_size = 0, gl_a_sep = 0, ls_a_sep = 0, ann_in = 0, ann_out = 0;
    char lensname_file[150], rannum_file[150], field_file[150], jack_file[150], field_manual[3], ran_manual[8], rgal_manual[8], rgal_manum[8], lgal_manum[8], lgal_manual[8], w_dir[20], color[9], fstr[30], b[5], clus_man[4], P_th[5];
    FILE *ranname, *fieldname, *jackname, *lensname;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 10) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    // Log radius bins, in arcsecs. From 0.003 deg to 3 deg.
    r_bincount = 20;
    r_min = 10.8;
    r_max = 10800;
    dlogr = (log10(r_max) - log10(r_min)) / r_bincount;
    
    double *r_width = malloc ((r_bincount + 1) * sizeof (double));
    double *r_cen = malloc ((r_bincount + 1) * sizeof (double)); // Center of the annulus
    
    // "Widths" of the annulus width.
    for (k = 0; k < (r_bincount + 1); k++) {
        r_width[k] = r_min * pow(10, k * dlogr);
        r_cen[k] = r_min * pow(10, (k - 0.5) * dlogr);   // In arcsecs
        //printf("%d\t%lf\n", k, r_cen[k] / 3600);
    }
    
    strcpy (b, argv[6]);
    strcpy (clus_man, argv[7]);
    strcpy (color, argv[8]);
    strcpy (P_th, argv[9]);
    
    strcpy (fstr, b);
    strcat (fstr, "_");
    strcat (fstr, P_th);
    strcat (fstr, "_");
    strcat (fstr, clus_man);
    strcat (fstr, "_lens");
    strcat (fstr, color);
    
    // Write in the lens files.
    strcpy (lensname_file, shear_gal_link);
    strcpy (w_dir, argv[5]);
    strcat (lensname_file, w_dir);
    strcat (lensname_file, "/");
    strcat (lensname_file, fstr);
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
    
    // Write in the random number files.
    strcpy (rannum_file, ran_link);
    strcat (rannum_file, "ran_num_iso.dat");
    
    ranname = fopen (rannum_file, "r");
    
    // Make sure the lens name file is present.
    if (!ranname) {
        printf("Cannot find the random number file!\n");
        exit(1);
    }
    
    // An array storing all file names.
    unsigned long *rannum = malloc (lensfile_num * sizeof (unsigned long));
    
    ranindex = 0;
    
    while (fscanf(ranname, "%*s %lu", &rannum[ranindex]) == 1) {
        ranindex++;
    }
    
    fclose (ranname);
    
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
    
    // Write in the jackknife region file of the wide fields.
    strcpy (jack_file, field_link);
    strcat (jack_file, "jack_region.dat");
    
    jackname = fopen (jack_file, "r");
    
    // Make sure the field jack-knife file is present.
    if (!jackname) {
        printf("Cannot find the file for jack-knife of the fields!\n");
        exit(1);
    }
    
    // An array storing all file names.
    jackreg *jack = malloc (field_num * sizeof *jack);
    
    fieldindex = 0;
    
    while (fscanf(jackname, "%*s %d %d %d", &jack[fieldindex].jack_num, &jack[fieldindex].jack_min, &jack[fieldindex].jack_max) == 3) {
        fieldindex++;
    }
    
    fclose (jackname);
    
    // For Field W1 to W4.
    strcpy (field_manual, argv[1]);
    field_int = atoi(field_manual);
    
    // Index of random realization.
    strcpy (ran_manual, argv[2]);
    
    for (i = field_int - 1; i < field_int; i++) {
        
        printf("\nField W%d\n", i + 1);
        
        FILE *lensfile, *ranfile;
        char lensfile_string[150], ranfile_string[150];
        
        strcpy (lensfile_string, shear_gal_link);
        strcat (lensfile_string, w_dir);
        strcat (lensfile_string, "/");
        strcat (lensfile_string, lens[i].name);
        
        lensfile = fopen (lensfile_string, "r");
        
        // Make sure the lens name file is present.
        if (!lensfile) {
            printf("Cannot find the lens file %s!\n", lens[i].name);
            exit(1);
        }
        
        // Open up the information of the lens galaxies.
        gallens *lgal = malloc (lens[i].num * sizeof * lgal);
        
        lgal_index = 0;
        
        // lens id, ra, dec.
        while (fscanf(lensfile, "%s %*s %lf %lf %*lf %*lf %*lf %*lf %d", &lgal[lgal_index].id, &lgal[lgal_index].alpha, &lgal[lgal_index].delta, &lgal[lgal_index].jack) == 4) {
        // while (fscanf(lensfile, "%s %*s %lf %lf %d", &lgal[lgal_index].id, &lgal[lgal_index].alpha, &lgal[lgal_index].delta, &lgal[lgal_index].jack) == 4) {
            lgal_index++;
            //if (lgal_index == 100) break;
        }
        
        fclose (lensfile);
        
        strcpy (ranfile_string, ran_link);
        strcat (ranfile_string, "W");
        strcat (ranfile_string, field_manual);
        strcat (ranfile_string, "/combined/ran_iso.dat");
        
        ranfile = fopen (ranfile_string, "r");
        
        // Make sure the lens name file is present.
        if (!ranfile) {
            printf("Cannot find the random realization file %s!\n", ran_manual);
            exit(1);
        }
        
        // Open up the information of the lens galaxies.
        rangal *rgal = malloc (rannum[i] * sizeof * rgal);
        
        rgal_index = 0;
        
        // Random galaxy ra, dec.
        while (fscanf(ranfile, "%d %lf %lf", &rgal[rgal_index].jack, &rgal[rgal_index].ra, &rgal[rgal_index].dec) == 3) {
            rgal_index++;
            // if (rgal_index == 10000) break;
        }
        
        fclose (ranfile);
        
        
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
        
        
        // Create a map to store up which mesh do the lens galaxies belong to, first loop over RA, then DEC. Another array is to store the index of the galaxies in the source list.
        int **lgal_count = malloc (decmesh_num * sizeof (int *));
        int **lgal_cum = malloc (decmesh_num * sizeof (int *));     // Cumulative number of lens galaxies.
        double ***lgal_coord = malloc (decmesh_num * sizeof (double **));     // Coordinates of the center of the grids.
        int **lgal_meshsort = malloc (lgal_index * sizeof (int *));
        
        for (k = 0; k < lgal_index; k++) {
            lgal_meshsort[k] = malloc (3 * sizeof (int));
            for (p = 0; p < 3; p++) {
                lgal_meshsort[k][p] = 0;    // [k][0] is mesh num in along RA, [k][1] is mesh num in along DEC, [k][2] is index of the galaxy.
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            lgal_count[k] = malloc (ramesh_num * sizeof (int));
            lgal_cum[k] = malloc (ramesh_num * sizeof (int));
            lgal_coord[k] = malloc (ramesh_num * sizeof (double *));
            for (p = 0; p < ramesh_num; p++) {
                lgal_count[k][p] = 0;
                lgal_cum[k][p] = 0;
                lgal_coord[k][p] = malloc (2 * sizeof (double));
                for (q = 0; q < 2; q++) {
                    lgal_coord[k][p][q] = 0;
                }
            }
        }
        
        // Store up the number of counts in the mesh grid
        for (k = 0; k < lgal_index; k++) {
            r_ind = min(ra_center, ramesh_num, lgal[k].alpha);
            d_ind = min(dec_center, decmesh_num, lgal[k].delta);
            lgal_meshsort[k][0] = r_ind;
            lgal_meshsort[k][1] = d_ind;
            lgal_meshsort[k][2] = k;
            lgal_count[d_ind][r_ind]++;
        }
        
        // Cumulative number of galaxies in each grid.
        for (k = 0; k < 1; k++) {
            lgal_cum[k][0] = lgal_count[0][0];
            for (p = 1; p < ramesh_num; p++) {
                lgal_cum[k][p] = lgal_cum[k][p - 1] + lgal_count[k][p];
            }
        }
        
        for (k = 1; k < decmesh_num; k++) {
            lgal_cum[k][0] = lgal_cum[k - 1][ramesh_num - 1] + lgal_count[k][0];
            for (p = 1; p < ramesh_num; p++) {
                lgal_cum[k][p] = lgal_cum[k][p - 1] + lgal_count[k][p];
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            for (p = 0; p < ramesh_num; p++) {
                for (q = 0; q < 2; q++) {
                    lgal_coord[k][p][0] = ra_center[p];
                    lgal_coord[k][p][1] = dec_center[k];
                }
            }
        }
        
        // Sort according to grid of RA; if the same, then sort according to grid of DEC; if still the same, then sort by galaxy index in the lens galaxy catalog.
        qsort(lgal_meshsort, lgal_index, sizeof lgal_meshsort[0], compare);
        
        
        // Create a map to store up which mesh do the random galaxies belong to, first loop over RA, then DEC. Another array is to store the index of the galaxies in the source list.
        int **rgal_count = malloc (decmesh_num * sizeof (int *));
        int **rgal_cum = malloc (decmesh_num * sizeof (int *));     // Cumulative number of lens galaxies.
        int **rgal_meshsort = malloc (rgal_index * sizeof (int *));
        
        for (k = 0; k < rgal_index; k++) {
            rgal_meshsort[k] = malloc (3 * sizeof (int));
            for (p = 0; p < 3; p++) {
                rgal_meshsort[k][p] = 0;    // [k][0] is mesh num in along RA, [k][1] is mesh num in along DEC, [k][2] is index of the galaxy.
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            rgal_count[k] = malloc (ramesh_num * sizeof (int));
            rgal_cum[k] = malloc (ramesh_num * sizeof (int));
            for (p = 0; p < ramesh_num; p++) {
                rgal_count[k][p] = 0;
                rgal_cum[k][p] = 0;
            }
        }
        
        // Store up the number of counts in the mesh grid
        for (k = 0; k < rgal_index; k++) {
            r_ind = min(ra_center, ramesh_num, rgal[k].ra);
            d_ind = min(dec_center, decmesh_num, rgal[k].dec);
            rgal_meshsort[k][0] = r_ind;
            rgal_meshsort[k][1] = d_ind;
            rgal_meshsort[k][2] = k;
            rgal_count[d_ind][r_ind]++;
        }
        
        // Cumulative number of galaxies in each grid.
        for (k = 0; k < 1; k++) {
            rgal_cum[k][0] = rgal_count[0][0];
            for (p = 1; p < ramesh_num; p++) {
                rgal_cum[k][p] = rgal_cum[k][p - 1] + rgal_count[k][p];
            }
        }
        
        for (k = 1; k < decmesh_num; k++) {
            rgal_cum[k][0] = rgal_cum[k - 1][ramesh_num - 1] + rgal_count[k][0];
            for (p = 1; p < ramesh_num; p++) {
                rgal_cum[k][p] = rgal_cum[k][p - 1] + rgal_count[k][p];
            }
        }
        
        
        // Sort according to grid of RA; if the same, then sort according to grid of DEC; if still the same, then sort by galaxy index in the random galaxy catalog.
        qsort(rgal_meshsort, rgal_index, sizeof rgal_meshsort[0], compare);
        
        // Multiple numbers of lens galaxies, starting from 0.
        strcpy (lgal_manual, argv[3]);
        strcpy (lgal_manum, argv[4]);
        lgal_eva = atoi(lgal_manum);
        lgal_num = atoi(lgal_manual) * lgal_eva;
        
        lgal_end = lgal_index;
        
        if ((atoi(lgal_manual) + 1) * lgal_eva < lgal_end) {
            lgal_end = (atoi(lgal_manual) + 1) * lgal_eva;
        }
        
        unsigned long **paircount = malloc ((lgal_end - lgal_num) * sizeof (unsigned long *));
        unsigned long ***paircount_jack = malloc ((lgal_end - lgal_num) * sizeof (unsigned long **));
        
        for (k = 0; k < (lgal_end - lgal_num); k++) {
            paircount[k] = malloc ((r_bincount + 1) * sizeof (unsigned long));
            paircount_jack[k] = malloc ((r_bincount + 1) * sizeof (unsigned long *));
            for (m = 0; m < (r_bincount + 1); m++) {
                paircount[k][m] = 0;
                paircount_jack[k][m] = malloc (jack[i].jack_num * sizeof (unsigned long));
                for (q = 0; q < jack[i].jack_num; q++) {
                    paircount_jack[k][m][q] = 0;
                }
            }
        }
        
        /*char twoptDR_file[150];
        FILE *twoptDRname;
        
        strcpy (twoptDR_file, twopt2_link);
        strcat (twoptDR_file, "W");
        strcat (twoptDR_file, field_manual);
        strcat (twoptDR_file, "/DR_");
        strcat (twoptDR_file, ran_manual);
        strcat (twoptDR_file, "_");
        strcat (twoptDR_file, lgal_manual);
        strcat (twoptDR_file, ".dat");
        
        twoptDRname = fopen (twoptDR_file, "w");
        
        char twoptDR_jackfile[150];
        FILE *twoptDRjackname;
        
        strcpy (twoptDR_jackfile, twopt2_link);
        strcat (twoptDR_jackfile, "W");
        strcat (twoptDR_jackfile, field_manual);
        strcat (twoptDR_jackfile, "/DR_");
        strcat (twoptDR_jackfile, ran_manual);
        strcat (twoptDR_jackfile, "_");
        strcat (twoptDR_jackfile, lgal_manual);
        strcat (twoptDR_jackfile, "jack.dat");
        
        twoptDRjackname = fopen (twoptDR_jackfile, "w");*/
        
        gal_count = 0;
        
        int *jack_count = malloc (jack[i].jack_num * sizeof (int));
        
        for (k = 0; k < jack[i].jack_num; k++) {
            jack_count[k] = 0;
        }
        
        // For a particular subset of lens
        for (k = lgal_num; k < lgal_end; k++) {
            
            // Count
            if ((lgal_end - k) % 1000 == 999) printf("%d / %d lens galaxies calculations done!\t#%d\n", (k - lgal_num), (lgal_end - lgal_num), k);
            
            // Chaining mesh.
            for (p = 0; p < decmesh_num; p++) { //decmesh_num
                for (q = 0; q < ramesh_num; q++) {  //ramesh_num
                    
                    // Angular separation between the lens and the grid.
                    gl_a_sep = ang_sep (lgal[k].alpha, lgal[k].delta, lgal_coord[p][q][0], lgal_coord[p][q][1]);
                    
                    for (m = 1; m < (r_bincount + 1); m++) {
                        
                        // Angular radius of the annulus.
                        ann_in = r_width[m - 1];
                        ann_out = r_width[m];
                        
                        // Considering also the edging effect. Make sure there is at least 1 random galaxy in the grid.
                        if ((gl_a_sep - ann_out <= 1.5 * ramesh_size * 3600) && (ann_in - gl_a_sep <= 1.5 * ramesh_size * 3600) && rgal_count[p][q] > 0) {
                            
                            // Create a dynamic array to store the index of the corresponding random galaxies.
                            int *rgrid_sort = malloc (rgal_count[p][q] * sizeof (int));
                            int *l_ind = malloc (rgal_count[p][q] * sizeof (int));
                            
                            for (s = 0; s < rgal_count[p][q]; s++) {    // lgal_count[p][q]
                                
                                rgrid_sort[s] = rgal_cum[p][q] - rgal_count[p][q] + s;  // Number of concerned random galaxy
                                l_ind[s] = rgal_meshsort[rgrid_sort[s]][2];     // The index of the concerned random galaxy in the lens catalog.
                                
                                
                                // Make sure the code is reading in the lens galaxies in the desired grid.
                                if (rgal_meshsort[rgrid_sort[s]][0] != q || rgal_meshsort[rgrid_sort[s]][1] != p) {
                                    printf("Error! A random galaxy is not in the concerned grid.\n");
                                    exit (1);
                                }
                                
                                // Lens-lens separation.
                                ls_a_sep = ang_sep (lgal[k].alpha, lgal[k].delta, rgal[l_ind[s]].ra, rgal[l_ind[s]].dec);
                                    
                                // For sources inside the annulus.
                                if (ls_a_sep <= ann_out && ls_a_sep >= ann_in) {
                                        paircount[k - lgal_num][m]++;
                                    //printf("%d\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\n", k, lgal[k].jack, m, rgrid_sort[s], rgal[l_ind[s]].ra, rgal[l_ind[s]].dec, rgal[l_ind[s]].jack, l_ind[s], paircount[k - lgal_num][m]);
                                        
                                    // Jack-knifed sample.
                                    for (t = jack[i].jack_min; t < jack[i].jack_max + 1; t++) {
                                        if (lgal[k].jack != t && rgal[l_ind[s]].jack != t) {
                                            paircount_jack[k - lgal_num][m][t - jack[i].jack_min]++;
                                            //printf("%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\t%d\n", t, k, lgal[k].jack, m, rgrid_sort[s], rgal[l_ind[s]].ra, rgal[l_ind[s]].dec, rgal[l_ind[s]].jack, l_ind[s], paircount_jack[k - lgal_num][m][t - jack[i].jack_min], t - jack[i].jack_min);
                                        }
                                    }
                                }
                                
                                
                            } // s loop, lens galaxies in that grid.
                            
                            free (rgrid_sort);
                            free (l_ind);
                        }
                    }  // m loop, for each annulus bin.
                } // q loop, ra loop
            }  // p loop, dec loop
            
            gal_count++;
            
            /*fprintf(twoptDRname, "%s\t", lgal[k].id);
            
            for (m = 1; m < r_bincount; m++) {
                fprintf(twoptDRname, "%lu\t", paircount[k - lgal_num][m]);
            }
            fprintf(twoptDRname, "%lu\n", paircount[k - lgal_num][r_bincount]);*/
            
            
            
            for (u = 0; u < jack[i].jack_num; u++) {
                if (lgal[k].jack - jack[i].jack_min != u) jack_count[u]++;
            }
            
            /*s = 48;
             
             if (lgal[k].jack != s) {
             fprintf(twoptDRjackname, "%d\t", k);
             
             for (m = 1; m < r_bincount; m++) {
             fprintf(twoptDRjackname, "%lu\t", paircount_jack[k - lgal_num][m][s - jack[i].jack_min]);
             }
             fprintf(twoptDRjackname, "%lu\n", paircount_jack[k - lgal_num][r_bincount][s - jack[i].jack_min]);
             }*/
            
        }  // k loop, lens galaxies.

        /*printf("Galaxy count = %d\n", gal_count);
        
        for (u = 0; u < jack[i].jack_num; u++) {
            printf("%d\t%d\n", u, jack_count[u]);
        }
        
        printf("%d\t%d\t%d\t%d\t%d\t%d\n", lgal_eva, lgal_num, lgal_end, jack[i].jack_num, jack[i].jack_min, jack[i].jack_max);*/
        
        char twoptDR_comfile[150];
        FILE *twoptDRcomname;
        
        strcpy (twoptDR_comfile, twopt_link);
        strcat (twoptDR_comfile, w_dir);
        strcat (twoptDR_comfile, "/W");
        strcat (twoptDR_comfile, field_manual);
        strcat (twoptDR_comfile, "/");
        strcat (twoptDR_comfile, fstr);
        strcat (twoptDR_comfile, "_DR_");
        strcat (twoptDR_comfile, lgal_manual);
        strcat (twoptDR_comfile, "_auto_all_wider.dat");
        
        twoptDRcomname = fopen (twoptDR_comfile, "w");
        
        unsigned long long *tot_paircount = malloc ((r_bincount + 1) * sizeof (unsigned long long));
        
        for (k = 0; k < (r_bincount + 1); k++) {
            tot_paircount[k] = 0;
        }
        
        for (k = 0; k < gal_count; k++) {
            for (m = 0; m < (r_bincount + 1); m++) {
                tot_paircount[m] += paircount[k][m];
            }
        }
        
        fprintf(twoptDRcomname, "%d\t", gal_count);
        
        for (m = 1; m < r_bincount; m++) {
            fprintf(twoptDRcomname, "%llu\t", tot_paircount[m]);
        }
        fprintf(twoptDRcomname, "%llu\n", tot_paircount[r_bincount]);
        
        
        char twoptDR_jackcomfile[150];
        FILE *twoptDRjackcomname;
        
        strcpy (twoptDR_jackcomfile, twopt_link);
        strcat (twoptDR_jackcomfile, w_dir);
        strcat (twoptDR_jackcomfile, "/W");
        strcat (twoptDR_jackcomfile, field_manual);
        strcat (twoptDR_jackcomfile, "/");
        strcat (twoptDR_jackcomfile, fstr);
        strcat (twoptDR_jackcomfile, "_DR_");
        strcat (twoptDR_jackcomfile, lgal_manual);
        strcat (twoptDR_jackcomfile, "_auto_jack_wider.dat");
        
        twoptDRjackcomname = fopen (twoptDR_jackcomfile, "w");
        
        unsigned long long **tot_jackpaircount = malloc ((r_bincount + 1) * sizeof (unsigned long long *));
        
        for (k = 0; k < (r_bincount + 1); k++) {
            tot_jackpaircount[k] = malloc (jack[i].jack_num * sizeof (unsigned long long));
            for (u = 0; u < jack[i].jack_num; u++) {
                tot_jackpaircount[k][u] = 0;
            }
        }
        
        for (k = 0; k < gal_count; k++) {
            for (m = 0; m < (r_bincount + 1); m++) {
                for (u = 0; u < jack[i].jack_num; u++) {
                    tot_jackpaircount[m][u] += paircount_jack[k][m][u];
                }
            }
        }
        
        for (u = 0; u < jack[i].jack_num; u++) {
            fprintf(twoptDRjackcomname, "%d\t%d\t", jack[i].jack_min + u, jack_count[u]);
            
            for (m = 1; m < r_bincount; m++) {
                fprintf(twoptDRjackcomname, "%llu\t", tot_jackpaircount[m][u]);
            }
            fprintf(twoptDRjackcomname, "%llu\n", tot_jackpaircount[r_bincount][u]);
        }
        
        
        
        
        for (k = 0; k < (lgal_end - lgal_num); k++) {
            for (m = 0; m < (r_bincount + 1); m++) {
                free (paircount_jack[k][m]);
            }
            free (paircount[k]);
            free (paircount_jack[k]);
        }
        for (m = 0; m < (r_bincount + 1); m++) {
            free (tot_jackpaircount[m]);
        }
        
        for (k = 0; k < decmesh_num; k++) {
            free(lgal_count[k]);
            free(lgal_cum[k]);
            free(rgal_count[k]);
            free(rgal_cum[k]);
            for (p = 0; p < ramesh_num; p++) {
                free(lgal_coord[k][p]);
            }
            free(lgal_coord[k]);
        }
        
        for (k = 0; k < lgal_index; k++) {
            free (lgal_meshsort[k]);
        }
        
        for (k = 0; k < rgal_index; k++) {
            free (rgal_meshsort[k]);
        }
        
        free (tot_jackpaircount);
        free (tot_paircount);
        free (paircount);
        free (paircount_jack);
        free (lgal_meshsort);
        free (lgal_cum);
        free (lgal_count);
        free (lgal_coord);
        free (jack_count);
        free (rgal_meshsort);
        free (rgal_cum);
        free (rgal_count);
        free (ra_center);
        free (dec_center);
        free (lgal);
        free (rgal);
        
        /*fclose (twoptDRname);
        fclose (twoptDRjackname);*/
        fclose (twoptDRcomname);
        fclose (twoptDRjackcomname);

    } // i-loop, field.
    
    
    free (fpos);
    free (lens);
    free (jack);
    free (r_width);
    free (r_cen);
    free (rannum);

    return 0;
}


