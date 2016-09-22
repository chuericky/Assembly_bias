/* This program calculates the 2pt auto-correlation function.
 Searching of random galaxies by chaining mesh technique.
 Use Landy-Szalay estimator. (DD - 2DR + RR) / RR.
 Including jack-knifed error estimation.
 
 Input: Field # (argv[1]), random realization (argv[2]), multiple of n lens (argv[3]), value of n (argv[4])
 
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
    char name[20];
    int num;
} namelens;

typedef struct {
    double ra, dec, ra_low, ra_up, dec_low, dec_up;
} fieldpos;

typedef struct {
    int jack_num, jack_min, jack_max;
} jackreg;

typedef struct {
    double ra, dec;
    int jack;
} rangal;


int main(int argc, char *argv[]) {
    int i = 0, k = 0, m = 0, p = 0, q = 0, s = 0, t = 0, u = 0, r_bincount = 0, ranindex = 0, lensfile_num = 4, spot = 0, file_count = 0, field_num = 0, fieldindex = 0, field_int = 0, rgal_index = 0, r_ind = 0, d_ind = 0, ramesh_num = 0, decmesh_num = 0, rgal_eva = 0, rgal_end = 0, rgal_num = 0, gal_count = 0;
    double r_min = 0, r_max = 0, dlogr = 0, d2r = M_PI / 180., ramap_size = 0, decmap_size = 0, ramesh_size = 0, decmesh_size = 0, gl_a_sep = 0, ls_a_sep = 0, ann_in = 0, ann_out = 0;
    char rannum_file[150], field_file[150], jack_file[150], field_manual[3], ran_manual[8], rgal_manual[8], rgal_manum[8], parameter_name[20], paracateg_name[20];
    FILE *ranname, *fieldname, *jackname;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 5) {
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
        printf("%d\t%lf\n", k, r_cen[k] / 3600);
    }
    
    getchar();
    
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
        
        FILE *ranfile;
        char ranfile_string[150];
        
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
            if (rgal_index == 1000) break;
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
        double ***lgal_coord = malloc (decmesh_num * sizeof (double **));     // Coordinates of the center of the grids.
        
        for (k = 0; k < decmesh_num; k++) {
            lgal_coord[k] = malloc (ramesh_num * sizeof (double *));
            for (p = 0; p < ramesh_num; p++) {
                lgal_coord[k][p] = malloc (2 * sizeof (double));
                for (q = 0; q < 2; q++) {
                    lgal_coord[k][p][q] = 0;
                }
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            for (p = 0; p < ramesh_num; p++) {
                lgal_coord[k][p][0] = ra_center[p];
                lgal_coord[k][p][1] = dec_center[k];
            }
        }
        
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
        
        // Multiple numbers of random galaxies, starting from 0.
        strcpy (rgal_manual, argv[3]);
        strcpy (rgal_manum, argv[4]);
        rgal_eva = atoi(rgal_manum);
        rgal_num = atoi(rgal_manual) * rgal_eva;
        
        rgal_end = rgal_index;
        
        if ((atoi(rgal_manual) + 1) * rgal_eva < rgal_end) {
            rgal_end = (atoi(rgal_manual) + 1) * rgal_eva;
        }
        
        unsigned long **paircount = malloc ((rgal_end - rgal_num) * sizeof (unsigned long *));
        unsigned long ***paircount_jack = malloc ((rgal_end - rgal_num) * sizeof (unsigned long **));
        
        for (k = 0; k < (rgal_end - rgal_num); k++) {
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
        
        /*char twoptRR_file[150];
        FILE *twoptRRname;
        
        strcpy (twoptRR_file, twopt2_link);
        strcat (twoptRR_file, "W");
        strcat (twoptRR_file, field_manual);
        strcat (twoptRR_file, "/RR_");
        strcat (twoptRR_file, ran_manual);
        strcat (twoptRR_file, "_");
        strcat (twoptRR_file, rgal_manual);
        strcat (twoptRR_file, ".dat");
        
        twoptRRname = fopen (twoptRR_file, "w");
        
        char twoptRR_jackfile[150];
        FILE *twoptRRjackname;
        
        strcpy (twoptRR_jackfile, twopt2_link);
        strcat (twoptRR_jackfile, "W");
        strcat (twoptRR_jackfile, field_manual);
        strcat (twoptRR_jackfile, "/RR_");
        strcat (twoptRR_jackfile, ran_manual);
        strcat (twoptRR_jackfile, "_");
        strcat (twoptRR_jackfile, rgal_manual);
        strcat (twoptRR_jackfile, "jack.dat");
        
        twoptRRjackname = fopen (twoptRR_jackfile, "w");*/
        
        gal_count = 0;
        
        int *jack_count = malloc (jack[i].jack_num * sizeof (int));
        
        for (k = 0; k < jack[i].jack_num; k++) {
            jack_count[k] = 0;
        }
        
        // For a particular subset of lens
        for (k = rgal_num; k < rgal_end; k++) {
            
            // Count
            if ((rgal_end - k) % 1000 == 999) printf("%d / %d random galaxies calculations done!\t#%d\n", (k - rgal_num), (rgal_end - rgal_num), k);
            
            // Chaining mesh.
            for (p = 0; p < decmesh_num; p++) { //decmesh_num
                for (q = 0; q < ramesh_num; q++) {  //ramesh_num
                    
                    // Angular separation between the lens and the grid.
                    gl_a_sep = ang_sep (rgal[k].ra, rgal[k].dec, lgal_coord[p][q][0], lgal_coord[p][q][1]);
                    
                    for (m = 1; m < (r_bincount + 1); m++) {
                        
                        // Angular radius of the annulus.
                        ann_in = r_width[m - 1];
                        ann_out = r_width[m];
                        
                        // Considering also the edging effect. Make sure there is at least 1 lens galaxy in the grid.
                        if ((gl_a_sep - ann_out <= 1.5 * ramesh_size * 3600) && (ann_in - gl_a_sep <= 1.5 * ramesh_size * 3600) && rgal_count[p][q] > 0) {
                            
                            // Create a dynamic array to store the index of the corresponding lens galaxies.
                            int *rgrid_sort = malloc (rgal_count[p][q] * sizeof (int));
                            int *r_ind = malloc (rgal_count[p][q] * sizeof (int));
                            
                            for (s = 0; s < rgal_count[p][q]; s++) {    // rgal_count[p][q]
                                
                                rgrid_sort[s] = rgal_cum[p][q] - rgal_count[p][q] + s;  // Number of concerned lens galaxy
                                r_ind[s] = rgal_meshsort[rgrid_sort[s]][2];     // The index of the concerned lens galaxy in the lens catalog.
                                
                                
                                // Make sure the code is reading in the lens galaxies in the desired grid.
                                if (rgal_meshsort[rgrid_sort[s]][0] != q || rgal_meshsort[rgrid_sort[s]][1] != p) {
                                    printf("Error! A random galaxy is not in the concerned grid.\n");
                                    exit (1);
                                }
                                
                                // Not the random galaxy itself
                                if (k != r_ind[s]) {
                                    
                                    // Lens-lens separation.
                                    ls_a_sep = ang_sep (rgal[k].ra, rgal[k].dec, rgal[r_ind[s]].ra, rgal[r_ind[s]].dec);
                                    
                                    // For sources inside the annulus.
                                    if (ls_a_sep <= ann_out && ls_a_sep >= ann_in) {
                                        paircount[k - rgal_num][m]++;
                                        //printf("%d\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\n", k, rgal[k].jack, m, rgrid_sort[s], rgal[r_ind[s]].ra, rgal[r_ind[s]].dec, rgal[r_ind[s]].jack, r_ind[s], paircount[k - rgal_num][m]);
                                        
                                        // Jack-knifed sample.
                                        for (t = jack[i].jack_min; t < jack[i].jack_max + 1; t++) {
                                            if (rgal[k].jack != t && rgal[r_ind[s]].jack != t) {
                                                paircount_jack[k - rgal_num][m][t - jack[i].jack_min]++;
                                                //printf("%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\t%d\n", t, k, rgal[k].jack, m, rgrid_sort[s], rgal[r_ind[s]].ra, rgal[r_ind[s]].dec, rgal[r_ind[s]].jack, r_ind[s], paircount_jack[k - rgal_num][m][t - jack[i].jack_min], t - jack[i].jack_min);
                                            }
                                        }
                                    }
                                }
                                
                            } // s loop, lens galaxies in that grid.
                            
                            free (rgrid_sort);
                            free (r_ind);
                        }
                    }  // m loop, for each annulus bin.
                } // q loop, ra loop
            }  // p loop, dec loop
            
            
            gal_count++;
            
            /*fprintf(twoptRRname, "%d\t", k);
             
             for (m = 1; m < r_bincount; m++) {
                 fprintf(twoptRRname, "%lu\t", paircount[k - rgal_num][m]);
             }
             fprintf(twoptRRname, "%lu\n", paircount[k - rgal_num][r_bincount]);*/
            
            
            for (u = 0; u < jack[i].jack_num; u++) {
                if (rgal[k].jack - jack[i].jack_min != u) jack_count[u]++;
            }
            
            /*s = 50;
             
             if (rgal[k].jack != s) {
                 fprintf(twoptRRjackname, "%d\t", k);
             
             for (m = 1; m < r_bincount; m++) {
                 fprintf(twoptRRjackname, "%lu\t", paircount_jack[k - rgal_num][m][s - jack[i].jack_min]);
             }
             fprintf(twoptRRjackname, "%lu\n", paircount_jack[k - rgal_num][r_bincount][s - jack[i].jack_min]);
             }*/
            
            //getchar();
            
        }  // k loop, lens galaxies.
        
        /*printf("Galaxy count = %d\n", gal_count);
        
        for (u = 0; u < jack[i].jack_num; u++) {
            printf("%d\t%d\n", u, jack_count[u]);
        }
        
        printf("%d\t%d\t%d\t%d\t%d\t%d\n", rgal_eva, rgal_num, rgal_end, jack[i].jack_num, jack[i].jack_min, jack[i].jack_max);*/
        
        char twoptRR_comfile[150];
        FILE *twoptRRcomname;
        
        strcpy (twoptRR_comfile, twopt_link);
        strcat (twoptRR_comfile, "without/W");
        strcat (twoptRR_comfile, field_manual);
        strcat (twoptRR_comfile, "/iso_RR_");
        strcat (twoptRR_comfile, rgal_manual);
        strcat (twoptRR_comfile, "_auto_all_wider.dat");
        
        twoptRRcomname = fopen (twoptRR_comfile, "w");
        
        unsigned long long *tot_paircount = malloc ((r_bincount + 1) * sizeof (unsigned long long));
        
        for (k = 0; k < (r_bincount + 1); k++) {
            tot_paircount[k] = 0;
        }
        
        for (k = 0; k < gal_count; k++) {
            for (m = 0; m < (r_bincount + 1); m++) {
                tot_paircount[m] += paircount[k][m];
            }
        }
        
        fprintf(twoptRRcomname, "%d\t", gal_count);
        
        for (m = 1; m < r_bincount; m++) {
            fprintf(twoptRRcomname, "%llu\t", tot_paircount[m]);
        }
        fprintf(twoptRRcomname, "%llu\n", tot_paircount[r_bincount]);
        
        
        char twoptRR_jackcomfile[150];
        FILE *twoptRRjackcomname;
        
        strcpy (twoptRR_jackcomfile, twopt_link);
        strcat (twoptRR_jackcomfile, "without/W");
        strcat (twoptRR_jackcomfile, field_manual);
        strcat (twoptRR_jackcomfile, "/iso_RR_");
        strcat (twoptRR_jackcomfile, rgal_manual);
        strcat (twoptRR_jackcomfile, "_auto_jack_wider.dat");
        
        twoptRRjackcomname = fopen (twoptRR_jackcomfile, "w");
        
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
            fprintf(twoptRRjackcomname, "%d\t%d\t", jack[i].jack_min + u, jack_count[u]);
            
            for (m = 1; m < r_bincount; m++) {
                fprintf(twoptRRjackcomname, "%llu\t", tot_jackpaircount[m][u]);
            }
            fprintf(twoptRRjackcomname, "%llu\n", tot_jackpaircount[r_bincount][u]);
        }


        for (k = 0; k < decmesh_num; k++) {
            free(rgal_count[k]);
            free(rgal_cum[k]);
            for (p = 0; p < ramesh_num; p++) {
                free(lgal_coord[k][p]);
            }
            free(lgal_coord[k]);
        }
        
        for (k = 0; k < rgal_index; k++) {
            free (rgal_meshsort[k]);
        }
        
        for (k = 0; k < (rgal_end - rgal_num); k++) {
            for (m = 0; m < (r_bincount + 1); m++) {
                free (paircount_jack[k][m]);
            }
            free (paircount[k]);
            free (paircount_jack[k]);
        }
        for (m = 0; m < (r_bincount + 1); m++) {
            free (tot_jackpaircount[m]);
        }
        
        free (tot_jackpaircount);
        free (tot_paircount);
        free (paircount);
        free (paircount_jack);
        free (jack_count);
        free (rgal);
        free (rgal_meshsort);
        free (rgal_cum);
        free (rgal_count);
        free (lgal_coord);
        free (ra_center);
        free (dec_center);
        
        /*fclose (twoptRRname);
        fclose (twoptRRjackname);*/
        fclose (twoptRRcomname);
        fclose (twoptRRjackcomname);
    } // i-loop, field.

    free (jack);
    free (fpos);
    free (rannum);
    free (r_width);
    free (r_cen);


    return 0;
}



