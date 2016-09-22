/* This program calculates the 2pt auto-correlation function.
 Searching of lens galaxies by chaining mesh technique.
 Use Landy-Szalay estimator. (DD - 2DR + RR) / RR.
 Including jack-knifed error estimation.
 
 Input: File index (argv[1]), multiple of n lens (argv[2]), value of n (argv[3])
 
 Updated: Aug 12, 2015
 */

// Struct for the galaxy catalog.
typedef struct {
    char name[20];
    int num;
} namelens;

typedef struct {
    char id[12];
    double x, y, z;
    int jack;
} gallens;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "tool.h"
#include "cosmopar.h"


int main(int argc, char *argv[]) {
    int i = 0, k = 0, p = 0, q = 0, m = 0, s = 0, t = 0, u = 0, r_bincount = 0, spot = 0, file_count = 0, lensfile_num = 0, lensindex = 0, cat_int = 0, lgal_index = 0, xmesh_num = 0, ymesh_num = 0, zmesh_num = 0, x_ind = 0, y_ind = 0, z_ind = 0, lgal_eva = 0, lgal_num = 0, lgal_end = 0, z_num = 0, gal_count = 0, jack_num = 0;
    double r_min = 0, r_max = 0, logr = 0, xmap_size = 0, ymap_size = 0, zmap_size = 0, xmesh_size = 0, ymesh_size = 0, zmesh_size = 0, z_max = 0, gl_a_sep = 0, ann_in = 0, ann_out = 0, ls_a_sep = 0, ls_z_sep = 0, z_ref = 0, Dc_ref = 0, d2r = M_PI / 180;
    char lensname_file[150], cat_manual[3], lgal_manual[5], lgal_manum[7], W_name[7];
    FILE *lensname;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 4) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    // Log radius bins, in arcsecs. From 0.003 deg to 3 deg.
    r_bincount = 30;
    r_min = 10.8;
    r_max = 10800;
    logr = (log10(r_max) - log10(r_min)) / r_bincount;
    
    double *r_width = malloc ((r_bincount + 1) * sizeof (double));
    double *r_cen = malloc ((r_bincount + 1) * sizeof (double));
    
    // "Widths" and centers of the annulus width.
    for (k = 0; k < (r_bincount + 1); k++) {
        r_width[k] = r_min * pow(10, k * logr);
        r_cen[k] = r_min * pow(10, (k - 0.5) * logr);
        //printf("%d\t%f\t%f\t%f\n", k, r_width[k], r_cen[k], r_cen[k]/3600);
    }
    
    // Redshift where simulation is located, and the comoving distance of that redshift.
    z_ref = 0.1;
    Dc_ref = Dc(z_ref);
    
    // Write in the lens files.
    strcpy (lensname_file, mock_lens_link);
    strcat (lensname_file, "all_lensname.txt");
    
    lensname = fopen (lensname_file, "r");
    
    // Make sure the lens name file is present.
    if (!lensname) {
        printf("Cannot find the lens name file!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count number of lines in the lens name file.
    lensfile_num = count_line (spot, lensname, file_count);
    
    // An array storing all file names.
    namelens *lens = malloc (lensfile_num * sizeof *lens);
    
    lensindex = 0;
    
    while (fscanf(lensname, "%d %s", &lens[lensindex].num, &lens[lensindex].name) == 2) {
        lensindex++;
    }
    
    fclose (lensname);
    
    
    // Which file to calculate.
    strcpy (cat_manual, argv[1]);
    cat_int = atoi (cat_manual);
    
    // Make sure the file index input does not exceed the number of files available.
    if (cat_int >= lensfile_num) {
        printf("File index input >= number of file available!\n");
        exit(1);
    }
    
    
    for (i = cat_int; i < cat_int + 1; i++) {
        
        printf("\nFile # %d\n", i);
        
        FILE *lensfile;
        char lensfile_string[150];
        
        strcpy (lensfile_string, mock_lens_link);
        strcat (lensfile_string, lens[i].name);
        
        lensfile = fopen (lensfile_string, "r");
        
        // Make sure the lens name file is present.
        if (!lensfile) {
            printf("Cannot find the lens file %s!\n", lens[i].name);
            exit(1);
        }
        
        // Open up the information of the lens galaxies.
        gallens *lgal = malloc (lens[i].num * sizeof *lgal);
        
        lgal_index = 0;
        
        // lens id, x, y, z, jackknife region.
        while (fscanf(lensfile, "%s %lf %lf %lf %*lf %*lf %*lf %*s %d", &lgal[lgal_index].id, &lgal[lgal_index].x, &lgal[lgal_index].y, &lgal[lgal_index].z, &lgal[lgal_index].jack) == 5) {
            //printf("%s\t%lf\t%lf\t%lf\t%d\n", lgal[lgal_index].id, lgal[lgal_index].x, lgal[lgal_index].y, lgal[lgal_index].z, lgal[lgal_index].jack);
            lgal_index++;
            if (lgal_index == 200) break;
        }
        
        fclose (lensfile);
        
        // Constructing the chaining mesh grids, in Mpc/h comoving.
        xmap_size = 250.;
        ymap_size = 250.;
        
        // Size of one mesh, and the number of meshes along x and y.
        xmesh_size = 5.;
        ymesh_size = 5.;
        
        xmesh_num = xmap_size / xmesh_size;
        ymesh_num = ymap_size / ymesh_size;
        
        // A map to store the particles in the corresponding mesh. Loop over x first, then y.
        double *x_cen = malloc (xmesh_num * sizeof (double));
        double *y_cen = malloc (ymesh_num * sizeof (double));
        
        for (k = 0; k < xmesh_num; k++) x_cen[k] = (k + 0.5) * xmesh_size;
        for (k = 0; k < ymesh_num; k++) y_cen[k] = (k + 0.5) * ymesh_size;
        
        // Create a map to store up which mesh do the lens galaxies belong to, first loop over x, then y. Another array is to store the index of the galaxies in the source list.
        int **lgal_count = malloc (ymesh_num * sizeof (int *));
        int **lgal_cum = malloc (ymesh_num * sizeof (int *));     // Cumulative number of lens galaxies.
        double ***lgal_coord = malloc (ymesh_num * sizeof (double **));     // Coordinates of the center of the grids.
        int **lgal_meshsort = malloc (lgal_index * sizeof (int *));
        
        for (k = 0; k < lgal_index; k++) {
            lgal_meshsort[k] = malloc (3 * sizeof (int));
            for (p = 0; p < 3; p++) {
                lgal_meshsort[k][p] = 0;    // [k][0] is mesh num in along x, [k][1] is mesh num in along y, [k][2] is index of the galaxy.
            }
        }
        
        for (k = 0; k < ymesh_num; k++) {
            lgal_count[k] = malloc (xmesh_num * sizeof (int));
            lgal_cum[k] = malloc (xmesh_num * sizeof (int));
            lgal_coord[k] = malloc (xmesh_num * sizeof (double *));
            for (p = 0; p < xmesh_num; p++) {
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
            x_ind = min(x_cen, xmesh_num, lgal[k].x);
            y_ind = min(y_cen, ymesh_num, lgal[k].y);

            lgal_meshsort[k][0] = x_ind;
            lgal_meshsort[k][1] = y_ind;
            lgal_meshsort[k][2] = k;
            lgal_count[y_ind][x_ind]++;
        }
        
        // Cumulative number of galaxies in each grid.
        for (k = 0; k < 1; k++) {
            lgal_cum[k][0] = lgal_count[0][0];
            for (p = 1; p < xmesh_num; p++) {
                lgal_cum[k][p] = lgal_cum[k][p - 1] + lgal_count[k][p];
            }
        }
        
        for (k = 1; k < ymesh_num; k++) {
            lgal_cum[k][0] = lgal_cum[k - 1][xmesh_num - 1] + lgal_count[k][0];
            for (p = 1; p < xmesh_num; p++) {
                lgal_cum[k][p] = lgal_cum[k][p - 1] + lgal_count[k][p];
            }
        }
        
        for (k = 0; k < ymesh_num; k++) {
            for (p = 0; p < xmesh_num; p++) {
                lgal_coord[k][p][0] = x_cen[p];
                lgal_coord[k][p][1] = y_cen[k];
            }
        }
        
        // Sort according to grid of x; if the same, then sort according to grid of y, and z; if still the same, then sort by galaxy index in the lens galaxy catalog.
        qsort(lgal_meshsort, lgal_index, sizeof lgal_meshsort[0], compare);
        
        /*for (k = 0; k < lgal_index; k++) {
            printf("%d\t%d\t%d\t%d\n", k, lgal_meshsort[k][0], lgal_meshsort[k][1], lgal_meshsort[k][2]);
        }*/
        
        // Multiple numbers of lens galaxies, starting from 0.
        strcpy (lgal_manual, argv[2]);
        strcpy (lgal_manum, argv[3]);
        lgal_eva = atoi(lgal_manum);
        lgal_num = atoi(lgal_manual) * lgal_eva;
        
        lgal_end = lgal_index;
        
        if ((atoi(lgal_manual) + 1) * lgal_eva < lgal_end) {
            lgal_end = (atoi(lgal_manual) + 1) * lgal_eva;
        }
        
        unsigned long **paircount = malloc ((lgal_end - lgal_num) * sizeof (unsigned long *));
        unsigned long ***paircount_jack = malloc ((lgal_end - lgal_num) * sizeof (unsigned long **));
        
        jack_num = 100;
        
        for (k = 0; k < (lgal_end - lgal_num); k++) {
            paircount[k] = malloc ((r_bincount + 1) * sizeof (unsigned long *));
            paircount_jack[k] = malloc ((r_bincount + 1) * sizeof (unsigned long **));
            for (m = 0; m < (r_bincount + 1); m++) {
                paircount[k][m] = 0;
                paircount_jack[k][m] = malloc (jack_num * sizeof (unsigned long));
                for (q = 0; q < jack_num; q++) {
                    paircount_jack[k][m][q] = 0;
                }
            }
        }
        
        gal_count = 0;

        int *jack_count = malloc (jack_num * sizeof (int));
                                                      
        for (k = 0; k < jack_num; k++) {
            jack_count[k] = 0;
        }
        
        // For a particular subset of lens
        for (k = lgal_num; k < lgal_end; k++) {
            
            // Count
            if ((lgal_end - k) % 1000 == 999) printf("%d / %d lens galaxies calculations done!\t#%d\n", (k - lgal_num), (lgal_end - lgal_num), k);
        
            // Chaining mesh.
            for (p = 0; p < ymesh_num; p++) { //ymesh_num
                for (q = 0; q < xmesh_num; q++) {  //xmesh_num
                    
                    // Distance between the lens and the grid.
                    gl_a_sep = pow((pow((lgal[k].x - lgal_coord[p][q][0]), 2) + pow((lgal[k].y - lgal_coord[p][q][1]), 2)), 0.5);
                    
                    for (m = 1; m < (r_bincount + 1); m++) {
                        
                        // Radius of the annulus, in comoving transverse distance (Mpc h^-1).
                        ann_in = r_width[m - 1] * Dc_ref * d2r / 3600;
                        ann_out = r_width[m] * Dc_ref * d2r / 3600;
                        
                        // Considering also the edging effect. Make sure there is at least 1 lens galaxy in the grid.
                        if ((gl_a_sep - ann_out <= 1.5 * xmesh_size) && (ann_in - gl_a_sep <= 1.5 * xmesh_size) && lgal_count[p][q] > 0) {
                            
                            // Create a dynamic array to store the index of the corresponding lens galaxies.
                            int *lgrid_sort = malloc (lgal_count[p][q] * sizeof (int));
                            int *l_ind = malloc (lgal_count[p][q] * sizeof (int));
                            
                            for (s = 0; s < lgal_count[p][q]; s++) {    // lgal_count[p][q]
                                
                                lgrid_sort[s] = lgal_cum[p][q] - lgal_count[p][q] + s;  // Number of concerned lens galaxy
                                l_ind[s] = lgal_meshsort[lgrid_sort[s]][2];     // The index of the concerned lens galaxy in the lens catalog.
                            
                                //printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", k, lgal_cum[p][q], lgal_count[p][q], s, lgrid_sort[s], l_ind[s], lgal_meshsort[lgrid_sort[s]][0], lgal_meshsort[lgrid_sort[s]][1]);
                                
                                // Make sure the code is reading in the lens galaxies in the desired grid.
                                if (lgal_meshsort[lgrid_sort[s]][0] != q || lgal_meshsort[lgrid_sort[s]][1] != p) {
                                    printf("Error! A lens galaxy is not in the concerned grid.\n");
                                    exit (1);
                                }
                                
                                // Not the lens galaxy itself
                                if (k != l_ind[s]) {
                                    // Lens-lens separation.
                                    ls_a_sep = pow((pow((lgal[k].x - lgal[l_ind[s]].x), 2) + pow((lgal[k].y - lgal[l_ind[s]].y), 2)), 0.5);
                                    
                                    // For lenses inside the annulus.
                                    if (ls_a_sep <= ann_out && ls_a_sep >= ann_in) {
                                        paircount[k - lgal_num][m]++;
                                        
                                        //printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\n", k, l_ind[s], lgal[k].x, lgal[k].y, lgal[l_ind[s]].x, lgal[l_ind[s]].y, ls_a_sep, m, ann_in, ann_out, r_width[m - 1] / 3600, r_width[m] / 3600);
                                        
                                        //getchar();

                                            
                                        //printf("%f\t%f\t%f\t%d\t%d\n", lgal[k].z, lgal[l_ind[s]].z, ls_z_sep, z_ind, paircount[k - lgal_num][m][z_ind]);
                                        //printf("%d\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%f\n", k, lgal[k].jack, m, lgrid_sort[s], lgal[l_ind[s]].x, lgal[l_ind[s]].y, lgal[l_ind[s]].jack, l_ind[s], ls_a_sep);
                                            
                                        // Jack-knifed sample.
                                        for (t = 0; t < jack_num; t++) {
                                            if (lgal[k].jack != t && lgal[l_ind[s]].jack != t) {
                                                paircount_jack[k - lgal_num][m][t]++;
                                                //printf("%d\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n\n", k, lgal[k].jack, m, lgrid_sort[s], lgal[l_ind[s]].x, lgal[l_ind[s]].y, lgal[l_ind[s]].jack, l_ind[s], paircount_jack[k - lgal_num][m][z_ind][t], t, z_ind);
                                                //getchar();
                                            }
                                        } // t loop, loop over the jack-knife regions.
                                            
                                    } // If the lens is inside the annulus.
                                } // If the lens galaxy does not count itself.
                            
                            } // s loop, loop over the galaxy counts in that mesh.
                            
                            // printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\n\n", m, gl_a_sep, ann_out, ann_in, gl_a_sep - ann_out, ann_in - gl_a_sep, 1.5 * xmesh_size);
                            
                            free (lgrid_sort);
                            free (l_ind);

                        } // if condition: if the mesh lies within the annulus.
                        
                    } // m loop, loop over radius bins.
                } // q loop, loop over x meshes
            } // p loop, loop over y meshes
            
            gal_count++;
            
            for (s = 0; s < jack_num; s++) {
                if (lgal[k].jack != s) jack_count[s]++;
            }

        } // k loop, loop over lens galaxies.
        
        
        
        
        unsigned long long *tot_paircount = malloc ((r_bincount + 1) * sizeof (unsigned long long));
        unsigned long long **tot_jackpaircount = malloc ((r_bincount + 1) * sizeof (unsigned long long *));
        
        for (k = 0; k < (r_bincount + 1); k++) {
            tot_paircount[k] = 0;
            tot_jackpaircount[k] = malloc (jack_num * sizeof (unsigned long long));
            for (p = 0; p < jack_num; p++) {
                tot_jackpaircount[k][p] = 0;
            }
        }
        
        for (k = 0; k < gal_count; k++) {
            for (m = 0; m < (r_bincount + 1); m++) {
                tot_paircount[m] += paircount[k][m];
                for (q = 0; q < jack_num; q++) {
                    tot_jackpaircount[m][q] += paircount_jack[k][m][q];
                }
            }
        }
        
        strncpy (W_name, lens[i].name, 7);
        W_name[7] = 0;
            
        char twoptDD_comfile[150];
        FILE *twoptDD_comname;
            
        strcpy (twoptDD_comfile, mock_2pt_link);
        strcat (twoptDD_comfile, W_name);
        strcat (twoptDD_comfile, "_DD_");
        strcat (twoptDD_comfile, lgal_manual);
        strcat (twoptDD_comfile, "_all.dat");
            
        twoptDD_comname = fopen (twoptDD_comfile, "w");
            
            
        fprintf(twoptDD_comname, "%d\t", gal_count);
            
        for (m = 1; m < r_bincount; m++) {
            fprintf(twoptDD_comname, "%llu\t", tot_paircount[m]);
        }
        fprintf(twoptDD_comname, "%llu\n", tot_paircount[r_bincount]);
            
        fclose (twoptDD_comname);
            
            
        // Jackknife files.
        char twoptDD_jackcomfile[150];
        FILE *twoptDD_jackcomname;
            
        strcpy (twoptDD_jackcomfile, mock_2pt_link);
        strcat (twoptDD_jackcomfile, W_name);
        strcat (twoptDD_jackcomfile, "_DD_");
        strcat (twoptDD_jackcomfile, lgal_manual);
        strcat (twoptDD_jackcomfile, "_jack.dat");
            
        twoptDD_jackcomname = fopen (twoptDD_jackcomfile, "w");
            
        for (u = 0; u < jack_num; u++) {
            fprintf(twoptDD_jackcomname, "%d\t%d\t", u, jack_count[u]);
            
            for (m = 1; m < r_bincount; m++) {
                fprintf(twoptDD_jackcomname, "%llu\t", tot_jackpaircount[m][u]);
            }
            fprintf(twoptDD_jackcomname, "%llu\n", tot_jackpaircount[r_bincount][u]);
        }
            
        fclose (twoptDD_jackcomname);
        
        //printf("\n");
        
        
        
        for (k = 0; k < ymesh_num; k++) {
            free(lgal_count[k]);
            free(lgal_cum[k]);
            for (p = 0; p < xmesh_num; p++) {
                free(lgal_coord[k][p]);
            }
            free(lgal_coord[k]);
        }
        
        for (k = 0; k < lgal_index; k++) {
            free (lgal_meshsort[k]);
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
        
        free (tot_jackpaircount);
        free (tot_paircount);
        free (jack_count);
        free (paircount);
        free (paircount_jack);
        free (lgal_meshsort);
        free (lgal_cum);
        free (lgal_count);
        free (lgal_coord);
        free (x_cen);
        free (y_cen);
        free (lgal);
    }
    
    

    free (lens);
    free (r_width);
    free (r_cen);
    
    return 0;
}