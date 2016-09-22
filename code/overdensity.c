/* This program calculates the overdensities of the lens galaxies at different scales, 3,5,7,9 and 11 h^-1 Mpc.
 Input: Field # (argv[1]),  multiple of n lens (argv[2]), value of n (argv[3])
 Updated: May 5, 2016
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
    char id[20];
    double ra, dec, z;
} gallens;



int main(int argc, char *argv[]) {
    int i = 0, j = 0, k = 0, p = 0, q = 0, m = 0, s = 0, z_bincount = 0, r_bincount = 0, spot = 0, file_count = 0, field_num = 0, fieldindex = 0, field_int = 0, lens_num = 0, nbr_num = 0, lensindex = 0, nbrindex = 0, lgal_index = 0, ngal_index = 0, ramesh_num = 0, decmesh_num = 0, r_ind = 0, d_ind = 0, lgal_eva = 0, lgal_num = 0, lgal_end = 0, z_range = 0;
    double z_min = 0, z_width = 0, r_min = 0, r_width = 0, ramap_size = 0, decmap_size = 0, ramesh_size = 0, decmesh_size = 0, gl_a_sep = 0, ls_a_sep = 0;
    char field_file[150], field_manual[150], lens_file[150], nbr_file[150], write_file[150], lgal_manual[10], lgal_manum[10];
    FILE *fieldname, *lensname, *nbrname, *writename;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 4) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    // Input which wide field
    strcpy (field_manual, argv[1]);
    field_int = atoi(field_manual);
    
    // Pre-compute the angular separations corresponding to different physical scales at different redshift slices.
    z_bincount = 21;
    z_min = 0.2;
    z_width = 0.01;
    
    double *z_bin = malloc ((z_bincount) * sizeof (double));    // Redshift bins.
    double *Da_bin = malloc ((z_bincount) * sizeof (double));   // Angular diameter distance bins.
    double **ang_bin = malloc ((z_bincount) * sizeof (double *));    // Angles corresponding to the scales at that redshift, in arcsecs.
    
    // Create an array to store up the physical radius of the circles for overdensity calculations.
    r_bincount = 5;
    r_min = 3.;
    r_width = 2.;
    
    double *rad_bin = malloc ((r_bincount) * sizeof (double));  // Radius bins.
    
    for (j = 0; j < r_bincount; j++) {
        rad_bin[j] = r_min + j * r_width;
    }
    
    for (i = 0; i < z_bincount; i++) {
        z_bin[i] = z_min + i * z_width;
        Da_bin[i] = Da(z_bin[i]);
        ang_bin[i] = malloc ((r_bincount) * sizeof (double));
        for (j = 0; j < r_bincount; j++) {
            ang_bin[i][j] = rad_bin[j] / Da_bin[i] * 180 * 3600 / M_PI;     // in arcsecs
        }
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
    
    strcpy (lens_file, gal_link);
    strcat (lens_file, "lensname.txt");
    strcpy (nbr_file, gal_link);
    strcat (nbr_file, "neighname.txt");
    
    lensname = fopen (lens_file, "r");
    nbrname = fopen (nbr_file, "r");
    
    // Make sure the lens name file is present.
    if (!lensname) {
        printf("Cannot find the lens name file!\n");
        exit(1);
    }
    
    // Make sure the neighbor name file is present.
    if (!nbrname) {
        printf("Cannot find the neighbor name file!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count numer of lines in the lens name file.
    lens_num = count_line (spot, lensname, file_count);
    
    // An array storing all file names.
    namelens *lens = malloc (lens_num * sizeof *lens);
    
    lensindex = 0;
    
    while (fscanf(lensname, "%d %s", &lens[lensindex].num, &lens[lensindex].name) == 2) {
        lensindex++;
    }
    
    fclose (lensname);
    
    spot = 0;
    file_count = -1;
    
    // Count numer of lines in the neighbor name file.
    nbr_num = count_line (spot, nbrname, file_count);
    
    // An array storing all file names.
    namelens *nbr = malloc (nbr_num * sizeof *nbr);
    
    nbrindex = 0;
    
    while (fscanf(nbrname, "%d %s", &nbr[nbrindex].num, &nbr[nbrindex].name) == 2) {
        nbrindex++;
    }
    
    fclose (nbrname);
    
    
    for (i = field_int - 1; i < field_int; i++) {
        
        printf("\nField W%d\n", i + 1);
        
        FILE *lensfile, *nbrfile;
        char lens_str[150], nbr_str[150];
        
        strcpy (lens_str, gal_link);
        strcat (lens_str, lens[i].name);
        
        lensfile = fopen (lens_str, "r");
        
        // Make sure the lens name file is present.
        if (!lensfile) {
            printf("Cannot find the lens file %s!\n", lens[i].name);
            exit(1);
        }
        
        // Open up the information of the lens galaxies.
        gallens *lgal = malloc (lens[i].num * sizeof * lgal);
        
        lgal_index = 0;
        
        // lens id, ra, dec and redshift
        while (fscanf(lensfile, "%s %*s %*lf %*lf %lf %lf %*lE %*lE %*lf %*lf %*d %*lE %*lE %*lf %*lE %*lE %*lf %*lf %*lf %*lf %*lE %*d %*lf %*d %lf %*lf %*lf %*lf %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*d %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lE %*lf %*lE %*lf %*lf %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE", &lgal[lgal_index].id, &lgal[lgal_index].ra, &lgal[lgal_index].dec, &lgal[lgal_index].z) == 4) {
            lgal_index++;
            //if (lgal_index == 3) break;
        }
        
        fclose (lensfile);
        
        
        strcpy (nbr_str, gal_link);
        strcat (nbr_str, nbr[i].name);
        
        nbrfile = fopen (nbr_str, "r");
        
        // Make sure the neighbor name file is present.
        if (!nbrfile) {
            printf("Cannot find the neighbor file %s!\n", nbr[i].name);
            exit(1);
        }
        
        // Open up the information of the neighbor galaxies.
        gallens *ngal = malloc (nbr[i].num * sizeof * ngal);
        
        ngal_index = 0;
        
        // Neighbor id, ra, dec.
        while (fscanf(nbrfile, "%s %*s %lf %lf %*lE %*lE %*lE %*lE %*lE %*lE %*lE", &ngal[ngal_index].id, &ngal[ngal_index].ra, &ngal[ngal_index].dec) == 3) {
            ngal_index++;
            //if (ngal_index == 100) break;
        }
        
        fclose (nbrfile);
        
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
        
        for (j = 0; j < ramesh_num; j++) {
            ra_center[j] = fpos[i].ra_low + (j + 0.5) * ramesh_size;
        }
        for (j = 0; j < decmesh_num; j++) {
            dec_center[j] = fpos[i].dec_low + (j + 0.5) * decmesh_size;
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
            r_ind = min(ra_center, ramesh_num, lgal[k].ra);
            d_ind = min(dec_center, decmesh_num, lgal[k].dec);
            lgal_meshsort[k][0] = r_ind;        // mesh num along RA
            lgal_meshsort[k][1] = d_ind;        // mesh num along DEC
            lgal_meshsort[k][2] = k;            // index of the galaxy
            lgal_count[d_ind][r_ind]++;         // Counts of lens in each mesh.
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
        
        // Create the mesh map.
        for (k = 0; k < decmesh_num; k++) {
            for (p = 0; p < ramesh_num; p++) {
                lgal_coord[k][p][0] = ra_center[p];
                lgal_coord[k][p][1] = dec_center[k];
            }
        }
        
        // Sort according to grid of RA; if the same, then sort according to grid of DEC; if still the same, then sort by galaxy index in the lens galaxy catalog.
        qsort(lgal_meshsort, lgal_index, sizeof lgal_meshsort[0], compare);
        
        // Create a map to store up which mesh do the neighbor galaxies belong to, first loop over RA, then DEC. Another array is to store the index of the galaxies in the source list.
        int **ngal_count = malloc (decmesh_num * sizeof (int *));
        int **ngal_cum = malloc (decmesh_num * sizeof (int *));     // Cumulative number of lens galaxies.
        int **ngal_meshsort = malloc (ngal_index * sizeof (int *));
        
        for (k = 0; k < ngal_index; k++) {
            ngal_meshsort[k] = malloc (3 * sizeof (int));
            for (p = 0; p < 3; p++) {
                ngal_meshsort[k][p] = 0;    // [k][0] is mesh num in along RA, [k][1] is mesh num in along DEC, [k][2] is index of the galaxy.
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            ngal_count[k] = malloc (ramesh_num * sizeof (int));
            ngal_cum[k] = malloc (ramesh_num * sizeof (int));
            for (p = 0; p < ramesh_num; p++) {
                ngal_count[k][p] = 0;
                ngal_cum[k][p] = 0;
            }
        }
        
        // Store up the number of counts in the mesh grid
        for (k = 0; k < ngal_index; k++) {
            r_ind = min(ra_center, ramesh_num, ngal[k].ra);
            d_ind = min(dec_center, decmesh_num, ngal[k].dec);
            ngal_meshsort[k][0] = r_ind;
            ngal_meshsort[k][1] = d_ind;
            ngal_meshsort[k][2] = k;
            ngal_count[d_ind][r_ind]++;
        }
        
        // Cumulative number of galaxies in each grid.
        for (k = 0; k < 1; k++) {
            ngal_cum[k][0] = ngal_count[0][0];
            for (p = 1; p < ramesh_num; p++) {
                ngal_cum[k][p] = ngal_cum[k][p - 1] + ngal_count[k][p];
            }
        }
        
        for (k = 1; k < decmesh_num; k++) {
            ngal_cum[k][0] = ngal_cum[k - 1][ramesh_num - 1] + ngal_count[k][0];
            for (p = 1; p < ramesh_num; p++) {
                ngal_cum[k][p] = ngal_cum[k][p - 1] + ngal_count[k][p];
            }
        }
        
        // Sort according to grid of RA; if the same, then sort according to grid of DEC; if still the same, then sort by galaxy index in the random galaxy catalog.
        qsort(ngal_meshsort, ngal_index, sizeof ngal_meshsort[0], compare);
        
        // Multiple numbers of lens galaxies, starting from 0.
        strcpy (lgal_manual, argv[2]);
        strcpy (lgal_manum, argv[3]);
        lgal_eva = atoi(lgal_manum);
        lgal_num = atoi(lgal_manual) * lgal_eva;
        
        lgal_end = lgal_index;
        
        if ((atoi(lgal_manual) + 1) * lgal_eva < lgal_end) {
            lgal_end = (atoi(lgal_manual) + 1) * lgal_eva;
        }
        
        // Neighbor count double arrays for each lens galaxy.
        unsigned long **nbr_count = malloc (lgal_index * sizeof (unsigned long *));
        
        for (k = 0; k < lgal_index; k++) {
            nbr_count[k] = malloc ((r_bincount) * sizeof (unsigned long));
            for (p = 0; p < r_bincount; p++) {
                nbr_count[k][p] = 0;
            }
        }
        
        for (k = lgal_num; k < lgal_end; k++) {
            
            // Count
            if ((lgal_end - k) % 1000 == 999) printf("%d / %d lens galaxies calculations done!\t#%d\n", (k - lgal_num), (lgal_end - lgal_num), k);
            
            // Redshift range of the lens galaxy.
            z_range = min(z_bin, z_bincount, lgal[k].z);
            
            // Boundary of the angular separations at the lens' redshift.
            double *ang_bound = malloc (r_bincount * sizeof (double));
            
            for (j = 0; j < r_bincount; j++) {
                ang_bound[j] = ang_bin[z_range][j];
            }
            
            
            // Chaining mesh.
            for (p = 0; p < decmesh_num; p++) { //decmesh_num
                for (q = 0; q < ramesh_num; q++) {  //ramesh_num
                    
                    // Angular separation between the lens and the grid.
                    gl_a_sep = ang_sep (lgal[k].ra, lgal[k].dec, lgal_coord[p][q][0], lgal_coord[p][q][1]);
                    
                    //printf("%d\t%lf\t%lf\t%lf\t%lf\t%lE\n", k, lgal[k].ra, lgal[k].dec, lgal_coord[p][q][0], lgal_coord[p][q][1], gl_a_sep);
                    
                    for (m = 0; m < r_bincount; m++) {  // m < r_bincount
                        // Considering also the edging effect. Make sure there is at least 1 random galaxy in the grid.
                        // 5400 = 1.5 * 3600 (3600 for arcsec in 1 deg)
                        if ((gl_a_sep - ang_bound[m] <= ramesh_size * 5400) && ngal_count[p][q] > 0) {
                            
                            // printf("%d\t%lf\t%lf\t%lf\t%lf\t%lE\t%lE\t%lE\n", k, lgal[k].ra, lgal[k].dec, lgal_coord[p][q][0], lgal_coord[p][q][1], gl_a_sep, ang_bound[m], gl_a_sep - ang_bound[m]);
                            
                            // Create a dynamic array to store the index of the corresponding neighboring galaxies.
                            int *ngrid_sort = malloc (ngal_count[p][q] * sizeof (int));
                            int *l_ind = malloc (ngal_count[p][q] * sizeof (int));
                            
                            for (s = 0; s < ngal_count[p][q]; s++) {    // lgal_count[p][q]
                                
                                ngrid_sort[s] = ngal_cum[p][q] - ngal_count[p][q] + s;  // Sorted index of concerned neighboring galaxy
                                l_ind[s] = ngal_meshsort[ngrid_sort[s]][2];     // The index of the concerned neighboring galaxy in the lens catalog.
                                
                                /*
                                printf("%d\t%s\t%g\t%g\t%g\t%g\t%g\n", k, lgal[k].id, lgal[k].ra, lgal[k].dec, lgal[k].z, Da(lgal[k].z), ang_bound[m]);
                                printf("%d\t%d\t%g\t%g\t%d\t%d\t%d\t%d\t%d\n\n", q, p, lgal_coord[p][q][0], lgal_coord[p][q][1], ngal_cum[p][q], ngal_count[p][q], s, ngrid_sort[s], ngal_meshsort[ngrid_sort[s]][2]);
                                */
                                
                                // Make sure the code is reading in the lens galaxies in the desired grid.
                                if (ngal_meshsort[ngrid_sort[s]][0] != q || ngal_meshsort[ngrid_sort[s]][1] != p) {
                                    printf("Error! A neighboring galaxy is not in the concerned grid.\n");
                                    exit (1);
                                }
                                
                                // Lens-neighbor separation.
                                ls_a_sep = ang_sep (lgal[k].ra, lgal[k].dec, ngal[l_ind[s]].ra, ngal[l_ind[s]].dec);
                                
                                // For the neighbor inside the circle.
                                if (ls_a_sep < ang_bound[m]) {
                                    
                                    // Avoid counting the same galaxy in the lens and neighbor catalogs.
                                    if (strcmp(lgal[k].id, ngal[l_ind[s]].id) != 0) {
                                    
                                        //printf("%d\t%s\t%g\t%g\t%g\t%g\n", k, lgal[k].id, lgal[k].ra, lgal[k].dec, lgal[k].z, Da(lgal[k].z));
                                        //printf("%d\t%s\t%g\t%g\t%g\t%g\n\n", l_ind[s], ngal[l_ind[s]].id, ngal[l_ind[s]].ra, ngal[l_ind[s]].dec, ls_a_sep, ang_bound[m]);
                                    
                                        nbr_count[k][m]++;
                                    }
                                }   // if loop, if the neighbor is inside the circle.
                                
                            } // s-loop, for the neighboring galaxies in that mesh.
                            
                            free (ngrid_sort);
                            free (l_ind);
                        
                        }   // if-loop,
                        
                    }   // m-loop, r_bincount
                    
                }   // q-loop, ramesh_num
            }       // p-loop, decmesh_num
            
            free (ang_bound);
            
        }   // k-loop, for each lens galaxy
        
        // Write the overdensity files.
        strcpy (write_file, overden_part_link);
        strcat (write_file, "W");
        strcat (write_file, field_manual);
        strcat (write_file, "/overdensity_");
        strcat (write_file, lgal_manual);
        strcat (write_file, ".dat");
        
        writename = fopen (write_file, "w");
        
        for (k = lgal_num; k < lgal_end; k++) {
            fprintf (writename, "%s\t", lgal[k].id);
            for (m = 0; m < r_bincount; m++) {
                fprintf (writename, "%lu\t", nbr_count[k][m]);
            }
            fprintf (writename, "\n");
        }
        
        fclose (writename);
        
        for (k = 0; k < lgal_index; k++) {  //(lgal_end - lgal_num)
            free (nbr_count[k]);
        }

        for (k = 0; k < decmesh_num; k++) {
            free(lgal_count[k]);
            free(lgal_cum[k]);
            free(ngal_count[k]);
            free(ngal_cum[k]);
            for (p = 0; p < ramesh_num; p++) {
                free(lgal_coord[k][p]);
            }
            free(lgal_coord[k]);
        }
        
        for (k = 0; k < lgal_index; k++) {
            free (lgal_meshsort[k]);
            free (ngal_meshsort[k]);
        }
        
        free (nbr_count);
        free (lgal_meshsort);
        free (lgal_cum);
        free (lgal_count);
        free (lgal_coord);
        free (ngal_meshsort);
        free (ngal_cum);
        free (ngal_count);
        free (ra_center);
        free (dec_center);
        free (lgal);
        free (ngal);
    }   // i-loop, wide fields
    
    
    
    
    
    for (i = 0; i < r_bincount; i++) {
        free (ang_bin[i]);
    }
    
    free (lens);
    free (nbr);
    free (fpos);
    free (ang_bin);
    free (z_bin);
    free (rad_bin);
    free (Da_bin);
    
    return 0;
}


