/* Calculates the overdensities of the mock galaxies.
 Updated: Jun 3, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "tool.h"
#include "cosmopar.h"

typedef struct {
    unsigned long id;
    double x, y, z;
} gal_m;

typedef struct {
    unsigned long id;
    double x, y, z;
} halo_m;

typedef struct {
    char id[20];
    unsigned long num;
} m_num;

int main(int argc, char *argv[]) {
    int i = 0, j = 0, k = 0, p = 0, q = 0, s = 0, m = 0, t = 0, xmesh_num = 0, ymesh_num = 0, gal_index = 0, halo_index = 0, spot = 0, file_count = 0, x_ind = 0, y_ind = 0;
    double xmesh_size = 0, ymesh_size = 0, lg_sep = 0, ln_sep = 0;
    const double xmap_size = 250.0, ymap_size = 250.0, cube_size = 250.0;
    const int r_bincount = 5;
    unsigned long halo_num = 0, gal_num = 0;
    char galname[150], haloname[150], halonumname[150], galnumname[150];
    FILE *galfile, *halofile, *halonumfile, *galnumfile;
    
    strcpy (galnumname, mock_lens_link);
    strcat (galnumname, "mock_lensname.txt");
    
    galnumfile = fopen (galnumname, "r");
    
    // Make sure the mock galaxy number file is present.
    if (!galnumfile) {
        printf("Cannot find the mock galaxy number file!\n");
        exit(1);
    }
    
    // Radius of the circles for overdensity calculations.
    double *den_r = malloc (5 * sizeof (double));
    
    for (i = 0; i < r_bincount; i++) {
        den_r[i] = (3. + (2 * i)) * h;
    }
    
    m_num *galcatnum = malloc (1 * sizeof *galcatnum);
    
    gal_index = 0;
    
    while (fscanf(galnumfile, "%lu %s", &galcatnum[gal_index].num, &galcatnum[gal_index].id) == 2) {
    }
    
    fclose (galnumfile);
    
    strcpy (galname, mock_lens_link);
    strcat (galname, galcatnum[0].id);
    
    galfile = fopen (galname, "r");
    
    // Make sure the mock galaxy catalog is present.
    if (!galfile) {
        printf("Cannot find the mock galaxy catalog!\n");
        exit(1);
    }
    
    // Count the number of mock galaxies.
    gal_num = galcatnum[0].num;
    
    // An array storing the mock galaxy info.
    gal_m *gal = malloc (gal_num * sizeof *gal);
    
    gal_index = 0;
    
    while (fscanf(galfile, "%lu %lf %lf %lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*d", &gal[gal_index].id, &gal[gal_index].x, &gal[gal_index].y, &gal[gal_index].z) == 4) {
        //printf("%lu\t%g\t%g\t%g\n", gal[gal_index].id, gal[gal_index].x, gal[gal_index].y, gal[gal_index].z);
        gal_index++;
        //if (gal_index == 100) break;
    }
    
    fclose (galfile);
    
    strcpy (halonumname, halonbr_link);
    strcat (halonumname, "halonum.dat");
    
    halonumfile = fopen (halonumname, "r");
    
    // Make sure the mock halo number file is present.
    if (!halonumfile) {
        printf("Cannot find the mock halo number file!\n");
        exit(1);
    }
    
    m_num *halocatnum = malloc (1 * sizeof *halocatnum);
    
    halo_index = 0;
    
    while (fscanf(halonumfile, "%lu %s", &halocatnum[halo_index].num, &halocatnum[halo_index].id) == 2) {
        //printf("%lu\t%s\n", halocatnum[halo_index].num, halocatnum[halo_index].id);
    }
    
    fclose (halonumfile);
    
    strcpy (haloname, halonbr_link);
    strcat (haloname, halocatnum[0].id);
    
    halofile = fopen (haloname, "r");
    
    // Make sure the mock halo catalog is present.
    if (!halofile) {
        printf("Cannot find the mock halo catalog!\n");
        exit(1);
    }
    
    // Count the number of mock halos.
    halo_num = halocatnum[0].num;
    
    // An array storing the mock halo info.
    halo_m *halo = malloc (halo_num * sizeof *halo);
    
    halo_index = 0;
    
    while (fscanf(galfile, "%lu %lf %lf %lf %*lf", &halo[halo_index].id, &halo[halo_index].x, &halo[halo_index].y, &halo[halo_index].z) == 4) {
        //printf("%s\t%g\t%g\t%g\n", halo[halo_index].id, halo[halo_index].x, halo[halo_index].y, halo[halo_index].z);
        halo_index++;
        //if (halo_index == 100) break;
    }
    
    fclose (halofile);

    
    // Chaining mesh grids, in Mpc h^-1.
    xmesh_size = 2.5;
    ymesh_size = 2.5;
    
    // Number of meshes along x and y.
    xmesh_num = xmap_size / xmesh_size;
    ymesh_num = ymap_size / ymesh_size;
    
    // A map to store the particles in the corresponding mesh. Loop over x first, then y.
    double *x_cen = malloc (xmesh_num * sizeof (double));
    double *y_cen = malloc (ymesh_num * sizeof (double));
    
    for (i = 0; i < xmesh_num; i++) x_cen[i] = (i + 0.5) * xmesh_size;
    for (i = 0; i < ymesh_num; i++) y_cen[i] = (i + 0.5) * ymesh_size;
    
    
    // Create a map to store up which mesh do the halos belong to, first loop over x, then y. Another array is to store the index of the halos in the source list.
    int **halo_count = malloc (ymesh_num * sizeof (int *));
    int **halo_cum = malloc (ymesh_num * sizeof (int *));     // Cumulative number of halos.
    double ***halo_coord = malloc (ymesh_num * sizeof (double **));     // Coordinates of the center of the grids.
    int **halo_meshsort = malloc (halo_index * sizeof (int *));
    
    for (k = 0; k < halo_index; k++) {
        halo_meshsort[k] = malloc (3 * sizeof (int));
        for (p = 0; p < 3; p++) {
            halo_meshsort[k][p] = 0;    // [k][0] is mesh num in along x, [k][1] is mesh num in along y, [k][2] is index of the halo.
        }
    }
    
    for (k = 0; k < ymesh_num; k++) {
        halo_count[k] = malloc (xmesh_num * sizeof (int));
        halo_cum[k] = malloc (xmesh_num * sizeof (int));
        halo_coord[k] = malloc (xmesh_num * sizeof (double *));
        for (p = 0; p < xmesh_num; p++) {
            halo_count[k][p] = 0;
            halo_cum[k][p] = 0;
            halo_coord[k][p] = malloc (2 * sizeof (double));
            for (q = 0; q < 2; q++) {
                halo_coord[k][p][q] = 0;
            }
        }
    }
    
    // Store up the number of counts in the mesh grid
    for (k = 0; k < halo_index; k++) {
        x_ind = min(x_cen, xmesh_num, halo[k].x);
        y_ind = min(y_cen, ymesh_num, halo[k].y);
        
        halo_meshsort[k][0] = x_ind;
        halo_meshsort[k][1] = y_ind;
        halo_meshsort[k][2] = k;
        halo_count[y_ind][x_ind]++;
        //printf("%d\t%lf\t%d\t%lf\t%d\t%d\n", halo_meshsort[k][2], halo[k].x, halo_meshsort[k][0], halo[k].y, halo_meshsort[k][1], halo_count[y_ind][x_ind]);
    }
    
    // Cumulative number of halos in each grid.
    for (k = 0; k < 1; k++) {
        halo_cum[k][0] = halo_count[0][0];
        for (p = 1; p < xmesh_num; p++) {
            halo_cum[k][p] = halo_cum[k][p - 1] + halo_count[k][p];
        }
    }
    
    for (k = 1; k < ymesh_num; k++) {
        halo_cum[k][0] = halo_cum[k - 1][xmesh_num - 1] + halo_count[k][0];
        for (p = 1; p < xmesh_num; p++) {
            halo_cum[k][p] = halo_cum[k][p - 1] + halo_count[k][p];
        }
    }
    
    for (k = 0; k < ymesh_num; k++) {
        for (p = 0; p < xmesh_num; p++) {
            halo_coord[k][p][0] = x_cen[p];
            halo_coord[k][p][1] = y_cen[k];
        }
    }
    
    // Sort according to grid of x; if the same, then sort according to grid of y, and z; if still the same, then sort by galaxy index in the halo catalog.
    qsort(halo_meshsort, halo_index, sizeof halo_meshsort[0], compare);
    
    unsigned long **numden = malloc (gal_index * sizeof (unsigned long *));
    
    for (k = 0; k < gal_index; k++) {
        numden[k] = malloc (r_bincount * sizeof (unsigned long));
        for (m = 0; m < r_bincount; m++) {
            numden[k][m] = 0;
        }
    }
    
    // Looping over the galaxies.
    for (p = 0; p < gal_index; p++) {       // gal_index
        
        if (p % 1000 == 0) printf("Galaxy #%d out of %d\n", p, gal_index);
        
        for (q = 0; q < ymesh_num; q++) {   // ymesh_num
            for (s = 0; s < xmesh_num; s++) {   //xmesh_num
                
                // Distance between the galaxy and the grid.
                lg_sep = dist_periodic(gal[p].x, gal[p].y, halo_coord[q][s][0], halo_coord[q][s][1], cube_size);
                //printf("%lf\t%lf\t%lf\t%lf\t%lf\n", gal[p].x, gal[p].y, halo_coord[q][s][0], halo_coord[q][s][1], lg_sep);
                
                for (m = 0; m < r_bincount; m++) {
                    
                    if (lg_sep - den_r[m] <= 1.5 * xmesh_size && halo_count[q][s] > 0) {
                        //printf("%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%lf\t%d\n", gal[p].x, gal[p].y, halo_coord[q][s][0], halo_coord[q][s][1], lg_sep, m, den_r[m], halo_count[q][s]);
                        // Create a dynamic array to store the index of the corresponding halos.
                        int *grid_sort = malloc (halo_count[q][s] * sizeof (int));
                        int *ind = malloc (halo_count[q][s] * sizeof (int));
                        
                        for (t = 0; t < halo_count[q][s]; t++) {
                            grid_sort[t] = halo_cum[q][s] - halo_count[q][s] + t;  // Number of concerned halo
                            ind[t] = halo_meshsort[grid_sort[t]][2];     // The index of the concerned halo in the halo catalog.
                            
                            // Make sure the code is reading in the halos in the desired grid.
                            if (halo_meshsort[grid_sort[t]][0] != s || halo_meshsort[grid_sort[t]][1] != q) {
                                printf("Error! A halo is not in the concerned grid.\n");
                                exit (1);
                            }
                            
                            // Not the mock galaxy itself.
                            if (gal[p].id != halo[ind[t]].id) {
                                // Mock galaxy and neighboring galaxy separation.
                                ln_sep = dist_periodic(gal[p].x, gal[p].y, halo[ind[t]].x, halo[ind[t]].y, cube_size);
                                
                                if (ln_sep <= den_r[m]) {
                                    numden[p][m]++;
                                    //printf("%lu\t%d\t%d\t%d\t%d\t%d\t%lu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lu\n", gal[p].id, t, halo_count[q][s], halo_cum[q][s], grid_sort[t], ind[t], halo[ind[t]].id, gal[p].x, gal[p].y, halo[ind[t]].x, halo[ind[t]].y, ln_sep, den_r[m], numden[p][m]);
                                }
                            
                                //printf("%lu\t%d\t%d\t%d\t%d\t%d\t%lu\t%lf\t%lf\t%lf\t%lf\t%lf\n", gal[p].id, t, halo_count[q][s], halo_cum[q][s], grid_sort[t], ind[t], halo[ind[t]].id, gal[p].x, gal[p].y, halo[ind[t]].x, halo[ind[t]].y, ln_sep);
                            }
                        }
                        
                        free (grid_sort);
                        free (ind);
                    }
                }   // m-loop, different radii.
            }   // s-loop, same y-coords
        }   // q-loop
    }   // p-loop, mock galaxies.
    
    char W_name[150];
    FILE *W_file;
    
    strcpy (W_name, extmock_link);
    strcat (W_name, "sm_overdensity.dat");
    
    W_file = fopen(W_name, "w");
    
    for (p = 0; p < gal_index; p++) {
        fprintf(W_file, "%lu\t", gal[p].id);
        for (t = 0; t < r_bincount; t++) {
            fprintf(W_file, "%lu\t", numden[p][t]);
        }
        fprintf(W_file, "\n");
    }
    
    fclose (W_file);
    
    for (k = 0; k < ymesh_num; k++) {
        free(halo_count[k]);
        free(halo_cum[k]);
        for (p = 0; p < xmesh_num; p++) {
            free (halo_coord[k][p]);
        }
        free (halo_coord[k]);
    }
    
    for (k = 0; k < halo_index; k++) {
        free (halo_meshsort[k]);
    }
    
    for (k = 0; k < gal_index; k++) {
        free (numden[k]);
    }
    
    free (numden);
    free (halo_meshsort);
    free (halo_cum);
    free (halo_count);
    free (halo_coord);
    free (x_cen);
    free (y_cen);
    free (gal);
    free (halo);
    free (den_r);
    
    return 0;
}


