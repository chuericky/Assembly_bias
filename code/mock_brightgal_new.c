/* This program removes satellites by finding if there is any more massive galaxies within some separation, by comparing with the FULL catalog.
    - Perpendicular to L.O.S.: D_0
    Input: D_0 (perp to L.O.S.) (argv[1]), color (argv[2]), star mass ratio (argv[3])

    Updated: Dec 15, 2015
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
    char id[20], host_id[20];
    double x, y, z, vx, vy, vz, sm, g_r, m_vir, vpeak, m_host;
    int jack, sat;
} gal_c;

typedef struct {
    char id[20];
    double x, y, z, sm;
} gal_f;

typedef struct {
    char id[120];
    int num;
} galcolor;

typedef struct {
    char id[120];
    int num;
} fullcat;


int main(int argc, char *argv[]) {
    int i = 0, j = 0, k = 0, p = 0, q = 0, s = 0, t = 0, u = 0, w = 0, par_negcount = 0, partner_count = 0, spot = 0, file_count = 0, color_num = 0, color_index = 0, lgal_index = 0, clus_id_num = 0, xmesh_num = 0, ymesh_num = 0, fd_count = 0, p_restore = 0, max_lum_ind = 0, x_ind = 0, y_ind = 0, par_counter = 0, color_int = 0, sat_count = 0, iso_count = 0, full_index = 0, fgal_index = 0;
    double D_0 = 0, xmap_size = 0, ymap_size = 0, xmesh_size = 0, ymesh_size = 0, gl_a_sep = 0, ll_a_sep = 0, cube_size = 250.0, sm = 0;
    char D_b_man[5], colorchecklink[150], color_man[5], fulllink[250], star_mass[5];
    FILE *colorcheckfile, *fullfile;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 4) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    // Transverse comoving distance in the projected plane (xy).
    strcpy (D_b_man, argv[1]);
    D_0 = atof (D_b_man);
    
    // Color
    strcpy (color_man, argv[2]);
    color_int = atoi (color_man);
    
    // More massive
    strcpy (star_mass, argv[3]);
    sm = log10(atof (star_mass));
    
    strcpy (colorchecklink, mock_lens_link);
    strcat (colorchecklink, "mock_lensname.txt");
    
    colorcheckfile = fopen (colorchecklink, "r");
    
    // Make sure the color list file is present.
    if (!colorcheckfile) {
        printf("Cannot find the color list!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count number of color list files.
    color_num = count_line (spot, colorcheckfile, file_count);
    
    // An array storing the color list info.
    galcolor *gcolor = malloc (color_num * sizeof *gcolor);
    
    color_index = 0;
    
    while (fscanf(colorcheckfile, "%d %s", &gcolor[color_index].num, &gcolor[color_index].id) == 2) {
        // printf("%d\t%s\n", gcolor[color_index].num, gcolor[color_index].id);
        color_index++;
    }
    
    fclose (colorcheckfile);
    
    strcpy (fulllink, mock_lens_link);
    strcat (fulllink, "mock_lensname.txt");
    
    fullfile = fopen (fulllink, "r");
    
    // Make sure the full list file is present.
    if (!fullfile) {
        printf("Cannot find the full list!\n");
        exit(1);
    }
    
    // An array storing the full list info.
    fullcat *full = malloc (color_num * sizeof *gcolor);
    
    full_index = 0;
    
    while (fscanf(fullfile, "%d %s", &full[full_index].num, &full[full_index].id) == 2) {
        //printf("%d\t%s\n", full[full_index].num, full[full_index].id);
        full_index++;
    }
    
    fclose (fullfile);
    
    // Chaining mesh grids, in Mpc h^-1.
    xmap_size = 250.;
    ymap_size = 250.;
    
    xmesh_size = 5.;
    ymesh_size = 5.;
    
    // Number of meshes along x and y.
    xmesh_num = xmap_size / xmesh_size;
    ymesh_num = ymap_size / ymesh_size;
    
    // A map to store the particles in the corresponding mesh. Loop over x first, then y.
    double *x_cen = malloc (xmesh_num * sizeof (double));
    double *y_cen = malloc (ymesh_num * sizeof (double));
    
    for (k = 0; k < xmesh_num; k++) x_cen[k] = (k + 0.5) * xmesh_size;
    for (k = 0; k < ymesh_num; k++) y_cen[k] = (k + 0.5) * ymesh_size;
    
    
    for (i = color_int; i < color_int + 1; i++) {   // color_index
        printf("%s begins!\n", gcolor[i].id);
        
        FILE *galfile, *f_file;
        char gal_string[150], full_string[150];
        
        strcpy (gal_string, mock_lens_link);
        strcat (gal_string, gcolor[i].id);
        
        galfile = fopen (gal_string, "r");
        
        // Make sure the galaxy color file exists.
        if (!galfile) {
            printf("Cannot find the galaxy color file %s!\n", gcolor[i].id);
            exit (1);
        }
        
        // Open up the information for the mock galaxies.
        gal_c *lgal = malloc (gcolor[i].num * sizeof *lgal);
        
        lgal_index = 0;
        
        // lens id, x, y, z, jackknife region.
        while (fscanf(galfile, "%s %lf %lf %lf %lf %lf %lf %lE %lf %lf %lf %lE %s", &lgal[lgal_index].id, &lgal[lgal_index].x, &lgal[lgal_index].y, &lgal[lgal_index].z, &lgal[lgal_index].vx, &lgal[lgal_index].vy, &lgal[lgal_index].vz, &lgal[lgal_index].m_vir, &lgal[lgal_index].vpeak, &lgal[lgal_index].sm, &lgal[lgal_index].g_r, &lgal[lgal_index].m_host, &lgal[lgal_index].host_id) == 13) {
            lgal[lgal_index].sat = 0;
            //printf("%d\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lE\t%lf\t%lf\t%lf\t%lE\t%s\n", lgal_index, lgal[lgal_index].id, lgal[lgal_index].x, lgal[lgal_index].y, lgal[lgal_index].z, lgal[lgal_index].vx, lgal[lgal_index].vy, lgal[lgal_index].vz, lgal[lgal_index].m_vir, lgal[lgal_index].vpeak, lgal[lgal_index].Mr, lgal[lgal_index].g_r, lgal[lgal_index].m_host, lgal[lgal_index].host_id);
            lgal_index++;
            //if (lgal_index == 200) break;
        }
        
        fclose (galfile);
        
        strcpy (full_string, mock_lens_link);
        strcat (full_string, full[i].id);
        
        f_file = fopen (full_string, "r");
        
        // Make sure the full galaxy file exists.
        if (!f_file) {
            printf("Cannot find the full galaxy file %s!\n", full[i].id);
            exit (1);
        }
        
        // Open up the information for the mock galaxies.
        gal_f *fgal = malloc (full[i].num * sizeof *fgal);
        
        fgal_index = 0;
        
        // full galaxy id, x, y, z, Mr.
        while (fscanf(f_file, "%s %lf %lf %lf %*lf %*lf %*lf %*lE %*lf %lf %*lf %*lE %*s", &fgal[fgal_index].id, &fgal[fgal_index].x, &fgal[fgal_index].y, &fgal[fgal_index].z, &fgal[fgal_index].sm) == 5) {
        //while (fscanf(f_file, "%s %*lE %lf %lf %lf %lf", &fgal[fgal_index].id, &fgal[fgal_index].x, &fgal[fgal_index].y, &fgal[fgal_index].z, &fgal[fgal_index].Mr) == 5) {
            //printf("%d\t%s\t%lf\t%lf\t%lf\t%lf\n", fgal_index, fgal[fgal_index].id, fgal[fgal_index].x, fgal[fgal_index].y, fgal[fgal_index].z, fgal[fgal_index].Mr);
            fgal_index++;
            //if (fgal_index == 200) break;
        }
        
        fclose (f_file);
        
        // Create a map to store up which mesh do the lens galaxies belong to, first loop over x, then y. Another array is to store the index of the galaxies in the source list.
        int **fgal_count = malloc (ymesh_num * sizeof (int *));
        int **fgal_cum = malloc (ymesh_num * sizeof (int *));     // Cumulative number of lens galaxies.
        double ***fgal_coord = malloc (ymesh_num * sizeof (double **));     // Coordinates of the center of the grids.
        int **fgal_meshsort = malloc (fgal_index * sizeof (int *));
        
        for (k = 0; k < fgal_index; k++) {
            fgal_meshsort[k] = malloc (3 * sizeof (int));
            for (p = 0; p < 3; p++) {
                fgal_meshsort[k][p] = 0;    // [k][0] is mesh num in along x, [k][1] is mesh num in along y, [k][2] is index of the galaxy.
            }
        }
        
        for (k = 0; k < ymesh_num; k++) {
            fgal_count[k] = malloc (xmesh_num * sizeof (int));
            fgal_cum[k] = malloc (xmesh_num * sizeof (int));
            fgal_coord[k] = malloc (xmesh_num * sizeof (double *));
            for (p = 0; p < xmesh_num; p++) {
                fgal_count[k][p] = 0;
                fgal_cum[k][p] = 0;
                fgal_coord[k][p] = malloc (2 * sizeof (double));
                for (q = 0; q < 2; q++) {
                    fgal_coord[k][p][q] = 0;
                }
            }
        }
        
        // Store up the number of counts in the mesh grid
        for (k = 0; k < fgal_index; k++) {
            x_ind = min(x_cen, xmesh_num, fgal[k].x);
            y_ind = min(y_cen, ymesh_num, fgal[k].y);
            
            fgal_meshsort[k][0] = x_ind;
            fgal_meshsort[k][1] = y_ind;
            fgal_meshsort[k][2] = k;
            fgal_count[y_ind][x_ind]++;
        }
        
        // Cumulative number of galaxies in each grid.
        for (k = 0; k < 1; k++) {
            fgal_cum[k][0] = fgal_count[0][0];
            for (p = 1; p < xmesh_num; p++) {
                fgal_cum[k][p] = fgal_cum[k][p - 1] + fgal_count[k][p];
            }
        }
        
        for (k = 1; k < ymesh_num; k++) {
            fgal_cum[k][0] = fgal_cum[k - 1][xmesh_num - 1] + fgal_count[k][0];
            for (p = 1; p < xmesh_num; p++) {
                fgal_cum[k][p] = fgal_cum[k][p - 1] + fgal_count[k][p];
            }
        }
        
        for (k = 0; k < ymesh_num; k++) {
            for (p = 0; p < xmesh_num; p++) {
                fgal_coord[k][p][0] = x_cen[p];
                fgal_coord[k][p][1] = y_cen[k];
            }
        }
        
        
        // Sort according to grid of x; if the same, then sort according to grid of y, and z; if still the same, then sort by galaxy index in the lens galaxy catalog.
        qsort(fgal_meshsort, fgal_index, sizeof fgal_meshsort[0], compare);
        
        // Looping over each seed lens over the grids.
        for (p = 0; p < lgal_index; p++) {
            
            if (p % 10000 == 0) {
                printf("Seed Galaxy, p = %d\n", p);
            }

            for (q = 0; q < ymesh_num; q++) { //decmesh_num
                for (s = 0; s < xmesh_num; s++) {  //ramesh_num

                    // Distance between the lens and the grid.
                    //gl_a_sep = pow((pow((lgal[p].x - fgal_coord[q][s][0]), 2) + pow((lgal[p].y - fgal_coord[q][s][1]), 2)), 0.5);
                    gl_a_sep = dist_periodic(lgal[p].x, lgal[p].y, fgal_coord[q][s][0], fgal_coord[q][s][1], cube_size);

                    // See which grid to look at.
                    if (gl_a_sep <= 1.5 * xmesh_size && fgal_count[q][s] > 0) {
                        // Create a dynamic array to store the index of the corresponding lens galaxies.
                        int *fgrid_sort = malloc (fgal_count[q][s] * sizeof (int));
                        int *f_ind = malloc (fgal_count[q][s] * sizeof (int));
                                
                        for (t = 0; t < fgal_count[q][s]; t++) {    // fgal_count[q][s]
                                    
                            fgrid_sort[t] = fgal_cum[q][s] - fgal_count[q][s] + t;  // Number of concerned lens galaxy
                            f_ind[t] = fgal_meshsort[fgrid_sort[t]][2];     // The index of the concerned lens galaxy in the lens catalog.
                
                            //ll_a_sep = pow((pow((lgal[p].x - fgal[f_ind[t]].x), 2) + pow((lgal[p].y - fgal[f_ind[t]].y), 2)), 0.5);
                            ll_a_sep = dist_periodic(lgal[p].x, lgal[p].y, fgal[f_ind[t]].x, fgal[f_ind[t]].y, cube_size);
                                        
                            // Those within transverse physical distance D_0.
                            if (ll_a_sep <= D_0 && ll_a_sep >= 1e-7) {
                                
                                //printf("%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", p, f_ind[t], lgal[p].x, fgal[f_ind[t]].x, lgal[p].y, fgal[f_ind[t]].y, ll_a_sep, lgal[p].Mr, fgal[f_ind[t]].Mr);
                                
                                // If the "stacked" sample galaxy is fainter than the full catalog galaxy.
                                //if (lgal[p].Mr - fgal[f_ind[t]].Mr >= mag_bright) {
                                if (fgal[f_ind[t]].sm - lgal[p].sm >= sm) {
                                    //printf("%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", p, f_ind[t], lgal[p].x, fgal[f_ind[t]].x, lgal[p].y, fgal[f_ind[t]].y, ll_a_sep, lgal[p].sm, fgal[f_ind[t]].sm);
                                    
                                    lgal[p].sat = 1;
                                    break;
                                }
                                // To make sure the galaxies concerned has the same cluster number as the principal lens.
                                //printf("%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\n", p, l_ind[t], lgal[p].x, lgal[l_ind[t]].x, lgal[p].y, lgal[l_ind[t]].y, ll_a_sep, lgal[p].clus_num, lgal[l_ind[t]].clus_num, partner_count);
                                            
                            } // Those lenses within physical distance D_0
                                        
                        } // t loop, for the lenses inside the grid.
                                
                    } // which grid to look at.
                            
                } // s loop, same dec.
            } // q loop, same ra.
        } // p-loop (lens catalog)
        
        
        char goodname_file[150], color_string[13];
        FILE *goodname;
        
        strncpy (color_string, gcolor[i].id, 13);
        color_string[13] = 0;
        
        
        strcpy (goodname_file, mock_lens_link);   // lens_without_link
        strcat (goodname_file, color_string);
        strcat (goodname_file, "_");
        strcat (goodname_file, D_b_man);
        strcat (goodname_file, "_");
        strcat (goodname_file, star_mass);
        strcat (goodname_file, ".dat");
        
        goodname = fopen (goodname_file, "w");
        
        // Satellite count
        sat_count = 0;
        iso_count = 0;
        
        for (k = 0; k < lgal_index; k++) {
            if (lgal[k].sat != 0) {
                //fprintf(othername, "%s\t%.*lf\t%.*lf\t%.*lf\t%.*lf\t%.*lf\t%.*lE\t%s\t%d\n", lgal[k].id, 3, lgal[k].x, 3, lgal[k].y, 3, lgal[k].z, 3, lgal[k].Mr, 3, lgal[k].g_r, 3, lgal[k].m_vir, lgal[k].host_id, lgal[k].jack);
                sat_count++;
            }
            if (lgal[k].sat == 0) {
                fprintf(goodname, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lE\t%lf\t%lf\t%lf\t%lE\t%s\n", lgal[k].id,  lgal[k].x,  lgal[k].y,  lgal[k].z,  lgal[k].vx,  lgal[k].vy,  lgal[k].vz,  lgal[k].m_vir,  lgal[k].vpeak,  lgal[k].sm,  lgal[k].g_r,  lgal[k].m_host, lgal[k].host_id);
                iso_count++;
            }
        }
        
        printf("Total # of galaxies = %d\n", gcolor[i].num);
        printf("# of central galaxies = %d (%.*f %%)\n", iso_count, 3, iso_count * 100. / gcolor[i].num);
        printf("# of satellite galaxies = %d (%.3f %%)\n", sat_count, 3, sat_count * 100. / gcolor[i].num);
        
        fclose (goodname);
        
        
        for (k = 0; k < ymesh_num; k++) {
            free(fgal_count[k]);
            free(fgal_cum[k]);
            for (p = 0; p < xmesh_num; p++) {
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
    }
    
    free (x_cen);
    free (y_cen);
    free (full);
    free (gcolor);
    
    return 0;
}