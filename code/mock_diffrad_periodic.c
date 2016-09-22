/* This program constructs the mock galaxy catalog by varying A, where d <= A * r_200c.
    Input: A (argv[1]), halo file (argv[2])
    Updated: Dec 14, 2015
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "tool.h"
#include "cosmopar.h"
#include "halo_model.h"

typedef struct {
    char halo_name[15], host_tag[15];
    double x, y, z, vx, vy, vz, m_vir, vpeak, M_r, g_r, m_host;
} namelens;

typedef struct {
    double x, y, z, r_200c;
} nameclus;

void writematrix(FILE *infile, namelens *w);


int main(int argc, char *argv[]) {
    char lens_file[150], clus_file[150], fac_char[5], lensout_string[150], mockhalo_str[7];
    FILE *lensname, *clusname, *lensout;
    int i = 0, j = 0, lens_num = 0, clus_num = 0, file_count = -1, spot = 0, lensindex = 0, clusindex = 0, counter = 0;
    double rho_0 = 0, lc_dist = 0, fac = 0, cube_size = 250.0;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 3) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    // A, the factor.
    strcpy (fac_char, argv[1]);
    fac = atof(fac_char);
    
    strcpy (mockhalo_str, argv[2]);
    
    // Critical density at z = 0
    rho_0 = rho_crit(0.0);
    
    // Read in the lens galaxies.
    strcpy (lens_file, mock_lens_link);
    strcat (lens_file, "sm_gr_fiducial_mock.dat");
    
    lensname = fopen (lens_file, "r");
    
    // Make sure the lens name file is present.
    if (!lensname) {
        printf("Cannot find the lens name file!\n");
        exit (1);
    }
    
    // Count number of lines in the lens name file.
    lens_num = count_line (spot, lensname, file_count);
    
    // An array to store all the info for the lens.
    namelens *lens = malloc (lens_num * sizeof *lens);
    
    lensindex = 0;
    
    while (fscanf(lensname, "%s %lf %lf %lf %lf %lf %lf %lE %lf %lf %lf %lE %s", &lens[lensindex].halo_name, &lens[lensindex].x, &lens[lensindex].y, &lens[lensindex].z, &lens[lensindex].vx, &lens[lensindex].vy, &lens[lensindex].vz, &lens[lensindex].m_vir, &lens[lensindex].vpeak, &lens[lensindex].M_r, &lens[lensindex].g_r, &lens[lensindex].m_host, &lens[lensindex].host_tag) == 13) {
        lensindex++;
    }
    
    // Read in the clusters.
    strcpy (clus_file, halocat_link);
    strcat (clus_file, "halo_");
    strcat (clus_file, mockhalo_str);
    strcat (clus_file, ".dat");
    
    clusname = fopen (clus_file, "r");
    
    // Make sure the cluster name file is present.
    if (!clusname) {
        printf("Cannot find the cluster name file!\n");
        exit (1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count number of lines in the cluster name file.
    clus_num = count_line (spot, clusname, file_count);
    
    // An array to store all the info for the cluster.
    nameclus *clus = malloc (clus_num * sizeof *clus);
    
    clusindex = 0;
    
    while (fscanf(clusname, "%*s %*d %*d %*lf %*lf %*lf %lf %lf %lf %*lf %lf", &clus[clusindex].x, &clus[clusindex].y, &clus[clusindex].z, &clus[clusindex].r_200c) == 4) {
        // Virial radius of the cluster, in Mpc
        // clus[clusindex].r_host = rad_200(rho_0, clus[clusindex].m_host) / Mpc_2_pc;
        //printf("%d\t%lE\t%lE\t%lE\t%lE\n", clusindex, clus[clusindex].x, clus[clusindex].y, clus[clusindex].z, clus[clusindex].r_200c);
        clusindex++;
        //if (clusindex == 100) break;
    }
    
    strcpy (lensout_string, mock_lens_link);
    strcat (lensout_string, "sm_");
    strcat (lensout_string, mockhalo_str);
    strcat (lensout_string, "_");
    strcat (lensout_string, fac_char);
    strcat (lensout_string, "_0.0.dat");
    
    lensout = fopen (lensout_string, "w");
    
    // Make sure the source name file is present.
    if (!lensout) {
        printf("Cannot open the lens file %s!\n", lensout_string);
        exit(1);
    }
    
    // Loop over the lens catalog.
    for (i = 0; i < lensindex; i++) {   // lensindex
        if (i % 5000 == 4999) printf("%d / %d lens completed!\n", i, lensindex);
        
        counter = 0;
        
        // Loop over the cluster catalog.
        for (j = 0; j < clus_num; j++) {
            // Projected distance between the lens and the cluster.
            //lc_dist = pow((pow((lens[i].x - clus[j].x), 2) + pow((lens[i].y - clus[j].y), 2)), 0.5);
            lc_dist = dist_periodic(lens[i].x, lens[i].y, clus[j].x, clus[j].y, cube_size);
            
            // Not belonging to this cluster.
            if (lc_dist > fac * clus[j].r_200c) {
                counter++;
                continue;   // Proceed to the next cluster.
            }
            
            // Belonging to this cluster.
            if (lc_dist <= fac * clus[j].r_200c) {
                //printf("%d\t%lf\t%lf\t%d\t%lf\t%lf\t%lE\t%lE\t%lE\n", i, lens[i].x, lens[i].y, j, clus[j].x, clus[j].y, lc_dist, clus[j].r_host, clus[j].r_host * fac);
                break;
            }
            
            //printf("%d\t%lf\t%lf\t%d\t%lf\t%lf\t%lE\t%lE\n", i, lens[i].x, lens[i].y, j, clus[j].x, clus[j].y, lc_dist, clus[j].r_host);
        }
        
        if (counter == clus_num)   writematrix (lensout, &lens[i]);
    }
    
    free (lens);
    free (clus);
    
    fclose (lensname);
    fclose (clusname);
    fclose (lensout);
    
    return 0;
}

// The function to write the entire matrix to a txt file.
void writematrix(FILE *infile, namelens *w) {
    fprintf(infile, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lE\t%lf\t%lf\t%lf\t%lE\t%s\n", w->halo_name, w->x, w->y, w->z, w->vx, w->vy, w->vz, w->m_vir, w->vpeak, w->M_r, w->g_r, w->m_host, w->host_tag);
}
