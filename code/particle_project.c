/* This program searches the cosmosim particles for the halos, where 2048^3 particles are involved.
    Particle lists are divided into 859 smaller lists. This program reads one sub-list, projects the particles, and iterate it by 859 times.
    Particle list in the halos are "appended", need to delete the halo files before using it.
    Search of particles is done by chaining mesh technique.
    
    Updated: Oct 1, 2015
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "tool.h"
#include "cosmopar.h"

typedef struct {
    double x, y, z;
} particle;

typedef struct {
    char id[15], color[5];
    double x, y, z, r_vir;
} gal_cat;

typedef struct {
    int num;
    char filename[10];
} halo_cat;

int main(int argc, char *argv[]) {
    
    int i = 0, j = 0, halo_cat_num = 0, spot = 0, file_count = -1, halo_cat_index = 0, halo_index = 0, tot_halo_num = 0;
    char cosmosim_link[150], halo_cat_link[150], color[5];
    FILE *halo_cat_file, *particle_file;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 1) {
        printf("Invalid # of passing arguments!\n");
        return 1;
    }
    
    // Halo catalog.
    strcpy (halo_cat_link, mock_lens_link);
    strcat (halo_cat_link, "lensname.txt");
    
    halo_cat_file = fopen (halo_cat_link, "r");
    
    if (!halo_cat_file) {
        printf("Cannot find the halo catalog file %s!\n", halo_cat_link);
        exit(1);
    }
    
    // Number of halo catalog files.
    halo_cat_num = count_line (spot, halo_cat_file, file_count);
    
    halo_cat *halo = malloc (halo_cat_num * sizeof *halo);
    
    halo_cat_index = 0;
    tot_halo_num = 0;
    
    // Get the names of the files, and the number of halos in that file.
    while (fscanf (halo_cat_file, "%d %s", &halo[halo_cat_index].num, &halo[halo_cat_index].filename) == 2) {
        tot_halo_num += halo[halo_cat_index].num;
        halo_cat_index++;
    }
    
    // Halo catalog, store up halo id, coordinates and R_vir.
    gal_cat *gal = malloc (tot_halo_num * sizeof *gal);

    halo_index = 0;
    
    for (i = 0; i < halo_cat_index; i++) {
        char halo_str[150];
        FILE *halo_indiv;
        
        strcpy (halo_str, mock_lens_link);
        strcat (halo_str, halo[i].filename);
        
        halo_indiv = fopen (halo_str, "r");
        
        if (!halo_indiv) {
            printf("No halo catalog file %s found!\n", halo[i].filename);
            exit(1);
        }
        
        while (fscanf (halo_indiv, "%s %lf %lf %lf %lE %*lf %*lf %*lf %*d %*d", &gal[halo_index].id, &gal[halo_index].x, &gal[halo_index].y, &gal[halo_index].z, &gal[halo_index].r_vir) == 5) {
            strncpy (gal[halo_index].color, halo[i].filename, 3);
            gal[halo_index].color[3] = 0;
            halo_index++;
        }
        
        fclose (halo_indiv);
    }
    
    for (i = 0; i < halo_index; i++) {
        printf("%d\t%s\t%s\t%lf\t%lf\t%lf\t%lE\n", i, gal[i].color, gal[i].id, gal[i].x, gal[i].y, gal[i].z, gal[i].r_vir);
    }
    
    /*
    strcpy (cosmosim_link, hard_drive_link);
    strcat (cosmosim_link, "000.csv");
    
    particle_file = fopen (cosmosim_link, "r");
    
    if (!particle_file) {
        printf("Cannot find the cosmosim file %s!\n", cosmosim_link);
        exit(1);
    }
    
    fclose (particle_file);
    */
    
    
    
    free (gal);
    free (halo);
    
    fclose (halo_cat_file);
    
    return 0;
}