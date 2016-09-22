/* This program selects the neighboring halos with mass >= 10^11 M_sun
 Updated: Jun 3, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tool.h"
#include "cosmopar.h"
#include "halo_model.h"

// Struct for the galaxy catalog.
typedef struct {
    unsigned long id;
    int pid, upid;
    double x, y, z, m_200c;
} halo;

int main() {
    int individual_range = 0, skip_line = 0;
    double d_vir = 0, cube_size = 0;
    unsigned long i = 0, halo_num = 0, halo_index = 0;
    char halo_file[150], halo_extract[150];
    FILE *halo_f, *halo_e;
    
    
    // Read the halo catalog.
    strcpy (halo_file, halocat_link);
    strcat (halo_file, "hlist_1.00030.list");
    
    halo_f = fopen (halo_file, "r");
    
    // Make sure the halo catalog is present.
    if (!halo_f) {
        printf("Cannot find the halo catalog!\n");
        exit(1);
    }
    
    halo_num = 11916318;
    cube_size = 250.0;        // Size of a side of the sim cube, 250 Mpc/h
    
    do {
        individual_range = fgetc(halo_f);
        if(individual_range == '\n')
        skip_line++;
        if (skip_line == 42)
            break;
    } while (individual_range != EOF);
    
    // A struct to store the information for the clusters.
    halo *clus = malloc (halo_num * sizeof *clus);
    
    halo_index = 0;
    
    while (fscanf(halo_f, "%*f %u %*f %*d %*d %d %d %*d %*d %*f %*lE %*lE %*lE %*f %*d %*f %*f %lf %lf %lf %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %lE %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f", &clus[halo_index].id, &clus[halo_index].pid, &clus[halo_index].upid, &clus[halo_index].x, &clus[halo_index].y, &clus[halo_index].z, &clus[halo_index].m_200c) == 7) {
        if (halo_index % 500000 == 0) printf("%u out of %u read!\n", halo_index, halo_num);
        //printf("%d\t%u\t%d\t%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", halo_index, clus[halo_index].id, clus[halo_index].pid, clus[halo_index].upid, clus[halo_index].m_vir, clus[halo_index].r_vir, clus[halo_index].r_s, clus[halo_index].x, clus[halo_index].y, clus[halo_index].z, clus[halo_index].m_200c);
        halo_index++;
        //if (halo_index == 10) break;
    }

    d_vir = Delta_vir(oL);

    fclose (halo_f);
    
    // Read the halo catalog.
    strcpy (halo_extract, halonbr_link);
    strcat (halo_extract, "halo_1.0e10.dat");
    
    halo_e = fopen (halo_extract, "w");
    
    // Make sure the halo catalog can be opened.
    if (!halo_e) {
        printf("Cannot open the halo catalog!\n");
        exit(1);
    }
    
    for (i = 0; i < halo_index; i++) {
        if (clus[i].m_200c >= 1e10) {
            
            fprintf(halo_e, "%u\t%lE\t%lE\t%lE\t%lE\n", clus[i].id, clus[i].x, clus[i].y, clus[i].z, clus[i].m_200c);
        }
    }
    
    fclose (halo_e);
    
    free (clus);
    
    return 0;
}