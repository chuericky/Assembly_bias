/* This program selects the host halos
 Updated: Dec 14, 2015
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
    double m_vir, r_vir, r_s, x, y, z, m_200c, r_200c;
} halo;

int main() {
    int individual_range = 0, skip_line = 0;
    double d_vir = 0, cube_size = 0;
    unsigned long i = 0, halo_num = 0, halo_index = 0;
    char halo_file[150], halo_extract[150], halo_ext_1[150], halo_ext_2[150], halo_ext_3[150], halo_ext_4[150], halo_ext_5[150];
    FILE *halo_f, *halo_e, *halo_1, *halo_2, *halo_3, *halo_4, *halo_5;
    
    
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
    
    while (fscanf(halo_f, "%*f %u %*f %*d %*d %d %d %*d %*d %*f %lE %lE %lE %*f %*d %*f %*f %lf %lf %lf %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %lE %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f", &clus[halo_index].id, &clus[halo_index].pid, &clus[halo_index].upid, &clus[halo_index].m_vir, &clus[halo_index].r_vir, &clus[halo_index].r_s, &clus[halo_index].x, &clus[halo_index].y, &clus[halo_index].z, &clus[halo_index].m_200c) == 10) {
        if (halo_index % 500000 == 0) printf("%u out of %u read!\n", halo_index, halo_num);
        //printf("%d\t%u\t%d\t%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", halo_index, clus[halo_index].id, clus[halo_index].pid, clus[halo_index].upid, clus[halo_index].m_vir, clus[halo_index].r_vir, clus[halo_index].r_s, clus[halo_index].x, clus[halo_index].y, clus[halo_index].z, clus[halo_index].m_200c);
        halo_index++;
        //if (halo_index == 10) break;
    }

    d_vir = Delta_vir(oL);

    fclose (halo_f);
    
    // Read the halo catalog.
    strcpy (halo_extract, halocat_link);
    strcat (halo_extract, "halo_1.0e11.dat");
    strcpy (halo_ext_1, halocat_link);
    strcat (halo_ext_1, "halo_1.5e13.dat");
    strcpy (halo_ext_2, halocat_link);
    strcat (halo_ext_2, "halo_5.0e13.dat");
    strcpy (halo_ext_3, halocat_link);
    strcat (halo_ext_3, "halo_1.0e14.dat");
    strcpy (halo_ext_4, halocat_link);
    strcat (halo_ext_4, "halo_2.5e14.dat");
    strcpy (halo_ext_5, halocat_link);
    strcat (halo_ext_5, "halo_3.0e14.dat");
    
    halo_e = fopen (halo_extract, "w");
    halo_1 = fopen (halo_ext_1, "w");
    halo_2 = fopen (halo_ext_2, "w");
    halo_3 = fopen (halo_ext_3, "w");
    halo_4 = fopen (halo_ext_4, "w");
    halo_5 = fopen (halo_ext_5, "w");
    
    // Make sure the halo catalog can be opened.
    if (!halo_e) {
        printf("Cannot open the halo catalog!\n");
        exit(1);
    }
    
    for (i = 0; i < halo_index; i++) {
        if (clus[i].m_200c >= 1e11 & clus[i].pid == -1 & clus[i].upid == -1) {
            clus[i].r_200c = pow((clus[i].m_200c/clus[i].m_vir) * (d_vir/200), 1./3) * clus[i].r_vir / 1000;
            
            fprintf(halo_e, "%u\t%d\t%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", clus[i].id, clus[i].pid, clus[i].upid, clus[i].m_vir, clus[i].r_vir / 1000, clus[i].r_s / 1000, clus[i].x, clus[i].y, clus[i].z, clus[i].m_200c, clus[i].r_200c);
            if (clus[i].m_200c >= 1.5e13 && clus[i].m_200c < 5.0e13) {
                fprintf(halo_1, "%u\t%d\t%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", clus[i].id, clus[i].pid, clus[i].upid, clus[i].m_vir, clus[i].r_vir / 1000, clus[i].r_s / 1000, clus[i].x, clus[i].y, clus[i].z, clus[i].m_200c, clus[i].r_200c);
            }
            if (clus[i].m_200c >= 5.0e13 && clus[i].m_200c < 1.0e14) {
                fprintf(halo_2, "%u\t%d\t%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", clus[i].id, clus[i].pid, clus[i].upid, clus[i].m_vir, clus[i].r_vir / 1000, clus[i].r_s / 1000, clus[i].x, clus[i].y, clus[i].z, clus[i].m_200c, clus[i].r_200c);
            }
            if (clus[i].m_200c >= 1.0e14 && clus[i].m_200c < 2.5e14) {
                fprintf(halo_3, "%u\t%d\t%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", clus[i].id, clus[i].pid, clus[i].upid, clus[i].m_vir, clus[i].r_vir / 1000, clus[i].r_s / 1000, clus[i].x, clus[i].y, clus[i].z, clus[i].m_200c, clus[i].r_200c);
            }
            if (clus[i].m_200c >= 2.5e14 && clus[i].m_200c < 3.0e14) {
                fprintf(halo_4, "%u\t%d\t%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", clus[i].id, clus[i].pid, clus[i].upid, clus[i].m_vir, clus[i].r_vir / 1000, clus[i].r_s / 1000, clus[i].x, clus[i].y, clus[i].z, clus[i].m_200c, clus[i].r_200c);
            }
            if (clus[i].m_200c >= 3.0e14) {
                fprintf(halo_5, "%u\t%d\t%d\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\t%lE\n", clus[i].id, clus[i].pid, clus[i].upid, clus[i].m_vir, clus[i].r_vir / 1000, clus[i].r_s / 1000, clus[i].x, clus[i].y, clus[i].z, clus[i].m_200c, clus[i].r_200c);
            }
        }
    }
    
    fclose (halo_e);
    fclose (halo_1);
    fclose (halo_2);
    fclose (halo_3);
    fclose (halo_4);
    fclose (halo_5);
    
    free (clus);
    
    return 0;
}