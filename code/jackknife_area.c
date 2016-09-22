/* This program assigns the jackknife index for the data galaxies.
    Input: with / without (argv[1]), b (argv[2]), clus_num (argv[3]), color (argv[4]), P_th (argv[5])
    Updated: Dec 6, 2015
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "tool.h"

// Struct for the galaxy catalog.
typedef struct {
    char name[20], field[7];
    double ra, dec, z, T_b, M_r, mass;
    int jackindex;
} namelens;

typedef struct {
    char field[7];
} nlens;

int main(int argc, char *argv[]) {
    int i = 0, j = 0, spot = 0, file_count = 0, field_num = 0, fieldindex = 0;
    FILE *datafile;
    char dataname[150], parameter_name[20], paracat_name[20], b[5], clus_man[4], color[5], fstr[30], P_th[10];
    
    // Make sure the exact number of arguments are passed.
    if (argc != 6) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    strcpy (paracat_name, argv[1]);
    strcpy (b, argv[2]);
    strcpy (clus_man, argv[3]);
    strcpy (color, argv[4]);
    strcpy (P_th, argv[5]);
    
    strcpy (fstr, b);
    strcat (fstr, "_");
    strcat (fstr, P_th);
    strcat (fstr, "_");
    strcat (fstr, clus_man);
    strcat (fstr, "_lens");
    strcat (fstr, color);
    
    strcpy (dataname, twopt_link);
    strcat (dataname, paracat_name);
    strcat (dataname, "/");
    strcat (dataname, fstr);
    strcat (dataname, "_data_count.dat");
    
    datafile = fopen (dataname, "w");
    
    long *data_count = malloc (54 * sizeof (long));
    
    for (i = 0; i < 54; i++) {
        data_count[i] = 0;
    }
    
    for (i = 0; i < 4; i++) {   // 199
        
        FILE *lensfile, *writelensfile;
        char lensname[150], writelensname[150], lens_index[4], lens_field[7];
        
        printf("%d\n", i);
        
        strcpy (lensname, shear_gal_link);
        strcat (lensname, paracat_name);
        strcat (lensname, "/W");
        sprintf (lens_index, "%d", i + 1);
        strcat (lensname, lens_index);
        strcat (lensname, "_");
        strcat (lensname, fstr);
        strcat (lensname, "_S1.tsv");
        strcpy (writelensname, shear_gal_link);
        strcat (writelensname, paracat_name);
        strcat (writelensname, "/W");
        strcat (writelensname, lens_index);
        strcat (writelensname, "_");
        strcat (writelensname, fstr);
        strcat (writelensname, "_S1_2.tsv");
        
        lensfile = fopen (lensname, "r");
        writelensfile = fopen (writelensname, "w");
        
        if (!lensfile) {
            printf("Cannot find the lens file %s!\n", lensfile);
            exit(1);
        }
        
        field_num = 0;
        spot = 0;
        file_count = -1;
        
        field_num = count_line (spot, lensfile, file_count);
        
        namelens *lens = malloc (field_num * sizeof *lens);
        
        fieldindex = 0;
        
        while (fscanf(lensfile, "%s %s %lf %lf %lf %lf %lf %lf", &lens[fieldindex].name, &lens[fieldindex].field, &lens[fieldindex].ra, &lens[fieldindex].dec, &lens[fieldindex].z, &lens[fieldindex].T_b, &lens[fieldindex].M_r, &lens[fieldindex].mass) == 8) {
            strcpy (lens_field, lens[fieldindex].field);
            lens[fieldindex].jackindex = jackarea(lens[fieldindex].field);
            data_count[lens[fieldindex].jackindex]++;
            fieldindex++;
        }
        
        for (j = 0; j < fieldindex; j++) {
            fprintf(writelensfile, "%s\t%s\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n", lens[j].name, lens[j].field, lens[j].ra, lens[j].dec, lens[j].z, lens[j].T_b, lens[j].M_r, lens[j].mass, lens[j].jackindex);
        }

        free (lens);
        
        fclose (lensfile);
        fclose (writelensfile);
        
    }
    
    for (i = 0; i < 54; i++) {
        fprintf(datafile, "%d\t%ld\n", i, data_count[i]);
    }
    
    long int total_count = 0;
    
    for (i = 0; i < 54; i++) {
        total_count += data_count[i];
    }
    
    printf("%ld\n", total_count);
    
    free (data_count);
    
    fclose (datafile);
    
    return 0;
}