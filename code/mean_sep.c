/* This program calculates the mean angular separation for CFHTLenS galaxies.
    Input: field # (argv[1])
    Updated: Oct 13, 2015
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
    char id[15];
    double alpha, delta;
} gallens;


int main(int argc, char *argv[]) {
    
    int i = 0, j = 0, k = 0, spot = 0, file_count = 0, lensfile_num = 0, lensindex = 0, f_int = 0;
    unsigned long long int lgal_index = 0, gal_pair = 0, pair_count = 0;
    double tot_sep = 0, mean_sep = 0;
    char field_file[150], field_man[3];
    FILE *fieldname;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 2) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    strcpy (field_file, lens_multi_link);
    strcat (field_file, "lensname.txt");
    
    fieldname = fopen (field_file, "r");
    
    // Make sure the lensname file is present.
    if (!fieldname) {
        printf("Cannot find the field number count file!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count number of lines in the lens name file.
    lensfile_num = count_line (spot, fieldname, file_count);
    
    // An array storing all the file names.
    namelens *lens = malloc (lensfile_num * sizeof * lens);
    
    lensindex = 0;
    
    while (fscanf(fieldname, "%d %s", &lens[lensindex].num, &lens[lensindex].name) == 2) {
        lensindex++;
    }
    
    fclose (fieldname);
    
    // Field W1 to W4.
    strcpy (field_man, argv[1]);
    f_int = atoi(field_man);
    
    for (i = f_int - 1; i < f_int; i++) {
        printf("\nField W%d\n", i + 1);
        
        FILE *lensfile;
        char lens_string[150];
        
        strcpy (lens_string, lens_multi_link);
        strcat (lens_string, lens[i].name);
        
        lensfile = fopen (lens_string, "r");
        
        // Make sure the lens name file is present.
        if (!lensfile) {
            printf("Cannot find the lens file %s!\n", lens[i].name);
            exit(1);
        }
        
        // Open up the information of the lens galaxies.
        gallens *lgal = malloc (lens[i].num * sizeof *lgal);
        
        lgal_index = 0;
        
        // Total 135 columns of the source catalog.
        while (fscanf(lensfile, "%s %*s %*lf %*lf %lf %lf %*lE %*lE %*lf %*lf %*d %*lE %*lE %*lf %*lE %*lE %*lf %*lf %*lf %*lf %*lE %*d %*lf %*d %*lf %*lf %*lf %*lf %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*d %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lE %*lf %*lE %*lf %*lf %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE", &lgal[lgal_index].id, &lgal[lgal_index].alpha, &lgal[lgal_index].delta) == 3) {
            
            //printf("%d\t%s\t%lf\t%lf\n", lgal_index, lgal[lgal_index].id, lgal[lgal_index].alpha, lgal[lgal_index].delta);
            
            lgal_index++;
            //if (lgal_index == 5) break;
        }
        fclose (lensfile);
        
        // Number of galaxy pairs.
        gal_pair = lgal_index * (lgal_index - 1) / 2;
        
        pair_count = 0;
        tot_sep = 0;
        mean_sep = 0;
        
        // Start from the 1st galaxy.
        for (j = 0; j < lgal_index - 1; j++) {
            if (j % 500 == 0) {
                printf("%d out of %d calculation finished!\n", j, lgal_index - 1);
            }
            for (k = j + 1; k < lgal_index; k++) {
                
                tot_sep += ang_sep (lgal[j].alpha, lgal[j].delta, lgal[k].alpha, lgal[k].delta);
                //printf("%d\t%s\t%lf\t%lf\t%d\t%s\t%lf\t%lf\t%lf\n", j, lgal[j].id, lgal[j].alpha, lgal[j].delta, k, lgal[k].id, lgal[k].alpha, lgal[k].delta, ang_sep (lgal[j].alpha, lgal[j].delta, lgal[k].alpha, lgal[k].delta));
                
                pair_count++;
                
            }   // k-loop
        }       // j-loop
        
        mean_sep = tot_sep / pair_count;
        
        printf("%lli\t%lf\t%lf\n", pair_count, tot_sep, mean_sep);
        
        
        if (gal_pair != pair_count) {
            printf("Galaxy count doesn't match up!\n");
            exit (1);
        }
        
        free (lgal);
    }
    
    return 0;
}




