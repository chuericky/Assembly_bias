/* This program generates the mock galaxy catalog with known correlation functions, w(theta) = A/theta, A = 1000.
    Field # (argv[1])
    Updated: Oct 19, 2016
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include <time.h>
#include "tool.h"
#include "cosmopar.h"

// Struct for galaxy catalog.
typedef struct {
    char name[20];
    int num;
} namelens;

typedef struct {
    double ra, dec, ra_rad, dec_rad;
} lensgal;


// Probability distribution function of the correlation function.
double corProb(double theta, double A) {
    return 0.5 * pow(theta,2) + A * pow(theta,1.2)/1.2;
}


int main(int argc, char *argv[]) {
    char lensname_file[150], field_manual[3];
    FILE *lensname;
    
    int i = 0, j = 0, k = 0, lensfile_num = 0, lensindex = 0, lgal_index = 0, spot = 0, file_count = -1, field_int = 0, nbr_count = 0, ln_sepint = 0;
    double prob_sep = 0, theta_sep = 0, phi = 0, ra_nbr = 0, dec_nbr = 0;
    const int PSIZE = 1000, RAN_NBR = 10;
    const double D2R = M_PI / 180, A = 100., TWOPI = 2 * M_PI;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 2) {
        printf("Invalid passing arguments!\n");
        exit(1);
    }
    
    srand (time(NULL));
    
    double *a_sep = malloc (PSIZE * sizeof(double));
    double *cumprob = malloc (PSIZE * sizeof(double));
    
    for (i = 0; i < PSIZE; i++) {
        a_sep[i] = 0.001 * (i+1) * D2R;
        // Assigning probability to assign locations of galaxies.
        cumprob[i] = corProb(a_sep[i],A) / corProb(D2R,A);
        //printf("%d\t%lf\t%lf\n", i, a_sep[i], cumprob[i]);
    }
    
    // Write in the lens files.
    strcpy (lensname_file, gal_link);
    strcat (lensname_file, "lensname.txt");
    
    lensname = fopen (lensname_file, "r");
    
    // Make sure the lens name file is present.
    if (!lensname) {
        printf("Cannot find the lens name file!\n");
        exit(1);
    }
    
    // Count number of lines in the lens name file.
    lensfile_num = count_line (spot, lensname, file_count);
    
    // An array storing all file names.
    namelens *lens = malloc (lensfile_num * sizeof *lens);
    
    lensindex = 0;
    
    while (fscanf(lensname, "%d %s", &lens[lensindex].num, &lens[lensindex].name) == 2) {
        //printf("%d\t%s\n", lens[lensindex].num, lens[lensindex].name);
        lensindex++;
    }
    
    // For Field W1 to W4
    strcpy (field_manual, argv[1]);
    field_int = atoi(field_manual);
    
    for (i = field_int - 1; i < field_int; i++) {
        
        printf("\nField W%d\n", i + 1);
        
        FILE *lensfile;
        char lensfile_string[150];
        
        strcpy (lensfile_string, gal_link);
        strcat (lensfile_string, lens[i].name);
        
        lensfile = fopen (lensfile_string, "r");
        
        // Make sure the source name file is present.
        if (!lensfile) {
            printf("Cannot find the lens file %s!\n", lens[i].name);
            exit(1);
        }
        
        // Open up the information of the source galaxies.
        lensgal *lgal = malloc (lens[i].num * sizeof *lgal);
        
        lgal_index = 0;
        // Total 135 columns of the source catalog.
        while (fscanf(lensfile, "%*s %*s %*lf %*lf %lf %lf %*lE %*lE %*lf %*lf %*d %*lE %*lE %*lf %*lE %*lE %*lf %*lf %*lf %*lf %*lE %*d %*lf %*d %*lf %*lf %*lf %*lf %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*d %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lf %*lE %*lf %*lE %*lf %*lf %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE %*lE", &lgal[lgal_index].ra, &lgal[lgal_index].dec) == 2) {
            
            lgal[lgal_index].ra_rad = lgal[lgal_index].ra * D2R;
            lgal[lgal_index].dec_rad = lgal[lgal_index].dec * D2R;
            
            lgal_index++;
            // if (lgal_index == 10) break;
        }
        
        fclose (lensfile);
        
        // Number of neighbors, RAN_NBR = number of neighbors for each lens.
        nbr_count = lens[i].num * RAN_NBR;
        
        // A dynamic array to store the positions of the neighbor galaxies.
        lensgal **ngal = malloc (lgal_index * sizeof **ngal);
        
        FILE *probfile;
        char probfile_string[150], f_ind[3];
        
        strcpy (probfile_string, ran_link);
        strcat (probfile_string, "W");
        strcpy (f_ind, field_manual);
        strcat (probfile_string, f_ind);
        strcat (probfile_string, "/combined/W");
        strcat (probfile_string, f_ind);
        strcat (probfile_string, "_probtesting.txt");
        
        probfile = fopen (probfile_string, "w");
        
        // For each lens
        for (j = 0; j < lgal_index; j++) {
            if (j % 5000 == 0) {
                printf("%d\n", j);
            }
            ngal[j] = malloc (RAN_NBR * sizeof *ngal);
            
            for (k = 0; k < RAN_NBR; k++) {
                prob_sep = fRand(0,1);
                ln_sepint = min(cumprob, PSIZE, prob_sep);
                theta_sep = a_sep[ln_sepint];
                
                // Assign the azimuthal angle for the neighbor.
                phi = fRand(0,TWOPI);

                ra_nbr = lgal[j].ra_rad + theta_sep * cos(phi);
                dec_nbr = lgal[j].dec_rad + theta_sep * sin(phi);
                
                //printf("%d\t%lf\t%d\t%lf\t%lf\n", k, prob_sep, ln_sepint, theta_sep, phi);
                fprintf(probfile, "%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\n", lgal[j].ra, lgal[j].dec, lgal[j].ra_rad, lgal[j].dec_rad, theta_sep, phi, ra_nbr, dec_nbr);
            }
        }
        
        fclose (probfile);
        
        
        for (j = 0; j < lgal_index; j++) {
            free (ngal[j]);
        }
        free (ngal);
        free (lgal);
    }

    free (lens);
    free (a_sep);
    free (cumprob);
    fclose (lensname);
    
    return 0;
}