// This program counts the number of pixels with masked value <= 1, and calculate the effective surface area, and the corresponding number of randoms being projected.
// Updated: Jun 11, 2015

#include <string.h>
#include <stdio.h>
#include <time.h>
#include "cfitsio/fitsio.h"

void printerror( int status);

typedef struct {
    char name[200], flag[200];
    long num;
    double eff_area;
    int ran_num;
} fitsin;

int main(int argc, char *argv[])
{
    fitsfile *fptr;//, *fptr2;       /* pointer to the FITS file, defined in fitsio.h */
    int status, status2,  nfound, anynull, anynull2, fitsindex = 0, fitsindex2 = 0, number;
    long naxes[2], fpixel, fpixel2, nbuffer, npixels, nbuffer_flag, npixels_flag, ii, jj, count = 0, count2 = 0, i = 0, j = 0, k = 0, tot_cou = 0, ran_check = 0;
    double pix_area = 5.166666789E-05 * 5.166666789E-05;
    unsigned long un_count = 0;
    int fitsnum[4] = {72, 25, 49, 25};
    double rannum[4] = {3240000., 1125000., 2205000., 1125000.};
    
    /* initialize random seed: */
    srand ( time(NULL) );
    
    for (i = 1; i < 5; i++) {
        FILE *lensname, *unmask_file, *flagname;
        char lensname_file[150], unmaskname_file[150], flag_file[150], str[2];
        
        strcpy (lensname_file, "/Users/rickyccy/Documents/Research_assemblybias/Masks/W");
        sprintf(str, "%d", i);
        strcat (lensname_file, str);
        strcat (lensname_file, "/fitsname.txt");
        
        lensname = fopen(lensname_file, "r");
        
        // Make sure the source name file is present.
        if (!lensname) {
            printf("Cannot find the fits name file!\n");
            exit(1);
        }
        
        /*strcpy (flag_file, "/Users/rickyccy/Documents/Research_assemblybias/Masks/W");
        sprintf(str, "%d", i);
        strcat (flag_file, str);
        strcat (flag_file, "/flag.txt");
        
        flagname = fopen(flag_file, "r");
        
        // Make sure the source name file is present.
        if (!flagname) {
            printf("Cannot find the flag name file!\n");
            exit(1);
        }*/
        
        strcpy (unmaskname_file, "/Users/rickyccy/Documents/Research_assemblybias/Masks/W");
        sprintf(str, "%d", i);
        strcat (unmaskname_file, str);
        strcat (unmaskname_file, "/unmask_num_all.dat");
        
        unmask_file = fopen(unmaskname_file, "w");
        
        fitsin *fits = malloc (fitsnum[i - 1] * sizeof *fits);
        
        fitsindex = 0;
        
        while (fscanf(lensname, "%s", &fits[fitsindex].name) == 1) {
            fitsindex++;
        }
        
        /*fitsindex2 = 0;
        
        while (fscanf(flagname, "%s", &fits[fitsindex2].flag) == 1) {
            fitsindex2++;
        }
        
        if (fitsindex != fitsindex2) {
            printf("Flag file and mask file count don't match!\n");
            exit(1);
        }*/
        
        tot_cou = 0;
        
        for (j = 0; j < fitsindex; j++) {  //fitsindex
            printf("%d\t%s\n", j, fits[j].name);
        
            #define buffsize 21000
            float nullval, nullval2;
            char filename[200], flagnamefile[200];
            
            strcpy (filename, "/Users/rickyccy/Documents/Research_assemblybias/Masks/W");
            strcat (filename, str);
            strcat (filename, "/");
            strcat (filename, fits[j].name);
            
            status = 0;
    
            if ( fits_open_file(&fptr, filename, READONLY, &status) )
                printerror( status );
            
            /*strcpy (flagnamefile, "/Users/rickyccy/Documents/Research_assemblybias/Masks/W");
            strcat (flagnamefile, str);
            strcat (flagnamefile, "/");
            strcat (flagnamefile, fits[j].flag);
            
            status2 = 0;
            
            if ( fits_open_file(&fptr2, flagnamefile, READONLY, &status2) )
                printerror( status2 );
            
            if ( fits_movnam_hdu(fptr2, 0, "COMPRESSED_IMAGE", 0, &status2) )
            printerror( status2 );*/
            
    
            // read the NAXIS1 and NAXIS2 keyword to get image size
            if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
                printerror( status );
    
            npixels  = naxes[0] * naxes[1];         // number of pixels in the image
            fpixel   = 1;
            fpixel2  = 1;
            nullval  = 0;                // don't check for null values in the image
            nullval2 = 0;
    
            count = 0;
            count2 = 0;
            
            short int **buffer = malloc (buffsize * sizeof (short int *));
            //short int **buffer_flag = malloc (buffsize * sizeof (short int *));
    
            for (jj = 0; jj < buffsize; jj++) {
                buffer[jj] = malloc (buffsize * sizeof (short int));
                //buffer_flag[jj] = malloc (buffsize * sizeof (short int));
                for (ii = 0; ii < buffsize; ii++) {
                    buffer[jj][ii] = 0;
                    //buffer_flag[jj][ii] = 0;
                }
            }
    
            un_count = 0;
            fits[j].num = 0;
    
            while (npixels > 0)
            {
                nbuffer = npixels;
                nbuffer_flag = npixels;
                if (npixels > buffsize)
                    nbuffer = buffsize;     // read as many pixels as will fit in buffer
                    nbuffer_flag = buffsize;     // read as many pixels as will fit in buffer
        
                // Note that even though the FITS images contains unsigned integer
                // pixel values (or more accurately, signed integer pixels with
                // a bias of 32768),  this routine is reading the values into a
                // float array.   Cfitsio automatically performs the datatype
                // conversion in cases like this.
        
                if ( fits_read_img(fptr, TSHORT, fpixel, nbuffer, &nullval,
                           buffer[count], &anynull, &status) )
                    printerror( status );
                
                /*if ( fits_read_img(fptr2, TSHORT, fpixel2, nbuffer_flag, &nullval2,
                                   buffer_flag[count2], &anynull2, &status2) )
                    printerror( status2 );
                
                
                if (count != count2) {
                    printf("Row count for mask file and flag file are not matching!\n");
                    exit(1);
                }
                
                if (fpixel != fpixel2) {
                    printf("Pixel count for mask file and flag file are not matching!\n");
                    exit(1);
                }*/
                
                for (ii = 0; ii < buffsize; ii++) {
                    if (buffer[count][ii] == 0) {// && buffer_flag[count2][ii] == 0) {
                        un_count++;
                        fits[j].num++;
                    }
                }
        
                npixels -= nbuffer;    // increment remaining number of pixels
                fpixel  += nbuffer;    // next pixel to be read in image
                fpixel2  += nbuffer;    // next pixel to be read in image
                
                count++;
                count2++;
            }
            
            tot_cou += fits[j].num;
            fits[j].eff_area = fits[j].num * pix_area;
    
            if ( fits_close_file(fptr, &status) )
                printerror( status );
            
            /*if ( fits_close_file(fptr2, &status2) )
                printerror( status2 );*/
         
            for (jj = 0; jj < buffsize; jj++) {
                free (buffer[jj]);
                //free (buffer_flag[jj]);
            }
    
            free (buffer);
            //free (buffer_flag);
            
        }
        
        // Create random maps.
        int **fits_cor_rannum = malloc (1 * sizeof (int *));
        
        for (k = 0; k < 1; k++) {
            fits_cor_rannum[k] = malloc (fitsindex * sizeof (int));
            for (j = 0; j < fitsindex; j++) {
                fits_cor_rannum[k][j] = 0;
            }
        }
        
        for (k = 0; k < 1; k++) {
            
            ran_check = 0;
        
            for (j = 0; j < fitsindex; j++) {
                fits[j].ran_num = (int)((rannum[i - 1] * fits[j].num) / tot_cou);
                ran_check += fits[j].ran_num;
            }
        
            if (ran_check > rannum[i - 1]) {
                do {
                    number = rand() % fitsindex;
                    fits[number].ran_num--;
                
                    ran_check = 0;
                    for (j = 0; j < fitsindex; j++) {
                        ran_check += fits[j].ran_num;
                    }
                } while (ran_check != rannum[i - 1]);
            }
        
            if (ran_check < rannum[i - 1]) {
                do {
                    number = rand() % fitsindex;
                    fits[number].ran_num++;
                
                    ran_check = 0;
                    for (j = 0; j < fitsindex; j++) {
                        ran_check += fits[j].ran_num;
                    }
                } while (ran_check != rannum[i - 1]);
            }
            
            for (j = 0; j < fitsindex; j++) {
                fits_cor_rannum[k][j] = fits[j].ran_num;
            }
        }
        
        for (j = 0; j < fitsindex; j++) {
            fprintf(unmask_file, "%s\t%lE\t", fits[j].name, fits[j].eff_area);
            fprintf(unmask_file, "%ld\n", fits_cor_rannum[0][j]);
        }
        
        
        fclose (lensname);
        fclose (unmask_file);
        fclose (flagname);
        free (fits);
        
        for (k = 0; k < 1; k++) {
            free (fits_cor_rannum[k]);
        }
        free (fits_cor_rannum);
    }
    
    return 0;
}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/
    
    
    if (status)
    {
        fits_report_error(stderr, status); /* print error report */
        
        exit( status );    /* terminate the program, returning error status */
    }
    return;
}