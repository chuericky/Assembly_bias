// This program projects random galaxies in the unmasked regions.
// Field num (argv[1]), CCD num (argv[2]), all or iso (argv[3])
// Updated: Sep 13, 2015

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "tool.h"
#include "cosmopar.h"
#include <time.h>
#include "cfitsio/fitsio.h"

typedef struct {
    char name[50], flag[50];
    int num[5];
} maskstruct;

void printerror( int status);

int main(int argc, char *argv[])
{
    fitsfile *fptr, *fptr2;       // pointer to the FITS file, defined in fitsio.h
    int status, status2, nfound, anynull, anynull2, buffsize = 21000, count, count2, dec_half = 0, ra_half = 0, ran_count = 0, ran_keepcheck = 0, dec_ind = 0, ra_ind = 0, field_num = 0, mask_num = 0, spot = 0, file_count = 0, maskfile_num = 0, maskindex = 0, mask_ind = 0;
    long i, j, naxes[2], fpixel, fpixel2, nbuffer, nbuffer_flag, npixels;
    double cen_pix_cord[2], pix_del[2], d2r = M_PI / 180, sdec = 0, cdec = 0, ra_ran = 0, dec_ran = 0, min_dec = 0, max_dec = 0, min_ra = 0, max_ra = 0;
    float nullval, nullval2;
    char filename[150], flagfile[150], maskfile[150], field_man[2], CCD_man[3], subCCD[7], s[2], gal_man[7];
    FILE *maskname, *flagname;
    
    // Make sure the exact number of arguments are passed.
    if (argc != 4) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    strcpy (maskfile, mask_link);
    strcat (maskfile, "W");
    strcpy (field_man, argv[1]);
    strcat (maskfile, field_man);
    strcat (maskfile, "/unmask_num_");
    strcpy (gal_man, argv[3]);
    strcat (maskfile, gal_man);
    strcat (maskfile, ".dat");
    
    maskname = fopen (maskfile, "r");
    
    // Make sure the mask name file is present.
    if (!maskname) {
        printf("Cannot find the mask name file!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    maskfile_num = count_line (spot, maskname, file_count);
    
    maskstruct *masks = malloc (maskfile_num * sizeof *masks);
    
    maskindex = 0;
    
    /*while (fscanf(maskname, "%s %s %*f %d %d %d %d %d", &masks[maskindex].name, &masks[maskindex].flag, &masks[maskindex].num[0], &masks[maskindex].num[1], &masks[maskindex].num[2], &masks[maskindex].num[3], &masks[maskindex].num[4]) == 7) {
        maskindex++;
    }*/
    while (fscanf(maskname, "%s %*f %d", &masks[maskindex].name, &masks[maskindex].num[0]) == 2) {
        maskindex++;
    }
    
    fclose (maskname);
    
    strcpy (CCD_man, argv[2]);
    mask_ind = atoi(CCD_man) - 1;
    
    if (mask_ind + 1 > maskindex) {
        printf("Invalid CCD index input!\n");
        exit (1);
    }
    
    printf("%s\n", masks[mask_ind].name);
    
    
    strcpy (filename, mask_link);
    strcat (filename, "W");
    strcat (filename, field_man);
    strcat (filename, "/");
    strcat (filename, masks[mask_ind].name);
    
    /*strcpy (flagfile, mask_link);
    strcat (flagfile, "W");
    strcat (flagfile, field_man);
    strcat (flagfile, "/");
    strcat (flagfile, masks[mask_ind].flag);*/
    
    strncpy (subCCD, masks[mask_ind].name, 6);
    
    /* initialize random seed: */
    srand ( time(NULL) );
    
    status = 0;
    status2 = 0;
    
    if ( fits_open_file(&fptr, filename, READONLY, &status) )
        printerror( status );
    
    /*if ( fits_open_file(&fptr2, flagfile, READONLY, &status2) )
        printerror( status2 );
    
    if ( fits_movnam_hdu(fptr2, 0, "COMPRESSED_IMAGE", 0, &status2) )
        printerror( status2 );*/
    
    // read the NAXIS1 and NAXIS2 keyword to get image size
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
        printerror( status );
    
    // Get the ra and dec of the central point.
    if ( fits_read_keys_dbl(fptr, "CRVAL", 1, 2, cen_pix_cord, &nfound, &status) )
        printerror( status );
    
    // Get the step in angular scale.
    if ( fits_read_keys_dbl(fptr, "CDELT", 1, 2, pix_del, &nfound, &status) )
        printerror( status );
    
    npixels  = naxes[0] * naxes[1];         // number of pixels in the image
    fpixel   = 1;
    fpixel2  = 1;
    nullval  = 0;                // don't check for null values in the image
    nullval2 = 0;
    
    short int **pix_val = malloc (buffsize * sizeof (short int *));
    //short int **flag_val = malloc (buffsize * sizeof (short int *));
    
    for (i = 0; i < buffsize; i++) {
        pix_val[i] = malloc (buffsize * sizeof (short int));
        //flag_val[i] = malloc (buffsize * sizeof (short int));
        for (j = 0; j < buffsize; j++) {
            pix_val[i][j] = 0;
            //flag_val[i][j] = 0;
        }
    }
    
    count = 0;
    count2 = 0;
    
    while (npixels > 0)
    {
        nbuffer = npixels;
        nbuffer_flag = npixels;
        if (npixels > buffsize)
            nbuffer = buffsize;     // read as many pixels as will fit in buffer
            nbuffer_flag = buffsize;
        
        if ( fits_read_img(fptr, TSHORT, fpixel, nbuffer, &nullval, pix_val[count], &anynull, &status) )
            printerror( status );
        
        /*if ( fits_read_img(fptr2, TSHORT, fpixel2, nbuffer_flag, &nullval2,
                           flag_val[count2], &anynull2, &status2) )
            printerror( status2 );
        
        if (count != count2) {
            printf("Row count for mask file and flag file are not matching!\n");
            exit(1);
        }
        
        if (fpixel != fpixel2) {
            printf("Pixel count for mask file and flag file are not matching!\n");
            exit(1);
        }*/
        
        
        npixels -= nbuffer;    // increment remaining number of pixels
        fpixel  += nbuffer;    // next pixel to be read in image
        fpixel2  += nbuffer;
        
        count++;
        count2++;
    }
    
    // Make sure the right number of pixels have been read.
    if (count != buffsize) {
        printf("Error in reading in pixel maps!\n");
        exit(1);
    }
    
    ra_half = naxes[0] / 2;
    dec_half = naxes[1] / 2;
    
    // Write in the RA and Dec of the maps.
    double **ra_map = malloc (buffsize * sizeof (double *));
    double *ra_temp = malloc (buffsize * sizeof (double));
    double *dec_temp = malloc (buffsize * sizeof (double));
    
    for (j = 0; j < buffsize; j++) {
        ra_temp[j] = cos(((j - ra_half + 0.5) * pix_del[0]) * d2r);
    }
    
    printf("Read fits files begin!\n");
    
    // rows (looping dec), in deg
    for (i = 0; i < buffsize; i++) {
        
        ra_map[i] = malloc (buffsize * sizeof (double));
        dec_temp[i] = cen_pix_cord[1] + (i - dec_half + 0.5) * pix_del[1];
        sdec = pow(sin(dec_temp[i] * d2r), 2);
        cdec = pow(cos(dec_temp[i] * d2r), 2);
        
        // columns (looping RA)
        for (j = 0; j < buffsize; j++) {
            if (j >= ra_half) {
                ra_map[i][j] = cen_pix_cord[0] - acos((ra_temp[j] - sdec) / cdec) / d2r;
            }
            else {
                ra_map[i][j] = cen_pix_cord[0] + acos((ra_temp[j] - sdec) / cdec) / d2r;
            }
        }
    }
    
    free (ra_temp);
    
    // To get to know the range where random numbers are drawn.
    double *dec_compare = malloc (2 * sizeof (double));
    double *ra_compare = malloc (4 * sizeof (double));
    
    dec_compare[0] = dec_temp[0];
    dec_compare[1] = dec_temp[buffsize - 1];
    ra_compare[0] = ra_map[0][0];
    ra_compare[1] = ra_map[0][buffsize - 1];
    ra_compare[2] = ra_map[buffsize - 1][0];
    ra_compare[3] = ra_map[buffsize - 1][buffsize - 1];
    
    min_dec = min_array(dec_compare, 2);
    max_dec = max_array(dec_compare, 2);
    min_ra = min_array(ra_compare, 4);
    max_ra = max_array(ra_compare, 4);
    
    free (dec_compare);
    free (ra_compare);
    
    for (i = 0; i < 1; i++) {
        
        FILE *ran_file;
        char random_link[150];
        
        strcpy (random_link, ran_link);;
        strcat (random_link, "W");
        strcat (random_link, field_man);
        strcat (random_link, "/");
        strcat (random_link, subCCD);
        strcat (random_link, "_ran_");
        strcat (random_link, gal_man);
        strcat (random_link, ".dat");
        
        ran_file = fopen(random_link, "w");
        
        // Number of random galaxies that must form in that CCD.
        ran_count = masks[mask_ind].num[i];

        ran_keepcheck = 0;
    
        // Project random galaxies.
        do {
            dec_ran = fRand(min_dec, max_dec);
            ra_ran = fRand(min_ra, max_ra);
        
            // Nearest dec grid.
            dec_ind = near_dec(buffsize, dec_temp, dec_ran);
        
            double *ra_temp2 = malloc (buffsize * sizeof (double));
            for (j = 0; j < buffsize; j++) {
                ra_temp2[j] = ra_map[dec_ind][j];
            }
            // Nearest ra grid.
            ra_ind = near_ra(buffsize, ra_temp2, ra_ran);
        
            // Those in unmasked regions, and flag = 0
            if (pix_val[dec_ind][ra_ind] == 0) {// && flag_val[dec_ind][ra_ind] == 0) {
                fprintf(ran_file, "%f\t%f\n", ra_ran, dec_ran);
                ran_keepcheck++;
            }
        
            free (ra_temp2);
        } while(ran_keepcheck != ran_count);
        
        fclose (ran_file);
    }
    
    if ( fits_close_file(fptr, &status) )
        printerror( status );
    
    /*if ( fits_close_file(fptr2, &status2) )
        printerror( status2 );*/
    
    for (i = 0; i < buffsize; i++) {
        free (pix_val[i]);
        //free (flag_val[i]);
        free (ra_map[i]);
    }
    
    free (pix_val);
    //free (flag_val);
    free (ra_map);
    free (dec_temp);
    free (masks);
    
    return 0;
}

void printerror(int status) {
    if (status) {
        fits_report_error(stderr, status); /* print error report */
        exit( status );    /* terminate the program, returning error status */
    }
    return;
}