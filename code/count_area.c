/* Counts the number of pixels inside the survey area and compare it to the circles with different radii, 3,5,7,9 and 11 h^-1 Mpc.
 Input: Field # (argv[1]), pointing # (argv[2]), multiple of n lens (argv[3]), value of n (argv[4])
 Updated: May 14, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "cfitsio/fitsio.h"
#include "tool.h"
#include "cosmopar.h"

// Struct for the galaxy catalog.
typedef struct {
    char name[50];
    int num;
} namelens;

typedef struct {
    double ra, dec, ra_low, ra_up, dec_low, dec_up;
} fieldpos;

typedef struct {
    char id[20];
    double ra, dec, z, ang_bound[5];
    signed long long int pix_num[5], unmask_num[5];
} gallens;

typedef struct {
    char name[200], flag[200];
    long num;
    double eff_area;
    int ran_num;
} fitsin;

typedef struct {
    char id[20];
    double cen_cord[2];
} pointing;

typedef struct {
    int pixel_ra[21000], pixel_dec;
} match_pixel;

void printerror( int status);

int main(int argc, char *argv[]) {
    int i = 0, j = 0, k = 0, p = 0, q = 0, m = 0, s = 0, z_bincount = 0, r_bincount = 0, spot = 0, file_count = 0, field_num = 0, pt_int = 0, fieldindex = 0, field_int = 0, f_ind = 0, lens_num = 0, lensindex = 0, lgal_index = 0, ramesh_num = 0, decmesh_num = 0, r_ind = 0, d_ind = 0, lgal_eva = 0, lgal_num = 0, lgal_end = 0, z_range = 0, status = 0, status2 = 0, nfound = 0, anynull = 0, anynull2 = 0, fitsindex = 0, fitsindex2 = 0, number = 0, ra_half = 0, dec_half = 0, pt_dec_num = 0, pt_dec = 0, cut_pixel_ind = 0, r_counter = 0;
    double z_min = 0, z_width = 0, r_min = 0, r_width = 0, ramap_size = 0, decmap_size = 0, ramesh_size = 0, decmesh_size = 0, gl_a_sep = 0, ls_a_sep = 0, sdec = 0, cdec = 0, theta_ang = 0, cos_theta = 0, l_bound_ra = 0, u_bound_dec = 0, cut_pixel = 0, ra_map = 0, percent = 0;
    long naxes[2], fpixel, fpixel2, nbuffer, npixels, nbuffer_flag, npixels_flag, ii, jj, kk, ll, mm, count = 0, count2 = 0, tot_cou = 0;
    double cen_pix_cord[2], pix_del[2];
    char field_file[150], field_manual[150], lens_file[150], write_file[150], lgal_manual[10], lgal_manum[10], pt_manual[10], point_name[7];
    FILE *fieldname, *lensname, *writename;
    fitsfile *fptr;//, *fptr2;       /* pointer to the FITS file, defined in fitsio.h */
    const double pix_area = 5.166666789E-05 * 5.166666789E-05, d2r = M_PI / 180, pix_size = 5.166666789E-05, pix_size_rad = pix_size * d2r;
    const int buffsize = 21000, nbr_pt_num = 9;
    unsigned long un_count = 0;
    int fitsnum[4] = {72, 25, 49, 25};
    
    // Make sure the exact number of arguments are passed.
    if (argc != 5) {
        printf("Invalid passing arguments!\n");
        return 1;
    }
    
    // Input which wide field
    strcpy (field_manual, argv[1]);
    field_int = atoi(field_manual);
    f_ind = field_int - 1;
    
    // Pre-compute the angular separations corresponding to different physical scales at different redshift slices.
    z_bincount = 21;
    z_min = 0.2;
    z_width = 0.01;
    
    double *z_bin = malloc ((z_bincount) * sizeof (double));    // Redshift bins.
    double *Da_bin = malloc ((z_bincount) * sizeof (double));   // Angular diameter distance bins.
    double **ang_bin = malloc ((z_bincount) * sizeof (double *));    // Angles corresponding to the scales at that redshift, in arcsecs.
    
    // Create an array to store up the physical radius of the circles for overdensity calculations.
    r_bincount = 5;
    r_counter = r_bincount - 1;
    r_min = 3.;
    r_width = 2.;
    
    double *rad_bin = malloc ((r_bincount) * sizeof (double));  // Radius bins.
    
    for (j = 0; j < r_bincount; j++) {
        rad_bin[j] = r_min + j * r_width;
    }
    
    for (i = 0; i < z_bincount; i++) {
        z_bin[i] = z_min + i * z_width;
        Da_bin[i] = Da(z_bin[i]);
        ang_bin[i] = malloc ((r_bincount) * sizeof (double));
        for (j = 0; j < r_bincount; j++) {
            ang_bin[i][j] = rad_bin[j] / Da_bin[i] * 180 * 3600 / M_PI;     // in arcsecs
        }
    }
    
    spot = 0;
    file_count = -1;
    
    // Write in the position file of the wide fields.
    strcpy (field_file, field_link);
    strcat (field_file, "position.dat");
    
    fieldname = fopen (field_file, "r");
    
    // Make sure the field position file is present.
    if (!fieldname) {
        printf("Cannot find the file for position of the fields!\n");
        exit(1);
    }
    
    // Count number of lines in the field position file.
    field_num = count_line (spot, fieldname, file_count);
    
    // An array storing all file names.
    fieldpos *fpos = malloc (field_num * sizeof *fpos);
    
    fieldindex = 0;
    
    while (fscanf(fieldname, "%*s %lf %lf %lf %lf %lf %lf", &fpos[fieldindex].ra, &fpos[fieldindex].dec, &fpos[fieldindex].ra_low, &fpos[fieldindex].ra_up, &fpos[fieldindex].dec_low, &fpos[fieldindex].dec_up) == 6) {
        fieldindex++;
    }
    
    fclose (fieldname);
    
    // Open the pointing name file.
    strcpy (lens_file, gal_link);
    strcat (lens_file, "W");
    strcat (lens_file, field_manual);
    strcat (lens_file, "/lenscount.txt");
    
    lensname = fopen (lens_file, "r");
    
    // Make sure the lens name file is present.
    if (!lensname) {
        printf("Cannot find the lens name file!\n");
        exit(1);
    }
    
    spot = 0;
    file_count = -1;
    
    // Count numer of lines in the lens name file.
    lens_num = count_line (spot, lensname, file_count);
    
    // An array storing all file names.
    namelens *lens = malloc (lens_num * sizeof *lens);
    
    lensindex = 0;
    
    while (fscanf(lensname, "%d %s", &lens[lensindex].num, &lens[lensindex].name) == 2) {
        lensindex++;
    }
    
    // Make sure the pointing count is correct.
    if (lensindex != fitsnum[f_ind]) {
        printf("Pointing counts are not correct!\n");
        exit (1);
    }
    
    fclose (lensname);
    
    // Input which pointing
    strcpy (pt_manual, argv[2]);
    pt_int = atoi(pt_manual);
    
    // Make sure the pointing input is valid.
    if (pt_int < 0 | pt_int >= fitsnum[f_ind]) {
        printf("Invalid pointing input!  Value should be between 0 and %d.\n", fitsnum[f_ind] - 1);
        exit (1);
    }
    
    
    for (i = pt_int; i < pt_int + 1; i++) {
        
        strncpy (point_name, lens[i].name, 6);
        point_name[6] = 0;
        
        printf("\nField W%s, pointing %s.\n", field_manual, point_name);
        
        FILE *lensfile;
        char lens_str[150];
        
        strcpy (lens_str, gal_link);
        strcat (lens_str, "W");
        strcat (lens_str, field_manual);
        strcat (lens_str, "/");
        strcat (lens_str, lens[i].name);
        
        lensfile = fopen (lens_str, "r");
        
        // Make sure the lens name file is present.
        if (!lensfile) {
            printf("Cannot find the lens file %s!\n", lens[i].name);
            exit(1);
        }
        
        // Open up the information of the lens galaxies.
        gallens *lgal = malloc (lens[i].num * sizeof * lgal);
        
        lgal_index = 0;
        
        // lens id, ra, dec and redshift
        while (fscanf(lensfile, "%s %*s %lf %lf %lf", &lgal[lgal_index].id, &lgal[lgal_index].ra, &lgal[lgal_index].dec, &lgal[lgal_index].z) == 4) {
            
            z_range = min(z_bin, z_bincount, lgal[lgal_index].z);
            
            for (j = 0; j < r_bincount; j++) {
                lgal[lgal_index].ang_bound[j] = ang_bin[z_range][j];
            }
            
            for (k = 0; k < r_bincount; k++) {
                lgal[lgal_index].pix_num[k] = 0;
                lgal[lgal_index].unmask_num[k] = 0;
            }
            
            lgal_index++;
            // if (lgal_index == 3) break;
        }
        
        
        fclose (lensfile);
        
        // Construct the chaining mesh grids, in arcsecs. Grid size along dec is defined as usual. Grid width along RA is defined according to the dec of the extremes of the fields.
        /*ramap_size = fabs(fpos[f_ind].ra_up - fpos[f_ind].ra_low) * 3600;
        decmap_size = fabs(fpos[f_ind].dec_up - fpos[f_ind].dec_low) * 3600;

        // Number of meshes along RA and DEC, followed by the size of each grid along RA and DEC (convert to arcsecs). // Sizes in deg.
        ramesh_size = 0.01;
        decmesh_size = 0.01;
        
        ramesh_num = ramap_size / (3600 * ramesh_size);
        decmesh_num = decmap_size / (3600 * decmesh_size);
        
        // A map to store the particles in the corresponding mesh. Loop over RA first, then DEC.
        double *ra_center = malloc (ramesh_num * sizeof (double));
        double *dec_center = malloc (decmesh_num * sizeof (double));
        
        for (j = 0; j < ramesh_num; j++) {
            ra_center[j] = fpos[f_ind].ra_low + (j + 0.5) * ramesh_size;
        }
        for (j = 0; j < decmesh_num; j++) {
            dec_center[j] = fpos[f_ind].dec_low + (j + 0.5) * decmesh_size;
        }
        
        // Create a map to store up which mesh do the lens galaxies belong to, first loop over RA, then DEC. Another array is to store the index of the galaxies in the source list.
        int **lgal_count = malloc (decmesh_num * sizeof (int *));
        int **lgal_cum = malloc (decmesh_num * sizeof (int *));     // Cumulative number of lens galaxies.
        double ***lgal_coord = malloc (decmesh_num * sizeof (double **));     // Coordinates of the center of the grids.
        int **lgal_meshsort = malloc (lgal_index * sizeof (int *));
        
        for (k = 0; k < lgal_index; k++) {
            lgal_meshsort[k] = malloc (3 * sizeof (int));
            for (p = 0; p < 3; p++) {
                lgal_meshsort[k][p] = 0;    // [k][0] is mesh num in along RA, [k][1] is mesh num in along DEC, [k][2] is index of the galaxy.
            }
        }
        
        for (k = 0; k < decmesh_num; k++) {
            lgal_count[k] = malloc (ramesh_num * sizeof (int));
            lgal_cum[k] = malloc (ramesh_num * sizeof (int));
            lgal_coord[k] = malloc (ramesh_num * sizeof (double *));
            for (p = 0; p < ramesh_num; p++) {
                lgal_count[k][p] = 0;
                lgal_cum[k][p] = 0;
                lgal_coord[k][p] = malloc (2 * sizeof (double));
                for (q = 0; q < 2; q++) {
                    lgal_coord[k][p][q] = 0;
                }
            }
        }
        
        // Store up the number of counts in the mesh grid
        for (k = 0; k < lgal_index; k++) {
            r_ind = min(ra_center, ramesh_num, lgal[k].ra);
            d_ind = min(dec_center, decmesh_num, lgal[k].dec);
            lgal_meshsort[k][0] = r_ind;        // mesh num along RA
            lgal_meshsort[k][1] = d_ind;        // mesh num along DEC
            lgal_meshsort[k][2] = k;            // index of the galaxy
            lgal_count[d_ind][r_ind]++;         // Counts of lens in each mesh.
        }
        
        // Cumulative number of galaxies in each grid.
        for (k = 0; k < 1; k++) {
            lgal_cum[k][0] = lgal_count[0][0];
            for (p = 1; p < ramesh_num; p++) {
                lgal_cum[k][p] = lgal_cum[k][p - 1] + lgal_count[k][p];
            }
        }
        
        for (k = 1; k < decmesh_num; k++) {
            lgal_cum[k][0] = lgal_cum[k - 1][ramesh_num - 1] + lgal_count[k][0];
            for (p = 1; p < ramesh_num; p++) {
                lgal_cum[k][p] = lgal_cum[k][p - 1] + lgal_count[k][p];
            }
        }
        
        // Create the mesh map.
        for (k = 0; k < decmesh_num; k++) {
            for (p = 0; p < ramesh_num; p++) {
                lgal_coord[k][p][0] = ra_center[p];
                lgal_coord[k][p][1] = dec_center[k];
            }
        }
        
        // Sort according to grid of RA; if the same, then sort according to grid of DEC; if still the same, then sort by galaxy index in the lens galaxy catalog.
        qsort(lgal_meshsort, lgal_index, sizeof lgal_meshsort[0], compare);*/
        
        FILE *fitsnamef;
        char fitsname_file[150];
        
        strcpy (fitsname_file, HD_mask_link);
        strcat (fitsname_file, "W");
        strcat (fitsname_file, field_manual);
        strcat (fitsname_file, "/");
        strcat (fitsname_file, point_name);
        strcat (fitsname_file, "_fitsname.txt");
        
        fitsnamef = fopen (fitsname_file, "r");
        
        // Make sure the FITS name file is present.
        if (!fitsnamef) {
            printf("Cannot find the FITS name file %s!\n", point_name);
            exit (1);
        }
        
        fitsin *fits = malloc (nbr_pt_num * sizeof *fits);
        
        fitsindex = 0;
        
        while (fscanf (fitsnamef, "%s", &fits[fitsindex].name) == 1) {
            fitsindex++;
        }
        
        fclose (fitsnamef);
        
        if (fitsindex != nbr_pt_num) {
            printf("Something wrong with reading neighboring pointing file!\n");
            exit (1);
        }
        
        pointing *pt = malloc (nbr_pt_num * sizeof *pt);
        
        // Read in the RA and Dec of the centers of the pointings.
        for (j = 0; j < nbr_pt_num; j++) {
            
            char fitsname[200];
            
            strcpy (fitsname, HD_mask_link);
            strcat (fitsname, "W");
            strcat (fitsname, field_manual);
            strcat (fitsname, "/");
            strcat (fitsname, fits[j].name);
            
            status = 0;
            
            if ( fits_open_file(&fptr, fitsname, READONLY, &status) )
                printerror( status );
            
            // read the NAXIS1 and NAXIS2 keyword to get image size
            if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
                printerror( status );
            
            // Get the ra and dec of the central point.
            if ( fits_read_keys_dbl(fptr, "CRVAL", 1, 2, pt[j].cen_cord, &nfound, &status) )
                printerror( status );
            
            // Get the step in angular scale.
            if ( fits_read_keys_dbl(fptr, "CDELT", 1, 2, pix_del, &nfound, &status) )
                printerror( status );
            
            if ( fits_close_file(fptr, &status) )
                printerror( status );
        }
        
        /*for (j = 0; j < nbr_pt_num; j++) {
            printf("%d\t%s\t%g\t%g\n", j, fits[j].name, pt[j].cen_cord[0], pt[j].cen_cord[1]);
        }*/
        
        pt_dec_num = nbr_pt_num/3;
        
        // Dec for each row. Pointings 0,1,2 share the same dec_temp, 3,4,5 and 6,7,8 share the other two.
        double **dec_temp = malloc (pt_dec_num * sizeof (double *));
        
        // Half the number of pixels in a row / column
        ra_half = naxes[0] / 2;
        dec_half = naxes[1] / 2;
        
        for (j = 0; j < pt_dec_num; j++) {
            dec_temp[j] = malloc (buffsize * sizeof (double));
            for (ii = 0; ii < buffsize; ii++) {
                dec_temp[j][ii] = pt[j * 3].cen_cord[1] + (ii - dec_half + 0.5) * pix_size;    // Along declination, the pixel size equals the dec difference.
            }
        }
        
        // Find the overlapping vertical pixels first, only need to do for pointings 0,1,3,4,6 and 7.
        theta_ang = (buffsize - dec_half - 0.5) * pix_size;
        cos_theta = cos(theta_ang * d2r);   // Angular separation of the boundary of the pointing to the center pixel with the same declination.
        
        match_pixel *match = malloc (nbr_pt_num * sizeof *match);
        
        // Mark the pixel number at which the overlapping occurs.
        for (j = 0; j < nbr_pt_num; j++) {  // nbr_pt_num
            
            // Find which dec row does the pointing belong to.
            pt_dec = j / 3;
            
            // Find the overlapping pixels in the vertical direction.
            if (j == 2 || j == 5 || j == 8) {
                for (k = 0; k < buffsize; k++) {
                    match[j].pixel_ra[k] = buffsize;
                }
            }
            
            if (j != 2 && j != 5 && j != 8) {
                // Loop over each dec row
                for (k = 0; k < buffsize; k++) {   // buffsize
                    
                    sdec = pow(sin(dec_temp[pt_dec][k] * d2r), 2);
                    cdec = pow(cos(dec_temp[pt_dec][k] * d2r), 2);
                    // RA of the left boundary of j + 1 pointing.
                    l_bound_ra = pt[j + 1].cen_cord[0] + acos((cos_theta - sdec) / cdec) / d2r;
                    
                    // Find out the index of the RA pixel in j pointing that should be cut away.
                    match[j].pixel_ra[k] = round(10500.5 + acos(cos((pt[j].cen_cord[0] - l_bound_ra) * d2r) * cdec + sdec) / pix_size_rad) - 1;     // Last ra index to be considered.

                    //printf("%d\t%.10lf\t%.10lf\t%.10lf\t%g\t%g\t%g\t%d\n", j, l_bound_ra, pt[j].cen_cord[0], dec_temp[pt_dec][k], sdec, cdec, pix_size_rad, cut_pixel_ind);
                }
            }
            
            
            // Find the overlapping pixels in the horizontal direction.
            if (j >= 6) {
                match[j].pixel_dec = -1;
            }
            
            if (j <= 5) {
                u_bound_dec = pt[j + 3].cen_cord[1] + theta_ang;
                match[j].pixel_dec = round(10500.5 - (pt[j].cen_cord[1] - u_bound_dec) / pix_size) - 1;     // This is the last dec index to be considered.
            }
            
        }  // j-loop, over 9 FITS files.

        
        for (j = 0; j < lgal_index; j++) {
            for (k = 0; k < r_bincount; k++) {
                if (lgal[j].pix_num[k] != 0 || lgal[j].unmask_num[k] != 0) {
                    printf("Pixel counter is not properly initialized!\n");
                    exit (1);
                }
            }
        }
        
        // Multiple numbers of lens galaxies, starting from 0.
        strcpy (lgal_manual, argv[3]);
        strcpy (lgal_manum, argv[4]);
        lgal_eva = atoi(lgal_manum);
        lgal_num = atoi(lgal_manual) * lgal_eva;
        
        lgal_end = lgal_index;
        
        if ((atoi(lgal_manual) + 1) * lgal_eva < lgal_end) {
            lgal_end = (atoi(lgal_manual) + 1) * lgal_eva;
        }
        
        if (lgal_num > lgal_end) {
            printf("Input pointing num is not correct!\n");
            exit (1);
        }
        
        // Read in the FITS files.
        for (j = 0; j < nbr_pt_num; j++) {  //nbr_pt_num
            
            printf("%d FITS file started.\n", j + 1);
            
            // Find which dec row does the pointing belong to.
            pt_dec = j / 3;
            
            float nullval = 0, nullval2 = 0;
            char fitsname[200];
            
            strcpy (fitsname, HD_mask_link);
            strcat (fitsname, "W");
            strcat (fitsname, field_manual);
            strcat (fitsname, "/");
            strcat (fitsname, fits[j].name);
            
            status = 0;
            
            if ( fits_open_file(&fptr, fitsname, READONLY, &status) )
                printerror( status );
            
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
            
            for (ii = 0; ii < buffsize; ii++) {
                pix_val[ii] = malloc (buffsize * sizeof (short int));
                for (jj = 0; jj < buffsize; jj++) {
                    pix_val[ii][jj] = 0;
                }
            }
            
            count = 0;
            count2 = 0;
            
            while (npixels > 0) {
                nbuffer = npixels;
                nbuffer_flag = npixels;
                if (npixels > buffsize)
                    nbuffer = buffsize;     // read as many pixels as will fit in buffer
                nbuffer_flag = buffsize;
                
                if ( fits_read_img(fptr, TSHORT, fpixel, nbuffer, &nullval, pix_val[count], &anynull, &status) )
                    printerror( status );
                
                npixels -= nbuffer;    // increment remaining number of pixels
                fpixel  += nbuffer;    // next pixel to be read in image
                fpixel2  += nbuffer;
                
                count++;
                count2++;
                
            }   // while-loop, loop over the pixels.
            
            // Make sure the right number of pixels have been read.
            if (count != buffsize) {
                printf("Error in reading in pixel maps!\n");
                exit(1);
            }
            
            // Write in the RA and Dec of the maps.
            // double **ra_map = malloc (buffsize * sizeof (double *));
            double *ra_temp = malloc (buffsize * sizeof (double));
            
            for (jj = 0; jj < buffsize; jj++) {
                ra_temp[jj] = cos(((jj - ra_half + 0.5) * pix_size) * d2r);
            }
            
            printf("Read fits files begin!\n");
            
            // rows (looping dec), in deg
            for (ii = 0; ii < buffsize; ii++) {    //buffsize
                
                if (ii % 5000 == 0)  printf("%d out of %d rows have been counted.\n", ii, buffsize);
                
                //ra_map[ii] = malloc (buffsize * sizeof (double));
                sdec = pow(sin(dec_temp[pt_dec][ii] * d2r), 2);
                cdec = pow(cos(dec_temp[pt_dec][ii] * d2r), 2);
                
                // columns (looping RA), RA for each pixel (not the same every column)
                for (jj = 0; jj < buffsize; jj++) {    //buffsize
                    if (jj >= ra_half) {
                        //ra_map[ii][jj] = cen_pix_cord[0] - acos((ra_temp[jj] - sdec) / cdec) / d2r;
                        ra_map = cen_pix_cord[0] - acos((ra_temp[jj] - sdec) / cdec) / d2r;
                    }
                    else {
                        //ra_map[ii][jj] = cen_pix_cord[0] + acos((ra_temp[jj] - sdec) / cdec) / d2r;
                        ra_map = cen_pix_cord[0] + acos((ra_temp[jj] - sdec) / cdec) / d2r;
                    }
                    
                    for (k = lgal_num; k < lgal_end; k++) {       // lgal_index
                        // Separation between the pixel and the lens, in arcsec
                        ls_a_sep = ang_sep (lgal[k].ra, lgal[k].dec, ra_map, dec_temp[pt_dec][ii]);
                        
                        // Check if the pixel is within the largest circle (to facilitate search.
                        if (ls_a_sep <= lgal[k].ang_bound[r_counter]) {
                            
                            // Count the number of pixels in the circle as well as the unmasked regions.
                            if (jj <= match[j].pixel_ra[ii] && ii > match[j].pixel_dec) {
                                lgal[k].pix_num[r_counter]++;
                                
                                // For other radii.
                                for (m = 0; m < r_counter; m++) {
                                    if (ls_a_sep <= lgal[k].ang_bound[m]) {
                                        lgal[k].pix_num[m]++;
                                    }
                                }
                            }
                            
                            // Count the unmasked pixels
                            if (pix_val[ii][jj] == 0) {
                                lgal[k].unmask_num[r_counter]++;
                                
                                for (m = 0; m < r_counter; m++) {
                                    if (ls_a_sep <= lgal[k].ang_bound[m]) {
                                        lgal[k].unmask_num[m]++;
                                    }
                                }
                            }
                            
                            /*printf("%d\t%d\t%d\t%d\t%lE\t%lE\t%lE\t%lE\t%d\t%d\t%lE\n", k, j, ii, jj, lgal[k].ra, lgal[k].dec, ra_map, dec_temp[pt_dec][ii], match[j].pixel_ra[ii], match[j].pixel_dec, ls_a_sep);
                            
                            for (m = 0; m < r_bincount; m++) {
                                printf("%d\t%lE\t%d\t%d\n", m, lgal[k].ang_bound[m], lgal[k].pix_num[m], lgal[k].unmask_num[m]);
                            }
                            printf("\n");*/
                        }
                    }
                }   // jj-loop,  looping RA
            }       // ii-loop,  looping DEC
            
            
            
            
            if ( fits_close_file(fptr, &status) )
                printerror( status );

            /*for (ii = 0; ii < buffsize; ii++) {
                free (ra_map[ii]);
            }*/
            
            free (ra_temp);
            free (pix_val);
            // free (ra_map);
            
        }   // j-loop, number of FITS files.
        
        FILE *outfile;
        char outname[150];
        
        strcpy (outname, gal_link);
        strcat (outname, "W");
        strcat (outname, field_manual);
        strcat (outname, "/area_count/");
        strcat (outname, point_name);
        strcat (outname, "_");
        strcat (outname, lgal_manual);
        strcat (outname, "_area_count.dat");
        
        outfile = fopen(outname, "w");
        
        for (k = lgal_num; k < lgal_end; k++) {
            fprintf(outfile, "%s\t", lgal[k].id);
            for (m = 0; m < r_bincount; m++) {
                percent = (lgal[k].unmask_num[m] + 0.) / lgal[k].pix_num[m];
                if (percent < 0 || percent > 1) {
                    printf("%d\t%s pixel count has some mistakes.\nPixel count = %d\nUnmasked pixel count = %d\n", k, lgal[k].id, lgal[k].pix_num[m], lgal[k].unmask_num[m]);
                    exit (1);
                }
                printf("%d\t%d\t%llu\t%llu\t%lE\n", k, m, lgal[k].unmask_num[m], lgal[k].pix_num[m], percent);
            }
            for (m = 0; m < r_bincount; m++) {
                fprintf(outfile, "%llu\t", lgal[k].unmask_num[m]);
            }
            for (m = 0; m < r_bincount; m++) {
                fprintf(outfile, "%llu\t", lgal[k].pix_num[m]);
            }
            for (m = 0; m < r_bincount; m++) {
                fprintf(outfile, "%lE\t", (lgal[k].unmask_num[m] + 0.) / lgal[k].pix_num[m]);
            }
            fprintf(outfile, "\n");
        }
        
        fclose (outfile);
        
        for (j = 0; j < pt_dec_num; j++) {
            free (dec_temp[j]);
        }
        
        /*for (k = 0; k < decmesh_num; k++) {
            free(lgal_count[k]);
            free(lgal_cum[k]);
            for (p = 0; p < ramesh_num; p++) {
                free(lgal_coord[k][p]);
            }
            free(lgal_coord[k]);
        }
        
        for (k = 0; k < lgal_index; k++) {
            free (lgal_meshsort[k]);
        }*/
        
        
        free (fits);
        /*free (ra_center);
        free (dec_center);
        free (lgal_meshsort);
        free (lgal_cum);
        free (lgal_count);
        free (lgal_coord);*/
        free (lgal);
        free (pt);
        free (dec_temp);
        free (match);
    }   // i-loop, pointings
    
    
    
    
    
    for (i = 0; i < r_bincount; i++) {
        free (ang_bin[i]);
    }
    
    free (lens);
    free (fpos);
    free (ang_bin);
    free (z_bin);
    free (rad_bin);
    free (Da_bin);
    
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
