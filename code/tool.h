// Updated: Aug 30, 2016

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef _TOOL_H_
#define _TOOL_H_

#define hard_drive_link "/Volumes/Seagate Expansion Drive/Bolshoi_sim/"
#define mock_shear_link "/Users/rickyccy/Documents/Research_assemblybias/mock/shear_profile/"
#define density_link "/Users/rickyccy/Documents/Research_assemblybias/mock/particles/density/"
#define halo_link "/Users/rickyccy/Documents/Research_assemblybias/mock/halo/"
#define trial_link "/Users/rickyccy/Documents/Research_assemblybias/mock/particles/all/"
#define catalog_link "/Users/rickyccy/Documents/Research_assemblybias/Data/"
#define halocat_link "/Users/rickyccy/Documents/Research_assemblybias/mock/halo_catalog/"
#define extquar_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/quartiles/"
#define ext2pt_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/quartiles/2pt_calculation/"
#define halonbr_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/mock/halo_catalog/"
#define extmock_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/mock/"
#define catalog_multi_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/"
#define gal_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/"
#define gal_shear_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/shear/"
#define HD_mask_link "/Volumes/RICKYCHUE/Masks/"
#define lens_primary_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/lens_primary/"
#define overden_part_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/overdensity_subpart/"
#define mask_link "/Users/rickyccy/Documents/Research_assemblybias/Masks/"
#define lens_link "/Users/rickyccy/Documents/Research_assemblybias/Data/lens/"
#define srcs_link "/Users/rickyccy/Documents/Research_assemblybias/Data/source/"
#define srcs2_link "/Users/rickyccy/Documents/Research_assemblybias/Data/source2/"
#define ran_link "/Users/rickyccy/Documents/Research_assemblybias/random_map/"
#define lens2_link "/Users/rickyccy/Documents/Research_assemblybias/Data/lens/"
#define ran2_link "/Users/rickyccy/Documents/Research_assemblybias/random_map_try/"
#define twopt2_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt_try/"
#define twopt_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/2pt/"
#define mock_lens_link "/Users/rickyccy/Documents/Research_assemblybias/mock/"
#define mock_lensnoab_link "/Users/rickyccy/Documents/Research_assemblybias/mock/erased_assembly_bias_mocks/"
#define mock_ran_link "/Users/rickyccy/Documents/Research_assemblybias/mock/ran_position/"
#define mock_2pt_link "/Users/rickyccy/Documents/Research_assemblybias/mock/2pt/"
#define mock_2ptnoab_link "/Users/rickyccy/Documents/Research_assemblybias/mock/2pt_noab/"
#define lens_multi_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens/"
#define lens_good_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_isolated/good/"
#define lens_FOF_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_FOF/"
#define lens_FOFiso_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_FOF_iso/"
#define lens_without_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_without/"
#define lens_onlybright_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_onlybright/"
#define lens_without_onlybright_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_without_onlybright/"
#define lens_with_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_with/"
#define srcs_multi_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/source/"
#define shear_multi_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/"
#define shear_gal_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/"
#define shear2_gal_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal_2/"
#define shear_gal2_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/gal/FOF/"
#define cross_gal_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/other_lens/"
#define shear_multi2_link "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/shear/try/"
#define lensmultisortlink "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_sorted2/"
#define lensmultinotsortlink "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_notremovesat/"
#define lenssortlink "/Users/rickyccy/Documents/Research_assemblybias/Data/lens_sorted2/combined/"
#define lensisosortlink "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens/"
#define lensnocluslink "/Users/rickyccy/Documents/Research_assemblybias/Data/lens_notgroup/"
#define lensisolink "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_isolated/"
#define lensfurtherlink "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/lens_further/"
#define lensnoclus2link "/Users/rickyccy/Documents/Research_assemblybias/Data/lens_notgroup2/"
#define lensnoclusmultifxlink "/Users/rickyccy/Documents/Research_assemblybias/Data/lens_notgroup/lum_size_gmass/"
#define lensassembiaslink "/Users/rickyccy/Documents/Research_assemblybias/Data/lens_assembias/"
#define lensnoclus_widerlink "/Users/rickyccy/Documents/Research_assemblybias/Data/lens_notgroup_widerbin/"
#define src2_link "/Users/rickyccy/Documents/Research_assemblybias/Data/source2/combined/"
#define src3_link "/Users/rickyccy/Documents/Research_assemblybias/Data/source2/combined2/"
#define field_link "/Users/rickyccy/Documents/Research_assemblybias/"
#define code_link "/Users/rickyccy/Documents/Research_assemblybias/code/"
#define shear_link "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/"
#define shearmultifx_link "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/lum_size_gmass/"
#define shearmultifx2_link "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/combined3/lum_size_gmass/"
#define shearassembias_link "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_assembias/"
#define shearassembias_link_combined "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_assembias/combined/"
#define shearmultifx_link_combined "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/lum_size_gmass/combined/"
#define shear_link_combined "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/combined/"
#define shear_link3_combined "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/combined2/"
#define shear_link_precombined "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/pre-combine/"
#define shear_link_assemcombined "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/assembly_bias/"
#define shear_link_assemcombined2 "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/assembly_bias/combined/"
#define shear_link_precombined_combine "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/pre-combine/combined/"
#define shear_link_cov "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/pre-combine/combined/cov/"
#define shear_assem_link_cov "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/assembly_bias/combined/cov/"
#define shear_link_precombined_lum "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/pre-combine/lum_size_gmass/"
#define shear_link4_combined "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile/combined3/"
#define shear_link2_combined "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile2/combined/"
#define shear_widerlink "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile_widerbin/"
#define shear_link2 "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile2/"
#define shear_widerlink_combined "/Users/rickyccy/Documents/Research_assemblybias/Data/shear_profile_widerbin/combined/"
#define fit_link "/Users/rickyccy/Documents/Research_assemblybias/Data/halo_fitting/"
#define cluster_link "/Users/rickyccy/Documents/Research_assemblybias/cfhtlens_3DMF_clusters/"
#define parameter_link "/Users/rickyccy/Documents/Research_assemblybias/Data/parameter_space/"

FILE *open_check(char *filename, char *mode);
int count_line(int spot, FILE *countfile, int line);
int min(double arr[], int size, double subtract);
int max(double arr[], int size, double subtract);
int ragrid_number(double r_search, double delta, double rasize);
int compare (const void *pa, const void *pb);
int compare2 (const void * a, const void * b);
int grid2consider (int **gridscon, int ra_grid, int ra_mesh, int dec_grid, int dec_mesh);
int lenslumclass (char* lenscolor, char* lum_class);
int lensmulticlass (char* lenscolor, int color_bin, char* lum_class, int lum_bin, char* size_class, int size_bin, char* mass_class, int mass_bin);
int lenslumclass_wider (char* lenscolor, char* lum_class);
int cmp(const void *x, const void *y);
int max_int(int arr[], int size);
int near_dec(int column, double *dec_temp, double gal_dec);
int near_ra(int column, double *ra_temp, double gal_ra);
int jackarea(char field[]);
int jackarea2(char field[]);
int jackarea3(char field[]);
double arr_mean(double arr[], int size);
double arr_var(double arr[], double avg, int size);
double max_compare(double a, double b);
double min_value(double arr[], int size);
double max_value(double arr[], int size);
double min_array(double arr[], int size);
double max_array(double arr[], int size);
double ang_sep (double ra_cen, double dec_cen, double ra, double dec);
double ang_sep_rad (double ra_cen, double dec_cen, double ra, double dec);
double ori_ang (double a1, double d1, double a2, double d2, double ang_sep);
double fRand(double fMin, double fMax);
double dist_periodic(double x1, double y1, double x2, double y2, double L);
void lens_filename (int lenslumclass, char name[]);
void lens_filename2 (int lenslumclass, char name[]);
void lens_filename_assembias (int lenslumclass, char name[]);
void lens_filenamemultifx (int lenslumclass, char name[], int color_bin, int lum_bin, int size_bin, int mass_bin);
void lens_filename2_wider (int lenslumclass_wider, char name[]);
void min_matrix(int row, int column, double **arr, int *row_ind, int *col_ind, double subtract);

#endif
