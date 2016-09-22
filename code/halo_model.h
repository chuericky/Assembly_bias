// Header file for halo fitting in Valender et. al. 2014.
// Updated: Nov 12, 2015

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nr.h"

#ifndef _HALO_MODEL_H_
#define _HALO_MODEL_H_

#define delta_crit  1.68647

double rad_200(double rho_cri, double M_200);
double Delta_vir(double x);
double con_NFW(double M_200, double z);
double scale_rad(double rho_cri, double M_200, double z);
double scale_rad_twopara(double rho_cri, double M_200, double z, double c_NFW);
double char_den(double c_NFW);
double bias_v(double v);
double f_v(double v);
double mass_lum_corr (char* bin_name, double m_opt);
double mass_lum_corr2 (char* bin_name, double m_opt);
double mass_lum_size_gmass_corr (char* bin_name, double m_opt);
void dsig_baryon(double star_mass, double r[], double sig_bar[], int N_size);
void dsig_cen1h(double M_200, double z, double x[], double sig_cen_1h[], int N_size);
void dsig_cen1h_twopara(double M_200, double z, double c_NFW, double r_phy[], double sig_cen_1h[], int N_size);
void dsig_sat1h(double rho_cri, double char_den, double scale_rad, double x[], double sig_sat_1h[], double M_200, double p_r[], int N_size);

#endif