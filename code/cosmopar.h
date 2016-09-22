// Header file for all the cosmological parameters and functions. Assume FLAT cosmology
// Updated: May 11, 2015

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nr.h"

#ifndef _COSMOPAR_H_
#define _COSMOPAR_H_

// Physical constants.
#define c          2.99792458e+10  // Speed of light (cm/s)
#define G_N        6.6738480e-8    // Newtonian Gravitational constant (cm^3 g^-1 s^-2)
#define Msun_2_g   1.9891e+33      // Mass of the Sun in gram
#define Mpc_2_cm   3.08567758e+24  // Mpc in cm
#define pc_2_cm    3.08567758e+18  // pc in cm
#define Mpc_2_pc   1e+6            // Mpc in pc
#define kpc_2_pc   1e+3            // kpc in pc

// Cosmological parameters.
#define oM         0.27            // Matter density
#define oL         0.73            // Lambda density
#define h          0.70            // Dimensionless Hubble constant
#define H0         7.0e+6          // Hubble constant (in cm s^-1 Mpc^-1)
#define sig_8      0.81            // Sigma-8
#define w_0        -1.0            // Equation of State

double evolution(double z);
double simpson(double z);
double Dc(double z);
double Da(double z);
double Dl(double z);
double Da12(double z1, double z2);
double exp_Da12(double z_hist[], double lPDF[], double sPDF[], int PDF_size);
double exp_invsig(double Da_hist[], double **Da12_hist, double lPDF[], double sPDF[], int PDF_size);
double exp_invsigsq(double Da_hist[], double **Da12_hist, double lPDF[], double sPDF[], int PDF_size);
double Sig_cri(double z1, double z2);
double H_z_sq(double z);
double rho_crit(double z);
double PFOF(double Dc_hist[], double lPDF_1[], double lPDF_2[], int PDF_size, double V_0);
void exp_invsig2(double Da_hist[], double **Da12_hist, double lPDF[], double sPDF[], int PDF_size, double phyr_width[], int phy_bin, double ls_a_sep_rad, double sum_eff_ls[]);
void exp_invsigsq2(double Da_hist[], double **Da12_hist, double lPDF[], double sPDF[], int PDF_size, double phyr_width[], int phy_bin, double ls_a_sep_rad, double sum_eff_weight[]);
void numsrcs(double Da_hist[], double lPDF[], int PDF_size, double phyr_width[], int phy_bin, double ls_a_sep_rad, double srcsnum[]);

#endif
