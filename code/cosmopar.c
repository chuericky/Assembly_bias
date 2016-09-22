// Updated: May 11, 2015

#include "cosmopar.h"
#include "nr.h"

double evolution(double z) {
    // E(z) in the integrand.
    return pow((oM * pow((1. + z), 3) + oL), -0.5);
}

double Dc(double z) {
    // Line of sight comoving distance from now to z. It is the same as traverse comoving distance if the cosmology is FLAT! Unit in Mpc.
    return (c / H0) * qromb(evolution, 0, z);
}

double Da(double z) {
    // Angular diameter distance from now to z. Unit in Mpc.
    return Dc(z) / (1. + z);
}

double Dl(double z) {
    // Luminosity distance from now to z. Unit in Mpc.
    return Dc(z) * (1. + z);
}

double Da12(double z1, double z2) {
    // Angular diameter distance between z1 and z2, (z2 > z1). Unit in Mpc.
    return (Dc(z2) - Dc(z1)) / (1. + z2);
}

double exp_Da12(double z_hist[], double lPDF[], double sPDF[], int PDF_size) {
    // A function to calculate the expected value of angular diameter distance Da12, given the PDFs of lens and source galaxies. Unit in Mpc
    int i = 0, j = 0;
    double sum = 0;
    
    for (i = 0; i < PDF_size; i++) {
        for (j = i + 1; j < PDF_size; j++) {
            sum += lPDF[i] * sPDF[j] * Da12(z_hist[i], z_hist[j]);
        }
    }
    return sum;
}

double PFOF(double Dc_hist[], double lPDF_1[], double lPDF_2[], int PDF_size, double V_0) {
    // Calculate the probability of lens 1 and lens 2 being within V_0 along l.o.s.
    int i = 0, j = 0;
    double sum = 0;
    // Special case: Redshift bins larger than V_0. Integral reduces to a single integral.
    for (i = 0; i < PDF_size; i++) {
        sum += lPDF_1[i] * lPDF_2[i];
    }
    return sum;
}


void numsrcs(double Da_hist[], double lPDF[], int PDF_size, double phyr_width[], int phy_bin, double ls_a_sep_rad, double srcsnum[]) {
    // A function to calculate the number of sources in each r bin. For a lens-source pair, due to the uncertainty in redshift of the lens, there is fraction of source galaxies in each r bin.
    int i = 0, j = 0;
    double R_bin[PDF_size];
    int x[PDF_size], y[PDF_size];
    
    for (i = 0; i < PDF_size; i++) {
        
        R_bin[i] = Da_hist[i] * ls_a_sep_rad;
        x[i] = 0;   // An array to store which R bin to go to.
        y[i] = -1;   // An array to store which redshift bin does that belong to.
        
        // If the physical R bin lie outside the concerned radius.
        if (R_bin[i] < phyr_width[0] || R_bin[i] > phyr_width[phy_bin - 1]) {
            continue;
        }
        
        else {
            // To determine which R bin to go to.
            for (j = 1; j < phy_bin; j++) {
                if (R_bin[i] < phyr_width[j] && R_bin[i] >= phyr_width[j - 1]) {
                    //printf("%d\t%lE\t%lE\n", i, Da_hist[i], ls_a_sep_rad);
                    x[i] = j;
                    y[i] = i;
                    break;
                }
            }
        }
    }
    
    // Number of source in each r bin.
    for (i = 0; i < PDF_size; i++) {
        if (y[i] > -1) {
            srcsnum[x[i]] += lPDF[y[i]];
        }
    }
    
}


void exp_invsig2(double Da_hist[], double **Da12_hist, double lPDF[], double sPDF[], int PDF_size, double phyr_width[], int phy_bin, double ls_a_sep_rad, double sum_eff_ls[]) {
    // A function to calculate the expected value of angular diameter distance ratios <1/Sigma_cr> = <D_lD_ls / D_s>, given the PDFs of lens and source galaxies. And assign the contribution of each redshift bin to the corresponding R bin. Unit in Mpc
    int i = 0, j = 0;
    double R_bin[PDF_size];
    int x[PDF_size];
    
    for (i = 0; i < PDF_size; i++) {
        
        R_bin[i] = Da_hist[i] * ls_a_sep_rad;
        x[i] = 0;   // An array to store which R bin to go to.
        
        // If the physical R bin lie outside the concerned radius.
        if (R_bin[i] < phyr_width[0] || R_bin[i] > phyr_width[phy_bin - 1]) {
            //printf("%d\t%d\t%d\t%f\t%f\t%f\tFalse\n", i, j, x[i], R_bin[i], phyr_width[0], phyr_width[phy_bin - 1]);
            continue;
        }
        
        else {
            // To determine which R bin to go to.
            for (j = 1; j < phy_bin; j++) {
                if (R_bin[i] < phyr_width[j] && R_bin[i] >= phyr_width[j - 1]) {
                    x[i] = j;
                    //printf("%d\t%d\t%d\t%f\t%f\t%f\tTrue\n", i, j, x[i], R_bin[i], phyr_width[j], phyr_width[j - 1]);
                    break;
                }
            }
        }
    }
    
    // i-loop is for lens, j-loop is for source
    for (i = 0; i < PDF_size; i++) {
        for (j = i + 1; j < PDF_size; j++) {
            sum_eff_ls[x[i]] += lPDF[i] * sPDF[j] * Da12_hist[i][j] * Da_hist[i] / Da_hist[j];
        }
    }
    
    /*double count = 0;
    
    for (i = 0; i < phy_bin; i++) {
        count += sum_eff_ls[i];
        printf("%lE\t", sum_eff_ls[i]);
    }
    printf("\n%lE\n", count);*/
}

void exp_invsigsq2(double Da_hist[], double **Da12_hist, double lPDF[], double sPDF[], int PDF_size, double phyr_width[], int phy_bin, double ls_a_sep_rad, double sum_eff_weight[]) {
    // A function to calculate the expected value of angular diameter distance ratios <1/Sigma_cr^2> = <(D_lD_ls / D_s)^2>, given the PDFs of lens and source galaxies. And assign the contribution of each redshift bin to the corresponding R bin. Unit in Mpc
    int i = 0, j = 0;
    double R_bin[PDF_size];
    int x[PDF_size];
    
    for (i = 0; i < PDF_size; i++) {
        R_bin[i] = Da_hist[i] * ls_a_sep_rad;
        x[i] = 0;   // An array to store which R bin to go to.
        
        if (R_bin[i] < phyr_width[0] || R_bin[i] > phyr_width[phy_bin - 1]) {
            //printf("%d\t%d\t%d\t%f\t%f\t%f\tFalse\n", i, j, x[i], R_bin[i], phyr_width[0], phyr_width[phy_bin - 1]);
            continue;
        }
        
        else {
            // To determine which R bin to go to.
            for (j = 1; j < phy_bin; j++) {
                if (R_bin[i] < phyr_width[j] && R_bin[i] >= phyr_width[j - 1]) {
                    x[i] = j;
                    //printf("%d\t%d\t%d\t%f\t%f\t%f\tTrue\n", i, j, x[i], R_bin[i], phyr_width[j], phyr_width[j - 1]);
                    break;
                }
            }
        }
    }
    
    // i-loop is for lens, j-loop is for source
    for (i = 0; i < PDF_size; i++) {
        for (j = i + 1; j < PDF_size; j++) {
            sum_eff_weight[x[i]] += lPDF[i] * sPDF[j] * pow((Da12_hist[i][j] * Da_hist[i] / Da_hist[j]), 2);
        }
    }
    
    /*double count = 0;
    
    for (i = 0; i < phy_bin; i++) {
        count += sum_eff_weight[i];
        printf("%lE\t", sum_eff_weight[i]);
    }
    printf("\n%lE\n", count);*/
}

double exp_invsig(double Da_hist[], double **Da12_hist, double lPDF[], double sPDF[], int PDF_size) {
    // A function to calculate the expected value of angular diameter distance ratios <1/Sigma_cr> = <D_lD_ls / D_s>, given the PDFs of lens and source galaxies. Unit in Mpc
    int i = 0, j = 0;
    double sum = 0;
    // i-loop is for lens, j-loop is for source
    for (i = 0; i < PDF_size; i++) {
        for (j = i + 1; j < PDF_size; j++) {
            sum += lPDF[i] * sPDF[j] * Da12_hist[i][j] * Da_hist[i] / Da_hist[j];
        }
    }
    return sum;
}

double exp_invsigsq(double Da_hist[], double **Da12_hist, double lPDF[], double sPDF[], int PDF_size) {
    // A function to calculate the expected value of angular diameter distance ratios <1/Sigma_cr^2> = <(D_lD_ls / D_s)^2>, given the PDFs of lens and source galaxies. Unit in Mpc^2
    int i = 0, j = 0;
    double sum = 0;
    // i-loop is for lens, j-loop is for source
    for (i = 0; i < PDF_size; i++) {
        for (j = i + 1; j < PDF_size; j++) {
            sum += lPDF[i] * sPDF[j] * pow((Da12_hist[i][j] * Da_hist[i] / Da_hist[j]), 2);
        }
    }
    return sum;
}

double Sig_cri(double z1, double z2) {
    // Critical Sigma. Unit in M_sun pc^-2.
    return (pow(c, 2) / (4 * M_PI * G_N)) * (Da(z2) / (Da(z1) * Da12(z1, z2))) * pc_2_cm / (Msun_2_g * Mpc_2_pc);
}

double H_z_sq(double z) {
    // Hubble parameter at redshift z. H(z) squared.
    return pow(H0 / evolution(z), 2);
}

double rho_crit(double z) {
    // Critical density of the Universe, in terms of M_sun / pc^-3.
    return 3 * H_z_sq(z) / (8 * M_PI * G_N) * pc_2_cm / (Msun_2_g * pow(Mpc_2_pc, 2));
}
