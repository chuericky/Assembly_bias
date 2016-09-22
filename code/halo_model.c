;// Updated: Nov 12, 2015

#include "cosmopar.h"
#include "halo_model.h"

double rad_200(double rho_cri, double M_200) {
    // r_200 of the halo, in units of pc.
    return pow((M_200 / (rho_cri * 800 * M_PI / 3)), 1./3);
}

double Delta_vir(double x) {
    // Characteristic density of virial.
    return 18 * pow(M_PI, 2) - 82 * x - 39 * pow(x, 2);
}

double con_NFW(double M_200, double z) {
    // Concentration of the halo, reference: VU11, Eqt. 7.
    return 5.71 * pow((M_200 / (2e12 / 0.7)), -0.084) * pow((1 + z), -0.47);
}

double scale_rad(double rho_cri, double M_200, double z) {
    // Scale radius in NFW profile.
    return rad_200(rho_cri, M_200) / con_NFW(M_200, z);
}

double scale_rad_twopara(double rho_cri, double M_200, double z, double c_NFW) {
    // Scale radius in NFW profile.
    return rad_200(rho_cri, M_200) / c_NFW;
}

double char_den(double c_NFW) {
    // Characteristics density.
    return (200./3) * pow(c_NFW, 3) / (log(1 + c_NFW) - (c_NFW / (1 + c_NFW)));
}

double bias_v(double v) {
    // Bias given by Sheth et al. (2001), incorportate the adjustments in Tinker et al. (2005).
    // The peak height v should be without the squared (take square root of the one shown in VU11).
    double a = 0.707, b = 0.35, c0 = 0.80;
    double pr = a * pow(v, 2.0);
    
    return 1 + (pow(a, 0.5) * pr + pow(a, 0.5) * b * pow(pr, (1-c0)) - pow(pr, c0)/(pow(pr, c0) + b * (1-c0) *(1-0.5*c0))) / (pow(a, 0.5) * delta_crit);
}

double f_v(double v) {
    // f(v) given by Eqt. (12) of VU11, normalization constant in that paper is wrong!
    double a = 0.707, p = 0.3, A = 0.108074459;
    double w = pow(v, 2.);
    
    return A * (1 + pow((a * w),-p)) * pow(w, -0.5) * exp(-a * w / 2);
}

double mass_lum_size_gmass_corr (char* bin_name, double m_opt) {
    // Correcting the fit halo mass because of photo-z correction (Velander et. al. 2014 Table B2). Hard-coded!
    double lum_fac = 0, mass_fac = 0;
    char color[3], lum[2], gmass[2];
    
    strncpy (color, bin_name, 3);
    color[3] = 0;
    strncpy (lum, bin_name+4, 2);
    lum[2] = 0;
    strncpy (gmass, bin_name+10, 2);
    gmass[2] = 0;
    
    printf("%s\t%s\t%s\n", color, lum, gmass);
    
    if ((strcmp(color, "blu") == 0) & (strcmp(lum, "L1") == 0)) lum_fac = 0.94;
    if ((strcmp(color, "blu") == 0) & (strcmp(lum, "L2") == 0)) lum_fac = 0.96;
    if ((strcmp(color, "blu") == 0) & (strcmp(lum, "L3") == 0)) lum_fac = 1.01;
    if ((strcmp(color, "blu") == 0) & (strcmp(lum, "L4") == 0)) lum_fac = 1.05;
    if ((strcmp(color, "blu") == 0) & (strcmp(lum, "L5") == 0)) lum_fac = 1.04;
    if ((strcmp(color, "blu") == 0) & (strcmp(lum, "L6") == 0)) lum_fac = 1.20;
    if ((strcmp(color, "blu") == 0) & (strcmp(lum, "L7") == 0)) lum_fac = 1;
    if ((strcmp(color, "blu") == 0) & (strcmp(lum, "L8") == 0)) lum_fac = 1;
    if ((strcmp(color, "red") == 0) & (strcmp(lum, "L1") == 0)) lum_fac = 0.86;
    if ((strcmp(color, "red") == 0) & (strcmp(lum, "L2") == 0)) lum_fac = 0.89;
    if ((strcmp(color, "red") == 0) & (strcmp(lum, "L3") == 0)) lum_fac = 0.96;
    if ((strcmp(color, "red") == 0) & (strcmp(lum, "L4") == 0)) lum_fac = 1.02;
    if ((strcmp(color, "red") == 0) & (strcmp(lum, "L5") == 0)) lum_fac = 1.09;
    if ((strcmp(color, "red") == 0) & (strcmp(lum, "L6") == 0)) lum_fac = 1.13;
    if ((strcmp(color, "red") == 0) & (strcmp(lum, "L7") == 0)) lum_fac = 1.16;
    if ((strcmp(color, "red") == 0) & (strcmp(lum, "L8") == 0)) lum_fac = 1.36;
    
    if ((strcmp(color, "blu") == 0) & (strcmp(gmass, "M1") == 0)) mass_fac = 1.18;
    if ((strcmp(color, "blu") == 0) & (strcmp(gmass, "M2") == 0)) mass_fac = 1.28;
    if ((strcmp(color, "blu") == 0) & (strcmp(gmass, "M3") == 0)) mass_fac = 1.50;
    if ((strcmp(color, "blu") == 0) & (strcmp(gmass, "M4") == 0)) mass_fac = 1.83;
    if ((strcmp(color, "blu") == 0) & (strcmp(gmass, "M5") == 0)) mass_fac = 1;
    if ((strcmp(color, "blu") == 0) & (strcmp(gmass, "M6") == 0)) mass_fac = 1;
    if ((strcmp(color, "blu") == 0) & (strcmp(gmass, "M7") == 0)) mass_fac = 1;
    if ((strcmp(color, "blu") == 0) & (strcmp(gmass, "M8") == 0)) mass_fac = 1;
    if ((strcmp(color, "red") == 0) & (strcmp(gmass, "M1") == 0)) mass_fac = 0.59;
    if ((strcmp(color, "red") == 0) & (strcmp(gmass, "M2") == 0)) mass_fac = 0.74;
    if ((strcmp(color, "red") == 0) & (strcmp(gmass, "M3") == 0)) mass_fac = 0.91;
    if ((strcmp(color, "red") == 0) & (strcmp(gmass, "M4") == 0)) mass_fac = 1.19;
    if ((strcmp(color, "red") == 0) & (strcmp(gmass, "M5") == 0)) mass_fac = 1.53;
    if ((strcmp(color, "red") == 0) & (strcmp(gmass, "M6") == 0)) mass_fac = 1.86;
    if ((strcmp(color, "red") == 0) & (strcmp(gmass, "M7") == 0)) mass_fac = 2.26;
    if ((strcmp(color, "red") == 0) & (strcmp(gmass, "M8") == 0)) mass_fac = 2.73;
    /*
    // For Velander bins.
    if (strcmp(bin_name, "blu_L1.dat") == 0) corr_fac = 0.94;
    if (strcmp(bin_name, "blu_L2.dat") == 0) corr_fac = 0.96;
    if (strcmp(bin_name, "blu_L3.dat") == 0) corr_fac = 1.01;
    if (strcmp(bin_name, "blu_L4.dat") == 0) corr_fac = 1.05;
    if (strcmp(bin_name, "blu_L5.dat") == 0) corr_fac = 1.04;
    if (strcmp(bin_name, "blu_L6.dat") == 0) corr_fac = 1.20;
    if (strcmp(bin_name, "blu_L7.dat") == 0) corr_fac = 1.20;
    if (strcmp(bin_name, "red_L1.dat") == 0) corr_fac = 0.86;
    if (strcmp(bin_name, "red_L2.dat") == 0) corr_fac = 0.89;
    if (strcmp(bin_name, "red_L3.dat") == 0) corr_fac = 0.96;
    if (strcmp(bin_name, "red_L4.dat") == 0) corr_fac = 1.02;
    if (strcmp(bin_name, "red_L5.dat") == 0) corr_fac = 1.09;
    if (strcmp(bin_name, "red_L6.dat") == 0) corr_fac = 1.13;
    if (strcmp(bin_name, "red_L7.dat") == 0) corr_fac = 1.16;
    if (strcmp(bin_name, "red_L8.dat") == 0) corr_fac = 1.36;*/
    
    return m_opt * lum_fac * mass_fac;
}

double mass_lum_corr (char* bin_name, double m_opt) {
    // Correcting the fit halo mass because of photo-z correction (Velander et. al. 2014 Table B2). Hard-coded!
    double corr_fac = 0;
    
    // For Velander bins.
    if (strcmp(bin_name, "blu_L1.dat") == 0) corr_fac = 0.94;
    if (strcmp(bin_name, "blu_L2.dat") == 0) corr_fac = 0.96;
    if (strcmp(bin_name, "blu_L3.dat") == 0) corr_fac = 1.01;
    if (strcmp(bin_name, "blu_L4.dat") == 0) corr_fac = 1.05;
    if (strcmp(bin_name, "blu_L5.dat") == 0) corr_fac = 1.04;
    if (strcmp(bin_name, "blu_L6.dat") == 0) corr_fac = 1.20;
    if (strcmp(bin_name, "blu_L7.dat") == 0) corr_fac = 1.20;
    if (strcmp(bin_name, "red_L1.dat") == 0) corr_fac = 0.86;
    if (strcmp(bin_name, "red_L2.dat") == 0) corr_fac = 0.89;
    if (strcmp(bin_name, "red_L3.dat") == 0) corr_fac = 0.96;
    if (strcmp(bin_name, "red_L4.dat") == 0) corr_fac = 1.02;
    if (strcmp(bin_name, "red_L5.dat") == 0) corr_fac = 1.09;
    if (strcmp(bin_name, "red_L6.dat") == 0) corr_fac = 1.13;
    if (strcmp(bin_name, "red_L7.dat") == 0) corr_fac = 1.16;
    if (strcmp(bin_name, "red_L8.dat") == 0) corr_fac = 1.36;
    
    printf("Luminosity correction factor\t=\t%f\n", corr_fac);
    return m_opt * corr_fac;
}

double mass_lum_corr2 (char* bin_name, double m_opt) {
    // Correcting the fit halo mass because of photo-z correction (Velander et. al. 2014 Table B2). Hard-coded!
    double corr_fac = 0;
    
    // For Velander bins.
    if (strcmp(bin_name, "blu_L1.dat") == 0) corr_fac = 1.05;
    if (strcmp(bin_name, "blu_L2.dat") == 0) corr_fac = 1.05;
    if (strcmp(bin_name, "blu_L3.dat") == 0) corr_fac = 1.20;
    if (strcmp(bin_name, "red_L1.dat") == 0) corr_fac = 0.8716;
    if (strcmp(bin_name, "red_L2.dat") == 0) corr_fac = 1.0996;
    if (strcmp(bin_name, "red_L3.dat") == 0) corr_fac = 1.16;
    
    printf("%f\n", corr_fac);
    return m_opt * corr_fac;
}

void dsig_baryon(double star_mass, double r[], double sig_bar[], int N_size) {
    // Calculates the differential Sigma of baryonic component.
    int i = 0;
    for (i = 0; i < N_size; i++) {
        sig_bar[i] = star_mass / (M_PI * pow(r[i], 2));
    }
}

void dsig_cen1h(double M_200, double z, double r_phy[], double sig_cen_1h[], int N_size) {
    // Calculates the differential Sigma of the Central 1-halo term.
    int i = 0;
    double a = 0, u = 0, rho_critical = 0, char_density = 0, scale_radius = 0, c_NFW = 0;
    
    rho_critical = rho_crit(z);
    c_NFW = con_NFW(M_200, z);
    char_density = char_den(c_NFW);
    scale_radius = scale_rad(rho_critical, M_200, z);
    
    for (i = 0; i < N_size; i++) {
        u = r_phy[i] / scale_radius;
        if (u < 1)
            a = (2/pow(u, 2) * (log(u/2) + acosh(1/u)/pow((1 - u*u), 0.5)) - (1 - acosh(1/u)/pow((1 - u*u), 0.5))/(u*u - 1));
        else if (u == 1)
            a = (2/pow(u, 2) * (1 + log(0.5)) - (1./3));
        else if (u > 1)
            a = (2/pow(u, 2) * (log(u/2) + acos(1/u)/pow((u*u - 1), 0.5)) - (1 - acos(1/u)/pow((u*u - 1), 0.5))/(u*u - 1));
        
        sig_cen_1h[i] = 2 * rho_critical * char_density * scale_radius * a;
    }
}

void dsig_cen1h_twopara(double M_200, double z, double c_NFW, double r_phy[], double sig_cen_1h[], int N_size) {
    // Calculates the differential Sigma of the Central 1-halo term. Concentration is a free parameter to input.
    int i = 0;
    double a = 0, u = 0, rho_critical = 0, char_density = 0, scale_radius = 0;
    
    rho_critical = rho_crit(z);
    char_density = char_den(c_NFW);
    scale_radius = scale_rad_twopara(rho_critical, M_200, z, c_NFW);
    
    for (i = 0; i < N_size; i++) {
        u = r_phy[i] / scale_radius;
        if (u < 1)
            a = (2/pow(u, 2) * (log(u/2) + acosh(1/u)/pow((1 - u*u), 0.5)) - (1 - acosh(1/u)/pow((1 - u*u), 0.5))/(u*u - 1));
        else if (u == 1)
            a = (2/pow(u, 2) * (1 + log(0.5)) - (1./3));
        else if (u > 1)
            a = (2/pow(u, 2) * (log(u/2) + acos(1/u)/pow((u*u - 1), 0.5)) - (1 - acos(1/u)/pow((u*u - 1), 0.5))/(u*u - 1));
        
        sig_cen_1h[i] = 2 * rho_critical * char_density * scale_radius * a;
    }
}

void dsig_sat1h(double rho_cri, double char_den, double scale_rad, double x[], double sig_sat_1h[], double M_200, double p_r[], int N_size) {
    // Calculates the differential Sigma of the Central 1-halo term.
    double r_200 = rad_200(rho_cri, M_200);
    double a = 0, b = 0, u = 0, x_rat = 0, r_cut = 0;
    int i = 0;
    
    x_rat = 0.4 * r_200 / scale_rad;
    r_cut = 0.4 * r_200;
    
    if (x_rat >= 1) {
        for (i = 0; i < N_size; i++) {
            u = x[i];
            
            if (u <= x_rat) {
                if (u < 1)
                    a = (2/pow(u, 2) * (log(u/2) + acosh(1/u)/pow((1 - u*u), 0.5)) - (1 - acosh(1/u)/pow((1 - u*u), 0.5))/(u*u - 1));
                else if (u == 1)
                    a = (2/pow(u, 2) * (1 + log(0.5)) - (1./3));
                else if (u > 1 && u <= x_rat)
                    a = (2/pow(u, 2) * (log(u/2) + acos(1/u)/pow((u*u - 1), 0.5)) - (1 - acos(1/u)/pow((u*u - 1), 0.5))/(u*u - 1));
                
                sig_sat_1h[i] = 2 * rho_cri * char_den * scale_rad * a;
            }
            
            else if (u > x_rat) {
                b = 2 * rho_cri * char_den * scale_rad * (2/pow(x_rat, 2) * (log(x_rat/2) + acos(1/x_rat)/pow((x_rat*x_rat - 1), 0.5)) - (1 - acos(1/x_rat)/pow((x_rat*x_rat - 1), 0.5))/(x_rat*x_rat - 1));
                
                sig_sat_1h[i] = b * pow((r_cut/p_r[i]), 2);
            }
        }
    }
    
    if (x_rat < 1) {
        for (i = 0; i < N_size; i++) {
            u = x[i];
            
            if (u <= x_rat) {
                a = (2/pow(u, 2) * (log(u/2) + acosh(1/u)/pow((1 - u*u), 0.5)) - (1 - acosh(1/u)/pow((1 - u*u), 0.5))/(u*u - 1));
                
                sig_sat_1h[i] = 2 * rho_cri * char_den * scale_rad * a;
            }
            
            else if (u > x_rat) {
                b = 2 * rho_cri * char_den * scale_rad * (2/pow(x_rat, 2) * (log(x_rat/2) + acos(1/x_rat)/pow((x_rat*x_rat - 1), 0.5)) - (1 - acos(1/x_rat)/pow((x_rat*x_rat - 1), 0.5))/(x_rat*x_rat - 1));
                
                sig_sat_1h[i] = b * pow((r_cut/p_r[i]), 2);
            }
        }
    }
}