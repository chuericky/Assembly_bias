CC = gcc
CFLAGS = -O3
CFLAGS2 = -g -O2  -Dg77Fortran -fPIC -fno-common
HOME_DIR = /Users/rickyccy/Documents/Research_assemblybias/code/
CFIT_DIR = /Users/rickyccy/Documents/Research_assemblybias/code/cfitsio/
LIBS = -lm 
CFITLIB = -lcfitsio
CFITPATH = -L./cfitsio/
CFLAGS1=-I/usr/local/include/gsl/
LDFLAGS=-lgsl -lgslcblas

all: convert category category2 category3 category_multi category_fore category_check category_check2 redred gg_lensing gg_lensing_parallel gg_lensing_bootstrap_parallel gg_lensing_bootstrap_parallel2 halo_fitting cluster_select agressive_sat_cut_extension lens_select lens_select_multi lens_select_extension gg_lensing_bootstrap_central gg_lensing_bootstrap_multifx gg_lensing_multifx_individualcor gg_lensing_multi gg_lensing_wider gg_lensing_assembias_individualcor gg_lensing_bootstrap_assembias gg_lensing_bootstrap_central_widerbin gg_lensing_bootstrap_central_adjusted halo_chisq_fit halo_chisq_fit_2 halo_chisq_fit_3 chisq_doublefit chisq_doublefit_multifx chisq_doublefit_assembias chisq_doublefit_assembias2 PDF_norm_exp double_Riemann gg_lensing_assembias_pdf_inexp gg_lensing_assembias_pdf_allexp gg_lensing_multifx_pdf_allexp chi2_assembly FOF halo_generate 2pt_corfun_DD 2pt_corfun_DR 2pt_corfun_RR 2pt_corfun_DD_try 2pt_corfun_DR_try 2pt_corfun_RR_try try_cfit try_read try try_create ran_generate jackknife_area jack_area_ext jackknife_area_new ran_jack_area 2pt_corfun_DD_try chi2_multi chi2_multi_wider mock_chi2 mock_chi2_old projected_FOF projected_FOF_new projected_PFOF mock_2pt_DD mock_2pt_DR mock_2pt_RR mock_2pt_DD_sim mock_2pt_DR_sim mock_2pt_RR_sim mock_2pt_DD_test mock_2pt_DR_test mock_2pt_cross_DD mock_2pt_cross_DD_noab mock_2pt_DD_noab mock_2pt_DR_noab mock_2pt_RR_test mock_lensing particle_project 2pt_auto_DD 2pt_auto_DR 2pt_auto_DD_wider 2pt_auto_DR_wider 2pt_auto_DD_ext 2pt_auto_DR_ext 2pt_auto_RR_ext 2pt_auto_DD_wider_2 2pt_auto_RR_wider_2 2pt_auto_RR_wider 2pt_auto_DR_wider_2 2pt_auto_DD_all 2pt_auto_DR_all 2pt_auto_DD_all_wider 2pt_auto_DR_all_wider 2pt_auto_RR_all_wider 2pt_cross_DD_wider 2pt_cross_DR_wider 2pt_cross_RD_wider 2pt_cross_RR_wider density_grid particle_to_halo mock_brightgal mock_brightgal_new read_halo_Sig PFOF mean_sep mock_diffrad mock_diffrad_periodic further_sat_remove mock_clus_halo mock_nbr_halo chi2fit sel_isogal get_iso count_area category_extension overdensity overdensity_uniden.c mock_overdensity select_good_redshift 

fit_nfw: fit_nfw.c Makefile
	$(CC) $(CFLAGS1) $(LDFLAGS) $< -o $@ $(LIBS)

fit_nfw_err: fit_nfw_err.c Makefile
	$(CC) $(CFLAGS1) $(LDFLAGS) $< -o $@ $(LIBS)

chi2fit: chi2fit.c tool.c tool.h Makefile
	$(CC) $(CFLAGS1) $(LDFLAGS) tool.c $< -o $@ $(LIBS)

particle_project: particle_project.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o particle_project particle_project.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_clus_halo: mock_clus_halo.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c halo_model.c halo_model.h Makefile
	$(CC) $(CFLAGS) -o mock_clus_halo mock_clus_halo.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c halo_model.c halo_model.h $(LIBS)

mock_nbr_halo: mock_nbr_halo.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c halo_model.c halo_model.h Makefile
	$(CC) $(CFLAGS) -o mock_nbr_halo mock_nbr_halo.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c halo_model.c halo_model.h $(LIBS)

mean_sep: mean_sep.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mean_sep mean_sep.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

read_halo_Sig: read_halo_Sig.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o read_halo_Sig read_halo_Sig.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

density_grid: density_grid.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o density_grid density_grid.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

particle_to_halo: particle_to_halo.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o particle_to_halo particle_to_halo.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

ran_jack_area: ran_jack_area.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o ran_jack_area ran_jack_area.c tool.c $(LIBS)

jackknife_area: jackknife_area.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o jackknife_area jackknife_area.c tool.c $(LIBS)

jack_area_ext: jack_area_ext.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o jack_area_ext jack_area_ext.c tool.c $(LIBS)

jackknife_area_new: jackknife_area_new.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o jackknife_area_new jackknife_area_new.c tool.c $(LIBS)

jackknife_area_cross: jackknife_area_cross.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o jackknife_area_cross jackknife_area_cross.c tool.c $(LIBS)

ran_generate: ran_generate.c $(CFIT_DIR)libcfitsio.a tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS2) -o ran_generate ran_generate.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(CFITPATH) $(CFITLIB) $(LIBS)

try: try.c $(CFIT_DIR)libcfitsio.a Makefile
	$(CC) $(CFLAGS2) -o try try.c $(CFITPATH) $(CFITLIB) $(LIBS)

try_create: try_create.c $(CFIT_DIR)libcfitsio.a Makefile
	$(CC) $(CFLAGS2) -o try_create try_create.c $(CFITPATH) $(CFITLIB) $(LIBS)

try_cfit: try_cfit.c $(CFIT_DIR)libcfitsio.a Makefile
	$(CC) $(CFLAGS2) -o try_cfit try_cfit.c $(CFITPATH) $(CFITLIB) $(LIBS)

try_read: try_read.c $(CFIT_DIR)libcfitsio.a Makefile
	$(CC) $(CFLAGS2) -o try_read try_read.c $(CFITPATH) $(CFITLIB) $(LIBS)

count_area: count_area.c $(CFIT_DIR)libcfitsio.a tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS2) -o count_area count_area.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(CFITPATH) $(CFITLIB) $(LIBS)

2pt_auto_DD: 2pt_auto_DD.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DD 2pt_auto_DD.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DD_all_wider: 2pt_auto_DD_all_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DD_all_wider 2pt_auto_DD_all_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DD_ext: 2pt_auto_DD_ext.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DD_ext 2pt_auto_DD_ext.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DR_ext: 2pt_auto_DR_ext.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DR_ext 2pt_auto_DR_ext.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_RR_ext: 2pt_auto_RR_ext.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_RR_ext 2pt_auto_RR_ext.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DR_all_wider: 2pt_auto_DR_all_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DR_all_wider 2pt_auto_DR_all_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_RR_all_wider: 2pt_auto_RR_all_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_RR_all_wider 2pt_auto_RR_all_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DD_wider_2: 2pt_auto_DD_wider_2.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DD_wider_2 2pt_auto_DD_wider_2.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DR_wider_2: 2pt_auto_DR_wider_2.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DR_wider_2 2pt_auto_DR_wider_2.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_RR_wider_2: 2pt_auto_RR_wider_2.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_RR_wider_2 2pt_auto_RR_wider_2.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_RR_wider: 2pt_auto_RR_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_RR_wider 2pt_auto_RR_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

overdensity: overdensity.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o overdensity overdensity.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

overdensity_uniden: overdensity_uniden.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o overdensity_uniden overdensity_uniden.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_overdensity: mock_overdensity.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_overdensity mock_overdensity.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_cross_DD_wider: 2pt_cross_DD_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_cross_DD_wider 2pt_cross_DD_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_cross_DR_wider: 2pt_cross_DR_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_cross_DR_wider 2pt_cross_DR_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_cross_RD_wider: 2pt_cross_RD_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_cross_RD_wider 2pt_cross_RD_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_cross_RR_wider: 2pt_cross_RR_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_cross_RR_wider 2pt_cross_RR_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DD_wider: 2pt_auto_DD_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DD_wider 2pt_auto_DD_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DD_all: 2pt_auto_DD_all.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DD_all 2pt_auto_DD_all.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DR: 2pt_auto_DR.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DR 2pt_auto_DR.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DR_wider: 2pt_auto_DR_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DR_wider 2pt_auto_DR_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_auto_DR_all: 2pt_auto_DR_all.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_auto_DR_all 2pt_auto_DR_all.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_corfun_DD: 2pt_corfun_DD.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_corfun_DD 2pt_corfun_DD.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_corfun_DR: 2pt_corfun_DR.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_corfun_DR 2pt_corfun_DR.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_corfun_RR: 2pt_corfun_RR.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_corfun_RR 2pt_corfun_RR.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_lensing: mock_lensing.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_lensing mock_lensing.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_cross_DD: mock_2pt_cross_DD.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_cross_DD mock_2pt_cross_DD.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_cross_DD_noab: mock_2pt_cross_DD_noab.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_cross_DD_noab mock_2pt_cross_DD_noab.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_DD_test: mock_2pt_DD_test.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_DD_test mock_2pt_DD_test.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_DR_test: mock_2pt_DR_test.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_DR_test mock_2pt_DR_test.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_DD_noab: mock_2pt_DD_noab.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_DD_noab mock_2pt_DD_noab.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_DR_noab: mock_2pt_DR_noab.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_DR_noab mock_2pt_DR_noab.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_RR_test: mock_2pt_RR_test.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_RR_test mock_2pt_RR_test.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_DD_sim: mock_2pt_DD_sim.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_DD_sim mock_2pt_DD_sim.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_DR_sim: mock_2pt_DR_sim.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_DR_sim mock_2pt_DR_sim.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_2pt_RR_sim: mock_2pt_RR_sim.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o mock_2pt_RR_sim mock_2pt_RR_sim.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_corfun_DD_try: 2pt_corfun_DD_try.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_corfun_DD_try 2pt_corfun_DD_try.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_corfun_DR_try: 2pt_corfun_DR_try.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_corfun_DR_try 2pt_corfun_DR_try.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_corfun_RR_try: 2pt_corfun_RR_try.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_corfun_RR_try 2pt_corfun_RR_try.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_corfun_DD_try_jack: 2pt_corfun_DD_try_jack.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_corfun_DD_try_jack 2pt_corfun_DD_try_jack.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_corfun_DR_try_jack: 2pt_corfun_DR_try_jack.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_corfun_DR_try_jack 2pt_corfun_DR_try_jack.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_corfun_RR_try_jack: 2pt_corfun_RR_try_jack.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_corfun_RR_try_jack 2pt_corfun_RR_try_jack.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_cross_DD: 2pt_cross_DD.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_cross_DD 2pt_cross_DD.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_cross_DR: 2pt_cross_DR.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_cross_DR 2pt_cross_DR.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_cross_RD: 2pt_cross_RD.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_cross_RD 2pt_cross_RD.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

2pt_cross_RR: 2pt_cross_RR.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o 2pt_cross_RR 2pt_cross_RR.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

category: category.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o category category.c tool.c $(LIBS)

convert: convert.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o convert convert.c tool.c $(LIBS)

redred: redred.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o redred redred.c tool.c $(LIBS)

category2: category2.c tool.c tool.h Makefile
	$(CC) $(CFLAGS) -o category2 category2.c tool.c $(LIBS)

category3: category3.c tool.c tool.h tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o category3 category3.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

category_multi: category_multi.c tool.c tool.h tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o category_multi category_multi.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

category_extension: category_extension.c tool.c tool.h tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o category_extension category_extension.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

select_good_redshift: select_good_redshift.c tool.c tool.h tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o select_good_redshift select_good_redshift.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

get_iso: get_iso.c tool.c tool.h tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o get_iso get_iso.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

sel_isogal: sel_isogal.c tool.c tool.h tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o sel_isogal sel_isogal.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

category_fore: category_fore.c tool.c tool.h tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o category_fore category_fore.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

category_check: category_check.c tool.c tool.h tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o category_check category_check.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

category_check2: category_check2.c tool.c tool.h tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o category_check2 category_check2.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing: gg_lensing.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing gg_lensing.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c $(LIBS)

gg_lensing_parallel: gg_lensing_parallel.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_parallel gg_lensing_parallel.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c $(LIBS)

gg_lensing_bootstrap_parallel: gg_lensing_bootstrap_parallel.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_bootstrap_parallel gg_lensing_bootstrap_parallel.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_bootstrap_central: gg_lensing_bootstrap_central.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_bootstrap_central gg_lensing_bootstrap_central.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_bootstrap_multifx: gg_lensing_bootstrap_multifx.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_bootstrap_multifx gg_lensing_bootstrap_multifx.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_bootstrap_assembias: gg_lensing_bootstrap_assembias.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_bootstrap_assembias gg_lensing_bootstrap_assembias.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_assembias_individualcor: gg_lensing_assembias_individualcor.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_assembias_individualcor gg_lensing_assembias_individualcor.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_multifx_pdf_allexp: gg_lensing_multifx_pdf_allexp.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_multifx_pdf_allexp gg_lensing_multifx_pdf_allexp.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_multi: gg_lensing_multi.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_multi gg_lensing_multi.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_wider: gg_lensing_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_wider gg_lensing_wider.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

FOF: FOF.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o FOF FOF.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

projected_FOF: projected_FOF.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o projected_FOF projected_FOF.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

projected_FOF_new: projected_FOF_new.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o projected_FOF_new projected_FOF_new.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

further_sat_remove: further_sat_remove.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o further_sat_remove further_sat_remove.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

agressive_sat_cut_extension: agressive_sat_cut_extension.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o agressive_sat_cut_extension agressive_sat_cut_extension.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_brightgal: mock_brightgal.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o mock_brightgal mock_brightgal.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

mock_brightgal_new: mock_brightgal_new.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o mock_brightgal_new mock_brightgal_new.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

projected_PFOF: projected_PFOF.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o projected_PFOF projected_PFOF.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

PFOF: PFOF.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o PFOF PFOF.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_assembias_pdf_inexp: gg_lensing_assembias_pdf_inexp.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_assembias_pdf_inexp gg_lensing_assembias_pdf_inexp.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_assembias_pdf_allexp: gg_lensing_assembias_pdf_allexp.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_assembias_pdf_allexp gg_lensing_assembias_pdf_allexp.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

PDF_norm_exp: PDF_norm_exp.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o PDF_norm_exp PDF_norm_exp.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_multifx_individualcor: gg_lensing_multifx_individualcor.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_multifx_individualcor gg_lensing_multifx_individualcor.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_bootstrap_central_widerbin: gg_lensing_bootstrap_central_widerbin.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_bootstrap_central_widerbin gg_lensing_bootstrap_central_widerbin.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_bootstrap_central_adjusted: gg_lensing_bootstrap_central_adjusted.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c shell.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_bootstrap_central_adjusted gg_lensing_bootstrap_central_adjusted.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c shell.c $(LIBS)

gg_lensing_bootstrap_parallel2: gg_lensing_bootstrap_parallel2.c erff.c gcf.c gammp.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c Makefile
	$(CC) $(CFLAGS) -o gg_lensing_bootstrap_parallel2 gg_lensing_bootstrap_parallel2.c tool.c cosmopar.c qromb.c nrutil.c polint.c trapzd.c erff.c gcf.c gammp.c $(LIBS)

halo_fitting: halo_fitting.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o halo_fitting halo_fitting.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

halo_chisq_fit: halo_chisq_fit.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o halo_chisq_fit halo_chisq_fit.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

halo_chisq_fit_2: halo_chisq_fit_2.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o halo_chisq_fit_2 halo_chisq_fit_2.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

halo_chisq_fit_3: halo_chisq_fit_3.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o halo_chisq_fit_3 halo_chisq_fit_3.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

chisq_doublefit: chisq_doublefit.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o chisq_doublefit chisq_doublefit.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

chisq_doublefit_multifx: chisq_doublefit_multifx.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o chisq_doublefit_multifx chisq_doublefit_multifx.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

chisq_doublefit_assembias: chisq_doublefit_assembias.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o chisq_doublefit_assembias chisq_doublefit_assembias.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

chisq_doublefit_assembias2: chisq_doublefit_assembias2.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o chisq_doublefit_assembias2 chisq_doublefit_assembias2.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

halo_generate: halo_generate.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o halo_generate halo_generate.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

chi2_assembly: chi2_assembly.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o chi2_assembly chi2_assembly.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

chi2_multi: chi2_multi.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o chi2_multi chi2_multi.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

chi2_multi_wider: chi2_multi_wider.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o chi2_multi_wider chi2_multi_wider.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

mock_chi2: mock_chi2.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o mock_chi2 mock_chi2.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

mock_chi2_old: mock_chi2_old.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o mock_chi2_old mock_chi2_old.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

double_Riemann: double_Riemann.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c gsl_math.h gsl_integration.h Makefile
	$(CC) $(CFLAGS) -o double_Riemann double_Riemann.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS) -fopenmp

cluster_select: cluster_select.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o cluster_select cluster_select.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

lens_select: lens_select.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o lens_select lens_select.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

lens_select_multi: lens_select_multi.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o lens_select_multi lens_select_multi.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

lens_select_extension: lens_select_extension.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c Makefile
	$(CC) $(CFLAGS) -o lens_select_extension lens_select_extension.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

mock_diffrad_periodic: mock_diffrad_periodic.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c halo_model.h Makefile
	$(CC) $(CFLAGS) -o mock_diffrad_periodic mock_diffrad_periodic.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)

mock_diffrad: mock_diffrad.c tool.c tool.h cosmopar.c cosmopar.h qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c halo_model.h Makefile
	$(CC) $(CFLAGS) -o mock_diffrad mock_diffrad.c tool.c cosmopar.c qromb.c nrutil.c nrutil.h polint.c trapzd.c halo_model.c $(LIBS)
