# Makefile for radiative transfer models
# Author: Alexey Shiklomanov, Dept of Earth and Environment / Boston University
# Email: ashiklom@bu.edu

# Compilers
FC = gfortran -ffree-form
FC_fixed = gfortran -ffixed-form

# Optimization
OPT = -O2

# Folder paths
bin = bin/
f_prospect = prospect/
f_sail = sail/
f_ps = prosail/
f_ds = dataSpec/

# Module shorthand
mod_f = $(f_ds)prosail.dataSpec_
mod_m = $(bin)mod_dataspec_

# RTM shared libraries
PROSPECT = $(bin)prospect.so
SAIL = $(bin)sail.so
PROSAIL = $(bin)prosail.so

# RTM dependencies
PROSPECT_depends = $(f_prospect)*.f90 $(mod_m)prospect5b.o $(mod_m)refractive.o
SAIL_depends = $(PROSPECT) $(f_sail)*.f90 $(mod_m)soil.o
PROSAIL_depends = $(PROSPECT) $(SAIL) $(f_ps)prosail.main.f90

$(PROSAIL) : $(RTM_depends)
	$(FC) -g -shared $(f_ps)prosail.main.f90 $(SAIL) $(PROSPECT) -o $(PROSAIL) -fPIC -I$(bin)

$(PROSPECT) : $(PROSPECT_depends)
	$(FC) -g -shared $(PROSPECT_depends) -o $(PROSPECT) -I$(bin) -fPIC

$(SAIL) : $(SAIL_depends)
	$(FC) -g -shared $(SAIL_depends) -o $(SAIL) -I$(bin) -fPIC

$(mod_m)prospect5b.o $(mod_m)refractive.o $(mod_m)soil.o : $(mod_f)wavelength.f90
	$(FC) -c $(mod_f)wavelength.f90 -J$(bin) -o $(mod_m)wavelength.o
	$(FC) -c $(mod_f)prospect5b.f90 -I$(bin) -J$(bin) -o $(mod_m)prospect5b.o 
	$(FC) -c $(mod_f)refractive.f90 -I$(bin) -J$(bin) -o $(mod_m)refractive.o 
	$(FC) -c $(mod_f)soil.f90 -I$(bin) -J$(bin) -o $(mod_m)soil.o 

clean :
	rm bin/*
