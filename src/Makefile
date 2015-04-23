# Makevars file for PEcAnRTM package
# Author: Alexey Shiklomanov
# Email: ashiklom@bu.edu

# Folder paths
f_pp = RTM/prospect/
f_sl = RTM/sail/
f_ps = RTM/prosail/
f_ds = RTM/dataSpec/

# Module shorthand
mod_f = $(f_ds)prosail.dataSpec_
mod_m = mod_dataspec_

# RTM shared libraries
PROSPECT = prospect.so
SAIL = sail.so
PROSAIL = prosail.so

# RTM dependencies
PROSPECT_depends = $(f_pp)*.f90 $(mod_m)prospect5b.o $(mod_m)refractive.o
SAIL_depends = $(PROSPECT) $(f_ss)*.f90 $(mod_m)soil.o
PROSAIL_depends = $(PROSPECT) $(SAIL) $(f_ps)prosail.main.f90

$(PROSAIL) : $(PROSAIL_depends)
	$(FC) -g -shared $(f_ps)prosail.main.f90 $(SAIL) $(PROSPECT) -o $(PROSAIL) -fPIC 

$(PROSPECT) : $(PROSPECT_depends)
	$(FC) -g -shared $(PROSPECT_depends) -o $(PROSPECT) -fPIC

$(SAIL) : $(SAIL_depends)
	$(FC) -g -shared $(SAIL_depends) -o $(SAIL) -fPIC

$(mod_m)prospect5b.o $(mod_m)refractive.o $(mod_m)soil.o : $(mod_f)wavelength.f90
	$(FC) -c $(mod_f)wavelength.f90 -o $(mod_m)wavelength.o
	$(FC) -c $(mod_f)prospect5b.f90 -o $(mod_m)prospect5b.o 
	$(FC) -c $(mod_f)refractive.f90 -o $(mod_m)refractive.o 
	$(FC) -c $(mod_f)soil.f90 -o $(mod_m)soil.o 

PKG_LIBS = -LRTM $(FLIBS)
$(SHLIB) : RTM/bin/pecan.rtm.so

clean :
	rm *.so *.o *.mod

