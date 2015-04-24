# Makevars file for PEcAnRTM package
# Author: Alexey Shiklomanov
# Email: ashiklom@bu.edu

# Module shorthand
mod_f = prosail.dataSpec_
mod_m = mod_dataspec_

$(mod_m)prospect5b.o $(mod_m)refractive.o $(mod_m)soil.o : $(mod_f)wavelength.f90
	gfortran -c $(mod_f)wavelength.f90 -o $(mod_m)wavelength.o
	gfortran -c $(mod_f)prospect5b.f90 -o $(mod_m)prospect5b.o 
	gfortran -c $(mod_f)refractive.f90 -o $(mod_m)refractive.o 
	gfortran -c $(mod_f)soil.f90 -o $(mod_m)soil.o 
