# Makefile for compiling Fortran optics codes

FC = gfortran
FFLAGS = -O2 -Wall -std=f2008
LDFLAGS = 

.PHONY: all clean

all: main_te_integrate

# Compile modules only (no main programs)
# We use -c to compile to object, which stops at the module boundary
multilayer_reflectance.o: multilayer_reflectance.f90
	$(FC) $(FFLAGS) -c -o $@ $<

te_greens.o: te_greens.f90 multilayer_reflectance.o
	$(FC) $(FFLAGS) -c -o $@ $<

tm_greens.o: tm_greens.f90 te_greens.o
	$(FC) $(FFLAGS) -c -o $@ $<

bessel_functions.o: bessel_functions.f90 multilayer_reflectance.o
	$(FC) $(FFLAGS) -c -o $@ $<

te_integrate.o: te_integrate.f90 te_greens.o tm_greens.o bessel_functions.o
	$(FC) $(FFLAGS) -c -o $@ $<

# Link main_te_integrate executable
# Note: order matters - dependent modules must come after the files that use them
main_te_integrate: te_integrate.o bessel_functions.o tm_greens.o te_greens.o multilayer_reflectance.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Clean up
clean:
	rm -f *.o *.mod main_te_integrate *.dat

# Clean intermediate files but keep executables
clean-obj:
	rm -f *.o *.mod
