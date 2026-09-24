FC = gfortran
FFLAGS = -O2

nbody: nbody.f90
	$(FC) $(FFLAGS) nbody.f90 -o nbody

clean:
	rm -f nbody *.mod *_osc.dat *_geo.dat *_osc.png *_geo.png *_comparacao.png
