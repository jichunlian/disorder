FC = ifort

all: disorder supercell

disorder: disorder.f90 stdout.o functions.o structure.o symmetry.o progress.o configurations.o outfiles.o
	$(FC) -o disorder disorder.f90 *.o
	cp disorder ../bin/

supercell: supercell.f90 structure.o functions.o
	$(FC) -o supercell supercell.f90 structure.o functions.o
	cp supercell ../bin/

functions.o: functions.f90
	$(FC) -c functions.f90

structure.o: structure.f90 functions.o
	$(FC) -c structure.f90

symmetry.o: symmetry.f90 functions.o stdout.o
	$(FC) -c symmetry.f90

progress.o: progress.f90
	$(FC) -c progress.f90

configurations.o: configurations.f90 progress.o stdout.o functions.o
	$(FC) -c configurations.f90

outfiles.o: outfiles.f90 functions.o
	$(FC) -c outfiles.f90

stdout.o: stdout.f90 functions.o
	$(FC) -c stdout.f90

clean:
	rm *.mod *.o
