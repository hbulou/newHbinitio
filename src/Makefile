INC=-I /usr/lib/x86_64-linux-gnu/openmpi/include
LIB=  #-L/usr/lib/x86_64-linux-gnu -lplplotd -lplplotf95d
INC=
FC=gfortran
FCPARA=mpif90
FCOPT= -fimplicit-none -finit-local-zero  -fopenmp -finit-local-zero # -Wunused -Wall -Wargument-mismatch -Wimplicit-interface
SRC=global.f90  FFT.f90 param.f90 tools.f90 IO.f90 poten.f90  numerov2.f90  mesh.f90 experiment.f90 davidson.f90   time_tracking.f90 ConjugateGradient.f90 tdse_mod.f90 operation.f90 Hbinitio.f90 

EXEC=Hbinitio.x
OBJ=time_tracking.o global.o FFT.o poten.o  IO.o tools.o numerov2.o   param.o   mesh.o  davidson.o ConjugateGradient.o experiment.o pseudopotential.o operation.o tdse_mod.o  

all: Hbinitio.f90 
	$(FC) $(FCOPT) global.f90 -c
	$(FC) -I/usr/include $(FCOPT) FFT.f90 -c	
	$(FC) $(FCOPT) time_tracking.f90 -c
	$(FC) $(FCOPT) poten.f90 -c
	$(FC) $(FCOPT) mesh.f90 -c
	$(FC) $(FCOPT) tools.f90 -c
	$(FC) $(FCOPT) IO.f90 -c
	$(FC) $(FCOPT) experiment.f90 -c
	$(FC) $(FCOPT) davidson.f90 -c	
	$(FC) $(FCOPT) tdse_mod.f90 -c
	$(FC) $(FCOPT) ConjugateGradient.f90 -c
	$(FC) $(FCOPT) numerov2.f90 -c
	$(FC) $(FCOPT) param.f90 -c
	$(FC) $(FCOPT) pseudopotential.f90 -c
	$(FC) $(FCOPT) operation.f90 -c
	$(FC) $(FCOPT)  $(OBJ) Hbinitio.f90  -o $(EXEC)  -lblas -llapack -lfftw3 

clean:
	rm *.mod *.o 
install:
	mv $(EXEC) ../bin

