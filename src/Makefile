INC=-I /usr/lib/x86_64-linux-gnu/openmpi/include
LIB=  #-L/usr/lib/x86_64-linux-gnu -lplplotd -lplplotf95d
INC=
FC=gfortran
FCPARA=mpif90
FCOPT= -fimplicit-none -finit-local-zero  -fopenmp -finit-local-zero # -Wunused -Wall -Wargument-mismatch -Wimplicit-interface
SRC=global.f90  param.f90 tools.f90 IO.f90 poten.f90  numerov2.f90  mesh.f90 davidson.f90   time_tracking.f90 ConjugateGradient.f90 tdse_mod.f90 operation.f90 Hbinitio.f90 

EXEC=Hbinitio.x
OBJ=time_tracking.o global.o poten.o IO.o param.o tools.o  mesh.o davidson.o ConjugateGradient.o numerov2.o  pseudopotential.o operation.o tdse_mod.o
all: $(OBJ) Hbinitio.f90
	$(FC) $(FCOPT)  $(OBJ) Hbinitio.f90  -o $(EXEC)  -lblas -llapack	
#serial:  $(SRC)
#	$(FC) $(FCOPT)  $(SRC) -o Hbinitio.x  -lblas -llapack

time_tracking.o: time_tracking.f90
	$(FC) $(FCOPT) $< -c
global.o: global.f90
	$(FC) $(FCOPT) $< -c
poten.o: poten.f90 global.o IO.o
	$(FC) $(FCOPT) $< -c
IO.o: IO.f90 global.o tools.o
	$(FC) $(FCOPT) $< -c
tools.o: tools.f90 param.o global.o mesh.o
	$(FC) $(FCOPT) $< -c
param.o: param.f90  global.o
	$(FC) $(FCOPT) $< -c
mesh.o: mesh.f90  global.o
	$(FC) $(FCOPT) $< -c
davidson.o: davidson.f90  global.o time_tracking.o mesh.o poten.o param.o
	$(FC) $(FCOPT) $< -c
numerov2.o: numerov2.f90  global.o mesh.o poten.o ConjugateGradient.o
	$(FC) $(FCOPT) $< -c
ConjugateGradient.o: ConjugateGradient.f90  global.o IO.o mesh.o
	$(FC) $(FCOPT) $< -c
operation.o: operation.f90 IO.o global.o mesh.o 
	$(FC) $(FCOPT) $< -c
pseudopotential.o: pseudopotential.f90
	$(FC) $(FCOPT) $< -c
tdse_mod.o: tdse_mod.f90 IO.o 
	$(FC) $(FCOPT) $< -c

clean:
	rm *.mod *.o 
install:
	mv $(EXEC) ../bin

#param.mod: param.f90 global.o
#	$(FC) $(FCOPT) $< -c



#IO.mod: IO.f90 global.o tools.o
#	$(FC) $(FCOPT) $< -c



#	$(FCPARA) $(FCOPT) Hbinitio.f90 -o Hbinitio.x $(INC) -lblas -llapack -lmpi
#pp: ppHbinitio.f90#
#	$(FC) ppHbinitio.f90 -o ppHbinitio.x $(LIB)
#parallel: mpi_Hbinitio.f90
#	$(FCPARA) $(FCOPT) mpi_Hbinitio.f90 -o mpi_Hbinitio.x $(INC) -lblas -llapack -lmpi
#dbg: dbg.f90
#	$(FCPARA) dbg.f90 -o dbg.x  $(INC) -lblas -llapack -lmpi
