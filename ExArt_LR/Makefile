# Makefile for the program ExArt 1D
# Questions or suggestions are always welcome: garesu@oa-cagliari.inaf.it

# Compiler
FC = gfortran
# Flags
FCFLAGS = -fbounds-check
FCFLAGS = -O2

PROGRAM = ExArt_LR

all: $(PROGRAM)

ExArt_LR.o: mod.o 
ExArt_LR: mod.o DISORT/DISORT.o DISORT/DISORTblas.o DISORT/DISORTlapack.o DISORT/ErrPack.o DISORT/RDI1MACH.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^

# Clean Up
	rm -f DISORT/*.o DISORT/*.mod 
	rm -f *.o *.mod
