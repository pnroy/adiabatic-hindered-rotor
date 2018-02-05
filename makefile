.SUFFIXES: .o .c .f .cxx
options= -O3 -DLINUX -DBLAS
CC=g++
cc=gcc
f77=gfortran
# *******objects for Hui's Lanczos code*******
# objects= matvec.o ran1.o four1.o prod.o gaulegf.o wangpot.o wanghehe.o\
# plgndr.o inter.o  gauleg.o sturmi.o inverr.o inverm.o isoev.o \
# lancbis.o bisec.o genran.o scalar.o hv.o trivec.o peckeris.o
# *******objects for direct diag code*******
objects= matvec.o prod.o  gaulegf.o four1.o \
plgndr.o inter.o  gauleg.o genran.o scalar.o cmdstuff.o 
H2Rb: H2Rb.o $(objects) matvec.o 
	$(CC) $(options) -o H2Rb H2Rb.o  $(objects) \
-lblas -llapack -lgfortran
clean:
	rm *.o
.cxx.o: 
	$(CC) -c  $(options) $<
.f.o:
	$(f77) -c  $<
.c.o:
	$(cc) -c $(options) $<
