#####################################################################
# Fortran compiler from Sun/Oracle Solaris Studio, relatively little
# support for F2003
#####################################################################

FC       =sunf95
FFLDCOMM =-openmp -traceback# -pg 
# Sun Performance Library is compiled with -dalign, program that calls it
# should use that option too (from SPL User's Guide)
OPT=-O4 -xtarget=native -m64 -dalign -w3 $(FFLDCOMM)# -xcheck=%all
#OPT=-g -ftrap=%all -dalign
FFLAGS  =$(OPT) -ansi
F77FLAGS=$(OPT)
LDFLAGS =$(FFLDCOMM)

# BLAS/LAPACK (Sun Performance Library)
#LIBS=-xlic_lib=sunperf

