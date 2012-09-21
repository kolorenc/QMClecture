##############################################################################
# Compiler definitions for pc190b.fzu.cz (Intel compiler)
##############################################################################

FC       = ifort
FFLDCOMM = -openmp -openmp-report=2# -check all
OPT      = -O2 -pc80 -ip -xHost $(FFLDCOMM)
# -stand f03 , -stand f95 , -stand f90 ... to check compliance with standards
# note that -warn by default also makes those stupid files *__genmod.f90
FFLAGS   = $(OPT) -warn -stand f03
# no warnings about obsolete Fortran constructs
F77FLAGS = $(OPT) -warn nousage
LDFLAGS  = $(FFLDCOMM)

#MKLPATH = /opt/intel/mkl/lib/intel64
#LIBS=-L$(MKLPATH) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

