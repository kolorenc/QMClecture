##############################################################################
# Compiler definitions for dorje.fzu.cz (Intel compiler 11.1)
##############################################################################

FC       = ifort
FFLDCOMM = -openmp -openmp-report=2
OPT      = -O2 -pc80 -ip $(FFLDCOMM)
# -stand f03 , -stand f95 , -stand f90 ... to check compliance with standards
# note that -warn by default also makes those stupid files *__genmod.f90
FFLAGS   = $(OPT) -warn -stand f03
# no warnings about obsolete Fortran constructs
F77FLAGS = $(OPT) -warn nousage
LDFLAGS  = $(FFLDCOMM)

#MKLPATH = /opt/intel/Compiler/11.1/064/mkl/lib/em64t
#LIBS    = -Wl,--start-group \$(MKLPATH)/libmkl_intel_lp64.a \$(MKLPATH)/libmkl_sequential.a \$(MKLPATH)/libmkl_core.a -Wl,--end-group -lguide

