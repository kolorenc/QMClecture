##############################################################################
# Compiler definitions for pc190b.fzu.cz (Intel compiler 11.1)
##############################################################################

#FC       = /opt/gcc-4.7.2/bin/gfortran-4.7
FC       = gfortran
FFLDCOMM = -fopenmp
# -flto -fwhole-program -fno-protect-parens  (not in 4.4.5)
# -ffast-math -funroll-loops
OPT      = -O3 -march=native -Wall -funroll-loops $(FFLDCOMM)
FFLAGS   = $(OPT) -std=f2003
# no warnings about obsolete Fortran constructs
F77FLAGS = $(OPT)
LDFLAGS  = $(FFLDCOMM)

