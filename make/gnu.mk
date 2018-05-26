##############################################################################
# Compiler definitions for gfortran
##############################################################################

FC       = gfortran
FFLDCOMM = -fopenmp
# -flto -fwhole-program -fno-protect-parens  (not in 4.4.5)
# -ffast-math -funroll-loops
OPT      = -O2 -march=native -Wall -funroll-loops $(FFLDCOMM)
FFLAGS   = $(OPT) -std=f2003
# no warnings about obsolete Fortran constructs
F77FLAGS = $(OPT)
LDFLAGS  = $(FFLDCOMM)

