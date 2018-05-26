##############################################################################
# Compiler definitions for MacBook (gfortran from MacPorts)
##############################################################################

FC       = gfortran-mp-7
FFLDCOMM = -fopenmp
OPT      = -O2 -Wall -funroll-loops $(FFLDCOMM)
# standards check/enforcement -std=f95 or -std=f2003, -fall-intrinsics
# allows use of GNU-extension functions/subroutines
FFLAGS   = $(OPT) -std=f2003
F77FLAGS = $(OPT)
LDFLAGS  = $(FFLDCOMM)

