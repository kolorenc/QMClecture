##############################################################################
# Compiler definitions for MacBook
##############################################################################

FC       = gfortran-mp-4.6
FFLDCOMM = -fopenmp
OPT      = -O2 -Wall $(FFLDCOMM)
# standards check/enforcement -std=f95 or -std=f2003, -fall-intrinsics
# allows use of GNU-extension functions/subroutines
FFLAGS   = $(OPT) -std=f2003
F77FLAGS = $(OPT)
LDFLAGS  = $(FFLDCOMM)

#LIBS=-framework veclib

