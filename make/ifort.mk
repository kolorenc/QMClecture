##############################################################################
# Compiler definitions for Intel compiler (ifort)
##############################################################################

FC       = ifort
FFLDCOMM = -openmp -openmp-report=2# -p -check all
OPT      = -O2 -pc80 -ip -xHost -traceback# -check all# -fpe0
OPT     += -fp-model strict -prec-div -prec-sqrt -no-ftz -assume protect_parens $(FFLDCOMM)
# -stand f03 , -stand f95 , -stand f90 ... to check compliance with standards
# note that -warn by default also makes those stupid files *__genmod.f90
FFLAGS   = $(OPT) -warn -stand f03
# no warnings about obsolete Fortran constructs
F77FLAGS = $(OPT) -warn nousage
LDFLAGS  = $(FFLDCOMM)

