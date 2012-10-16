#####################################################################
# Fortran compiler from AMD 
#####################################################################

FC=/opt/open64/bin/openf95

OPT=-O3 -mp -m64# -xcheck=%all
#OPT=-g -ftrap=%all -dalign
FFLAGS  =$(OPT)
F77FLAGS=$(OPT)
LDFLAGS =-mp -m64

# ACML 
#LIBS=-L/opt/acml/open64_64/lib -lacml -Wl,-rpath=/opt/acml/open64_64/lib 

