PLATFORM := default
include make/$(PLATFORM).mk

SYS    := paraHe orthoHe H2
SYSRUN := $(addsuffix _run, $(SYS))
SYSOBJ := $(addsuffix _input.o, $(SYS))

all: $(SYS)

paraHe: cleansys
	make -f makefile.sys SYS=$@

orthoHe: cleansys
	rm -f qmc_run.o qmc.o $(SYSOBJ)
	make -f makefile.sys SYS=$@

H2: cleansys
	rm -f qmc_run.o qmc.o $(SYSOBJ)
	make -f makefile.sys SYS=$@

cleansys:
	rm -f qmc_run.o qmc.o $(SYSOBJ)

clean: cleansys
	make -f makefile.sys SYS=paraHe clean

mrproper: clean
	rm -f $(SYSRUN)
	make -f makefile.sys SYS=paraHe mrproper
