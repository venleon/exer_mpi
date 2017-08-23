# ============================================================================
# Name        : Makefile
# Author      : Venleon
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello MPI World in Fortran
# ============================================================================

USEMPI=yes
ifdef USEMPI
	program = exer8.mpi
else
	program = exer8.serial
endif
all: $(program)

OBJS:=  param.o precision_module.o timing_module.o \
	xmpi.o readkey.o params.o mnemonic.o \
	general_mpi.o spaceparams.o \
	sub.o main.o

MODS = precision_module.mod timing_module.mod \
       xmpi.mod readkey.mod mnemonic.mod \
       general_mpi.mod spaceparams.mod

SRCS = $(patsubst %.o,%.f95,$(OBJS))

$(program): $(OBJS)
	$(F95) $(F95FLAGS) -o $@ $(OBJS)

install: all
	mkdir -p ./bin
	cp $(program) ./bin

%.o: %.f95
	@echo compiling $<
	$(F95) -c $(F95FLAGS) $<

$(MODS): $(patsubst %.mod,%.f95,$(MODS))
	@echo compiling $< to create $@
	$(F95) -c $(F95FLAGS) $@

F95FLAGS += -DHAVE_MPI_WTIME
F95FLAGS += -g -cpp
F95FLAGS += -fbacktrace -Wall -Wextra


F95:= gfortran
ifdef USEMPI
  ifdef USEXLF
	F95:=mpifort
  else
	F95:=mpif90
  endif
  F95FLAGS += -DUSEMPI
endif


clean:
	rm -f *.o *.mod $(program)
