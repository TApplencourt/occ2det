IRPF90 = irpf90  #-a -d
FC     = ifort
FCFLAGS= -O2 -xHost -ffree-line-length-none -I .
NINJA  = ninja
AR     = ar
RANLIB = ranlib

SRC=
OBJ=
LIB=

include irpf90.make
export

irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile 
	$(IRPF90)


run: occ2det
	./occ2det 0 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1
