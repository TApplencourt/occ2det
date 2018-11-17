IRPF90 = irpf90 --codelet single_index:100000000 --codelet gen_dets:500 #-a -d
FC     = ifort 
FCFLAGS= -O3 -g -xHost -I .
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


run: occ2det codelet_single_index
	./occ2det 0 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1

bench: codelet_single_index
	./codelet_single_index 0 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1

bench2: codelet_gen_dets
	./codelet_gen_dets 0 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1
