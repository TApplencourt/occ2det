#IRPF90 = ~/irpf90/bin/irpf90 --codelet single_index:100000000 --codelet gen_dets:500000 #--codelet bench_v:5000 #-a -d
IRPF90 = ~/irpf90/bin/irpf90 --codelet single_index:100000000 --codelet gen_dets:5000 #--codelet bench_v:5000 #-a -d
FC     = gfortran
FCFLAGS= -O2 -g -march=native -ffree-line-length-none -I .
#FC     = ifort
#FCFLAGS= -O2 -g -xHost -I .
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

bench_single: codelet_single_index
	./codelet_single_index 0 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1

bench: codelet_gen_dets check
	./codelet_gen_dets 0 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1

check: occ2det
	./occ2det 1 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1 | head -100 > out
	diff out REF

bench_v: codelet_bench_v
	./codelet_bench_v 0 36 15 1 1 2 2 1 1 1 1 1 1 1 0 2 1 0 1 2 1 1 1 0 2 1 1 2 1 1 2 1 1 1 1 2 1 0 1

