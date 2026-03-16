# Fortran SLICOT benchmark driver.
# Run from project root: make -C scripts/slicot-fortran/bench -f Makefile SLICOT_REF=../../SLICOT-Reference
# Or from this dir: make SLICOT_REF=../../../SLICOT-Reference
# Needs: slicot.a, lpkaux.a. Links with -lblas -llapack (or -lopenblas if BLASLIB set).

SHELL = /bin/sh
# When run via make -C scripts/slicot-fortran/bench, CURDIR is bench dir
BENCH_DIR := $(CURDIR)
TOP       := $(BENCH_DIR)/../../..
SLICOT_REF ?= $(TOP)/SLICOT-Reference
FORTRAN = gfortran
OPTS    = -O2 -fPIC -g
# Prefer OpenBLAS if available; fallback to ref BLAS/LAPACK
BLASLIB   ?= -lopenblas
LAPACKLIB ?= -lopenblas
# Use ref BLAS/LAPACK: make BLASLIB=-lblas LAPACKLIB=-llapack

BENCH_SRC = bench_slicot.f90
BENCH_EXE = bench_slicot

.PHONY: all run clean

all: $(BENCH_EXE)

$(BENCH_EXE): $(BENCH_SRC)
	$(FORTRAN) $(OPTS) -o $(BENCH_EXE) $(BENCH_SRC) \
	  $(SLICOT_REF)/slicot.a $(SLICOT_REF)/lpkaux.a $(LAPACKLIB) $(BLASLIB)

run: $(BENCH_EXE)
	./$(BENCH_EXE)

clean:
	rm -f $(BENCH_EXE) *.o *.mod
