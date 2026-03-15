# Top-level SLICOT Fortran build: lpkaux -> slicot -> example drivers
# Invoke from project root: make -C SLICOT-Reference -f scripts/slicot-fortran/Makefile
# Prerequisite: gfortran, system OpenBLAS. See docs/FORTRAN_BUILD.md

SHELL = /bin/sh
FORTRAN = gfortran
OPTS = -O2 -fPIC -g
ARCH = ar
ARCHFLAGS = cr

LPKAUX_DIR := src/lapack_aux
LPKAUX_SRCS := $(wildcard $(LPKAUX_DIR)/*.f)
LPKAUX_OBJS := $(patsubst $(LPKAUX_DIR)/%.f,build/lpkaux/%.o,$(LPKAUX_SRCS))

.PHONY: all lpkaux.a slicot examples clean

all: lpkaux.a slicot examples

lpkaux.a: $(LPKAUX_OBJS)
	@mkdir -p .
	$(ARCH) $(ARCHFLAGS) $@ $^

build/lpkaux/%.o: $(LPKAUX_DIR)/%.f
	@mkdir -p build/lpkaux
	$(FORTRAN) $(OPTS) -c $< -o $@

slicot: lpkaux.a
	$(MAKE) -C src -f makefile_Unix

examples: slicot
	$(MAKE) -C examples -f makefile_Unix

clean:
	rm -f lpkaux.a slicot.a
	rm -rf build
	$(MAKE) -C src -f makefile_Unix clean
	$(MAKE) -C examples -f makefile_Unix clean
