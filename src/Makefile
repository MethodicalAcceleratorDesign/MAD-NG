#
# o----------------------------------------------------------------------------o
# |
# | MAD Makefile
# |
# | Methodical Accelerator Design (Copyleft 2015+)
# | Authors: L. Deniau, laurent.deniau at cern.ch
# | Contrib: -
# |
# o----------------------------------------------------------------------------o
# | You can redistribute this file and/or modify it under the terms of the GNU
# | General Public License GPLv3 (or later), as published by the Free Software
# | Foundation. This file is distributed in the hope that it will be useful, but
# | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
# o----------------------------------------------------------------------------o
#

CC      := gcc-mp-4.8
LJ      := ../bin/lj2-macosx

CFLAGS  := -std=c11 -W -Wall -Wextra -pedantic
CFLAGS  += -O3 -flto -ffast-math -ftree-vectorize -fPIC $(OMP)
CFLAGS  += -I../lib/luajit/OSX/include/luajit-2.1 

LJFLAGS :=

LDFLAGS := -pagezero_size 10000 -image_base 100000000
LDFLAGS += -static-libgfortran # -static-libgcc 
LDFLAGS += -Wl,-macosx_version_min,10.8
LDFLAGS += ../lib/luajit/src/libluajit.a \
           ../lib/lapack/liblapack.a \
           ../lib/lapack/librefblas.a \
           /opt/local/lib/gcc48/libgfortran.a \
           /opt/local/lib/gcc48/libquadmath.a

# OMP := -fopenmp

CSRC    := $(filter-out %_test.c,$(wildcard *.c))
LSRC    := # $(filter-out %_test.lua,$(wildcard *.lua))
OBJ     := $(patsubst %.c,%.o,$(CSRC))
OBJ     += $(patsubst %.lua,%.o,$(LSRC))

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.lua
	$(LJ) $(LJFLAGS) -b $< $@

mad: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

dep:
	otool -L mad

clean:
	rm -f mad

cleanall: clean
	rm -f *.o