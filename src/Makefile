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

# setup
CC  := gcc
LJ  := ../bin/luajit
DIR := build

# depend
DFLAGS  := -std=c11 -MM 

# compiler
CFLAGS  := -std=c11 -W -Wall -Wextra -pedantic
CFLAGS  += -O3 -ffast-math -ftree-vectorize -fPIC
CFLAGS  += # -flto
CFLAGS  += # -fopenmp

# lua
LFLAGS  :=

# linker
LDFLAGS := -pagezero_size 10000 -image_base 100000000
LDFLAGS += -static-libgfortran # -static-libgcc -> crash LuaJIT
LDFLAGS += -Wl,-macosx_version_min,10.8
LDFLAGS += ../lib/luajit/src/libluajit.a \
           ../lib/lapack/liblapack.a \
           ../lib/lapack/librefblas.a \
           -lgfortran $(shell gcc -print-file-name=libquadmath.a)

# files setup
CSRC    := $(wildcard *.c sse/*.c)
LSRC    := # $(wildcard *.lua) # compile also lua modules

DEP     := $(patsubst %.c,%.d,$(CSRC))
DEP     := $(addprefix $(DIR)/,$(DEP))

OBJ     := $(patsubst %.c,%.o,$(CSRC))
OBJ     += $(patsubst %.lua,%.o,$(LSRC))
OBJ     := $(addprefix $(DIR)/,$(OBJ))

DOIT    := $(shell mkdir -p $(DIR) $(DIR)/sse)

# rules
mad: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

$(DIR)/%.d: %.c
	$(CC) $(DFLAGS) $< > $@

$(DIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(DIR)/%.o: %.lua
	$(LJ) $(LFLAGS) -b $< $@

dep:
	otool -L mad

clean:
	rm -f mad

cleanall: clean
	rm -rf $(DIR)

# files specific setup
$(DIR)/mad_main.o: CFLAGS += -I../lib/luajit/OSX/include/luajit-2.1

# include dependencies
ifeq ($(filter-out mad,$(MAKECMDGOALS)),)
-include $(DEP)
endif
