LIBNAME = libcmad.so

SRC = $(wildcard mad_*.c)
# mad_cnum.c mad_vec.c mad_mat.c

COPT = -O3 -ffast-math -ftree-vectorize

CWARNS = -W -Wall -Wextra -pedantic \
	-Waggregate-return \
	-Wcast-align \
	-Wcast-qual \
	-Wdisabled-optimization \
	-Wpointer-arith \
	-Wshadow \
	-Wsign-compare \
	-Wundef \
	-Wwrite-strings \
	-Wbad-function-cast \
	-Wmissing-prototypes \
	-Wnested-externs \
	-Wstrict-prototypes \
    -Wunreachable-code \


CFLAGS = $(CWARNS) $(COPT) -fPIC
CC = gcc -std=c99

OBJ = $(SRC:.c=.o)

# For Linux
linux:
	make libcmad.so "LFLAGS = -shared -static-libgcc"

# For Mac OS
macosx:
	make libcmad.so "LFLAGS = -shared -static-libgcc"
#	make libcmad.so "DLLFLAGS = -bundle -undefined dynamic_lookup"

windows:
	echo "NYI, bye"

libcmad.so: $(OBJ)
	env $(CC) $(LFLAGS) $(OBJ) -o $(LIBNAME)

clean:
	rm -f $(OBJ) libcmad*.so

mad_cnum.o: mad_cnum.c mad_cnum.h
mad_vec.o:  mad_vec.c mad_vec.h
mad_mat.o:  mad_mat.c mad_mat.h

