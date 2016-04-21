/*
 o----------------------------------------------------------------------------o
 |
 | Memory module test suite
 |
 | Methodical Accelerator Design (Copyleft 2015+)
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
*/

/*
  single-threaded versions:
  gcc -std=c11 -Wall -W -pedantic -O3 -flto -ffast-math -ftree-vectorize ../mad_mem.c ../mad_log.c mad_mem_test.c -o mad_mem_test
  // icc -std=c11 -Wall -O3 -static-intel -static-libstdc++ ../mad_mem.c ../mad_log.c mad_mem_test.c -o mad_mem_test

  multi-threaded versions:
  gcc -std=c11 -Wall -W -pedantic -O3 -flto -ffast-math -ftree-vectorize -fopenmp ../mad_mem.c ../mad_log.c mad_mem_test.c -o mad_mem_test
  // icc -std=c11 -Wall -O3 -openmp -openmp-link -static-intel -static-libstdc++ ../mad_mem.c ../mad_log.c mad_mem_test.c -o mad_mem_test
  // export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/usr/bin/icc-14.0-base/compiler/lib

  usage:
  mad_mem_test object_size total_loops [inner_loops] 

  example:
  ./mad_mem_test 5 200000000 10

  performance (gcc 4.8.2):
  ~ 200M mad_malloc+mad_free/second           (single-threaded) on 2.3Ghz i5
  ~  90M mad_malloc+mad_free/second/thread    ( multi-threaded) on 2.3Ghz i5
  -DMAD_MEM_AUTOCOLLECT slows down by about 10-20%
  
  performance (gcc 4.8.2):
  ~  31M malloc+free/second                   (single-threaded) on 2.3Ghz i5
  ~   3M malloc+free/second/thread            ( multi-threaded) on 2.3Ghz i5
*/

#include <stdio.h>
#include <stdlib.h>

#include "../mad_mem.h"

static long object_size = 5;
static long total_loops = 200000000;
static long inner_loops = 10;

void
mad_mem_test(void)
{
  long outer_loops = total_loops/inner_loops;
  fprintf(stderr, "[%d] object_size = %ld bytes, outer_loops = %ld, inner_loops = %ld, total_loops = %ld\n",
          omp_get_thread_num(), object_size, outer_loops, inner_loops, outer_loops*inner_loops);
  fprintf(stderr, "[%d] cached = %5zu bytes (initial)\n", omp_get_thread_num(), mad_mem_cached());

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int n=0; n < omp_get_num_threads(); n++) {
    void *s[inner_loops];

    if (omp_get_num_threads() > 1)
      printf("thread no %d\n", omp_get_thread_num());

    for (long i=0; i < outer_loops; i++) {
      for (long j=0; j < inner_loops; j++)
        s[j] = malloc(object_size);
  
      for (long j=0; j < inner_loops; j++)
        free(s[j]);
    }

    fprintf(stderr, "[%d] cached = %5zu bytes (before collect)\n", omp_get_thread_num(), mad_mem_cached());
    mad_mem_collect();
    fprintf(stderr, "[%d] cached = %5zu bytes (after  collect)\n", omp_get_thread_num(), mad_mem_cached());
  }
}

#ifndef MAD_MEM_TEST_NOMAIN
int
main(int argc, char *argv[])
{
  if (argc<3) {
    printf( "usage:\n  mem_alloc object_size total_loops [inner_loops] \n"
            "example:\n  mem_alloc 5 200000000 10\n");
    return EXIT_FAILURE;
  }

  // retrieve arguments
  object_size = strtol(argv[1],0,0);
  inner_loops = argc > 3 ? strtol(argv[3],0,0) : 100;
  total_loops = strtol(argv[2],0,0);

  if (object_size < 0) object_size = 1;
  if (total_loops < 1) total_loops = 1;

  mad_mem_test();
}
#endif
