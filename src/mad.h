#ifndef MAD_H
#define MAD_H

/*
 o----------------------------------------------------------------------------o
 |
 | MAD definitions
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
  
  Purpose:
  - provide some global and portable definitions commonly used in MAD

  Information:
  - formal parameters ending with an underscope can be null (optional).

 o----------------------------------------------------------------------------o
 */

// --- includes --------------------------------------------------------------o

#include <stddef.h>

// --- types -----------------------------------------------------------------o

typedef const char*      str_t;
typedef int              idx_t;
typedef double           num_t;
typedef double _Complex cnum_t;

// --- constants -------------------------------------------------------------o

#define TRUE  1
#define FALSE 0

// --- macros ----------------------------------------------------------------o

#define MIN(a,b)    ((a)<(b) ? (a):(b))
#define MAX(a,b)    ((a)>(b) ? (a):(b))
#define MIN3(a,b,c)  MIN(a,MIN(b,c))
#define MAX3(a,b,c)  MAX(a,MAX(b,c))
#define SWAP(a,b,t) ((t)=(a), (a)=(b), (b)=(t))

#define MKSTR(a)     MKSTR_OP_(a)
#define MKSTR_OP_(a) #a

#define MKNAME(a,b)      MKNAME_OP_(a,b)
#define MKNAME3(a,b,c)   MKNAME(a,MKNAME(b,c))
#define MKNAME4(a,b,c,d) MKNAME(a,MKNAME(b,MKNAME(c,d)))
#define MKNAME_OP_(a,b)  a##b

// --- GNU C -----------------------------------------------------------------o

#ifndef __GNUC__
#define __attribute__(a)
#endif

// --- Open Multi-Processing -------------------------------------------------o

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_procs()   1
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#define omp_get_thread_num()  0
#endif

// --- POSIX & WIN -----------------------------------------------------------o

#if !defined(_WIN32) && !defined(_WIN64) && (defined(__unix__) || \
     defined(__unix) || (defined(__APPLE__) && defined(__MACH__)))

#define POSIX_VERSION _POSIX_VERSION
/* 198808L for     POSIX.1-1988
   199009L for     POSIX.1-1990
   199506L for ISO POSIX.1-1996
   200112L for ISO POSIX.1-2001
   200809L for ISO POSIX.1-2008 */

#elif defined(_WIN64)
#define WIN_VERSION 64

#elif defined(_WIN32)
#define WIN_VERSION 32

#else
#error "unsupported platform"
#endif

// ---------------------------------------------------------------------------o

#endif // MAD_H
