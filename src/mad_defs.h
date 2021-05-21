#ifndef MAD_H
#define MAD_H

/*
 o-----------------------------------------------------------------------------o
 |
 | MAD definitions
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  Purpose:
  - provide some global and portable definitions commonly used in MAD

  Information:
  - formal parameters ending with an underscope can be null (optional).

 o-----------------------------------------------------------------------------o
 */

// --- includes --------------------------------------------------------------o

#include <stddef.h>
#include <stdint.h>

// --- types -----------------------------------------------------------------o

typedef const char*      str_t;
typedef const void*      ptr_t;
typedef _Bool            log_t;
typedef int32_t          idx_t;
typedef int32_t          ssz_t;
typedef uint32_t         u32_t;
typedef uint64_t         u64_t;
typedef double           num_t;
typedef double _Complex cnum_t;

// --- constants -------------------------------------------------------------o

#define TRUE  1
#define FALSE 0

// --- macros ----------------------------------------------------------------o

#define SQR(a)        ((a)*(a))
#define CUB(a)        ((a)*(a)*(a))
#define MIN(a,b)      ((b)<(a) ? (b):(a))
#define MAX(a,b)      ((b)>(a) ? (b):(a))
#define MIN3(a,b,c)   ((b)<(a) ? MIN(b,c):MIN(a,c))
#define MAX3(a,b,c)   ((b)>(a) ? MAX(b,c):MAX(a,c))
#define MIN4(a,b,c,d) ((b)<(a) ? MIN3(b,c,d):MIN3(a,c,d))
#define MAX4(a,b,c,d) ((b)>(a) ? MAX3(b,c,d):MAX3(a,c,d))
#define SWAP(a,b,t)   ((t)=(a), (a)=(b), (b)=(t))
#define SIGN(a)       (((a) >  0) - ((a) < 0)) // -1, 0, 1
#define SIGN1(a)      (((a) >= 0) - ((a) < 0)) // -1, 1

#define MKSTR(...)     MKSTR_OP_(__VA_ARGS__)
#define MKSTR_OP_(...) #__VA_ARGS__

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
#define __thread
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

// --- special assert (loop) --------------------------------------------------o

// #define DBGASSERT

#ifdef DBGASSERT
#undef  assert
#define assert(c) \
        ((void)( (c) || (__assert_fail(#c, __FILE__, __LINE__, __func__),1) ))

void __assert_fail(const char *assertion, const char *file, int line,
                   const char *function) __attribute__((noreturn));
#endif

// ---------------------------------------------------------------------------o

#endif // MAD_H
