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

#define MIN(a,b)    ((b)<(a) ? (b):(a))
#define MAX(a,b)    ((b)>(a) ? (b):(a))
#define MIN3(a,b,c)  MIN(a,MIN(b,c))
#define MAX3(a,b,c)  MAX(a,MAX(b,c))
#define SWAP(a,b,t) ((t)=(a), (a)=(b), (b)=(t))

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

// --- MATH Constants --------------------------------------------------------o

#ifndef M_E
#define M_E           2.7182818284590452353602874713526625  /* e */
#define M_LOG2E       1.4426950408889634073599246810018921  /* log_2 e */
#define M_LOG10E      0.4342944819032518276511289189166051  /* log_10 e */
#define M_LN2         0.6931471805599453094172321214581766  /* log_e 2 */
#define M_LN10        2.3025850929940456840179914546843642  /* log_e 10 */
#define M_PI          3.1415926535897932384626433832795029  /* pi */
#define M_PI_2        1.5707963267948966192313216916397514  /* pi/2 */
#define M_PI_4        0.7853981633974483096156608458198757  /* pi/4 */
#define M_1_PI        0.3183098861837906715377675267450287  /* 1/pi */
#define M_2_PI        0.6366197723675813430755350534900574  /* 2/pi */
#define M_2_SQRTPI    1.1283791670955125738961589031215452  /* 2/sqrt(pi) */
#define M_SQRT2       1.4142135623730950488016887242096981  /* sqrt(2) */
#define M_SQRT1_2     0.7071067811865475244008443621048490  /* 1/sqrt(2) */
#endif

// ---------------------------------------------------------------------------o

#endif // MAD_H
