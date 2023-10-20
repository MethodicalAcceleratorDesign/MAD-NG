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
#include <stdbool.h>

// --- types -----------------------------------------------------------------o

typedef bool             log_t;
typedef int32_t          idx_t;
typedef int32_t          ssz_t;
typedef uint32_t         u32_t;
typedef uint64_t         u64_t;
typedef double           num_t;
typedef double _Complex  cpx_t;
typedef const char*      str_t;
typedef const void*      ptr_t;

// --- constants -------------------------------------------------------------o

#define TRUE  1
#define FALSE 0

// --- macros ----------------------------------------------------------------o

#define SQR(a)           ((a)*(a))
#define CUB(a)           ((a)*(a)*(a))
#define SWAP(a,b,t)      ((t)=(a), (a)=(b), (b)=(t))
#define SIGN(a)          (((a) >  0) - ((a) < 0)) // -1, 0, 1
#define SIGN1(a)         (((a) >= 0) - ((a) < 0)) // -1, 1

#define MIN(a,...)       MKNAME(MIN_,NARG(__VA_ARGS__))(a,__VA_ARGS__)
#define MAX(a,...)       MKNAME(MAX_,NARG(__VA_ARGS__))(a,__VA_ARGS__)
#define MIN_1(a,b)       ((b)<(a) ? b : a)
#define MAX_1(a,b)       ((b)>(a) ? b : a)
#define MIN_2(a,b,c)     ((b)<(a) ? MIN_1(b,c) : MIN_1(a,c))
#define MAX_2(a,b,c)     ((b)>(a) ? MAX_1(b,c) : MAX_1(a,c))
#define MIN_3(a,b,c,d)   ((b)<(a) ? MIN_2(b,c,d) : MIN_2(a,c,d))
#define MAX_3(a,b,c,d)   ((b)>(a) ? MAX_2(b,c,d) : MAX_2(a,c,d))

#define FOR(i,...)       MKNAME(FOR_,NARG(__VA_ARGS__))(i,__VA_ARGS__)
#define FOR_1(i,n)       for (idx_t i=  0 ; i<(n); i++)
#define FOR_2(i,i0,n)    for (idx_t i=(i0); i<(n); i++)
#define FOR_3(i,i0,n,s)  for (idx_t i=(i0); i<(n); i+=(s))

#define RFOR(i,...)      MKNAME(RFOR_,NARG(__VA_ARGS__))(i,__VA_ARGS__)
#define RFOR_1(i,n)      for (idx_t i=(n)-1; i>=  0 ; i--)
#define RFOR_2(i,n,i0)   for (idx_t i=(n)-1; i>=(i0); i--)
#define RFOR_3(i,n,i0,s) for (idx_t i=(n)-1; i>=(i0); i-=(s))

#define MKSTR(...)       MKSTR_OP_(__VA_ARGS__)
#define MKSTR_OP_(...)   #__VA_ARGS__

#define MKNAME(a,b)      MKNAME_OP_(a,b)
#define MKNAME3(a,b,c)   MKNAME(a,MKNAME(b,c))
#define MKNAME4(a,b,c,d) MKNAME(a,MKNAME(b,MKNAME(c,d)))
#define MKNAME_OP_(a,b)  a##b

// my old (2006) famous PP_NARG... for more features consider P99.
// https://gitlab.inria.fr/gustedt/p99
// https://gustedt.gitlabpages.inria.fr/p99/p99-html/index.html

#define NARG(...)       NARG_VCAT_(__VA_ARGS__,NARG_RSEQ_())
#define NARG_VCAT_(...) NARG_NSEQ_(__VA_ARGS__)
#define NARG_NSEQ_( \
     _1, _2, _3, _4, _5, _6, _7, _8, _9,_10, \
    _11,_12,_13,_14,_15,_16,_17,_18,_19,_20, \
    _21,_22,_23,_24,_25,_26,_27,_28,_29,_30, \
    _31,_32,_33,_34,_35,_36,_37,_38,_39,_40, \
    _41,_42,_43,_44,_45,_46,_47,_48,_49,_50, \
    _51,_52,_53,_54,_55,_56,_57,_58,_59,_60, \
    _61,_62,_63,  N, ...) N
#define NARG_RSEQ_() \
    63,62,61,60,                   \
    59,58,57,56,55,54,53,52,51,50, \
    49,48,47,46,45,44,43,42,41,40, \
    39,38,37,36,35,34,33,32,31,30, \
    29,28,27,26,25,24,23,22,21,20, \
    19,18,17,16,15,14,13,12,11,10, \
     9, 8, 7, 6, 5, 4, 3, 2, 1, 0

// --- GNU C -----------------------------------------------------------------o

#ifndef __GNUC__
#define __attribute__(a)
// #define expect(a,v)     (a)
// #define expect_p(a,v,p) (a)
// #define unreachable()
// #else
// #define unreachable()   __builtin_unreachable()
// #define expect(a,v)     __builtin_expect((a),v) // questionnable effect...
// #if __GNUC__ >= 9
// #define expect_p(a,v,p) __builtin_expect_with_probability((a),v,p)
// #else
// #define expect_p(a,v,p) ((p) >= 0.9 ? __builtin_expect((a),v) : (a))
// #endif
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
