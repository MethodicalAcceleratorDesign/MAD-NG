#ifndef MAD_SSE_H
#define MAD_SSE_H

/*
 o----------------------------------------------------------------------------o
 |
 | SSE2 & AVX2 & AVX512 interface
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
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
 */

#include <immintrin.h>

// --- macros ----------------------------------------------------------------o

// -- SSE2 ---

// char
#define MAD_SSE2_CSIZ 16
#define MAD_SSE2_CMSK            (MAD_SSE2_CSIZ-1)
#define MAD_SSE2_CRND(n) ((n) & ~(MAD_SSE2_CSIZ-1))
#define MAD_SSE2_CMOD(n) ((n) &  (MAD_SSE2_CSIZ-1))

// int
#define MAD_SSE2_ISIZ 4
#define MAD_SSE2_IMSK            (MAD_SSE2_ISIZ-1)
#define MAD_SSE2_IRND(n) ((n) & ~(MAD_SSE2_ISIZ-1))
#define MAD_SSE2_IMOD(n) ((n) &  (MAD_SSE2_ISIZ-1))

// double
#define MAD_SSE2_DSIZ 2
#define MAD_SSE2_DMSK            (MAD_SSE2_DSIZ-1)
#define MAD_SSE2_DRND(n) ((n) & ~(MAD_SSE2_DSIZ-1))
#define MAD_SSE2_DMOD(n) ((n) &  (MAD_SSE2_DSIZ-1))

#ifdef __SSE3__
#define _mm_loadu_si128(a) _mm_lddqu_si128(a)
#endif

// -- AVX2 ---

// char
#define MAD_AVX2_CSIZ 32
#define MAD_AVX2_CMSK            (MAD_AVX2_CSIZ-1)
#define MAD_AVX2_CRND(n) ((n) & ~(MAD_AVX2_CSIZ-1))
#define MAD_AVX2_CMOD(n) ((n) &  (MAD_AVX2_CSIZ-1))

// int
#define MAD_AVX2_ISIZ 8
#define MAD_AVX2_IMSK            (MAD_AVX2_ISIZ-1)
#define MAD_AVX2_IRND(n) ((n) & ~(MAD_AVX2_ISIZ-1))
#define MAD_AVX2_IMOD(n) ((n) &  (MAD_AVX2_ISIZ-1))

// double
#define MAD_AVX2_DSIZ 4
#define MAD_AVX2_DMSK            (MAD_AVX2_DSIZ-1)
#define MAD_AVX2_DRND(n) ((n) & ~(MAD_AVX2_DSIZ-1))
#define MAD_AVX2_DMOD(n) ((n) &  (MAD_AVX2_DSIZ-1))

// -- AVX512 ---

// char
#define MAD_AVX512_CSIZ 64
#define MAD_AVX512_CMSK            (MAD_AVX512_CSIZ-1)
#define MAD_AVX512_CRND(n) ((n) & ~(MAD_AVX512_CSIZ-1))
#define MAD_AVX512_CMOD(n) ((n) &  (MAD_AVX512_CSIZ-1))

// int
#define MAD_AVX512_ISIZ 16
#define MAD_AVX512_IMSK            (MAD_AVX512_ISIZ-1)
#define MAD_AVX512_IRND(n) ((n) & ~(MAD_AVX512_ISIZ-1))
#define MAD_AVX512_IMOD(n) ((n) &  (MAD_AVX512_ISIZ-1))

// double
#define MAD_AVX512_DSIZ 8
#define MAD_AVX512_DMSK            (MAD_AVX512_DSIZ-1)
#define MAD_AVX512_DRND(n) ((n) & ~(MAD_AVX512_DSIZ-1))
#define MAD_AVX512_DMOD(n) ((n) &  (MAD_AVX512_DSIZ-1))

// --- globals ---------------------------------------------------------------o

extern const unsigned char mad_sse2_msk1[16][16];
extern const unsigned char mad_sse2_msk2[16][16];

extern const unsigned char mad_avx2_msk1[32][32];
extern const unsigned char mad_avx2_msk2[32][32];

// ---------------------------------------------------------------------------o

#endif // MAD_SSE_H
