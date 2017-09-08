#ifndef MAD_SSE_H
#define MAD_SSE_H

/*
 o----------------------------------------------------------------------------o
 |
 | SSE2 interface
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
 */

#include <immintrin.h>

// --- macros ----------------------------------------------------------------o

// char
#define MAD_SSE_CSIZ 16
#define MAD_SSE_CMSK            (MAD_SSE_CSIZ-1)
#define MAD_SSE_CRND(n) ((n) & ~(MAD_SSE_CSIZ-1))
#define MAD_SSE_CMOD(n) ((n) &  (MAD_SSE_CSIZ-1))

// int
#define MAD_SSE_ISIZ 4
#define MAD_SSE_IMSK            (MAD_SSE_ISIZ-1)
#define MAD_SSE_IRND(n) ((n) & ~(MAD_SSE_ISIZ-1))
#define MAD_SSE_IMOD(n) ((n) &  (MAD_SSE_ISIZ-1))

// double
#define MAD_SSE_DSIZ 2
#define MAD_SSE_DMSK            (MAD_SSE_DSIZ-1)
#define MAD_SSE_DRND(n) ((n) & ~(MAD_SSE_DSIZ-1))
#define MAD_SSE_DMOD(n) ((n) &  (MAD_SSE_DSIZ-1))

#ifdef __SSE3__
#define _mm_loadu_si128(a) _mm_lddqu_si128(a)
#endif

// --- globals ---------------------------------------------------------------o

extern const unsigned char mad_sse_msk1[16][16];
extern const unsigned char mad_sse_msk2[16][16];

// ---------------------------------------------------------------------------o

#endif // MAD_SSE_H
